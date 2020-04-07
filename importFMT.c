/* The MIT License

   Copyright (C) 2018-2020 Giulio Genovese

   Author: Giulio Genovese <giulio.genovese@gmail.com>

   Permission is hereby granted, free of charge, to any person obtaining a copy
   of this software and associated documentation files (the "Software"), to deal
   in the Software without restriction, including without limitation the rights
   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
   copies of the Software, and to permit persons to whom the Software is
   furnished to do so, subject to the following conditions:

   The above copyright notice and this permission notice shall be included in
   all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
   THE SOFTWARE.

 */

#include <getopt.h>
#include <errno.h>
#include <htslib/vcf.h>
#include <htslib/synced_bcf_reader.h>
#include "bcftools.h"

#define IMPORTFMT_VERSION "2020-04-07"

const char *about(void)
{
	return "Import format fields from a source VCF file into a target VCF file\n";
}

static const char *usage_text(void)
{
	return "\n"
	       "About: Import format fields from a source VCF into a target VCF. (version " IMPORTFMT_VERSION
	       " https://github.com/freeseek/mocha)\n"
	       "Usage: bcftools +importFMT [options] --formats <list> <target.vcf.gz> <source.vcf.gz>\n"
	       "\n"
	       "Plugin options:\n"
	       "    -f, --formats <list>               list of formats to be imported from the source VCF\n"
	       "        --no-version                   do not append version and command line to the header\n"
	       "    -o, --output <file>                write output to a file [standard output]\n"
	       "    -O, --output-type b|u|z|v          b: compressed BCF, u: uncompressed BCF, z: compressed VCF, v: uncompressed VCF [v]\n"
	       "    -r, --regions <region>             restrict to comma-separated list of regions\n"
	       "    -R, --regions-file <file>          restrict to regions listed in a file\n"
	       "    -t, --targets [^]<region>          similar to -r but streams rather than index-jumps. Exclude regions with \"^\" prefix\n"
	       "    -T, --targets-file [^]<file>       similar to -R but streams rather than index-jumps. Exclude regions with \"^\" prefix\n"
	       "        --threads <int>                number of extra output compression threads [0]\n"
	       "\n"
	       "Example:\n"
	       "    bcftools +importFMT --formats Ldev,Bdev,Bdev_Phase target.dose.bcf source.mocha.bcf\n"
	       "\n";
}

// code from bcf_hdr_register_hrec() in vcf.c as I did not understand a better way to get the
// format type
static uint32_t bcf_hdr_get_type(bcf_hrec_t *hrec)
{
	if (!strcmp(hrec->key, "INFO") && !strcmp(hrec->key, "FILTER")
	    && !strcmp(hrec->key, "FORMAT"))
		error("Header record %s=%s is not INFO/FILTER/FORMAT\n", hrec->key,
		      hrec->value);

	for (int i = 0; i < hrec->nkeys; i++) {
		if (!strcmp(hrec->keys[i], "Type")) {
			if (!strcmp(hrec->vals[i], "Integer"))
				return BCF_HT_INT;
			else if (!strcmp(hrec->vals[i], "Float"))
				return BCF_HT_REAL;
			else if (!strcmp(hrec->vals[i], "String"))
				return BCF_HT_STR;
			else if (!strcmp(hrec->vals[i], "Character"))
				return BCF_HT_STR;
			else if (!strcmp(hrec->vals[i], "Flag"))
				return BCF_HT_FLAG;
			else {
				hts_log_warning(
					"The type \"%s\" is not supported, assuming \"String\"",
					hrec->vals[i]);
				return BCF_HT_STR;
			}
		}
	}
	error("Header record is missing the Type key\n");
}

int run(int argc, char *argv[])
{
	char *format_list = NULL;
	char *output_fname = NULL;
	int output_type = FT_VCF;
	int n_threads = 0;
	int record_cmd_line = 1;
	char *targets_list = NULL;
	int targets_is_file = 0;
	char *regions_list = NULL;
	int regions_is_file = 0;

	static struct option loptions[] = {{"formats", required_argument, NULL, 'f'},
					   {"output", required_argument, NULL, 'o'},
					   {"output-type", required_argument, NULL, 'O'},
					   {"no-version", no_argument, NULL, 8},
					   {"targets", required_argument, NULL, 't'},
					   {"targets-file", required_argument, NULL, 'T'},
					   {"regions", required_argument, NULL, 'r'},
					   {"regions-file", required_argument, NULL, 'R'},
					   {"threads", required_argument, NULL, 9},
					   {NULL, 0, NULL, 0}};
	int c;
	while ((c = getopt_long(argc, argv, "f:t:T:r:R:h?o:O:89", loptions, NULL)) >= 0) {
		switch (c) {
		case 'f':
			format_list = optarg;
			break;
		case 'o':
			output_fname = optarg;
			break;
		case 'O':
			switch (optarg[0]) {
			case 'b':
				output_type = FT_BCF_GZ;
				break;
			case 'u':
				output_type = FT_BCF;
				break;
			case 'z':
				output_type = FT_VCF_GZ;
				break;
			case 'v':
				output_type = FT_VCF;
				break;
			default:
				error("The output type \"%s\" not recognised\n", optarg);
			}
			break;
		case 9:
			n_threads = strtol(optarg, 0, 0);
			break;
		case 8:
			record_cmd_line = 0;
			break;

		case 't':
			targets_list = optarg;
			break;
		case 'T':
			targets_list = optarg;
			targets_is_file = 1;
			break;
		case 'r':
			regions_list = optarg;
			break;
		case 'R':
			regions_list = optarg;
			regions_is_file = 1;
			break;

		case 'h':
		case '?':
		default:
			error("%s", usage_text());
		}
	}
	if (argc - optind != 2)
		error("%s", usage_text());
	if (!format_list)
		error("Missing list of formats to be imported with --formats option\n%s",
		      usage_text());

	bcf_srs_t *srs = bcf_sr_init();
	bcf_sr_set_opt(srs, BCF_SR_REQUIRE_IDX);
	if (regions_list) {
		if (bcf_sr_set_regions(srs, regions_list, regions_is_file) < 0)
			error("Failed to read the regions: %s\n", regions_list);
	}
	if (targets_list) {
		if (bcf_sr_set_targets(srs, targets_list, targets_is_file, 0) < 0)
			error("Failed to read the targets: %s\n", targets_list);
	}
	if (bcf_sr_set_threads(srs, n_threads) < 0)
		error("Failed to create threads\n");
	while (optind < argc) {
		if (!bcf_sr_add_reader(srs, argv[optind]))
			error("Failed to open %s: %s\n", argv[optind],
			      bcf_sr_strerror(srs->errnum));
		optind++;
	}

	htsFile *out_fh =
		hts_open(output_fname ? output_fname : "-", hts_bcf_wmode(output_type));
	if (!out_fh)
		error("Can't write to \"%s\": %s\n", output_fname, strerror(errno));
	if (n_threads)
		hts_set_opt(out_fh, HTS_OPT_THREAD_POOL, srs->p);

	bcf_hdr_t *trgt_hdr = bcf_sr_get_header(srs, 0);
	bcf_hdr_t *src_hdr = bcf_sr_get_header(srs, 1);
	bcf_hdr_t *out_hdr = bcf_hdr_dup(trgt_hdr);

	// import source VCF format fields into target VCF header
	int moff = 0, *off = NULL;
	int ncols = ksplit_core(format_list, ',', &moff, &off);
	int *type = (int *)malloc(ncols * sizeof(int));
	for (int i = 0; i < ncols; i++) {
		bcf_hrec_t *hrec =
			bcf_hdr_get_hrec(src_hdr, BCF_HL_FMT, NULL, &format_list[off[i]], NULL);
		if (!hrec)
			error("Format %s not present in source VCF %s\n", &format_list[off[i]],
			      argv[optind - 1]);
		bcf_hdr_add_hrec(out_hdr, bcf_hrec_dup(hrec));
		type[i] = bcf_hdr_get_type(hrec);
		if (type[i] != BCF_HT_INT && type[i] != BCF_HT_REAL)
			error("TODO: %s:%d .. type=%d\n", __FILE__, __LINE__, type[i]);
	}
	if (record_cmd_line)
		bcf_hdr_append_version(out_hdr, argc, argv, "bcftools_+importFMT");
	if (bcf_hdr_write(out_fh, out_hdr) < 0)
		error("Unable to write to output VCF file\n");

	// generate a sample map from source header to target header
	int *imap = (int *)malloc(bcf_hdr_nsamples(trgt_hdr) * sizeof(int));
	for (int i = 0; i < bcf_hdr_nsamples(trgt_hdr); i++) {
		imap[i] = bcf_hdr_id2int(src_hdr, BCF_DT_SAMPLE,
					 bcf_hdr_int2id(trgt_hdr, BCF_DT_SAMPLE, i));
	}

	int msrc_arr = 0, mtrgt_arr = 0;
	void *src_arr = NULL, *trgt_arr = NULL;

	while (bcf_sr_next_line(srs)) {
		if (!bcf_sr_has_line(srs, 0))
			continue;
		bcf1_t *trgt_line = bcf_sr_get_line(srs, 0);
		if (bcf_sr_has_line(srs, 1)) {
			bcf1_t *src_line = bcf_sr_get_line(srs, 1);
			for (int i = 0; i < ncols; i++) {
				switch (type[i]) {
				case BCF_HT_INT: {
					int nsrc_arr = bcf_get_format_values(
						src_hdr, src_line, &format_list[off[i]],
						&src_arr, &msrc_arr, BCF_HT_INT);
					if (nsrc_arr > 0) {
						int32_t *src_int32_arr = (int32_t *)src_arr;
						int32_t *trgt_int32_arr = (int32_t *)trgt_arr;
						int num = nsrc_arr / bcf_hdr_nsamples(src_hdr);
						hts_expand(int32_t,
							   num * bcf_hdr_nsamples(trgt_hdr),
							   mtrgt_arr, trgt_int32_arr);
						for (int j = 0; j < bcf_hdr_nsamples(trgt_hdr);
						     j++) {
							for (int k = 0; k < num; k++) {
								trgt_int32_arr[num * j + k] =
									imap[j] < 0
										? bcf_int32_missing
										: src_int32_arr
											[num * imap[j]
											 + k];
							}
						}
						if (bcf_update_format_int32(
							    out_hdr, trgt_line,
							    &format_list[off[i]],
							    trgt_int32_arr,
							    num * bcf_hdr_nsamples(out_hdr))
						    < 0)
							error("Could not update %s format field\n",
							      &format_list[off[i]]);
						trgt_arr = (void *)trgt_int32_arr;
					}
					break;
				}
				case BCF_HT_REAL: {
					int nsrc_arr = bcf_get_format_values(
						src_hdr, src_line, &format_list[off[i]],
						&src_arr, &msrc_arr, BCF_HT_REAL);
					if (nsrc_arr > 0) {
						float *src_float_arr = (float *)src_arr;
						float *trgt_float_arr = (float *)trgt_arr;
						int num = nsrc_arr / bcf_hdr_nsamples(src_hdr);
						hts_expand(float,
							   num *bcf_hdr_nsamples(trgt_hdr),
							   mtrgt_arr, trgt_float_arr);
						for (int j = 0; j < bcf_hdr_nsamples(trgt_hdr);
						     j++) {
							for (int k = 0; k < num; k++) {
								trgt_float_arr[num * j + k] =
									imap[j] < 0
										? bcf_float_missing
										: src_float_arr
											[num * imap[j]
											 + k];
							}
						}
						if (bcf_update_format_float(
							    out_hdr, trgt_line,
							    &format_list[off[i]],
							    trgt_float_arr,
							    num * bcf_hdr_nsamples(out_hdr))
						    < 0)
							error("Could not update %s format field\n",
							      &format_list[off[i]]);
						trgt_arr = (void *)trgt_float_arr;
					}
					break;
				}
				case BCF_HT_STR:
				default:
					error("TODO: %s:%d .. type=%d\n", __FILE__, __LINE__,
					      type[i]);
				}
			}
		}
		if (bcf_write(out_fh, out_hdr, trgt_line) < 0)
			error("Unable to write to output VCF file\n");
	}

	free(src_arr);
	free(trgt_arr);
	free(imap);
	free(type);
	free(off);

	bcf_hdr_destroy(out_hdr);
	hts_close(out_fh);
	bcf_sr_destroy(srs);

	return 0;
}
