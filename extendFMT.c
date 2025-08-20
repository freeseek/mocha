/* The MIT License

   Copyright (C) 2018-2025 Giulio Genovese

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

#include <stdio.h>
#include <unistd.h>
#include <getopt.h>
#include <errno.h>
#include <htslib/vcf.h>
#include <htslib/synced_bcf_reader.h>
#include "mocha.h"
#include "bcftools.h"
#include "rbuf.h"

#define EXTENDFMT_VERSION "2025-08-19"

/******************************************
 * CIRCULAR BUFFER                        *
 ******************************************/

// most of this code was inspired by Petr Danecek's code in bcftools/vcfbuf.c

typedef struct {
    bcf1_t *line; // BCF record
    void *vals;   // auxiliary array for values from the format to be extended
    int *dist;    // distances to closest predetermined value
} data_t;

typedef struct {
    int win;
    int nsmpl; // number of samples
    int gt_id, fmt_id;
    int type; // whether BCF_HT_INT or BCF_HT_REAL values
    int size; // size of elements in the auxiliary vals array
    int8_t *phase_arr;
    data_t *data;
    rbuf_t rbuf;
} auxbuf_t;

static auxbuf_t *auxbuf_init(int win, int nsmpl, int type, int fmt_id, int gt_id) {
    auxbuf_t *buf = (auxbuf_t *)calloc(1, sizeof(auxbuf_t));
    buf->win = win;
    buf->nsmpl = nsmpl;
    buf->type = type;
    switch (type) {
    case BCF_HT_INT:
        buf->size = sizeof(int32_t);
        break;
    case BCF_HT_REAL:
        buf->size = sizeof(float);
        break;
    default:
        error("Unexpected type %d", type);
    }
    buf->fmt_id = fmt_id;
    buf->gt_id = gt_id;
    buf->phase_arr = (int8_t *)malloc(buf->nsmpl * sizeof(int8_t));
    rbuf_init(&buf->rbuf, 0);
    return buf;
}

static void auxbuf_destroy(auxbuf_t *buf) {
    int i;
    for (i = 0; i < buf->rbuf.m; i++) {
        if (buf->data[i].line) bcf_destroy(buf->data[i].line);
        if (buf->data[i].vals) free(buf->data[i].vals);
        if (buf->data[i].dist) free(buf->data[i].dist);
    }
    free(buf->data);
    free(buf->phase_arr);
    free(buf);
}

static inline int bcf_int8_is_missing(int8_t value) { return value == bcf_int8_missing; }
static inline int bcf_int16_is_missing(int16_t value) { return value == bcf_int16_missing; }
static inline int bcf_int32_is_missing(int32_t value) { return value == bcf_int32_missing; }
#define bcf_int8_set_missing(x) x = bcf_int8_missing
#define bcf_int16_set_missing(x) x = bcf_int16_missing
#define bcf_int32_set_missing(x) x = bcf_int32_missing

// push a new record into the buffer
static void auxbuf_push(auxbuf_t *buf, bcf1_t *line) {
    rbuf_expand0(&buf->rbuf, data_t, buf->rbuf.n + 1, buf->data);
    int k, curr = rbuf_append(&buf->rbuf);

    // if a record is already present in the buffer, destroy it
    if (buf->data[curr].line) bcf_destroy(buf->data[curr].line);
    buf->data[curr].line = bcf_dup(line);
    // allocate auxiliary arrays if they have not been previously allocated (this minimizes
    // necessary allocations)
    if (!buf->data[curr].vals) buf->data[curr].vals = malloc(buf->nsmpl * buf->size);
    if (!buf->data[curr].dist) buf->data[curr].dist = (int *)malloc(buf->nsmpl * sizeof(int));
    memset(buf->data[curr].dist, 0, buf->nsmpl * sizeof(int));

    int phase =
        (buf->gt_id < 0) ? 0 : bcf_get_genotype_phase(bcf_get_fmt_id(line, buf->gt_id), buf->phase_arr, buf->nsmpl);

    bcf_fmt_t *fmt = bcf_get_fmt_id(line, buf->fmt_id);
    if (fmt && fmt->n != 1) error("Format vector has incorrect number of values\n");

    // extract information from BCF record into an auxiliary array
    if (fmt) {
#define BRANCH(ht_type_t, bt_type_t, is_missing, set_missing)                                                          \
    {                                                                                                                  \
        ht_type_t *vals = (ht_type_t *)buf->data[curr].vals;                                                           \
        bt_type_t *p = (bt_type_t *)fmt->p;                                                                            \
        for (k = 0; k < buf->nsmpl; k++) {                                                                             \
            if (is_missing(p[k])) {                                                                                    \
                set_missing(vals[k]);                                                                                  \
            } else {                                                                                                   \
                vals[k] = p[k];                                                                                        \
                if (vals[k] && phase) {                                                                                \
                    if (buf->phase_arr[k] == bcf_int8_missing || buf->phase_arr[k] == bcf_int8_vector_end)             \
                        vals[k] = (ht_type_t)0;                                                                        \
                    else                                                                                               \
                        vals[k] *= buf->phase_arr[k];                                                                  \
                }                                                                                                      \
            }                                                                                                          \
        }                                                                                                              \
    }
        if (buf->type == BCF_HT_INT && fmt->type == BCF_BT_INT8) {
            BRANCH(int32_t, int8_t, bcf_int8_is_missing, bcf_int32_set_missing);
        } else if (buf->type == BCF_HT_INT && fmt->type == BCF_BT_INT16) {
            BRANCH(int32_t, int16_t, bcf_int16_is_missing, bcf_int32_set_missing);
        } else if (buf->type == BCF_HT_INT && fmt->type == BCF_BT_INT32) {
            BRANCH(int32_t, int32_t, bcf_int32_is_missing, bcf_int32_set_missing);
        } else if (buf->type == BCF_HT_REAL && fmt->type == BCF_BT_FLOAT) {
            BRANCH(float, float, bcf_float_is_missing, bcf_float_set_missing);
        } else {
            error("Unexpected type combination %d %d\n", buf->type, fmt->type);
        }
#undef BRANCH
    } else {
        memset(buf->data[curr].vals, 0, buf->nsmpl * buf->size);
    }

    int prev = curr;
    rbuf_prev(&buf->rbuf, &prev);
    int prev_dist = buf->data[curr].line->pos - buf->data[prev].line->pos;

// propagate information backwards
#define BRANCH(ht_type_t, is_missing)                                                                                  \
    {                                                                                                                  \
        for (k = 0; k < buf->nsmpl; k++) {                                                                             \
            ht_type_t *curr_vals = (ht_type_t *)buf->data[curr].vals;                                                  \
            if (curr_vals[k] && !is_missing(curr_vals[k])) {                                                           \
                int i = curr;                                                                                          \
                while (rbuf_prev(&buf->rbuf, &i)) {                                                                    \
                    ht_type_t *vals = (ht_type_t *)buf->data[i].vals;                                                  \
                    int i_dist = buf->data[curr].dist[k] + buf->data[curr].line->pos - buf->data[i].line->pos;         \
                    if (i_dist <= buf->win                                                                             \
                        && (vals[k] == (ht_type_t)0 || is_missing(vals[k]) || i_dist < buf->data[i].dist[k])) {        \
                        vals[k] = curr_vals[k];                                                                        \
                        buf->data[i].dist[k] = i_dist;                                                                 \
                    }                                                                                                  \
                }                                                                                                      \
            } else if (prev != curr) {                                                                                 \
                ht_type_t *prev_vals = (ht_type_t *)buf->data[prev].vals;                                              \
                if (prev_vals[k] && !is_missing(prev_vals[k]) && buf->data[prev].dist[k] + prev_dist <= buf->win) {    \
                    curr_vals[k] = prev_vals[k];                                                                       \
                    buf->data[curr].dist[k] = buf->data[prev].dist[k] + prev_dist;                                     \
                }                                                                                                      \
            }                                                                                                          \
        }                                                                                                              \
    }
    if (buf->type == BCF_HT_INT) {
        BRANCH(int32_t, bcf_int32_is_missing);
    } else if (buf->type == BCF_HT_REAL) {
        BRANCH(float, bcf_float_is_missing);
    } else {
        error("Unexpected type %d\n", buf->type);
    }
#undef BRANCH
}

static data_t *auxbuf_flush(auxbuf_t *buf, int flush_all) {
    if (buf->rbuf.n == 0) return NULL;
    int first = rbuf_kth(&buf->rbuf, 0);
    int last = rbuf_last(&buf->rbuf);
    if (!flush_all && buf->data[last].line->pos - buf->data[first].line->pos <= buf->win) return NULL;

    int i = rbuf_shift(&buf->rbuf);
    bcf1_t *line = buf->data[i].line;
    int phase =
        (buf->gt_id < 0) ? 0 : bcf_get_genotype_phase(bcf_get_fmt_id(line, buf->gt_id), buf->phase_arr, buf->nsmpl);

    // fix the phase before returning the VCF record
    if (phase) {
        int k;
#define BRANCH(type_t, is_missing)                                                                                     \
    {                                                                                                                  \
        type_t *vals = (type_t *)buf->data[i].vals;                                                                    \
        for (k = 0; k < buf->nsmpl; k++) {                                                                             \
            if (is_missing(vals[k])) continue;                                                                         \
            if (buf->phase_arr[k] == bcf_int8_missing || buf->phase_arr[k] == bcf_int8_vector_end)                     \
                vals[k] = (type_t)0;                                                                                   \
            else                                                                                                       \
                vals[k] *= buf->phase_arr[k];                                                                          \
        }                                                                                                              \
    }
        if (buf->type == BCF_HT_INT) {
            BRANCH(int32_t, bcf_int32_is_missing);
        } else if (buf->type == BCF_HT_REAL) {
            BRANCH(float, bcf_float_is_missing);
        } else {
            error("Unexpected type %d\n", buf->type);
        }
#undef BRANCH
    }

    return &buf->data[i];
}

/******************************************
 * PLUGIN                                 *
 ******************************************/

const char *about(void) { return "Extend format fields to nearby variants.\n"; }

static const char *usage_text(void) {
    return "\n"
           "About: Extend format fields to nearby variants. (version " EXTENDFMT_VERSION
           " http://github.com/freeseek/mocha)\n"
           "Usage: bcftools +extendFMT [options] --format <ID> <in.vcf.gz>\n"
           "\n"
           "Plugin options:\n"
           "    -f, --format <tag>              FORMAT tag to be extended\n"
           "    -p, --phase                     whether the format to be extended is for phased heterozygotes\n"
           "    -d, --dist <int>                maximum distance used to extend the calls [1e6]\n"
           "        --no-version                do not append version and command line to the header\n"
           "    -o, --output <file>             write output to a file [standard output]\n"
           "    -O, --output-type u|b|v|z[0-9]  u/b: un/compressed BCF, v/z: un/compressed VCF, 0-9: compression level "
           "[v]\n"
           "    -r, --regions <region>          restrict to comma-separated list of regions\n"
           "    -R, --regions-file <file>       restrict to regions listed in a file\n"
           "        --regions-overlap 0|1|2     Include if POS in the region (0), record overlaps (1), variant "
           "overlaps (2) [1]\n"
           "    -t, --targets [^]<region>       similar to -r but streams rather than index-jumps. Exclude regions "
           "with \"^\" prefix\n"
           "    -T, --targets-file [^]<file>    similar to -R but streams rather than index-jumps. Exclude regions "
           "with \"^\" prefix\n"
           "        --targets-overlap 0|1|2     Include if POS in the region (0), record overlaps (1), variant "
           "overlaps (2) [0]\n"
           "        --threads <int>             number of extra output compression threads [0]\n"
           "    -s, --samples [^]<list>         comma separated list of samples to include (or exclude with \"^\" "
           "prefix)\n"
           "    -S, --samples-file [^]<file>    file of samples to include (or exclude with \"^\" prefix)\n"
           "        --force-samples             only warn about unknown subset samples\n"
           "    -W, --write-index[=FMT]         Automatically index the output files [off]\n"
           "\n"
           "Example:\n"
           "    bcftools +extendFMT --format AS --phase --dist 500000 file.bcf\n"
           "\n";
}

// code from bcf_hdr_register_hrec() in vcf.c as I did not understand a better way to get the
// format type
static uint32_t bcf_hdr_get_type(bcf_hrec_t *hrec) {
    if (!strcmp(hrec->key, "INFO") && !strcmp(hrec->key, "FILTER") && !strcmp(hrec->key, "FORMAT"))
        error("Header record %s=%s is not INFO/FILTER/FORMAT\n", hrec->key, hrec->value);

    int i;
    for (i = 0; i < hrec->nkeys; i++) {
        if (!strcmp(hrec->keys[i], "Type")) {
            if (!strcmp(hrec->vals[i], "Integer")) {
                return BCF_HT_INT;
            } else if (!strcmp(hrec->vals[i], "Float")) {
                return BCF_HT_REAL;
            } else if (!strcmp(hrec->vals[i], "String")) {
                return BCF_HT_STR;
            } else if (!strcmp(hrec->vals[i], "Character")) {
                return BCF_HT_STR;
            } else if (!strcmp(hrec->vals[i], "Flag")) {
                return BCF_HT_FLAG;
            } else {
                hts_log_warning("The type \"%s\" is not supported, assuming \"String\"", hrec->vals[i]);
                return BCF_HT_STR;
            }
        }
    }
    error("Header record is missing the Type key\n");
}

static void flush(auxbuf_t *buf, htsFile *out_fh, bcf_hdr_t *hdr, char *format, int flush_all) {
    data_t *data;
    while ((data = auxbuf_flush(buf, flush_all))) {
        bcf_update_format(hdr, data->line, format, data->vals, buf->nsmpl, buf->type);
        if (bcf_write(out_fh, hdr, data->line) < 0) error("Unable to write to output VCF file\n");
    }
}

int run(int argc, char **argv) {
    char *format = NULL;
    int phase_format = 0;
    int dist = 1e6;
    char *output_fname = NULL;
    char *index_fname;
    int output_type = FT_VCF;
    int regions_overlap = 1;
    int targets_overlap = 0;
    int clevel = -1;
    int n_threads = 0;
    int record_cmd_line = 1;
    int write_index = 0;
    char *targets_list = NULL;
    int targets_is_file = 0;
    char *regions_list = NULL;
    int regions_is_file = 0;
    char *sample_names = NULL;
    int sample_is_file = 0;
    int force_samples = 0;

    static struct option loptions[] = {
        {"format", required_argument, NULL, 'f'},        {"phase", no_argument, NULL, 'p'},
        {"dist", required_argument, NULL, 'd'},          {"output", required_argument, NULL, 'o'},
        {"output-type", required_argument, NULL, 'O'},   {"threads", required_argument, NULL, 9},
        {"regions", required_argument, NULL, 'r'},       {"regions-file", required_argument, NULL, 'R'},
        {"regions-overlap", required_argument, NULL, 2}, {"targets", required_argument, NULL, 't'},
        {"targets-file", required_argument, NULL, 'T'},  {"targets-overlap", required_argument, NULL, 3},
        {"samples", required_argument, NULL, 's'},       {"samples-file", required_argument, NULL, 'S'},
        {"force-samples", no_argument, NULL, 1},         {"no-version", no_argument, NULL, 8},
        {"write-index", optional_argument, NULL, 'W'},   {0, 0, 0, 0}};
    int c;
    char *tmp;
    while ((c = getopt_long(argc, argv, "h?f:pd:o:O:r:R:t:T:s:S:W::", loptions, NULL)) >= 0) {
        switch (c) {
        case 'f':
            format = optarg;
            break;
        case 'p':
            phase_format = 1;
            break;
        case 'd':
            dist = (int)strtol(optarg, &tmp, 0);
            if (*tmp) error("Could not parse: --dist %s\n", optarg);
            if (dist <= 0) error("Distance to extend calls needs to be positive: --dist %s\n", optarg);
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
            default: {
                clevel = strtol(optarg, &tmp, 10);
                if (*tmp || clevel < 0 || clevel > 9) error("The output type \"%s\" not recognised\n", optarg);
            }
            }
            if (optarg[1]) {
                clevel = strtol(optarg + 1, &tmp, 10);
                if (*tmp || clevel < 0 || clevel > 9)
                    error("Could not parse argument: --compression-level %s\n", optarg + 1);
            }
            break;
        case 9:
            n_threads = strtol(optarg, &tmp, 0);
            if (*tmp) error("Could not parse argument: --threads %s\n", optarg);
            break;
        case 8:
            record_cmd_line = 0;
            break;
        case 'r':
            regions_list = optarg;
            break;
        case 'R':
            regions_list = optarg;
            regions_is_file = 1;
            break;
        case 2:
            if (!strcasecmp(optarg, "0"))
                regions_overlap = 0;
            else if (!strcasecmp(optarg, "1"))
                regions_overlap = 1;
            else if (!strcasecmp(optarg, "2"))
                regions_overlap = 2;
            else
                error("Could not parse: --regions-overlap %s\n", optarg);
            break;
        case 't':
            targets_list = optarg;
            break;
        case 'T':
            targets_list = optarg;
            targets_is_file = 1;
            break;
        case 3:
            if (!strcasecmp(optarg, "0"))
                targets_overlap = 0;
            else if (!strcasecmp(optarg, "1"))
                targets_overlap = 1;
            else if (!strcasecmp(optarg, "2"))
                targets_overlap = 2;
            else
                error("Could not parse: --targets-overlap %s\n", optarg);
            break;
        case 's':
            sample_names = optarg;
            break;
        case 'S':
            sample_names = optarg;
            sample_is_file = 1;
            break;
        case 1:
            force_samples = 1;
            break;
        case 'W':
            if (!(write_index = write_index_parse(optarg))) error("Unsupported index format '%s'\n", optarg);
            break;
        case 'h':
        case '?':
        default:
            error("%s", usage_text());
            break;
        }
    }
    if (!format) error("Format ID needs to be specified with --format\n%s", usage_text());

    char *input_fname = NULL;
    if (optind == argc) {
        if (!isatty(fileno((FILE *)stdin))) {
            input_fname = "-"; // reading from stdin
        } else {
            error("%s", usage_text());
        }
    } else if (optind + 1 != argc) {
        error("%s", usage_text());
    } else {
        input_fname = argv[optind];
    }

    bcf_srs_t *srs = bcf_sr_init();
    bcf_sr_set_opt(srs, BCF_SR_PAIR_LOGIC, BCF_SR_PAIR_EXACT);
    if (regions_list) {
        bcf_sr_set_opt(srs, BCF_SR_REGIONS_OVERLAP, regions_overlap);
        if (bcf_sr_set_regions(srs, regions_list, regions_is_file) < 0)
            error("Failed to read the regions: %s\n", regions_list);
    }
    if (targets_list) {
        bcf_sr_set_opt(srs, BCF_SR_TARGETS_OVERLAP, targets_overlap);
        if (bcf_sr_set_targets(srs, targets_list, targets_is_file, 0) < 0)
            error("Failed to read the targets: %s\n", targets_list);
    }
    if (bcf_sr_set_threads(srs, n_threads) < 0) error("Failed to create threads\n");
    if (!bcf_sr_add_reader(srs, input_fname))
        error("Failed to open %s: %s\n", input_fname, bcf_sr_strerror(srs->errnum));

    bcf_hdr_t *hdr = bcf_sr_get_header(srs, 0);
    if (record_cmd_line) bcf_hdr_append_version(hdr, argc, argv, "bcftools_plugin");

    if (sample_names) {
        int ret = bcf_hdr_set_samples(hdr, sample_names, sample_is_file);
        if (ret < 0)
            error("Error parsing the list of samples: %s\n", sample_names);
        else if (force_samples && ret > 0)
            error("Sample name mismatch: sample #%d not found in the header\n", ret);
    }

    // get the type of format BCF_HT_INT/BCF_HT_REAL
    bcf_hrec_t *hrec = bcf_hdr_get_hrec(hdr, BCF_HL_FMT, NULL, format, NULL);
    int format_type = bcf_hdr_get_type(hrec);
    if (format_type != BCF_HT_INT && format_type != BCF_HT_REAL)
        error("TODO: %s:%d .. type=%d\n", __FILE__, __LINE__, format_type);

    int gt_id = phase_format ? bcf_hdr_id2int(hdr, BCF_DT_ID, "GT") : -1;
    if (phase_format && gt_id < 0) error("Format GT was not found in the input header\n");
    int fmt_id = bcf_hdr_id2int(hdr, BCF_DT_ID, format);
    if (fmt_id < 0) error("Format %s was not found in the input header\n", format);

    char wmode[8];
    set_wmode(wmode, output_type, (char *)output_fname, clevel);
    htsFile *out_fh = hts_open(output_fname ? output_fname : "-", wmode);
    if (out_fh == NULL) error("[%s] Error: cannot write to \"%s\": %s\n", __func__, output_fname, strerror(errno));
    if (n_threads) hts_set_opt(out_fh, HTS_OPT_THREAD_POOL, srs->p);
    if (bcf_hdr_write(out_fh, hdr) < 0) error("Unable to write to output VCF file\n");
    if (init_index2(out_fh, hdr, output_fname, &index_fname, write_index) < 0)
        error("Error: failed to initialise index for %s\n", output_fname);

    auxbuf_t *buf = auxbuf_init(dist, bcf_hdr_nsamples(hdr), format_type, fmt_id, gt_id);
    int prev_rid = -1;
    while (bcf_sr_next_line(srs)) {
        bcf1_t *line = bcf_sr_get_line(srs, 0);
        if (prev_rid != line->rid) flush(buf, out_fh, hdr, format, 1);
        auxbuf_push(buf, line);
        flush(buf, out_fh, hdr, format, 0);
        prev_rid = line->rid;
    }
    flush(buf, out_fh, hdr, format, 1);
    auxbuf_destroy(buf);

    if (write_index) {
        if (bcf_idx_save(out_fh) < 0) {
            if (hts_close(out_fh) != 0) error("Close failed %s\n", strcmp(output_fname, "-") ? output_fname : "stdout");
            error("Error: cannot write to index %s\n", index_fname);
        }
        free(index_fname);
    }
    hts_close(out_fh);
    bcf_sr_destroy(srs);

    return 0;
}
