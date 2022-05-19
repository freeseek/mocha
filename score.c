/* The MIT License

   Copyright (C) 2021-2022 Giulio Genovese

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
#include <stdlib.h>
#include <getopt.h>
#include <dirent.h>
#include <sys/stat.h>
#include <htslib/vcf.h>
#include <htslib/synced_bcf_reader.h>
#include <htslib/kseq.h>
#include <htslib/khash_str2int.h>
#include "mocha.h" // tsv_init_delimiter() tsv_read_string() tsv_read_float() tsv_parse_delimiter()
#include "bcftools.h"
#include "tsv2vcf.h"
#include "filter.h"

#define SCORE_VERSION "2022-05-18"

#define FLT_INCLUDE (1 << 0)
#define FLT_EXCLUDE (1 << 1)

// ##FORMAT=<ID=GT,Number=1,Type=String,Description="Phased genotypes">
#define SCORE_GT 1
// ##FORMAT=<ID=DS,Number=A,Type=Float,Description="Genotype dosage">
#define SCORE_DS 2 // DS = AP1 + AP2
// ##FORMAT=<ID=HDS,Number=2,Type=Float,Description="Estimated Haploid Alternate Allele Dosage ">
#define SCORE_HDS 3
// ##FORMAT=<ID=AP1,Number=A,Type=Float,Description="ALT allele probability of first haplotype">
// ##FORMAT=<ID=AP2,Number=A,Type=Float,Description="ALT allele probability of second haplotype">
#define SCORE_AP 4
// ##FORMAT=<ID=GP,Number=G,Type=Float,Description="Estimated Genotype Probability">
#define SCORE_GP 5
// ##FORMAT=<ID=AS,Number=1,Type=Integer,Description="Allelic shift (1/-1 if the alternate allele is over/under
// represented)">
#define SCORE_AS 6

/****************************************
 * HELPER FUNCTIONS                     *
 ****************************************/

static inline char **get_file_list(const char *pathname, int *nfiles) {
    char **filenames = NULL;
    struct stat statbuf;
    if (stat(pathname, &statbuf) < 0) error("Can't open \"%s\": %s\n", pathname, strerror(errno));
    if (S_ISDIR(statbuf.st_mode)) {
        DIR *d = opendir(pathname);
        struct dirent *dir;
        int mfiles = 0;
        int p = strlen(pathname);
        while ((dir = readdir(d))) {
            if (strcmp(dir->d_name, ".") == 0 || strcmp(dir->d_name, "..") == 0) continue;
            hts_expand0(char *, *nfiles + 1, mfiles, filenames);
            int q = strlen(dir->d_name);
            filenames[*nfiles] = (char *)malloc((p + q + 2) * sizeof(char));
            memcpy(filenames[*nfiles], pathname, p);
            filenames[*nfiles][p] = '/';
            memcpy(filenames[*nfiles] + p + 1, dir->d_name, q + 1);
            (*nfiles)++;
        }
        closedir(d);
    } else {
        filenames = hts_readlines(pathname, nfiles);
        if (!filenames) error("Failed to read from file %s\n", pathname);
    }
    if (*nfiles == 0) error("No files found in %s\n", pathname);
    return filenames;
}

static inline int bcf_hdr_name2id_flexible(const bcf_hdr_t *hdr, char *chr) {
    if (!chr) return -1;
    char buf[] = {'c', 'h', 'r', '\0', '\0', '\0'};
    int rid = bcf_hdr_name2id(hdr, chr);
    if (rid >= 0) return rid;
    if (strncmp(chr, "chr", 3) == 0) rid = bcf_hdr_name2id(hdr, chr + 3);
    if (rid >= 0) return rid;
    if (strlen(chr) > 2) return -1;
    strcpy(buf + 3, chr);
    rid = bcf_hdr_name2id(hdr, buf);
    if (rid >= 0) return rid;
    if (strcmp(chr, "23") == 0 || strcmp(chr, "25") == 0 || strcmp(chr, "XY") == 0 || strcmp(chr, "XX") == 0) {
        rid = bcf_hdr_name2id(hdr, "X");
        if (rid >= 0) return rid;
        rid = bcf_hdr_name2id(hdr, "chrX");
    } else if (strcmp(chr, "24") == 0) {
        rid = bcf_hdr_name2id(hdr, "Y");
        if (rid >= 0) return rid;
        rid = bcf_hdr_name2id(hdr, "chrY");
    } else if (strcmp(chr, "26") == 0 || strcmp(chr, "MT") == 0 || strcmp(chr, "chrM") == 0) {
        rid = bcf_hdr_name2id(hdr, "MT");
        if (rid >= 0) return rid;
        rid = bcf_hdr_name2id(hdr, "chrM");
    }
    return rid;
}

static int tsv_setter_chrom_flexible(tsv_t *tsv, bcf1_t *rec, void *usr) {
    char tmp = *tsv->se;
    *tsv->se = 0;
    rec->rid = bcf_hdr_name2id_flexible((bcf_hdr_t *)usr, tsv->ss);
    *tsv->se = tmp;
    return 0;
}

/****************************************
 * PGS FILE IMPLEMENTATION              *
 ****************************************/

typedef struct {
    int32_t rid;
    hts_pos_t pos;
    char *a1;
    char *a2;
    float beta;
    float or ;
    float p;
} marker_t;

typedef struct {
    int use_snp;
    int snp_tsv;
    int chr_tsv;
    int bp_tsv;
    int a1_tsv;
    int a2_tsv;
    int beta_tsv;
    int or_tsv;
    int p_tsv;
    void *str2id;
    marker_t *markers;
    int n_markers;
    int m_markers;
} summary_t;

// check also on https://github.com/bulik/ldsc/blob/master/munge_sumstats.py
static const char *snp_hdr_str[] = {"snp", "snpid", "SNPID", "MarkerName", "Marker",
                                    "SNP", "rsID",  "Name",  "variant_id", "ID"};
static const char *chr_hdr_str[] = {"CHR", "Chromosome", "Chrom", "Chr", "chr_name", "chromosome", "#CHROM"};
static const char *bp_hdr_str[] = {"bp", "POS", "Position", "Pos", "BP", "chr_position", "base_pair_location"};
static const char *a1_hdr_str[] = {"effect_allele", "a1", "A1", "Effect_allele", "Allele1", "allele1"};
static const char *a2_hdr_str[] = {"other_allele",    "a2", "A2", "Non_Effect_allele", "Allele2", "allele2",
                                   "reference_allele"};
static const char *beta_hdr_str[] = {"effect", "BETA", "Beta", "Effect", "beta", "effect_weight", "A1Effect"};
static const char *p_hdr_str[] = {"pvalue", "pval", "Pvalue", "P.value", "P", "PValue", "p_value"};
static const char *or_hdr_str[] = {"OR"};

static summary_t *summary_init(const char *fn, bcf_hdr_t *hdr, int snp_id_mode) {
    summary_t *summary = (summary_t *)calloc(1, sizeof(summary_t));
    summary->str2id = khash_str2int_init();

    htsFile *fp = hts_open(fn, "r");
    if (fp == NULL) error("Could not open %s: %s\n", fn, strerror(errno));

    kstring_t str = {0, 0, NULL};
    if (hts_getline(fp, KS_SEP_LINE, &str) <= 0) error("Error reading from file: %s\n", fn);
    while (str.s[0] == '#' && strncmp(str.s, "#CHROM", 6) != 0) hts_getline(fp, KS_SEP_LINE, &str);
    // some formats are tab-delimited and some formats (e.g. PLINK and SBayesR) are not
    // here we make a determination based on the first header row
    char delimiter = strchr(str.s, '\t') ? '\t' : '\0';
    marker_t *marker = (marker_t *)calloc(1, sizeof(marker_t));
    tsv_t *tsv = tsv_init_delimiter(str.s, delimiter);

    for (int i = 0; i < sizeof(snp_hdr_str) / sizeof(char *); i++)
        if (tsv_register(tsv, snp_hdr_str[i], tsv_setter_id, NULL) == 0) summary->snp_tsv = 1;
    if (snp_id_mode && !summary->snp_tsv)
        error("Column for marker name is not provided in file %s but required with option --snp\n", fn);

    for (int i = 0; i < sizeof(chr_hdr_str) / sizeof(char *); i++)
        if (tsv_register(tsv, chr_hdr_str[i], tsv_setter_chrom_flexible, (void *)hdr) == 0) summary->chr_tsv = 1;

    for (int i = 0; i < sizeof(bp_hdr_str) / sizeof(char *); i++)
        if (tsv_register(tsv, bp_hdr_str[i], tsv_setter_pos, NULL) == 0) summary->bp_tsv = 1;

    if (snp_id_mode || !summary->chr_tsv || !summary->bp_tsv) summary->use_snp = 1;
    if (!summary->snp_tsv && summary->use_snp)
        error("Columns for chromosome and position required if column for marker name is not provided in file %s\n",
              fn);

    for (int i = 0; i < sizeof(a1_hdr_str) / sizeof(char *); i++)
        if (tsv_register(tsv, a1_hdr_str[i], tsv_read_string, (void *)&marker->a1) == 0) summary->a1_tsv = 1;

    if (!summary->a1_tsv) error("Column for effect allele required but missing from file %s\n", fn);

    for (int i = 0; i < sizeof(a2_hdr_str) / sizeof(char *); i++)
        if (tsv_register(tsv, a2_hdr_str[i], tsv_read_string, (void *)&marker->a2) == 0) summary->a2_tsv = 1;

    for (int i = 0; i < sizeof(beta_hdr_str) / sizeof(char *); i++)
        if (tsv_register(tsv, beta_hdr_str[i], tsv_read_float, (void *)&marker->beta) == 0) summary->beta_tsv = 1;

    for (int i = 0; i < sizeof(or_hdr_str) / sizeof(char *); i++)
        if (tsv_register(tsv, or_hdr_str[i], tsv_read_float, (void *)&marker->or) == 0) summary->or_tsv = 1;

    if (!summary->beta_tsv && !summary->or_tsv)
        error("Column for either effect weight or odds ratio required but missing from file %s\n", fn);

    for (int i = 0; i < sizeof(p_hdr_str) / sizeof(char *); i++)
        if (tsv_register(tsv, p_hdr_str[i], tsv_read_float, (void *)&marker->p) == 0) summary->p_tsv = 1;

    int force_warned = 0;
    bcf1_t *rec = bcf_init();
    char chr_bp_str[2 * (sizeof(int32_t) + sizeof(hts_pos_t)) + 1];
    while (hts_getline(fp, KS_SEP_LINE, &str) > 0) {
        if (str.s[0] == '#') continue; // skip comments
        rec->rid = -1;
        rec->pos = -1;
        hts_expand(marker_t, summary->n_markers + 1, summary->m_markers, summary->markers);
        if (!tsv_parse_delimiter(tsv, rec, str.s, delimiter)) {
            if (rec->rid < 0 && summary->chr_tsv && !force_warned) {
                fprintf(stderr,
                        "Warning: could not recognize chromosome in line:\n%s\n(This warning is printed only once.)\n",
                        str.s);
                force_warned = 1;
                continue;
            }
            marker->rid = rec->rid;
            marker->pos = rec->pos;
            if (!summary->beta_tsv) marker->beta = logf(marker->or);
            char *key;
            if (summary->use_snp) {
                key = strdup(rec->d.id);
            } else {
                sprintf(chr_bp_str, "%08" PRIx32 "%016" PRIx64, marker->rid, marker->pos);
                key = strdup(chr_bp_str);
            }
            int size = khash_str2int_size(summary->str2id);
            if (khash_str2int_inc(summary->str2id, key) < size) {
                if (!force_warned) {
                    if (summary->use_snp)
                        fprintf(stderr,
                                "Warning: could not include marker name %s as present multiple times\n(This warning is "
                                "printed only once.)\n",
                                rec->d.id);
                    else
                        fprintf(stderr,
                                "Warning: could not include chromosome position %s %" PRId64
                                " as present multiple times\n(This warning is printed only once.)\n",
                                bcf_hdr_id2name(hdr, rec->rid), rec->pos + 1);
                    force_warned = 1;
                }
                free(key);
                free(marker->a1);
                free(marker->a2);
            } else {
                memcpy((void *)&summary->markers[summary->n_markers], (const void *)marker, sizeof(marker_t));
                summary->n_markers++;
            }
        } else {
            error("Could not parse line: %s\n", str.s);
        }
    }
    bcf_destroy(rec);
    tsv_destroy(tsv);
    free(marker);
    free(str.s);
    hts_close(fp);

    return summary;
}

static void summary_destroy(summary_t *summary) {
    khash_str2int_destroy_free(summary->str2id);
    for (int i = 0; i < summary->n_markers; i++) {
        free(summary->markers[i].a1);
        free(summary->markers[i].a2);
    }
    free(summary->markers);
    free(summary);
}

/****************************************
 * PLUGIN                               *
 ****************************************/

const char *about(void) { return "Compute polygenic scores for an input cohort.\n"; }

static const char *usage_text(void) {
    return "\n"
           "About: Compute polygenic scores for an input cohort. (version " SCORE_VERSION
           " https://github.com/freeseek/mocha)\n"
           "\n"
           "Usage: bcftools +score [options] <in.vcf.gz> [<score1.gz> <score2.gz> ...]\n"
           "Plugin options:\n"
           "       --use <tag>               FORMAT tag to use to compute allele dosages: GP, AP, HDS, DS, GT, AS\n"
           "       --summaries <dir|file>    summary statistics files from directory or list from file\n"
           "       --snp-id                  use SNP ID to match variants\n"
           "   -v, --vcf                     summary statistics in VCF Summary Statistics format\n"
           "       --q-score-thr LIST        comma separated list of p-value thresholds\n"
           "       --counts                  include SNP counts in the output table\n"
           "   -o, --output <file.tsv>       write output to a file [standard output]\n"
           "       --sample-header           header for sample ID column [SAMPLE]\n"
           "   -e, --exclude <expr>          exclude sites for which the expression is true\n"
           "   -f, --apply-filters <list>    require at least one of the listed FILTER strings (e.g. \"PASS,.\")\n"
           "   -i, --include <expr>          select sites for which the expression is true\n"
           "   -r, --regions <region>        restrict to comma-separated list of regions\n"
           "   -R, --regions-file <file>     restrict to regions listed in a file\n"
           "   -t, --targets <region>        similar to -r but streams rather than index-jumps\n"
           "   -T, --targets-file <file>     similar to -R but streams rather than index-jumps\n"
           "   -s, --samples [^]<list>       comma separated list of samples to include (or exclude with \"^\" "
           "prefix)\n"
           "   -S, --samples-file [^]<file>  file of samples to include (or exclude with \"^\" prefix)\n"
           "       --force-samples           only warn about unknown subset samples\n"
           "\n"
           "Examples:\n"
           "   bcftools +score --use GT -o scores.tsv --q-score-thr 1e-8,1e-7,1e-6,1e-5,1e-4,0.001,0.05 input.bcf "
           "PGS000001.txt.gz\n"
           "   bcftools +score --use DS -o scores.tsv -i 'INFO>0.8 && AF>0.01 && AF<0.99' input.bcf PGS000001.txt.gz "
           "PGS000002.txt.gz\n"
           "\n";
}

static double *parse_list(const char *str, int *n) {
    char *endptr;
    char **s = hts_readlist(str, 0, n);
    double *v = (double *)malloc(*n * sizeof(double));
    for (int i = 0; i < *n; i++) {
        v[i] = strtof(s[i], &endptr);
        if (*endptr) error("Could not parse element: %s\n", s[i]);
        free(s[i]);
    }
    free(s);
    return v;
}

int run(int argc, char **argv) {
    int vcf_mode = 0;
    int snp_id_mode = 0;
    int use_tag = 0;
    int display_cnts = 0;
    int filter_logic = 0;
    int regions_is_file = 0;
    int targets_is_file = 0;
    int sample_is_file = 0;
    int force_samples = 0;
    const char *q_score_thr_str = NULL;
    const char *pathname = NULL;
    const char *output_fname = "-";
    const char *sample_header = "SAMPLE";
    const char *filter_str = NULL;
    const char *regions_list = NULL;
    const char *targets_list = NULL;
    const char *sample_names = NULL;
    filter_t *filter = NULL;
    bcf_srs_t *sr = bcf_sr_init();

    static struct option loptions[] = {{"use", required_argument, NULL, 1},
                                       {"summaries", required_argument, NULL, 2},
                                       {"snp-id", no_argument, NULL, 3},
                                       {"vcf", no_argument, NULL, 'v'},
                                       {"q-score-thr", required_argument, NULL, 4},
                                       {"counts", no_argument, NULL, 5},
                                       {"output", required_argument, NULL, 'o'},
                                       {"sample-header", required_argument, NULL, 6},
                                       {"exclude", required_argument, NULL, 'e'},
                                       {"apply-filters", required_argument, NULL, 'f'},
                                       {"include", required_argument, NULL, 'i'},
                                       {"regions", required_argument, NULL, 'r'},
                                       {"regions-file", required_argument, NULL, 'R'},
                                       {"targets", required_argument, NULL, 't'},
                                       {"targets-file", required_argument, NULL, 'T'},
                                       {"samples", required_argument, NULL, 's'},
                                       {"samples-file", required_argument, NULL, 'S'},
                                       {"force-samples", no_argument, NULL, 7},
                                       {NULL, 0, NULL, 0}};
    int c;
    while ((c = getopt_long(argc, argv, "h?vw:p:o:e:f:i:r:R:t:T:s:S:", loptions, NULL)) >= 0) {
        switch (c) {
        case 1:
            if (!strcasecmp(optarg, "GT"))
                use_tag = SCORE_GT;
            else if (!strcasecmp(optarg, "DS"))
                use_tag = SCORE_DS;
            else if (!strcasecmp(optarg, "HDS"))
                use_tag = SCORE_HDS;
            else if (!strcasecmp(optarg, "AP"))
                use_tag = SCORE_AP;
            else if (!strcasecmp(optarg, "GP"))
                use_tag = SCORE_GP;
            else if (!strcasecmp(optarg, "AS"))
                use_tag = SCORE_AS;
            else
                error("The argument not recognised, expected --use GT, DS, HDS, AP, GP, or AS: %s\n", optarg);
            break;
        case 2:
            pathname = optarg;
            break;
        case 3:
            snp_id_mode = 1;
            break;
        case 'v':
            vcf_mode = 1;
            break;
        case 4:
            q_score_thr_str = optarg;
            break;
        case 5:
            display_cnts = 1;
            break;
        case 'o':
            output_fname = optarg;
            break;
        case 6:
            sample_header = optarg;
            break;
        case 'e':
            filter_str = optarg;
            filter_logic |= FLT_EXCLUDE;
            break;
        case 'f':
            sr->apply_filters = optarg;
            break;
            break;
        case 'i':
            filter_str = optarg;
            filter_logic |= FLT_INCLUDE;
            break;
        case 'r':
            regions_list = optarg;
            break;
        case 'R':
            regions_list = optarg;
            regions_is_file = 1;
            break;
        case 't':
            targets_list = optarg;
            break;
        case 'T':
            targets_list = optarg;
            targets_is_file = 1;
            break;
        case 's':
            sample_names = optarg;
            break;
        case 'S':
            sample_names = optarg;
            sample_is_file = 1;
            break;
        case 7:
            force_samples = 1;
            break;
        case 'h':
        case '?':
        default:
            error("%s", usage_text());
            break;
        }
    }

    if ((pathname && optind + 1 != argc) || (!pathname && optind + 2 > argc)) error("%s", usage_text());

    sr->require_index = 1;
    if (filter_logic == (FLT_EXCLUDE | FLT_INCLUDE)) error("Only one of --include or --exclude can be given.\n");
    if (regions_list) {
        if (bcf_sr_set_regions(sr, regions_list, regions_is_file) < 0)
            error("Failed to read the regions: %s\n", regions_list);
    }
    if (targets_list) {
        if (bcf_sr_set_targets(sr, targets_list, targets_is_file, 0) < 0)
            error("Failed to read the targets: %s\n", targets_list);
        sr->collapse |= COLLAPSE_BOTH;
    }

    if (!bcf_sr_add_reader(sr, argv[optind]))
        error("Error opening %s: %s\n", argv[optind], bcf_sr_strerror(sr->errnum));

    bcf_hdr_t *hdr = bcf_sr_get_header(sr, 0);
    if (filter_str) filter = filter_init(hdr, filter_str);
    int n_q_score_thr = 1;
    double *q_score_thr = q_score_thr_str ? parse_list(q_score_thr_str, &n_q_score_thr) : NULL;
    if (q_score_thr)
        for (int i = 0; i < n_q_score_thr; i++) q_score_thr[i] = -log10(q_score_thr[i]);

    // subset VCF file
    if (sample_names) {
        int ret = bcf_hdr_set_samples(hdr, sample_names, sample_is_file);
        if (ret < 0)
            error("Error parsing the sample list\n");
        else if (ret > 0) {
            if (force_samples)
                fprintf(stderr, "Warn: sample #%d not found in the header... skipping\n", ret);
            else
                error(
                    "Error: sample #%d not found in the header. Use \"--force-samples\" to "
                    "ignore this error\n",
                    ret);
        }
        if (bcf_hdr_nsamples(hdr) == 0) error("Error: subsetting has removed all samples\n");
    }
    int n_smpls = bcf_hdr_nsamples(hdr);

    int n_files, n_prs, *prs2vcf = NULL, m_prs2vcf = 0, *prs2idx = NULL, m_prs2idx = 0;
    char **filenames = NULL;
    if (pathname) {
        filenames = get_file_list(pathname, &n_files);
    } else {
        n_files = argc - optind - 1;
        filenames = argv + optind + 1;
    }
    summary_t **summaries = NULL;
    char **prs_names = NULL;
    int m_prs_names = 0;
    if (vcf_mode) {
        n_prs = 0;
        for (int i = 0; i < n_files; i++) {
            if (!bcf_sr_add_reader(sr, filenames[i]))
                error("Error opening %s: %s\n", filenames[i], bcf_sr_strerror(sr->errnum));
            hdr = bcf_sr_get_header(sr, i + 1);
            hts_expand(int, n_prs + bcf_hdr_nsamples(hdr), m_prs2vcf, prs2vcf);
            hts_expand(int, n_prs + bcf_hdr_nsamples(hdr), m_prs2idx, prs2idx);
            hts_expand(char *, n_prs + bcf_hdr_nsamples(hdr), m_prs_names, prs_names);
            for (int j = 0; j < bcf_hdr_nsamples(hdr); j++) {
                prs2vcf[n_prs + j] = i + 1;
                prs2idx[n_prs + j] = j;
                prs_names[n_prs + j] = strdup(hdr->samples[j]);
            }
            n_prs += bcf_hdr_nsamples(hdr);
            int es_id = bcf_hdr_id2int(hdr, BCF_DT_ID, "ES");
            if (!bcf_hdr_idinfo_exists(hdr, BCF_HL_FMT, es_id))
                error("VCF summary statistics file %s does not include the ES FORMAT field\n", filenames[i]);
            if (q_score_thr) {
                int lp_id = bcf_hdr_id2int(hdr, BCF_DT_ID, "LP");
                if (!bcf_hdr_idinfo_exists(hdr, BCF_HL_FMT, lp_id))
                    error("VCF summary statistics file %s does not include the LP FORMAT field\n", filenames[i]);
            }
        }
    } else {
        n_prs = n_files;
        summaries = (summary_t **)malloc((n_prs) * sizeof(summary_t *));
        prs_names = (char **)malloc((n_prs) * sizeof(char *));
        hdr = bcf_sr_get_header(sr, 0);
        for (int i = 0; i < n_prs; i++) {
            summaries[i] = summary_init(filenames[i], hdr, snp_id_mode);
            fprintf(stderr, "Read %d markers from file %s and matching by %s\n", summaries[i]->n_markers, filenames[i],
                    summaries[i]->use_snp ? "marker name" : "chromosome position");
            char *ptr, *ext_str[] = {"gz", "txt", "tsv", "vcf", "bcf"};
            int j = 0;
            while (j < sizeof(ext_str) / sizeof(char *) && (ptr = strrchr(filenames[i], '.')))
                for (j = 0; j < sizeof(ext_str) / sizeof(char *); j++)
                    if (strcmp(ptr + 1, ext_str[j]) == 0) {
                        *ptr = '\0';
                        break;
                    }
            prs_names[i] = strdup(strrchr(filenames[i], '/') ? strrchr(filenames[i], '/') + 1 : filenames[i]);
        }
    }

    hdr = bcf_sr_get_header(sr, 0);
    int gt_id = bcf_hdr_id2int(hdr, BCF_DT_ID, "GT");
    int ds_id = bcf_hdr_id2int(hdr, BCF_DT_ID, "DS");
    int hds_id = bcf_hdr_id2int(hdr, BCF_DT_ID, "HDS");
    int ap1_id = bcf_hdr_id2int(hdr, BCF_DT_ID, "AP1");
    int ap2_id = bcf_hdr_id2int(hdr, BCF_DT_ID, "AP2");
    int gp_id = bcf_hdr_id2int(hdr, BCF_DT_ID, "GP");
    int as_id = bcf_hdr_id2int(hdr, BCF_DT_ID, "AS");

    if (!use_tag) {
        if (bcf_hdr_idinfo_exists(hdr, BCF_HL_FMT, gt_id)) use_tag = SCORE_GT;
        if (bcf_hdr_idinfo_exists(hdr, BCF_HL_FMT, ds_id)) use_tag = SCORE_DS;
        if (bcf_hdr_idinfo_exists(hdr, BCF_HL_FMT, hds_id)) use_tag = SCORE_HDS;
        if (bcf_hdr_idinfo_exists(hdr, BCF_HL_FMT, ap1_id) && bcf_hdr_idinfo_exists(hdr, BCF_HL_FMT, ap2_id))
            use_tag = SCORE_AP;
        if (bcf_hdr_idinfo_exists(hdr, BCF_HL_FMT, gp_id)) use_tag = SCORE_GP;
        if (!use_tag) error("VCF file %s does not include any of the GT, GP, or DS FORMAT fields\n", argv[optind]);
    } else {
        switch (use_tag) {
        case SCORE_GT:
            if (!bcf_hdr_idinfo_exists(hdr, BCF_HL_FMT, gt_id))
                error("VCF file %s does not include the GT FORMAT field\n", argv[optind]);
            break;
        case SCORE_DS:
            if (!bcf_hdr_idinfo_exists(hdr, BCF_HL_FMT, ds_id))
                error("VCF file %s does not include the DS FORMAT field\n", argv[optind]);
            break;
        case SCORE_HDS: // only for Minimac4 VCFs
            if (!bcf_hdr_idinfo_exists(hdr, BCF_HL_FMT, hds_id))
                error("VCF file %s does not include the HDS FORMAT field\n", argv[optind]);
            break;
        case SCORE_AP:
            if (!bcf_hdr_idinfo_exists(hdr, BCF_HL_FMT, ap1_id) || !bcf_hdr_idinfo_exists(hdr, BCF_HL_FMT, ap2_id))
                error("VCF file %s does not include either the AP1 or the AP2 FORMAT fields\n", argv[optind]);
            break;
        case SCORE_GP:
            if (!bcf_hdr_idinfo_exists(hdr, BCF_HL_FMT, gp_id))
                error("VCF file %s does not include the GP FORMAT field\n", argv[optind]);
            break;
        case SCORE_AS:
            if (!bcf_hdr_idinfo_exists(hdr, BCF_HL_FMT, as_id))
                error("VCF file %s does not include the AS FORMAT field\n", argv[optind]);
            break;
        }
    }

    fprintf(stderr, "Using %s to compute polygenic scores\n",
            use_tag == SCORE_GT    ? "genotypes (GT)"
            : use_tag == SCORE_GP  ? "genotype probabilities (GP)"
            : use_tag == SCORE_AP  ? "ALT haplotype probabilities (AP)"
            : use_tag == SCORE_HDS ? "haploid alternate allele dosage (HDS)"
            : use_tag == SCORE_DS  ? "genotype dosages (DS)"
                                   : "allelic shifts (AS)");

    FILE *out_fh = strcmp("-", output_fname) ? fopen(output_fname, "w") : stdout;
    if (!out_fh) error("Error: cannot write to %s\n", output_fname);

    int m_int32 = 0, m_float = 0, n_float, m_alleles = 2 * n_smpls;
    int32_t *int32_arr = NULL;
    float *float_arr = NULL;
    char *str = NULL;
    float *alleles = (float *)malloc(m_alleles * sizeof(float));
    float *missing = (float *)malloc(n_smpls * sizeof(float));
    int *idxs = vcf_mode ? NULL : (int *)malloc(n_prs * sizeof(int));
    float *scores = (float *)calloc(n_prs * n_q_score_thr * n_smpls, sizeof(float));
    int *cnts = (int *)calloc(n_prs * n_q_score_thr * n_smpls, sizeof(int));

    while (bcf_sr_next_line(sr)) {
        if (!bcf_sr_has_line(sr, 0)) continue;
        bcf1_t *line = bcf_sr_get_line(sr, 0);
        if (filter) {
            int ret = filter_test(filter, line, NULL);
            if ((filter_logic == FLT_INCLUDE && !ret) || ret) continue;
        }

        int skip_line = 1;
        if (vcf_mode) {
            for (int i = 0; i < n_files; i++) {
                if (!bcf_sr_has_line(sr, i + 1)) continue;
                skip_line = 0;
                break;
            }
        } else {
            char chr_bp_str[2 * (sizeof(int32_t) + sizeof(hts_pos_t)) + 1];
            sprintf(chr_bp_str, "%08" PRIx32 "%016" PRIx64, line->rid, line->pos);
            memset((void *)idxs, -1, n_prs * sizeof(int));
            for (int i = 0; i < n_prs; i++) {
                if (khash_str2int_get(summaries[i]->str2id, summaries[i]->use_snp ? line->d.id : chr_bp_str, &idxs[i])
                    < 0)
                    continue;
                skip_line = 0;
            }
        }
        if (skip_line) continue;

        hdr = bcf_sr_get_header(sr, 0);
        hts_expand(float, line->n_allele *n_smpls, m_alleles, alleles);
        memset((void *)alleles, 0, line->n_allele * n_smpls * sizeof(float));
        memset((void *)missing, 0, n_smpls * sizeof(int));
        int number;
        char *ap_str[] = {"AP1", "AP2"};
        switch (use_tag) {
        case SCORE_GT:
            number = bcf_get_genotypes(hdr, line, &int32_arr, &m_int32);
            number /= bcf_hdr_nsamples(hdr);
            assert(number == 2);
            for (int k = 0; k < n_smpls; k++) {
                int32_t *ptr = int32_arr + (number * k);
                if (bcf_gt_is_missing(ptr[0]) || bcf_gt_is_missing(ptr[1])) {
                    missing[k] = 1;
                } else {
                    alleles[bcf_gt_allele(ptr[0]) * n_smpls + k]++;
                    alleles[bcf_gt_allele(ptr[1]) * n_smpls + k]++;
                }
            }
            break;
        case SCORE_DS:
            number = bcf_get_format_float(hdr, line, "DS", &float_arr, &m_float);
            number /= bcf_hdr_nsamples(hdr);
            assert(number == line->n_allele - 1);
            if (number == 1) { // line->n_allele == 2
                for (int k = 0; k < n_smpls; k++) {
                    if (bcf_float_is_missing(float_arr[k])) {
                        missing[k] = 1;
                    } else {
                        alleles[k] += 2.0f - float_arr[k]; // check whether this should be replaced with 2.0f
                        alleles[n_smpls + k] += float_arr[k];
                    }
                }
            } else {
                for (int k = 0; k < n_smpls; k++) {
                    float *ptr = float_arr + (number * k);
                    alleles[k] += 2.0f; // check whether this should be replaced with 2.0f
                    for (int idx = 0; idx < number; idx++) {
                        if (bcf_float_is_missing(ptr[idx])) {
                            missing[k] = 1;
                        } else {
                            alleles[k] -= ptr[idx];
                            alleles[(idx + 1) * n_smpls + k] += ptr[idx];
                        }
                    }
                }
            }
            break;
        case SCORE_HDS: // only for Minimac4 VCFs
            number = bcf_get_format_float(hdr, line, "HDS", &float_arr, &m_float);
            number /= bcf_hdr_nsamples(hdr);
            assert(number == 2 && line->n_allele == 2);
            for (int k = 0; k < n_smpls; k++) {
                if (bcf_float_is_missing(float_arr[2 * k]) || bcf_float_is_missing(float_arr[2 * k + 1])) {
                    missing[k] = 1;
                } else {
                    alleles[k] += 1.0f - float_arr[2 * k] - float_arr[2 * k + 1];
                    alleles[n_smpls + k] += float_arr[2 * k] + float_arr[2 * k + 1];
                }
            }
            break;
        case SCORE_AP:
            for (int ap = 0; ap < sizeof(ap_str) / sizeof(char *); ap++) {
                number = bcf_get_format_float(hdr, line, ap_str[ap], &float_arr, &m_float);
                number /= bcf_hdr_nsamples(hdr);
                assert(number == line->n_allele - 1);
                if (number == 1) { // line->n_allele == 2
                    for (int k = 0; k < n_smpls; k++) {
                        if (bcf_float_is_missing(float_arr[k])) {
                            missing[k] = 1;
                        } else {
                            alleles[k] += 1.0f - float_arr[k];
                            alleles[n_smpls + k] += float_arr[k];
                        }
                    }
                } else {
                    for (int k = 0; k < n_smpls; k++) {
                        float *ptr = float_arr + (number * k);
                        alleles[k] += 1.0f;
                        for (int idx = 0; idx < number; idx++) {
                            if (bcf_float_is_missing(ptr[idx])) {
                                missing[k] = 1;
                            } else {
                                alleles[k] -= ptr[idx];
                                alleles[(idx + 1) * n_smpls + k] += ptr[idx];
                            }
                        }
                    }
                }
            }
            break;
        case SCORE_GP:
            number = bcf_get_format_float(hdr, line, "GP", &float_arr, &m_float);
            number /= bcf_hdr_nsamples(hdr);
            assert(number == (line->n_allele) * (line->n_allele + 1) / 2);
            if (number == 3) { // line->n_allele == 2
                for (int k = 0; k < n_smpls; k++) {
                    float *ptr = float_arr + (number * k);
                    if (bcf_float_is_missing(ptr[0]) || bcf_float_is_missing(ptr[1]) || bcf_float_is_missing(ptr[2])) {
                        missing[k] = 1;
                    } else {
                        alleles[k] += 2.0f * ptr[0] + ptr[1];
                        alleles[n_smpls + k] += ptr[1] + 2.0f * ptr[2];
                    }
                }
            } else {
                for (int k = 0; k < n_smpls; k++) {
                    float *ptr = float_arr + (number * k);
                    // The Variant Call Format Specification
                    // for P=2 and N=2, the ordering is 00,01,11,02,12,22
                    // for P=2, the index of the genotype “a/b”, where a≤b, is b(b+ 1)/2 +a
                    for (int b = 0; b < line->n_allele; b++) {
                        for (int a = 0; a <= b; a++) {
                            int idx = b * (b + 1) / 2 + a;
                            if (bcf_float_is_missing(ptr[idx])) {
                                missing[k] = 1;
                            } else {
                                alleles[a * n_smpls + k] += ptr[idx];
                                alleles[b * n_smpls + k] += ptr[idx];
                            }
                        }
                    }
                }
            }
            break;
        case SCORE_AS:
            number = bcf_get_format_int32(hdr, line, "AS", &int32_arr, &m_int32);
            number /= bcf_hdr_nsamples(hdr);
            assert(number == 1 && line->n_allele == 2);
            for (int k = 0; k < n_smpls; k++) {
                if (int32_arr[k] == 0 || int32_arr[k] == bcf_int32_missing) {
                    missing[k] = 1;
                } else {
                    alleles[k] -= (float)int32_arr[k];
                    alleles[n_smpls + k] += (float)int32_arr[k];
                }
            }
            break;
        }

        float beta, minus_log10p = 0.0f;
        for (int i = 0; i < n_prs; i++) {
            char *a1;
            if (vcf_mode) {
                if (!bcf_sr_has_line(sr, prs2vcf[i])) continue;
                hdr = bcf_sr_get_header(sr, prs2vcf[i]);
                line = bcf_sr_get_line(sr, prs2vcf[i]);
                a1 = line->d.allele[1];
                n_float = bcf_get_format_float(hdr, line, "ES", &float_arr, &m_float);
                if (n_float < 0) continue;
                if (n_float != bcf_hdr_nsamples(hdr))
                    error("VCF file %s has incorrect number of ES fields at position %" PRId64 "\n",
                          filenames[prs2vcf[i] - 1], line->pos + 1);
                beta = float_arr[prs2idx[i]];
                if (q_score_thr) {
                    n_float = bcf_get_format_float(hdr, line, "LP", &float_arr, &m_float);
                    if (n_float < 0) continue;
                    if (n_float != bcf_hdr_nsamples(hdr))
                        error("VCF file %s has incorrect number of LP fields at position %" PRId64 "\n",
                              filenames[prs2vcf[i] - 1], line->pos + 1);
                    minus_log10p = float_arr[prs2idx[i]];
                }
            } else {
                if (idxs[i] < 0) continue;
                marker_t *marker = &summaries[i]->markers[idxs[i]];
                a1 = marker->a1;
                beta = marker->beta;
                minus_log10p = -log10(marker->p);
            }
            // find effect allele
            int idx_allele;
            for (idx_allele = 0; idx_allele < line->n_allele; idx_allele++)
                if (strcmp(a1, line->d.allele[idx_allele]) == 0) break;
            if (idx_allele == line->n_allele) continue;
            for (int j = 0; j < n_q_score_thr; j++) {
                if (q_score_thr && minus_log10p < q_score_thr[j]) continue;
                float *ptr = scores + (i * n_q_score_thr + j) * n_smpls;
                float *ptr2 = alleles + idx_allele * n_smpls;
                int *ptr3 = cnts + (i * n_q_score_thr + j) * n_smpls;
                for (int k = 0; k < n_smpls; k++)
                    if (missing[k] == 0) {
                        ptr[k] += beta * ptr2[k];
                        ptr3[k]++;
                    }
            }
        }
    }

    hdr = bcf_sr_get_header(sr, 0);
    fprintf(out_fh, "%s", sample_header);
    for (int i = 0; i < n_prs; i++) {
        for (int j = 0; j < n_q_score_thr; j++) {
            fprintf(out_fh, "\t%s", prs_names[i]);
            if (q_score_thr) fprintf(out_fh, "_p%.4g", exp(-M_LN10 * q_score_thr[j]));
            if (display_cnts) {
                fprintf(out_fh, "\t%s_CNT", prs_names[i]);
                if (q_score_thr) fprintf(out_fh, "_p%.4g", exp(-M_LN10 * q_score_thr[j]));
            }
        }
    }
    fprintf(out_fh, "\n");
    for (int k = 0; k < n_smpls; k++) {
        fprintf(out_fh, "%s", hdr->samples[k]);
        for (int i = 0; i < n_prs; i++)
            for (int j = 0; j < n_q_score_thr; j++) {
                fprintf(out_fh, "\t%f", scores[(i * n_q_score_thr + j) * n_smpls + k]);
                if (display_cnts) fprintf(out_fh, "\t%d", cnts[(i * n_q_score_thr + j) * n_smpls + k]);
            }
        fprintf(out_fh, "\n");
    }

    if (filter) filter_destroy(filter);
    if (pathname) {
        for (int i = 0; i < n_files; i++) free(filenames[i]);
        free(filenames);
    }
    for (int i = 0; i < n_prs; i++) free(prs_names[i]);
    free(prs_names);
    if (out_fh != stdout) fclose(out_fh);
    free(scores);
    free(cnts);
    free(idxs);
    free(alleles);
    free(missing);
    free(str);
    free(int32_arr);
    free(float_arr);
    free(q_score_thr);
    free(prs2vcf);
    free(prs2idx);
    if (summaries) {
        for (int i = 0; i < n_prs; i++) summary_destroy(summaries[i]);
        free(summaries);
    }
    bcf_sr_destroy(sr);

    return 0;
}
