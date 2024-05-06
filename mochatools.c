/* The MIT License

   Copyright (C) 2017-2024 Giulio Genovese

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
#include <htslib/vcf.h>
#include <htslib/faidx.h>
#include <htslib/kfunc.h>
#include "mocha.h"
#include "bcftools.h"

#define MOCHATOOLS_VERSION "2024-05-05"

#define TAG_LIST_DFLT "none"
#define GC_WIN_DFLT "200"

#define INFO_GC (1 << 0)
#define INFO_CPG (1 << 1)
#define INFO_ALLELE_A (1 << 2)
#define INFO_ALLELE_B (1 << 3)
#define INFO_AC_SEX (1 << 4)
#define INFO_MISSING_SEX (1 << 5)
#define INFO_pAC_HET (1 << 6)
#define INFO_AD_HET (1 << 7)
#define INFO_pBAF_STATS (1 << 8)
#define INFO_COR_BAF_LRR (1 << 9)
#define INFO_MACH (1 << 10)

// see Marchini, J., Howie, B. Genotype imputation for genome-wide association studies. Nat Rev Genet 11, 499–511
// (2010). https://doi.org/10.1038/nrg2796
// ##FORMAT=<ID=DS,Number=A,Type=Float,Description="Genotype dosage">
#define SCORE_DS 1 // DS = AP1 + AP2
// ##FORMAT=<ID=HDS,Number=2,Type=Float,Description="Estimated Haploid Alternate Allele Dosage ">
#define SCORE_HDS 2
// ##FORMAT=<ID=AP1,Number=A,Type=Float,Description="ALT allele probability of first haplotype">
// ##FORMAT=<ID=AP2,Number=A,Type=Float,Description="ALT allele probability of second haplotype">
#define SCORE_AP 3
// ##FORMAT=<ID=GP,Number=G,Type=Float,Description="Estimated Genotype Probability">
#define SCORE_GP 4

static inline double sq(double x) { return x * x; }

typedef struct {
    int flags;
    int adjust; // whether to adjust BAF and LRR
    int gc_win; // length of window to compute GC and CpG content
    int *gender;
    char *as_str;   // format field ID to summarize for allelic shift count summary
    int phase;      // whether phase information should be included in the allelic shift
    char *info_str; // info field for binomial or Fisher's exact test
    int phred;      // whether the test result should be phred scaled
    int odds_ratio; // whether the odds ratio should be computed
    int log_odds;   // whether the log odds ratio should be computed
    int dosage_tag; // which dosage tag should be used
    int nsmpl, allele_a_id, allele_b_id, adjust_id, gt_id, ad_id, baf_id, lrr_id, info_id, as_id, ds_id, hds_id, ap1_id,
        ap2_id, gp_id;
    int8_t *gt_phase_arr, *as_arr;
    int16_t *gt0_arr, *gt1_arr, *ad0_arr, *ad1_arr;
    int32_t *info_arr;
    int m_info;
    float *adjust_arr;
    int m_adjust;
    float *baf_arr[2];
    int *imap_arr;
    faidx_t *fai;
    bcf_hdr_t *in_hdr, *out_hdr;
    kstring_t summary_str, test_str, odds_ratio_str, log_odds_str;
} args_t;

args_t *args;

static int *read_computed_gender(bcf_hdr_t *hdr, char *fname) {
    htsFile *fp = hts_open(fname, "r");
    if (fp == NULL) error("Could not open %s: %s\n", fname, strerror(errno));

    kstring_t str = {0, 0, NULL};
    if (hts_getline(fp, KS_SEP_LINE, &str) <= 0) error("Empty file: %s\n", fname);
    tsv_t *tsv = tsv_init_delimiter(str.s, '\t');
    int idx, computed_gender;
    if (tsv_register(tsv, "sample_id", tsv_read_sample_id, (void *)&idx) < 0)
        error("File %s is missing the sample_id column\n", fname);
    if (tsv_register(tsv, "computed_gender", tsv_read_computed_gender, (void *)&computed_gender) < 0)
        error("File %s is missing the computed_gender column\n", fname);

    int i = 0;
    int *gender = (int *)calloc((size_t)bcf_hdr_nsamples(hdr), sizeof(int));
    while (hts_getline(fp, KS_SEP_LINE, &str) > 0) {
        if (!tsv_parse_delimiter(tsv, (bcf1_t *)hdr, str.s, '\t')) {
            if (idx < 0) continue;
            gender[idx] = computed_gender;
            i++;
        } else {
            error("Could not parse line: %s\n", str.s);
        }
    }

    tsv_destroy(tsv);
    free(str.s);
    hts_close(fp);
    return gender;
}

const char *about(void) { return "MOsaic CHromosomal Alterations tools.\n"; }

const char *usage(void) {
    return "\n"
           "About: tools for the MOsaic CHromosomal Alterations pipeline. "
           "(version " MOCHATOOLS_VERSION
           " https://github.com/freeseek/mocha)\n"
           "Usage: bcftools +mochatools [General Options] -- [Plugin Options]\n"
           "Options:\n"
           "   run \"bcftools plugin\" for a list of common options\n"
           "\n"
           "Plugin options:\n"
           "   -l, --list-tags               list available INFO tags with description for VCF output\n"
           "   -t, --tags LIST               list of output INFO tags [" TAG_LIST_DFLT
           "]\n"
           "       --adjust                  adjust BAF and LRR using INFO/ADJ_COEFF\n"
           "   -f, --fasta-ref <file>        reference sequence to compute GC and CpG content\n"
           "       --gc-window-size <int>    window size in bp used to compute the GC and CpG content [" GC_WIN_DFLT
           "]\n"
           "   -x, --sex <file>              file including information about the gender of the samples\n"
           "       --dosage-tag              select dosage tag to use for dosage computations (DS, HDS, AP, GP)\n"
           "   -s, --samples [^]<list>       comma separated list of samples to include (or exclude with \"^\" "
           "prefix)\n"
           "   -S, --samples-file [^]<file>  file of samples to include (or exclude with \"^\" prefix)\n"
           "       --force-samples           only warn about unknown subset samples\n"
           "   -G, --drop-genotypes          drop individual genotype information (after running statistical tests)\n"
           "\n"
           "       --summary <tag>           allelic shift FORMAT tag to summarize\n"
           "       --phase                   whether phase information should be included in the allelic shift\n"
           "       --test <tag>              performs binomial or Fisher's exact test for INFO tag\n"
           "       --phred                   reports p-values as phred scaled\n"
           "       --odds-ratio              includes odds ratios for the INFO tag\n"
           "       --log-odds                includes log odds ratios for the INFO tag\n"
           "\n"
           "Examples:\n"
           "    bcftools +mochatools -- -l\n"
           "    bcftools +mochatools file.bcf -Ob -o output.bcf -- -t GC -f human_g1k_v37.fasta\n"
           "    bcftools +mochatools file.bcf -Ob -o output.bcf -- -t ALLELE_A,ALLELE_B\n"
           "    bcftools +mochatools file.bcf -- --summary AS --test AS --drop-genotypes\n"
           "    bcftools +mochatools file.bcf -- --summary AS --phase --test pAS --drop-genotypes\n"
           "\n";
}

static int parse_tags(const char *str) {
    int i, flags = 0, n_tags;
    if (!strcasecmp(str, "none")) return flags;
    char **tags = hts_readlist(str, 0, &n_tags);
    for (i = 0; i < n_tags; i++) {
        if (!strcasecmp(tags[i], "GC"))
            flags |= INFO_GC;
        else if (!strcasecmp(tags[i], "CpG"))
            flags |= INFO_CPG;
        else if (!strcasecmp(tags[i], "ALLELE_A"))
            flags |= INFO_ALLELE_A;
        else if (!strcasecmp(tags[i], "ALLELE_B"))
            flags |= INFO_ALLELE_B;
        else if (!strcasecmp(tags[i], "AC_Sex"))
            flags |= INFO_AC_SEX;
        else if (!strcasecmp(tags[i], "MISSING_Sex"))
            flags |= INFO_MISSING_SEX;
        else if (!strcasecmp(tags[i], "pAC_Het"))
            flags |= INFO_pAC_HET;
        else if (!strcasecmp(tags[i], "AD_Het"))
            flags |= INFO_AD_HET;
        else if (!strcasecmp(tags[i], "pBAF_Stats"))
            flags |= INFO_pBAF_STATS;
        else if (!strcasecmp(tags[i], "Cor_BAF_LRR"))
            flags |= INFO_COR_BAF_LRR;
        else if (!strcasecmp(tags[i], "MACH"))
            flags |= INFO_MACH;
        else
            error("Error parsing \"--tags %s\": the tag \"%s\" is not supported\n", str, tags[i]);
        free(tags[i]);
    }
    if (n_tags) free(tags);
    return flags;
}

static void list_tags(void) {
    error(
        "INFO/GC           Number:1  Type:Float    ..  GC ratio content around the marker\n"
        "INFO/CpG          Number:1  Type:Float    ..  CpG ratio content around the marker\n"
        "INFO/ALLELE_A     Number:1  Type:Integer  ..  A allele\n"
        "INFO/ALLELE_B     Number:1  Type:Integer  ..  B allele\n"
        "INFO/AC_Sex       Number:4  Type:Integer  ..  Number of homozygous and heterozygous genotypes by gender\n"
        "INFO/MISSING_Sex  Number:4  Type:Integer  ..  Number of non-missing and missing genotypes by gender\n"
        "INFO/pAC_Het      Number:2  Type:Integer  ..  Number of heterozygous genotypes by transmission type\n"
        "INFO/AD_Het       Number:2  Type:Integer  ..  Allelic depths for the reference and alternate alleles across "
        "heterozygous genotypes\n"
        "INFO/pBAF_Stats   Number:4  Type:Float    ..  Welch's t-test and Mann-Whitney U test for BAF by transmission "
        "type across heterozygous genotypes\n"
        "INFO/Cor_BAF_LRR  Number:3  Type=Float    ..  Pearson correlation for BAF and LRR at AA, AB, and BB "
        "genotypes\n"
        "INFO/MACH         Number:3  Type=Float    ..  Statistics to infer the MACH r2 imputation measure\n");
}

int init(int argc, char **argv, bcf_hdr_t *in, bcf_hdr_t *out) {
    args = (args_t *)calloc(1, sizeof(args_t));
    const char *tag_list = TAG_LIST_DFLT;
    char *tmp;
    args->gc_win = (int)strtol(GC_WIN_DFLT, NULL, 0);
    args->in_hdr = in;
    args->out_hdr = out;
    args->info_str = NULL;
    args->as_str = NULL;
    int i;
    int sample_is_file = 0;
    int force_samples = 0;
    int sites_only = 0;
    char *sample_names = NULL;
    char *gender_fname = NULL;
    char *ref_fname = NULL;
    kstring_t str = {0, 0, NULL};

    int c;
    static struct option loptions[] = {{"list-tags", no_argument, NULL, 'l'},
                                       {"tags", required_argument, NULL, 't'},
                                       {"adjust", no_argument, NULL, 1},
                                       {"fasta-ref", required_argument, NULL, 'f'},
                                       {"gc-window-size", required_argument, NULL, 2},
                                       {"sex", required_argument, NULL, 'x'},
                                       {"dosage-tag", required_argument, NULL, 3},
                                       {"samples", required_argument, NULL, 's'},
                                       {"samples-file", required_argument, NULL, 'S'},
                                       {"force-samples", no_argument, NULL, 4},
                                       {"drop-genotypes", no_argument, NULL, 'G'},
                                       {"summary", required_argument, NULL, 5},
                                       {"phase", no_argument, NULL, 6},
                                       {"test", required_argument, NULL, 7},
                                       {"phred", no_argument, NULL, 8},
                                       {"odds-ratio", no_argument, NULL, 9},
                                       {"log-odds", no_argument, NULL, 10},
                                       {NULL, 0, NULL, 0}};

    while ((c = getopt_long(argc, argv, "h?lt:f:x:s:S:G", loptions, NULL)) >= 0) {
        switch (c) {
        case 'l':
            list_tags();
            break;
        case 't':
            tag_list = optarg;
            break;
        case 1:
            args->adjust = 1;
            break;
        case 'f':
            ref_fname = optarg;
            break;
        case 2:
            args->gc_win = (int)strtol(optarg, &tmp, 0);
            if (*tmp) error("Could not parse: -w %s\n", optarg);
            if (args->gc_win <= 0) error("Window size is not positive: -w %s\n", optarg);
            break;
        case 'x':
            gender_fname = optarg;
            break;
        case 's':
            sample_names = optarg;
            break;
        case 'S':
            sample_names = optarg;
            sample_is_file = 1;
            break;
        case 3:
            if (!strcasecmp(optarg, "DS"))
                args->dosage_tag = SCORE_DS;
            else if (!strcasecmp(optarg, "HDS"))
                args->dosage_tag = SCORE_HDS;
            else if (!strcasecmp(optarg, "AP"))
                args->dosage_tag = SCORE_AP;
            else if (!strcasecmp(optarg, "GP"))
                args->dosage_tag = SCORE_GP;
            else
                error("The argument not recognised, expected --dosage-tag DS, HDS, AP, or GP: %s\n", optarg);
            break;
        case 4:
            force_samples = 1;
            break;
        case 'G':
            sites_only = 1;
            break;
        case 5:
            args->as_str = optarg;
            break;
        case 6:
            args->phase = 1;
            break;
        case 7:
            args->info_str = optarg;
            break;
        case 8:
            args->phred = 1;
            break;
        case 9:
            args->odds_ratio = 1;
            break;
        case 10:
            args->log_odds = 1;
            break;
        case 'h':
        case '?':
        default:
            error("%s", usage());
            break;
        }
    }

    fprintf(stderr,
            "MOCHATOOLS_VERSION " MOCHATOOLS_VERSION " https://github.com/freeseek/mocha BCFtools %s HTSlib %s\n",
            bcftools_version(), hts_version());
    if (!in || !out) error("Expected input VCF\n%s", usage());
    args->flags |= parse_tags(tag_list);

    // this ugly workaround is required to make sure we can set samples on both headers even
    // when sample_is_file is true and sample_names is stdin
    if (sample_names) {
        int nsmpl;
        char **smpl = hts_readlist(sample_names[0] == '^' ? sample_names + 1 : sample_names, sample_is_file, &nsmpl);
        if (!smpl) error("Failed to parse %s", sample_names);
        kstring_t tmp = {0, 0, NULL};
        if (nsmpl) {
            if (sample_names[0] == '^')
                ksprintf(&tmp, "^%s", smpl[0]);
            else
                ksprintf(&tmp, "%s", smpl[0]);
            for (i = 1; i < nsmpl; i++) ksprintf(&tmp, ",%s", smpl[i]);
        }
        int ret = bcf_hdr_set_samples(args->in_hdr, tmp.s, 0);
        if (ret < 0)
            error("Error parsing the sample list\n");
        else if (ret > 0) {
            if (force_samples)
                fprintf(stderr,
                        "Warn: subset called for sample that does not exist in header: "
                        "\"%s\"... skipping\n",
                        smpl[ret - 1]);
            else
                error(
                    "Error: subset called for sample that does not exist in header: \"%s\". "
                    "Use \"--force-samples\" to ignore this error.\n",
                    smpl[ret - 1]);
        }
        if (bcf_hdr_nsamples(args->in_hdr) == 0) fprintf(stderr, "Warn: subsetting has removed all samples\n");
        if (bcf_hdr_set_samples(args->out_hdr, tmp.s, 0) < 0) error("Error parsing the sample list\n");
        free(tmp.s);
        for (i = 0; i < nsmpl; i++) free(smpl[i]);
        free(smpl);
    }

    if (gender_fname) args->gender = read_computed_gender(args->in_hdr, gender_fname);

    if (args->flags & (INFO_GC | INFO_CPG)) {
        if (!ref_fname) error("Reference sequence not provided, cannot infer GC or CpG content\n");
        args->fai = fai_load(ref_fname);
        if (!args->fai) error("Failed to load the fai index: %s\n", ref_fname);
        bcf_hdr_append(args->out_hdr,
                       "##INFO=<ID=GC,Number=1,Type=Float,Description=\"GC ratio content around the marker\">");
        bcf_hdr_append(args->out_hdr,
                       "##INFO=<ID=CpG,Number=1,Type=Float,Description=\"CpG ratio content around the marker\">");
    }

    args->allele_a_id = bcf_hdr_id2int(args->in_hdr, BCF_DT_ID, "ALLELE_A");
    if (!bcf_hdr_idinfo_exists(args->in_hdr, BCF_HL_INFO, args->allele_a_id)) args->allele_a_id = -1;
    if (args->allele_a_id >= 0 && bcf_hdr_id2type(args->in_hdr, BCF_HL_INFO, args->allele_a_id) != BCF_HT_INT)
        error("Error: input VCF file ALLELE_A info field is not of integer type\n");

    args->allele_b_id = bcf_hdr_id2int(args->in_hdr, BCF_DT_ID, "ALLELE_B");
    if (!bcf_hdr_idinfo_exists(args->in_hdr, BCF_HL_INFO, args->allele_a_id)) args->allele_b_id = -1;
    if (args->allele_b_id >= 0 && bcf_hdr_id2type(args->in_hdr, BCF_HL_INFO, args->allele_b_id) != BCF_HT_INT)
        error("Error: input VCF file ALLELE_B info field is not of float type\n");

    args->adjust_id = bcf_hdr_id2int(args->in_hdr, BCF_DT_ID, "ADJ_COEFF");
    if (!bcf_hdr_idinfo_exists(args->in_hdr, BCF_HL_INFO, args->adjust_id)) args->adjust_id = -1;
    if (args->adjust_id >= 0 && bcf_hdr_id2type(args->in_hdr, BCF_HL_INFO, args->adjust_id) != BCF_HT_REAL)
        error("Error: input VCF file ADJ_COEFF info field is not of integer type\n");
    if (args->adjust_id >= 0
        && (bcf_hdr_id2length(args->in_hdr, BCF_HL_INFO, args->adjust_id) != BCF_VL_FIXED
            || bcf_hdr_id2number(args->in_hdr, BCF_HL_INFO, args->adjust_id) != 9))
        error("Error: input VCF file ADJ_COEFF info field has wrong number of values\n");

    if (args->adjust && args->adjust_id < 0)
        error("Error: ADJ_COEFF field is not present, cannot perform --adjust correction\n");

    args->nsmpl = bcf_hdr_nsamples(args->in_hdr);

    args->gt_id = bcf_hdr_id2int(args->in_hdr, BCF_DT_ID, "GT");
    if (!bcf_hdr_idinfo_exists(args->in_hdr, BCF_HL_FMT, args->gt_id)) args->gt_id = -1;

    args->ad_id = bcf_hdr_id2int(args->in_hdr, BCF_DT_ID, "AD");
    if (!bcf_hdr_idinfo_exists(args->in_hdr, BCF_HL_FMT, args->ad_id)) args->ad_id = -1;
    if (args->ad_id >= 0 && bcf_hdr_id2type(args->in_hdr, BCF_HL_FMT, args->ad_id) != BCF_HT_INT)
        error("Error: input VCF file AD format field is not of integer type\n");
    if (args->ad_id >= 0 && bcf_hdr_id2length(args->in_hdr, BCF_HL_FMT, args->ad_id) != BCF_VL_R)
        error("Error: input VCF file AD format field has wrong number of values\n");

    args->baf_id = bcf_hdr_id2int(args->in_hdr, BCF_DT_ID, "BAF");
    if (!bcf_hdr_idinfo_exists(args->in_hdr, BCF_HL_FMT, args->baf_id)) args->baf_id = -1;
    if (args->baf_id >= 0 && bcf_hdr_id2type(args->in_hdr, BCF_HL_FMT, args->baf_id) != BCF_HT_REAL)
        error("Error: input VCF file BAF format field is not of float type\n");

    args->lrr_id = bcf_hdr_id2int(args->in_hdr, BCF_DT_ID, "LRR");
    if (!bcf_hdr_idinfo_exists(args->in_hdr, BCF_HL_FMT, args->lrr_id)) args->lrr_id = -1;
    if (args->lrr_id >= 0 && bcf_hdr_id2type(args->in_hdr, BCF_HL_FMT, args->lrr_id) != BCF_HT_REAL)
        error("Error: input VCF file LRR format field is not of float type\n");

    if (args->flags & (INFO_ALLELE_A | INFO_ALLELE_B)) {
        if (!bcf_hdr_idinfo_exists(args->in_hdr, BCF_HL_FMT, args->gt_id))
            error("Error: GT format field is not present, cannot infer ALLELE_A or ALLELE_B\n");
        if (!bcf_hdr_idinfo_exists(args->in_hdr, BCF_HL_FMT, args->baf_id))
            error("Error: BAF format field is not present, cannot infer ALLELE_A or ALLELE_B\n");
        if (args->flags & INFO_ALLELE_A) {
            if (bcf_hdr_id2int(args->in_hdr, BCF_DT_ID, "ALLELE_A") >= 0)
                error("Field ALLELE_A already present in the VCF.\n");
            bcf_hdr_append(args->out_hdr, "##INFO=<ID=ALLELE_A,Number=1,Type=Integer,Description=\"A allele\">");
        }
        if (args->flags & INFO_ALLELE_B) {
            if (bcf_hdr_id2int(args->in_hdr, BCF_DT_ID, "ALLELE_B") >= 0)
                error("Field ALLELE_B already present in the VCF.\n");
            if (args->flags & INFO_ALLELE_B)
                bcf_hdr_append(args->out_hdr, "##INFO=<ID=ALLELE_B,Number=1,Type=Integer,Description=\"B allele\">");
        }
    }

    if (args->flags & (INFO_AC_SEX)) {
        if (!args->gender) error("Error: AC_Sex require --sex\n");
        bcf_hdr_append(args->out_hdr,
                       "##INFO=<ID=AC_Sex,Number=4,Type=Integer,Description=\"Number of homozygous and heterozygous "
                       "genotypes by gender\">");
    }

    if (args->flags & (INFO_MISSING_SEX)) {
        if (!args->gender) error("Error: MISSING_Sex require --sex\n");
        bcf_hdr_append(args->out_hdr,
                       "##INFO=<ID=MISSING_Sex,Number=4,Type=Integer,Description=\"Number of non-missing and missing "
                       "genotypes by gender\">");
    }

    if (args->flags & INFO_pAC_HET) {
        if (args->gt_id < 0) error("Error: GT format field is not present, cannot compute pAC_Het\n");
        bcf_hdr_append(args->out_hdr,
                       "##INFO=<ID=pAC_Het,Number=2,Type=Integer,Description=\"Number of heterozygous genotypes "
                       "by transmission type\">");
    }

    if (args->flags & INFO_AD_HET) {
        if (args->gt_id < 0) error("Error: GT format field is not present, cannot compute AD_Het\n");
        if (args->ad_id < 0) error("Error: AD format field is not present, cannot compute AD_Het\n");
        bcf_hdr_append(args->out_hdr,
                       "##INFO=<ID=AD_Het,Number=2,Type=Integer,Description=\"Allelic depths for the reference and "
                       "alternate alleles across heterozygous genotypes\">");
    }

    if (args->flags & INFO_pBAF_STATS) {
        if (args->gt_id < 0) error("Error: GT format field is not present, cannot compute pBAF_Stats\n");
        if (args->ad_id < 0 && args->baf_id < 0)
            error("Error: Neither AD nor BAF format fields are present, cannot compute pBAF_Stats\n");
        bcf_hdr_append(
            args->out_hdr,
            "##INFO=<ID=pBAF_Stats,Number=4,Type=Float,Description=\"Welch'"
            "s t-test and Mann-Whitney U test for allelic transmission ratios across heterozygous genotypes\">");
    }

    if (args->flags & INFO_COR_BAF_LRR) {
        if (args->allele_a_id < 0)
            error("Error: ALLELE_A field is not present, cannot perform --cor-BAF-LRR analysis\n");
        if (args->allele_b_id < 0)
            error("Error: ALLELE_B field is not present, cannot perform --cor-BAF-LRR analysis\n");
        if (args->baf_id < 0) error("Error: BAF format is not present, cannot perform --cor-BAF-LRR analysis\n");
        if (args->lrr_id < 0) error("Error: LRR format is not present, cannot perform --cor-BAF-LRR analysis\n");
        bcf_hdr_append(args->out_hdr,
                       "##INFO=<ID=Cor_BAF_LRR,Number=3,Type=Float,Description=\"Pearson "
                       "correlation for BAF and LRR at AA, AB, and BB genotypes\">");
    }

    args->ds_id = bcf_hdr_id2int(args->in_hdr, BCF_DT_ID, "DS");
    args->hds_id = bcf_hdr_id2int(args->in_hdr, BCF_DT_ID, "HDS");
    args->ap1_id = bcf_hdr_id2int(args->in_hdr, BCF_DT_ID, "AP1");
    args->ap2_id = bcf_hdr_id2int(args->in_hdr, BCF_DT_ID, "AP2");
    args->gp_id = bcf_hdr_id2int(args->in_hdr, BCF_DT_ID, "GP");
    if (args->flags & INFO_MACH) {
        if (!args->dosage_tag) {
            if (bcf_hdr_idinfo_exists(args->in_hdr, BCF_HL_FMT, args->hds_id)) args->dosage_tag = SCORE_HDS;
            if (bcf_hdr_idinfo_exists(args->in_hdr, BCF_HL_FMT, args->ap1_id)
                && bcf_hdr_idinfo_exists(args->in_hdr, BCF_HL_FMT, args->ap2_id))
                args->dosage_tag = SCORE_AP;
            if (bcf_hdr_idinfo_exists(args->in_hdr, BCF_HL_FMT, args->gp_id)) args->dosage_tag = SCORE_GP;
            if (bcf_hdr_idinfo_exists(args->in_hdr, BCF_HL_FMT, args->ds_id)) args->dosage_tag = SCORE_DS;
            if (!args->dosage_tag)
                error("VCF file %s does not include any of the DS, HDS, AP1/AP2, or DS FORMAT fields\n", argv[optind]);
        } else {
            switch (args->dosage_tag) {
            case SCORE_DS:
                if (!bcf_hdr_idinfo_exists(args->in_hdr, BCF_HL_FMT, args->ds_id))
                    error("VCF file %s does not include the DS FORMAT field\n", argv[optind]);
                break;
            case SCORE_HDS: // only for Minimac4 VCFs
                if (!bcf_hdr_idinfo_exists(args->in_hdr, BCF_HL_FMT, args->hds_id))
                    error("VCF file %s does not include the HDS FORMAT field\n", argv[optind]);
                break;
            case SCORE_AP:
                if (!bcf_hdr_idinfo_exists(args->in_hdr, BCF_HL_FMT, args->ap1_id)
                    || !bcf_hdr_idinfo_exists(args->in_hdr, BCF_HL_FMT, args->ap2_id))
                    error("VCF file %s does not include either the AP1 or the AP2 FORMAT fields\n", argv[optind]);
                break;
            case SCORE_GP:
                if (!bcf_hdr_idinfo_exists(args->in_hdr, BCF_HL_FMT, args->gp_id))
                    error("VCF file %s does not include the GP FORMAT field\n", argv[optind]);
                break;
            }
        }
        fprintf(stderr, "Using %s to compute genotype dosage count\n",
                args->dosage_tag == SCORE_GP    ? "genotype probabilities (GP)"
                : args->dosage_tag == SCORE_AP  ? "ALT haplotype probabilities (AP)"
                : args->dosage_tag == SCORE_HDS ? "haploid alternate allele dosage (HDS)"
                                                : "genotype dosages (DS)");
        bcf_hdr_append(
            args->out_hdr,
            "##INFO=<ID=MACH,Number=3,Type=Float,Description=\"Statistics to infer the MACH r2 imputation measure\">");
    }

    if (args->as_str) {
        args->as_id = bcf_hdr_id2int(args->in_hdr, BCF_DT_ID, args->as_str);
        if (!bcf_hdr_idinfo_exists(args->in_hdr, BCF_HL_FMT, args->as_id))
            error("Error: %s format field is not present, cannot perform --summary\n", args->as_str);
        if (args->phase) kputc('p', &args->summary_str);
        kputs(args->as_str, &args->summary_str);
        ksprintf(&str, "##INFO=<ID=%s,Number=2,Type=Integer,Description=\"%sAllelic shift counts for %s\">",
                 args->summary_str.s, args->phase ? "phased " : "", args->as_str);
        bcf_hdr_append(args->out_hdr, str.s);
    }

    // this analysis is performed using a field that can be added by the previous analyses
    if (args->info_str) {
        args->info_id = bcf_hdr_id2int(args->out_hdr, BCF_DT_ID, args->info_str);
        if (bcf_hdr_sync(args->out_hdr) < 0) error_errno("[%s] Failed to update header", __func__);
        if (!bcf_hdr_idinfo_exists(args->out_hdr, BCF_HL_INFO, args->info_id))
            error("Error: %s info field is not present, cannot perform --test\n", args->info_str);
        int number = bcf_hdr_id2number(args->out_hdr, BCF_HL_INFO, args->info_id);
        if (args->phred) kputc('p', &args->test_str);
        if (number == 2)
            kputs("binom_", &args->test_str);
        else if (number == 4)
            kputs("fisher_", &args->test_str);
        else
            error("Error: %s info field must contain 2 or 4 elements but it contains %d\n", args->info_str, number);
        kputs(args->info_str, &args->test_str);
        str.l = 0;
        ksprintf(&str, "##INFO=<ID=%s,Number=1,Type=Float,Description=\"%s%s test for %s\">", args->test_str.s,
                 args->phred ? "Phred scaled " : "", number == 2 ? "Binomial" : "Fisher's exact", args->info_str);
        bcf_hdr_append(args->out_hdr, str.s);
        if (args->odds_ratio) {
            kputs("odds_ratio_", &args->odds_ratio_str);
            kputs(args->info_str, &args->odds_ratio_str);
            str.l = 0;
            ksprintf(&str, "##INFO=<ID=%s,Number=1,Type=Float,Description=\"Odds ratio for %s\">",
                     args->odds_ratio_str.s, args->info_str);
            bcf_hdr_append(args->out_hdr, str.s);
        }
        if (args->log_odds) {
            kputs("log_odds_", &args->log_odds_str);
            kputs(args->info_str, &args->log_odds_str);
            str.l = 0;
            ksprintf(&str, "##INFO=<ID=%s,Number=1,Type=Float,Description=\"Log odds ratio for %s\">",
                     args->log_odds_str.s, args->info_str);
            bcf_hdr_append(args->out_hdr, str.s);
        }
    }

    if (sites_only) {
        if (bcf_hdr_set_samples(args->out_hdr, NULL, 0) < 0) error("Error parsing the sample list\n");
        bcf_hdr_remove(args->out_hdr, BCF_HL_FMT, NULL);
    }

    args->gt_phase_arr = (int8_t *)malloc(args->nsmpl * sizeof(int8_t));
    args->as_arr = (int8_t *)malloc(args->nsmpl * sizeof(int8_t));
    args->gt0_arr = (int16_t *)malloc(args->nsmpl * sizeof(int16_t));
    args->gt1_arr = (int16_t *)malloc(args->nsmpl * sizeof(int16_t));
    args->ad0_arr = (int16_t *)malloc(args->nsmpl * sizeof(int16_t));
    args->ad1_arr = (int16_t *)malloc(args->nsmpl * sizeof(int16_t));
    args->baf_arr[0] = (float *)malloc(args->nsmpl * sizeof(float));
    args->baf_arr[1] = (float *)malloc(args->nsmpl * sizeof(float));
    args->imap_arr = (int *)malloc(args->nsmpl * sizeof(int));

    free(str.s);
    return 0;
}

// Heng Li's implementation in htslib/kfunc.c
#define KF_GAMMA_EPS 1e-14
#define KF_TINY 1e-290

static double log_kf_betai_aux(double a, double b, double x) {
    double C, D, f;
    int j;
    if (x == 0.) return 0.;
    if (x == 1.) return 1.;
    f = 1.;
    C = f;
    D = 0.;
    // Modified Lentz's algorithm for computing continued fraction
    for (j = 1; j < 200; ++j) {
        double aa, d;
        int m = j >> 1;
        aa = (j & 1) ? -(a + m) * (a + b + m) * x / ((a + 2 * m) * (a + 2 * m + 1))
                     : m * (b - m) * x / ((a + 2 * m - 1) * (a + 2 * m));
        D = 1. + aa * D;
        if (D < KF_TINY) D = KF_TINY;
        C = 1. + aa / C;
        if (C < KF_TINY) C = KF_TINY;
        D = 1. / D;
        d = C * D;
        f *= d;
        if (fabs(d - 1.) < KF_GAMMA_EPS) break;
    }
    return (kf_lgamma(a + b) - kf_lgamma(a) - kf_lgamma(b) + a * log(x) + b * log(1. - x)) - log(a * f);
}

// returns -10 * (log(2) + pbinom(min(k, n - k), n, 1/2, log.p = TRUE)) / log(10)
static double phred_pbinom(int k, int n) {
    int i, j;
    static double *dbinom = NULL, *pbinom = NULL;
    static size_t n_size = 0, m_dbinom = 0, m_pbinom = 0;

    if (n < 0 && k < 0) {
        free(dbinom);
        free(pbinom);
        return NAN;
    }

    if (n < 0 || k < 0 || k > n) return NAN;

    if (k == n >> 1) return 0.0;

    if (k << 1 > n) k = n - k;

    if (n > 1000) return -10.0 * M_LOG10E * (M_LN2 + log_kf_betai_aux(n - k, k + 1, .5));

    if (n >= n_size) {
        size_t len = (size_t)(1 + (1 + (n >> 1)) * ((n + 1) >> 1));
        hts_expand(double, len, m_dbinom, dbinom);
        hts_expand(double, len, m_pbinom, pbinom);
        dbinom[0] = 1.0;
        for (i = n_size ? (int)n_size : 1; i < n + 1; i++) {
            int prev_idx = i - 1 ? 1 + ((i - 1) >> 1) * (i >> 1) : 0;
            int curr_idx = 1 + (i >> 1) * ((i + 1) >> 1);
            dbinom[curr_idx] = dbinom[prev_idx] * 0.5;
            pbinom[curr_idx] = dbinom[curr_idx];
            for (j = 1; j < ((i + 1) >> 1); j++) {
                curr_idx++;
                dbinom[curr_idx] = (double)i / (double)j * dbinom[prev_idx] * 0.5;
                pbinom[curr_idx] = pbinom[curr_idx - 1] + dbinom[curr_idx];
                prev_idx++;
            }
        }
        n_size = (size_t)(n + 1);
    }

    int idx = 1 + (n >> 1) * ((n + 1) >> 1) + k;
    return -10.0 * M_LOG10E * (M_LN2 + log(pbinom[idx]));
}

static int sample_mean_var(const float *x, int n, double *xm, double *xss) {
    if (n < 2) return -1;
    *xm = 0;
    *xss = 0;
    int i, j = 0;
    for (i = 0; i < n; i++) {
        if (!isnan(x[i])) {
            *xm += (double)x[i];
            *xss += sq((double)x[i]);
            j++;
        }
    }
    if (j <= 1) return -1;
    *xm /= (double)j;
    *xss -= sq(*xm) * (double)j;
    *xss /= (double)(j - 1);
    return 0;
}

static double welch_t_test(float *a, float *b, int na, int nb) {
    double mua, mub, sa2, sb2, t, v;
    if (na < 2 || nb < 2) return HUGE_VAL;
    sample_mean_var(a, na, &mua, &sa2);
    sample_mean_var(b, nb, &mub, &sb2);
    t = (mua - mub) / sqrt(sa2 / na + sb2 / nb);
    v = (sa2 / na + sb2 / nb);
    v *= v;
    v /= sq(sa2) / na / na / (na - 1) + sq(sb2) / nb / nb / (nb - 1);
    return kf_betai(v / 2.0f, 0.5, v / (v + sq(t)));
}

// Petr Danecek's and James Bonfield's implementation in bcftools/bam2bcf.c
double mann_whitney_1947_cdf(int n, int m, int U);

static int cmpfunc(const void *a, const void *b) { return (*(float *)a > *(float *)b) - (*(float *)a < *(float *)b); }

// it currently does not handle nans
// adapted from Petr Danecek's implementation of calc_mwu_bias_cdf() in bcftools/bam2bcf.c
static double mann_whitney_u(float *a, float *b, int na, int nb) {
    qsort(a, (size_t)na, sizeof(float), cmpfunc);
    qsort(b, (size_t)nb, sizeof(float), cmpfunc);

    int i = 0, j = 0, ca, cb;
    double U = 0, ties = 0;
    while (i < na || j < nb) {
        double curr = (j == nb || (i < na && a[i] < b[j])) ? a[i] : b[j];
        for (ca = 0; i < na && a[i] == curr; i++) ca++;
        for (cb = 0; j < nb && b[j] == curr; j++) cb++;
        U += ca * (j - cb * 0.5);
        if (ca && cb) {
            double tie = ca + cb;
            ties += (sq(tie) - 1.0) * tie;
        }
    }
    if (!na || !nb) return HUGE_VAL;

    double U_min = ((double)na * nb) - U;
    if (U < U_min) U_min = U;

    if (na == 1) return 2.0 * (floor(U_min) + 1.0) / (double)(nb + 1);
    if (nb == 1) return 2.0 * (floor(U_min) + 1.0) / (double)(na + 1);

    // Normal approximation, very good for na>=8 && nb>=8 and reasonable if na<8 or nb<8
    if (na >= 8 || nb >= 8) {
        double mean = ((double)na * nb) * 0.5;
        // Correction for ties:
        double N = na + nb;
        double var2 = (sq(N) - 1) * N - ties;
        if (var2 == 0) return 1.0;
        var2 *= ((double)na * nb) / N / (N - 1) / 12.0;
        // No correction for ties:
        // float var2 = ((double)na*nb) * (na + nb + 1) / 12.0;
        double z = (U_min - mean) / sqrt(2.0 * var2); // z is N(0,1)
        // return 2.0 - kf_erfc(z);  // which is 1 + erf(z)
        return kf_erfc(-z); // which is 1 - erf(-z)
    }

    // Exact calculation
    double pval = 2.0 * mann_whitney_1947_cdf(na, nb, (int)U_min);
    return pval > 1.0 ? 1.0 : pval;
}

static inline int bcf_int8_is_vector_end(int8_t value) { return value == bcf_int8_vector_end; }
static inline int bcf_int16_is_vector_end(int16_t value) { return value == bcf_int16_vector_end; }
static inline int bcf_int32_is_vector_end(int32_t value) { return value == bcf_int32_vector_end; }
static inline int bcf_int8_is_missing(int8_t value) { return value == bcf_int8_missing; }
static inline int bcf_int16_is_missing(int16_t value) { return value == bcf_int16_missing; }
static inline int bcf_int32_is_missing(int32_t value) { return value == bcf_int32_missing; }

// retrieve phase information from BCF record
// assumes little endian architecture
static int bcf_get_format_sign(bcf_fmt_t *fmt, int8_t *as_arr, int nsmpl) {
    if (!fmt || fmt->n != 1) return 0;
    int i;

#define BRANCH(type_t, is_vector_end, is_missing)                                                                      \
    {                                                                                                                  \
        type_t *p = (type_t *)fmt->p;                                                                                  \
        for (i = 0; i < nsmpl; i++) {                                                                                  \
            if (is_vector_end(p[i]) || is_missing(p[i]))                                                               \
                as_arr[i] = bcf_int8_missing;                                                                          \
            else if (p[i] == (type_t)0)                                                                                \
                as_arr[i] = (int8_t)0;                                                                                 \
            else if (p[i] > (type_t)0)                                                                                 \
                as_arr[i] = (int8_t)1;                                                                                 \
            else                                                                                                       \
                as_arr[i] = (int8_t) - 1;                                                                              \
        }                                                                                                              \
    }
    switch (fmt->type) {
    case BCF_BT_INT8:
        BRANCH(int8_t, bcf_int8_is_vector_end, bcf_int8_is_missing);
        break;
    case BCF_BT_INT16:
        BRANCH(int16_t, bcf_int16_is_vector_end, bcf_int16_is_missing);
        break;
    case BCF_BT_INT32:
        BRANCH(int32_t, bcf_int32_is_vector_end, bcf_int32_is_missing);
        break;
    case BCF_BT_FLOAT:
        BRANCH(int32_t, bcf_float_is_vector_end, bcf_float_is_missing);
        break;
    default:
        error("Unexpected type %d\n", fmt->type);
    }
#undef BRANCH

    return 1;
}

bcf1_t *process(bcf1_t *rec) {
    int i, j, k;
    // compute GC and CpG content for each site
    if (args->flags & (INFO_GC | INFO_CPG)) {
        int fa_len;
        int at_cnt = 0, cg_cnt = 0, cpg_cnt = 0;
        const char *ref = rec->d.allele[0];
        char *fa = faidx_fetch_seq(args->fai, bcf_seqname(args->in_hdr, rec), rec->pos - args->gc_win,
                                   rec->pos + (int)strlen(ref) - 1 + args->gc_win, &fa_len);
        if (!fa)
            error("fai_fetch_seq failed at %s:%" PRId64 "\n", bcf_hdr_id2name(args->in_hdr, rec->rid), rec->pos + 1);
        for (i = 0; i < fa_len; i++) {
            if (fa[i] > 96) fa[i] = (char)(fa[i] - 32);
            if (fa[i] == 'A' || fa[i] == 'T') at_cnt++;
            if (fa[i] == 'C' || fa[i] == 'G') cg_cnt++;
            if (i > 0)
                if (fa[i - 1] == 'C' && fa[i] == 'G') cpg_cnt += 2;
        }
        free(fa);
        float ratio;
        if (args->flags & INFO_GC) {
            ratio = (float)(cg_cnt) / (float)(at_cnt + cg_cnt);
            bcf_update_info_float(args->out_hdr, rec, "GC", &ratio, 1);
        }
        if (args->flags & INFO_CPG) {
            ratio = (float)cpg_cnt / (float)(fa_len);
            bcf_update_info_float(args->out_hdr, rec, "CpG", &ratio, 1);
        }
    }

    // if no samples present in the

    int ac_sex[] = {0, 0, 0, 0};
    int missing_sex[] = {0, 0, 0, 0};
    int pac_het[] = {0, 0};
    int ad_het[] = {0, 0};
    int summary[] = {0, 0};
    int psummary[] = {0, 0};
    float ret[4];

    // if samples are not present, skip to the end
    if (args->nsmpl == 0) goto ret;

    // extract format information from VCF format records
    bcf_fmt_t *gt_fmt = bcf_get_fmt_id(rec, args->gt_id);
    int gt_phase = bcf_get_genotype_phase(gt_fmt, args->gt_phase_arr, args->nsmpl);

    // if genotypes are not present, skip to the end
    if (!bcf_get_unphased_genotype_alleles(gt_fmt, args->gt0_arr, args->gt1_arr, args->nsmpl)) goto ret;

    bcf_fmt_t *as_fmt = bcf_get_fmt_id(rec, args->as_id);
    int as_sign = (args->as_str && as_fmt) ? bcf_get_format_sign(as_fmt, args->as_arr, args->nsmpl) : 0;

    bcf_fmt_t *ad_fmt = bcf_get_fmt_id(rec, args->ad_id);
    int ad = ad_fmt && ad_fmt->n == rec->n_allele ? bcf_get_allelic_depth(ad_fmt, args->gt0_arr, args->gt1_arr,
                                                                          args->ad0_arr, args->ad1_arr, args->nsmpl)
                                                  : 0;

    bcf_fmt_t *baf_fmt = bcf_get_fmt_id(rec, args->baf_id);
    int baf = baf_fmt && baf_fmt->n == 1 && baf_fmt->type == BCF_BT_FLOAT ? 1 : 0;

    bcf_fmt_t *lrr_fmt = bcf_get_fmt_id(rec, args->lrr_id);
    int lrr = lrr_fmt && lrr_fmt->n == 1 && lrr_fmt->type == BCF_BT_FLOAT ? 1 : 0;

    if ((args->flags & (INFO_ALLELE_A | INFO_ALLELE_B)) && baf) {
        int alleles[2];
        switch (rec->n_allele) {
        case 1:
            alleles[0] = -1;
            alleles[1] = -1;
            break;
        case 2:
            alleles[0] = 0;
            alleles[1] = 1;
            break;
        case 3:
            alleles[0] = 1;
            alleles[1] = 2;
            break;
        default:
            error("Observed wrong number of alleles at %s:%" PRId64 "\n", bcf_hdr_id2name(args->in_hdr, rec->rid),
                  rec->pos + 1);
        }
        float median[2] = {NAN, NAN};
        int alleles_idx[2] = {-1, -1};
        for (i = 0; i < 2; i++) {
            int n = 0;
            for (j = 0; j < args->nsmpl; j++)
                if (args->gt0_arr[j] == alleles[i] && args->gt1_arr[j] == alleles[i]) args->imap_arr[n++] = j;
            median[i] = get_median((float *)baf_fmt->p, n, args->imap_arr);
            if (median[i] < .5)
                alleles_idx[i] = alleles[0];
            else if (median[i] > .5)
                alleles_idx[i] = alleles[1];
        }
        if (alleles_idx[0] == alleles_idx[1]) {
            alleles_idx[0] = -1;
            alleles_idx[1] = -1;
            fprintf(stderr, "Unable to infer the A and B alleles while parsing the site %s:%" PRId64 "\n",
                    bcf_hdr_id2name(args->in_hdr, rec->rid), rec->pos + 1);
        } else if (alleles_idx[0] == -1) {
            alleles_idx[0] = alleles_idx[1] == alleles[0] ? alleles[1] : alleles[0];
        } else if (alleles_idx[1] == -1) {
            alleles_idx[1] = alleles_idx[0] == alleles[0] ? alleles[1] : alleles[0];
        }
        if (args->flags & INFO_ALLELE_A) bcf_update_info_int32(args->out_hdr, rec, "ALLELE_A", &alleles_idx[0], 1);
        if (args->flags & INFO_ALLELE_B) bcf_update_info_int32(args->out_hdr, rec, "ALLELE_B", &alleles_idx[1], 1);
    }

    // determine ALLELE_A and ALLELE_B from VCF record
    bcf_info_t *allele_a_info = bcf_get_info_id(rec, args->allele_a_id);
    int8_t allele_a = allele_a_info ? ((int8_t *)allele_a_info->vptr)[0] : -1;
    bcf_info_t *allele_b_info = bcf_get_info_id(rec, args->allele_b_id);
    int8_t allele_b = allele_b_info ? ((int8_t *)allele_b_info->vptr)[0] : -1;

    // adjust BAF and LRR values
    if (args->adjust && baf && lrr) {
        if (bcf_get_info_float(args->in_hdr, rec, "ADJ_COEFF", &args->adjust_arr, &args->m_adjust) < 0)
            error("Could not retrieve adjusting coefficients at %s:%" PRId64 "\n",
                  bcf_hdr_id2name(args->in_hdr, rec->rid), rec->pos + 1);
        if (args->m_adjust < 9)
            error("Not enough adjusting coefficients at %s:%" PRId64 "\n", bcf_hdr_id2name(args->in_hdr, rec->rid),
                  rec->pos + 1);
        for (i = 0; i < args->nsmpl; i++) {
            int n_a = (args->gt0_arr[i] == allele_a) + (args->gt1_arr[i] == allele_a);
            int n_b = (args->gt0_arr[i] == allele_b) + (args->gt1_arr[i] == allele_b);
            if (n_a + n_b != 2) continue;
            ((float *)baf_fmt->p)[i] -=
                args->adjust_arr[3 * n_b + 1] * ((float *)lrr_fmt->p)[i] - args->adjust_arr[3 * n_b];
            ((float *)lrr_fmt->p)[i] -= args->adjust_arr[3 * n_b + 2];
        }
    }

    for (i = 0; i < args->nsmpl; i++) {
        float curr_baf = NAN;

        int is_missing = args->gt0_arr[i] == bcf_int16_missing || args->gt1_arr[i] == bcf_int16_missing;
        if (args->gender) {
            switch (args->gender[i]) {
            case GENDER_MALE: // male
                missing_sex[is_missing]++;
                break;
            case GENDER_FEMALE: // female
                missing_sex[2 + is_missing]++;
                break;
            default:
                break;
            }
        }

        // if genotype is missing, skip
        if (is_missing) continue;

        int idx_as_sign =
            (as_sign && args->as_arr[i] != bcf_int8_missing && args->as_arr[i] != 0) ? (1 - args->as_arr[i]) / 2 : -1;
        if (idx_as_sign >= 0) summary[idx_as_sign]++;

        int is_het = args->gt0_arr[i] != args->gt1_arr[i];
        if (args->gender) {
            switch (args->gender[i]) {
            case GENDER_MALE: // male
                ac_sex[is_het]++;
                break;
            case GENDER_FEMALE: // female
                ac_sex[2 + is_het]++;
                break;
            default:
                break;
            }
        }

        // if genotype is not heterozygous reference / alternate, skip
        if (!is_het || (args->gt0_arr[i] != 0 && args->gt1_arr[i] != 0)) continue;

        int idx_gt_phase = gt_phase && (args->gt_phase_arr[i] == -1 || args->gt_phase_arr[i] == 1)
                               ? (1 - args->gt_phase_arr[i]) / 2
                               : -1;
        if (idx_gt_phase >= 0) pac_het[idx_gt_phase]++;

        int idx_fmt_phase =
            (idx_gt_phase >= 0 && idx_as_sign >= 0) ? (1 - args->as_arr[i] * args->gt_phase_arr[i]) / 2 : -1;
        if (idx_fmt_phase >= 0) psummary[idx_fmt_phase]++;

        if (ad) {
            int ref_cnt = args->ad0_arr[i];
            int alt_cnt = args->ad1_arr[i];
            ad_het[0] += ref_cnt;
            ad_het[1] += alt_cnt;
            curr_baf = ((float)alt_cnt + 0.5f) / ((float)ref_cnt + (float)alt_cnt + 1.0f);
        }
        if (baf) curr_baf = ((float *)baf_fmt->p)[i];
        if (idx_gt_phase >= 0 && !isnan(curr_baf)) {
            args->baf_arr[idx_gt_phase][pac_het[idx_gt_phase] - 1] = curr_baf;
        }
    }

    if (args->flags & INFO_AC_SEX) bcf_update_info_int32(args->out_hdr, rec, "AC_Sex", &ac_sex, 4);
    if (args->flags & INFO_MISSING_SEX) bcf_update_info_int32(args->out_hdr, rec, "MISSING_Sex", &missing_sex, 4);
    if (args->flags & INFO_pAC_HET) bcf_update_info_int32(args->out_hdr, rec, "pAC_Het", &pac_het, 2);
    if ((args->flags & INFO_AD_HET) && ad) bcf_update_info_int32(args->out_hdr, rec, "AD_Het", &ad_het, 2);

    if ((args->flags & INFO_pBAF_STATS) && (ad || baf)) {
        ret[0] = get_median(args->baf_arr[0], pac_het[0], NULL);
        ret[1] = get_median(args->baf_arr[1], pac_het[1], NULL);
        ret[2] = welch_t_test(args->baf_arr[0], args->baf_arr[1], pac_het[0], pac_het[1]);
        ret[3] = mann_whitney_u(args->baf_arr[0], args->baf_arr[1], pac_het[0], pac_het[1]);
        bcf_update_info_float(args->out_hdr, rec, "pBAF_Stats", &ret, 4);
    }

    if ((args->flags & INFO_COR_BAF_LRR) && baf && lrr) {
        float rho[3];
        for (i = 0; i < 3; i++) {
            int n = 0;
            for (j = 0; j < args->nsmpl; j++) {
                int n_a = (args->gt0_arr[j] == allele_a) + (args->gt1_arr[j] == allele_a);
                int n_b = (args->gt0_arr[j] == allele_b) + (args->gt1_arr[j] == allele_b);
                if (n_a == 2 - i && n_b == i) args->imap_arr[n++] = j;
            }
            // compute the Pearson correlation
            float xss = 0.0f, yss = 0.0f, xyss = 0.0f;
            get_cov((float *)baf_fmt->p, (float *)lrr_fmt->p, n, args->imap_arr, &xss, &yss, &xyss);
            rho[i] = xyss / sqrtf(xss * yss);
        }
        bcf_update_info_float(args->out_hdr, rec, "Cor_BAF_LRR", &rho, 3);
    }

    if ((args->flags & INFO_MACH)) {
        float p1, p2, sum_ds[3] = {0.0f, 0.0f, 0.0f};
        bcf_fmt_t *ds_fmt, *hds_fmt, *ap1_fmt, *ap2_fmt, *gp_fmt;
        switch (args->dosage_tag) {
        case SCORE_DS:
            ds_fmt = bcf_get_fmt_id(rec, args->ds_id);
            if (ds_fmt) {
                for (k = 0; k < args->nsmpl; k++) {
                    p1 = ((float *)ds_fmt->p)[k];
                    if (bcf_float_is_missing(p1)) continue;
                    sum_ds[0]++;
                    sum_ds[1] += p1;
                    sum_ds[2] += p1 * p1;
                }
                bcf_update_info_float(args->out_hdr, rec, "MACH", &sum_ds, 3);
            }
            break;
        case SCORE_HDS: // only for Minimac4 VCFs
            hds_fmt = bcf_get_fmt_id(rec, args->hds_id);
            if (hds_fmt) {
                for (k = 0; k < args->nsmpl; k++) {
                    p1 = ((float *)hds_fmt->p)[2 * k];
                    p2 = ((float *)hds_fmt->p)[2 * k + 1];
                    if (bcf_float_is_missing(p1) || bcf_float_is_missing(p2)) continue;
                    sum_ds[0]++;
                    sum_ds[1] += p1 + p2;
                    sum_ds[2] += (p1 + p2) * (p1 + p2);
                }
                bcf_update_info_float(args->out_hdr, rec, "MACH", &sum_ds, 3);
            }
            break;
        case SCORE_AP:
            ap1_fmt = bcf_get_fmt_id(rec, args->ap1_id);
            ap2_fmt = bcf_get_fmt_id(rec, args->ap2_id);
            if (ap1_fmt && ap2_fmt) {
                for (k = 0; k < args->nsmpl; k++) {
                    p1 = ((float *)ap1_fmt->p)[k];
                    p2 = ((float *)ap2_fmt->p)[k];
                    if (bcf_float_is_missing(p1) || bcf_float_is_missing(p2)) continue;
                    sum_ds[0]++;
                    sum_ds[1] += p1 + p2;
                    sum_ds[2] += (p1 + p2) * (p1 + p2);
                }
                bcf_update_info_float(args->out_hdr, rec, "MACH", &sum_ds, 3);
            }
            break;
        case SCORE_GP:
            gp_fmt = bcf_get_fmt_id(rec, args->gp_id);
            if (gp_fmt) {
                for (k = 0; k < args->nsmpl; k++) {
                    p1 = ((float *)gp_fmt->p)[3 * k + 1];
                    p2 = ((float *)gp_fmt->p)[3 * k + 2];
                    if (bcf_float_is_missing(p1) || bcf_float_is_missing(p2)) continue;
                    sum_ds[0]++;
                    sum_ds[1] += p1 + 2.0f * p2;
                    sum_ds[2] += (p1 + 2.0f * p2) * (p1 + 2.0f * p2);
                }
                bcf_update_info_float(args->out_hdr, rec, "MACH", &sum_ds, 3);
            }
            break;
        }
    }

ret:
    if (args->as_str)
        bcf_update_info_int32(args->out_hdr, rec, args->summary_str.s, args->phase ? &psummary : &summary, 2);

    // perform binomial or Fisher's exact test on existing INFO field
    if (args->info_str) {
        bcf_info_t *info = bcf_get_info_id(rec, args->info_id);
        if (!(info->type & (BCF_BT_INT8 | BCF_BT_INT16 | BCF_BT_INT32)) || (info->len != 2 && info->len != 4))
            error("INFO field %s at %s:%" PRId64 " should contain two or four integers\n", args->info_str,
                  bcf_hdr_id2name(args->in_hdr, rec->rid), rec->pos + 1);
        bcf_get_info_int32(args->out_hdr, rec, args->info_str, &args->info_arr, &args->m_info);
        float value, odds_ratio;
        if (info->len == 2) {
            double phred_pval = phred_pbinom(args->info_arr[0], args->info_arr[0] + args->info_arr[1]);
            value = (float)(args->phred ? phred_pval : exp(-0.1 * M_LN10 * phred_pval));
            odds_ratio = (float)args->info_arr[1] / (float)args->info_arr[0];
        } else {
            double left, right, fisher;
            double pval = kt_fisher_exact(args->info_arr[0], args->info_arr[1], args->info_arr[2], args->info_arr[3],
                                          &left, &right, &fisher);
            value = (float)(args->phred ? -10.0 * M_LOG10E * log(pval) : pval);
            odds_ratio = ((float)args->info_arr[1] * (float)args->info_arr[2])
                         / ((float)args->info_arr[0] * (float)args->info_arr[3]);
        }
        bcf_update_info_float(args->out_hdr, rec, args->test_str.s, &value, 1);
        if (args->odds_ratio) bcf_update_info_float(args->out_hdr, rec, args->odds_ratio_str.s, &odds_ratio, 1);
        if (args->log_odds) {
            odds_ratio = logf(odds_ratio);
            bcf_update_info_float(args->out_hdr, rec, args->log_odds_str.s, &odds_ratio, 1);
        }
    }

    // remove all samples if sites_only was selected
    if (bcf_hdr_nsamples(args->out_hdr) == 0) bcf_subset(args->out_hdr, rec, 0, NULL);
    return rec;
}

void destroy(void) {
    phred_pbinom(-1, -1);
    free(args->gender);
    free(args->gt_phase_arr);
    free(args->as_arr);
    free(args->gt0_arr);
    free(args->gt1_arr);
    free(args->ad0_arr);
    free(args->ad1_arr);
    free(args->info_arr);
    free(args->adjust_arr);
    free(args->baf_arr[0]);
    free(args->baf_arr[1]);
    free(args->imap_arr);
    free(args->summary_str.s);
    free(args->test_str.s);
    free(args->odds_ratio_str.s);
    free(args->log_odds_str.s);
    free(args);
}
