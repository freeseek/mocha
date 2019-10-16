/* The MIT License

   Copyright (C) 2017-2018 Giulio Genovese

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
#include <htslib/kseq.h>
#include <htslib/faidx.h>
#include <htslib/kfunc.h>
#include "bcftools.h"

#define VERSION "2018-10-18"

static inline double sq(double x) { return x*x; }

typedef struct
{
    int phase; // whether to include genotype phase in the asymmetry test
    int ad; // whether to perform the allelic depth test
    int ws; // length of window to compute GC and CpG content
    char *format; // format field for sign balance test
    int nsmpl, gt_id, ad_id, baf_id, fmt_id;
    int *sex;
    int8_t *gt_phase_arr, *fmt_sign_arr;
    int16_t *gt0_arr, *gt1_arr, *ad0_arr, *ad1_arr;
    float *baf_arr[2];
    faidx_t *fai;
    bcf_hdr_t *in_hdr, *out_hdr;
}
args_t;

args_t *args;

const char *about(void)
{
    return "MOsaic CHromosomal Alterations tools.\n";
}

const char *usage(void)
{
    return
"\n"
"About: tools for the MOsaic CHromosomal Alterations pipeline. ("VERSION")\n"
"\n"
"Usage: bcftools +mochatools [General Options] -- [Plugin Options]\n"
"Options:\n"
"   run \"bcftools plugin\" for a list of common options\n"
"\n"
"Plugin options:\n"
"   -b, --balance <ID>            performs binomial test for sign balance of format field ID\n"
"   -p, --phase                   integrates genotype phase in the balance tests\n"
"   -a, --ad-het                  performs binomial test for reference / alternate allelic depth (AD)\n"
"   -x  --sex <file>              file including information about sex of samples\n"
"   -f, --fasta-ref <file>        reference sequence to compute GC and CpG content\n"
"   -w, --window-size <int>       Window size in bp used to compute the GC and CpG content [200]\n"
"   -s, --samples [^]<list>       comma separated list of samples to include (or exclude with \"^\" prefix)\n"
"   -S, --samples-file [^]<file>  file of samples to include (or exclude with \"^\" prefix)\n"
"       --force-samples           only warn about unknown subset samples\n"
"   -G, --drop-genotypes          drop individual genotype information (after running statistical tests)\n"
"\n"
"Example:\n"
"    bcftools +mochatools file.bcf -- --balance Bdev_Phase --drop-genotypes\n";
}

static void parse_sex(bcf_hdr_t *hdr,
                      char *fname,
                      int *sex)
{
    htsFile *fp = hts_open(fname, "r");
    if ( !fp ) error("Could not read: %s\n", fname);

    kstring_t str = {0, 0, NULL};
    if ( hts_getline(fp, KS_SEP_LINE, &str) <= 0 ) error("Empty file: %s\n", fname);

    // args->sex = (int *)calloc((size_t)bcf_hdr_nsamples(args->in_hdr), sizeof(int));

    int moff = 0, *off = NULL;
    char *tmp = NULL;
    do
    {
        int ncols = ksplit_core(str.s, 0, &moff, &off);
        if ( ncols<2 ) error("Could not parse the sex file: %s\n", str.s);

        int sample = bcf_hdr_id2int(hdr, BCF_DT_SAMPLE, &str.s[off[0]]);
        if ( sample>=0 )
        {
            if (str.s[off[1]]=='M') sex[sample] = 1;
            else if (str.s[off[1]]=='F') sex[sample] = 2;
            else sex[sample] = (int)strtol(&str.s[off[1]], &tmp, 0);
        }
    } while ( hts_getline(fp, KS_SEP_LINE, &str)>=0 );

    free(str.s);
    free(off);
    hts_close(fp);
}

int init(int argc,
         char **argv,
         bcf_hdr_t *in,
         bcf_hdr_t *out)
{
    args = (args_t *)calloc(1, sizeof(args_t));
    args->ws = 200;
    args->in_hdr = in;
    args->out_hdr = out;
    args->format = NULL;
    int sample_is_file = 0;
    int force_samples = 0;
    int sites_only = 0;
    char *sample_names = NULL;
    char *sex_fname = NULL;
    char *ref_fname = NULL;
    char *tmp = NULL;

    int c;
    static struct option loptions[] =
    {
        {"balance", required_argument, NULL, 'b'},
        {"ad-het", no_argument, NULL, 'a'},
        {"sex", required_argument, NULL, 'x'},
        {"phase", no_argument, NULL, 'p'},
        {"fasta-ref", required_argument, NULL, 'f'},
        {"window-size", required_argument, NULL, 'w'},
        {"samples", required_argument, NULL, 's'},
        {"samples-file", required_argument, NULL, 'S'},
        {"force-samples", no_argument, NULL, 1},
        {"drop-genotypes", no_argument,NULL, 'G'},
        {NULL,0,NULL,0}
    };

    while ((c = getopt_long(argc, argv, "h?b:ax:pf:w:s:S:G",loptions,NULL)) >= 0)
    {
        switch (c)
        {
            case 'b': args->format = optarg; break;
            case 'a': args->ad = 1; break;
            case 'x': sex_fname = optarg; break;
            case 'p': args->phase = 1; break;
            case 'f': ref_fname = optarg; break;
            case 'w':
                args->ws = (int)strtol(optarg, &tmp, 10);
                if ( *tmp ) error("Could not parse: -w %s\n", optarg);
                if( args->ws <= 0 ) error("Window size is not positive: -w %s\n", optarg);
                break;
            case 's': sample_names = optarg; break;
            case 'S': sample_names = optarg; sample_is_file = 1; break;
            case  1 : force_samples = 1; break;
            case 'G': sites_only = 1; break;
            case 'h':
            case '?':
            default: error("%s", usage()); break;
        }
    }

    // this ugly workaround is required to make sure we can set samples on both headers even when sample_is_file is true and sample_names is stdin
    if ( sample_names )
    {
        int nsmpl;
        char **smpl = hts_readlist(sample_names[0]=='^' ? sample_names + 1 : sample_names, sample_is_file, &nsmpl);
        kstring_t tmp = {0, 0, NULL};
        if (sample_names[0]=='^') ksprintf(&tmp, "^%s", smpl[0]);
        else ksprintf(&tmp, "%s", smpl[0]);
        for (int i=1; i<nsmpl; i++) ksprintf(&tmp, ",%s", smpl[i]);
        int ret = bcf_hdr_set_samples(args->in_hdr, tmp.s, 0);
        if ( ret < 0 ) error("Error parsing the sample list\n");
        else if ( ret > 0 )
        {
            if ( force_samples )
                fprintf(stderr, "Warn: subset called for sample that does not exist in header: \"%s\"... skipping\n", smpl[ret-1]);
            else
                error("Error: subset called for sample that does not exist in header: \"%s\". Use \"--force-samples\" to ignore this error.\n", smpl[ret-1]);
        }
        if ( bcf_hdr_nsamples(args->in_hdr) == 0 )
            error("Error: subsetting has removed all samples\n");
        if ( bcf_hdr_set_samples(args->out_hdr, tmp.s, 0) < 0 ) error("Error parsing the sample list\n");
        free(tmp.s);
        for (int i=0; i<nsmpl; i++) free(smpl[i]);
        free(smpl);
    }

    if ( sex_fname )
    {
        args->sex = (int *)calloc((size_t)bcf_hdr_nsamples(args->in_hdr), sizeof(int));
        parse_sex(args->in_hdr, sex_fname, args->sex);
    }

    if ( ref_fname )
    {
        args->fai = fai_load(ref_fname);
        if ( !args->fai ) error("Failed to load the fai index: %s\n", ref_fname);
        bcf_hdr_append(args->out_hdr, "##INFO=<ID=GC,Number=1,Type=Float,Description=\"GC ratio content around the variant\">");
        bcf_hdr_append(args->out_hdr, "##INFO=<ID=CpG,Number=1,Type=Float,Description=\"CpG ratio content around the variant\">");
    }

    args->nsmpl = bcf_hdr_nsamples(args->in_hdr);
    if ( args->nsmpl == 0 ) error("Error: input VCF file has no samples\n");
    args->gt_id = bcf_hdr_id2int(args->in_hdr, BCF_DT_ID, "GT");
    args->ad_id = bcf_hdr_id2int(args->in_hdr, BCF_DT_ID, "AD");
    args->baf_id = bcf_hdr_id2int(args->in_hdr, BCF_DT_ID, "BAF");
    args->fmt_id = args->format ? bcf_hdr_id2int(args->in_hdr, BCF_DT_ID, args->format) : -1;

    if ( args->format && args->fmt_id < 0 ) error("Error: %s format field is not present, cannot perform --balance analysis\n", args->format);
    if ( args->ad && ( args->gt_id < 0 || args->ad_id < 0 ) ) error("Error: Either GT or AD format fields are not present, cannot perform --ad-het analysis\n");
    if ( args->phase && ( args->gt_id < 0 || (args->ad_id < 0 && args->baf_id < 0 && args->fmt_id < 0 ) ) ) error("Error: Either GT or AD/BAF/%s format fields are not present, cannot perform --phase analysis\n", args->format);

    if ( args->format )
    {
        bcf_hdr_append(args->out_hdr, "##INFO=<ID=Bal,Number=2,Type=Integer,Description=\"Reference alternate allelic shift counts\">");
        bcf_hdr_append(args->out_hdr, "##INFO=<ID=Bal_Test,Number=1,Type=Float,Description=\"Reference alternate allelic shift binomial test -log10(P)\">");
        if ( args->phase )
        {
            bcf_hdr_append(args->out_hdr, "##INFO=<ID=Bal_Phase,Number=2,Type=Integer,Description=\"Paternal maternal allelic shift counts\">");
            bcf_hdr_append(args->out_hdr, "##INFO=<ID=Bal_Phase_Test,Number=1,Type=Float,Description=\"Paternal maternal allelic shift binomial test -log10(P)\">");
        }
    }

    bcf_hdr_append(args->out_hdr, "##INFO=<ID=AC_Het,Number=1,Type=Integer,Description=\"Number of heterozygous genotypes\">");
    if ( args->sex )
    {
        bcf_hdr_append(args->out_hdr, "##INFO=<ID=AC_Het_Sex,Number=2,Type=Integer,Description=\"Number of heterozygous genotypes by sex\">");
        bcf_hdr_append(args->out_hdr, "##INFO=<ID=AC_Sex_Test,Number=1,Type=Float,Description=\"Fisher's exact test for alternate alleles and sex\">");
    }

    if ( args->ad && args->ad_id >= 0 )
    {
        bcf_hdr_append(args->out_hdr, "##INFO=<ID=AD_Het,Number=2,Type=Integer,Description=\"Allelic depths for the reference and alternate alleles across heterozygous genotypes\">");
        bcf_hdr_append(args->out_hdr, "##INFO=<ID=AD_Het_Test,Number=1,Type=Float,Description=\"Binomial test for reference and alternate allelic depth across heterozygous genotypes -log10(P)\">");
    }
    if ( args->phase )
    {
        bcf_hdr_append(args->out_hdr, "##INFO=<ID=AC_Het_Phase,Number=2,Type=Integer,Description=\"Number of heterozygous genotypes by transmission type\">");
        bcf_hdr_append(args->out_hdr, "##INFO=<ID=AC_Het_Phase_Test,Number=1,Type=Float,Description=\"Binomial test for allelic transmission bias across heterozygous genotypes -log10(P)\">");
        if ( args->ad_id >= 0 || args->baf_id >= 0 )
            bcf_hdr_append(args->out_hdr, "##INFO=<ID=BAF_Phase_Test,Number=4,Type=Float,Description=\"Welch's t-test and Mann-Whitney U test for allelic transmission ratios across heterozygous genotypes\">");
    }

    if ( sites_only ) if ( bcf_hdr_set_samples(args->out_hdr, NULL, 0)  < 0 ) error("Error parsing the sample list\n");

    args->gt_phase_arr = (int8_t *)malloc(args->nsmpl * sizeof(int8_t));
    args->fmt_sign_arr = (int8_t *)malloc(args->nsmpl * sizeof(int8_t));
    args->gt0_arr = (int16_t *)malloc(args->nsmpl * sizeof(int16_t));
    args->gt1_arr = (int16_t *)malloc(args->nsmpl * sizeof(int16_t));
    args->ad0_arr = (int16_t *)malloc(args->nsmpl * sizeof(int16_t));
    args->ad1_arr = (int16_t *)malloc(args->nsmpl * sizeof(int16_t));
    args->baf_arr[0] = (float *)malloc(args->nsmpl * sizeof(float));
    args->baf_arr[1] = (float *)malloc(args->nsmpl * sizeof(float));

    return 0;
}

// Petr Danecek's implementation in bcftools/mcall.c
double binom_dist(int N, double p, int k);

// returns 2*pbinom(k,n,1/2) if k<n/2 by precomputing values in a table
static double binom_exact(int k, int n)
{
    static double *dbinom = NULL, *pbinom = NULL;
    static size_t n_size = 0, m_dbinom = 0, m_pbinom = 0;

    if ( n < 0 && k < 0 )
    {
        free(dbinom);
        free(pbinom);
        return NAN;
    }

    if ( n < 0 || k < 0 || k > n ) return NAN;

    if ( n > 1000 ) return binom_dist(n, 0.5, k);

    if ( k == n>>1 ) return 1.0;

    if ( k<<1 > n ) k = n - k;

    if ( n >= n_size )
    {
        size_t len = (size_t)(1+(1+(n>>1))*((n+1)>>1));
        hts_expand(double, len, m_dbinom, dbinom);
        hts_expand(double, len, m_pbinom, pbinom);
        dbinom[0] = 1.0;
        for (int i=n_size?(int)n_size:1; i<n+1; i++)
        {
            int prev_idx = i-1?1+((i-1)>>1)*(i>>1):0;
            int curr_idx = 1+(i>>1)*((i+1)>>1);
            dbinom[curr_idx] = dbinom[prev_idx] * 0.5;
            pbinom[curr_idx] = dbinom[curr_idx];
            for (int j=1; j<(i+1)>>1; j++)
            {
                curr_idx++;
                dbinom[curr_idx] = (double)i / (double)j * dbinom[prev_idx] * 0.5;
                pbinom[curr_idx] = pbinom[curr_idx-1] + dbinom[curr_idx];
                prev_idx++;
            }
        }
        n_size = (size_t)(n + 1);
    }

    int idx = 1+(n>>1)*((n+1)>>1)+k;
    return 2.0 * pbinom[idx];
}

// Giulio Genovese's implementation in bcftools/vcfmocha.c
static int sample_mean_var(const float *x,
                           int n,
                           double *mu,
                           double *s2)
{
    if (n<2) return -1;
    *mu = 0;
    *s2 = 0;
    int j = 0;
    for (int i=0; i<n; i++)
    {
        if ( !isnan(x[i]) )
        {
            *mu += (double)x[i];
            *s2 += sq((double)x[i]);
            j++;
        }
    }
    if ( j <= 1 ) return -1;
    *mu /= (double)j;
    *s2 -= sq(*mu)*(double)j;
    *s2 /= (double)(j-1);
    return 0;
}

static double welch_t_test(float *a,
                           float *b,
                           int na,
                           int nb)
{
    double mua, mub, sa2, sb2, t, v;
    if (na < 2 || nb < 2) return HUGE_VAL;
    sample_mean_var(a, na, &mua, &sa2);
    sample_mean_var(b, nb, &mub, &sb2);
    t = ( mua - mub ) / sqrt( sa2 / na + sb2 / nb );
    v = ( sa2 / na + sb2 / nb );
    v *= v;
    v /= sq(sa2) / na / na / (na - 1) + sq(sb2) / nb / nb / (nb - 1);
    return kf_betai(v/2.0f, 0.5, v/(v+sq(t)));
}

// Petr Danecek's and James Bonfield's implementation in bcftools/bam2bcf.c
double mann_whitney_1947_cdf(int n, int m, int U);

static int cmpfunc(const void * a, const void * b) {
   return ( *(float*)a > *(float*)b ) - ( *(float*)a < *(float*)b );
}

// it currently does not handle nans
// adapted from Petr Danecek's implementation of calc_mwu_bias_cdf() in bcftools/bam2bcf.c
static double mann_whitney_u(float *a,
                             float *b,
                             int na,
                             int nb)
{
    qsort (a, (size_t)na, sizeof(float), cmpfunc);
    qsort (b, (size_t)nb, sizeof(float), cmpfunc);

    int i = 0, j = 0, ca, cb;
    double U = 0, ties = 0;
    while ( i<na || j<nb )
    {
        double curr = (j==nb || (i<na && a[i]<b[j])) ? a[i] : b[j];
        for (ca=0; i<na && a[i]==curr; i++) ca++;
        for (cb=0; j<nb && b[j]==curr; j++) cb++;
        U += ca * (j - cb * 0.5);
        if ( ca && cb )
        {
            double tie = ca + cb;
            ties += (sq(tie) - 1.0) * tie;
        }
    }
    if ( !na || !nb ) return HUGE_VAL;

    double U_min = ((double)na * nb) - U;
    if ( U < U_min ) U_min = U;

    if ( na==1 ) return 2.0 * (floor(U_min) + 1.0) / (double)(nb + 1);
    if ( nb==1 ) return 2.0 * (floor(U_min) + 1.0) / (double)(na + 1);

    // Normal approximation, very good for na>=8 && nb>=8 and reasonable if na<8 or nb<8
    if ( na>=8 || nb>=8 )
    {
        double mean = ((double)na*nb)*0.5;
        // Correction for ties:
        double N = na + nb;
        double var2 = (sq(N) - 1) * N - ties;
        if ( var2==0 ) return 1.0;
        var2 *= ((double)na*nb) / N / (N - 1) / 12.0;
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

// Giulio Genovese's implementation in bcftools/vcfmocha.c
float get_median(const float *v, int n, const int *imap);
int bcf_get_genotype_phase(bcf_fmt_t *fmt, int8_t *gt_phase_arr, int nsmpl);
int bcf_get_genotype_alleles(const bcf_fmt_t *fmt, int16_t *gt0_arr, int16_t *gt1_arr, int nsmpl);
int bcf_get_allelic_depth(const bcf_fmt_t *fmt, const int16_t *gt0_arr, const int16_t *gt1_arr, int16_t *ad0_arr, int16_t *ad1_arr, int nsmpl);

// retrieve phase information from BCF record
// assumes little endian architecture
static int bcf_get_format_sign(bcf_fmt_t *fmt,
                               int8_t *fmt_sign_arr,
                               int nsmpl)
{
    if ( !fmt || fmt->n != 1 ) return 0;

    #define BRANCH(type_t, bcf_type_vector_end, bcf_type_missing) { \
        type_t *p = (type_t *)fmt->p; \
        for (int i=0; i<nsmpl; i++) \
        { \
            if ( p[i] == bcf_type_vector_end || p[i] == bcf_type_missing ) fmt_sign_arr[i] = bcf_int8_missing; \
            else if ( p[i] == (type_t)0 ) fmt_sign_arr[i] = (int8_t)0; \
            else if ( p[i] > (type_t)0 ) fmt_sign_arr[i] = (int8_t)1; \
            else fmt_sign_arr[i] = (int8_t)-1; \
        } \
    }
    switch (fmt->type) {
        case BCF_BT_INT8:  BRANCH(int8_t, bcf_int8_vector_end, bcf_int8_missing); break;
        case BCF_BT_INT16: BRANCH(int16_t, bcf_int16_vector_end, bcf_int16_missing); break;
        case BCF_BT_INT32: BRANCH(int32_t, bcf_int32_vector_end, bcf_int32_missing); break;
        case BCF_BT_FLOAT: BRANCH(int32_t, bcf_float_vector_end, bcf_float_missing); break;
        default: error("Unexpected type %d", fmt->type);
    }
    #undef BRANCH

    return 1;
}

bcf1_t *process(bcf1_t *rec)
{
    // compute GC and CpG content for each site
    if (args->fai)
    {
        int fa_len;
        int at_cnt = 0, cg_cnt = 0, cpg_cnt = 0;
        const char *ref = rec->d.allele[0];
        char *fa = faidx_fetch_seq(args->fai, bcf_seqname(args->in_hdr,rec), rec->pos-args->ws, rec->pos+(int)strlen(ref)-1+args->ws, &fa_len);
        if ( !fa ) error("fai_fetch_seq failed at %s:%d\n", bcf_hdr_id2name(args->in_hdr,rec->rid), rec->pos+1);
        for (int i=0; i<fa_len; i++)
        {
            if ( fa[i] > 96 ) fa[i] = (char)(fa[i] - 32);
            if (fa[i]=='A' || fa[i]=='T') at_cnt++;
            if (fa[i]=='C' || fa[i]=='G') cg_cnt++;
            if (i>0) if (fa[i-1]=='C' && fa[i]=='G') cpg_cnt+=2;
        }
        free(fa);
        float ratio = (float)(cg_cnt) / (float)(at_cnt + cg_cnt);
        bcf_update_info_float(args->out_hdr, rec, "GC", &ratio, 1);
        ratio = (float)cpg_cnt / (float)(fa_len);
        bcf_update_info_float(args->out_hdr, rec, "CpG", &ratio, 1);
    }

    // extract format information from VCF format records
    bcf_fmt_t *gt_fmt = bcf_get_fmt_id(rec, args->gt_id);
    int gt_phase = bcf_get_genotype_phase(gt_fmt, args->gt_phase_arr, args->nsmpl);
    if ( !bcf_get_genotype_alleles(gt_fmt, args->gt0_arr, args->gt1_arr, args->nsmpl) ) goto ret;
    bcf_fmt_t *fmt = bcf_get_fmt_id(rec, args->fmt_id);
    int fmt_sign = ( args->format && fmt ) ? bcf_get_format_sign(fmt, args->fmt_sign_arr, args->nsmpl) : 0;
    bcf_fmt_t *ad_fmt = bcf_get_fmt_id(rec, args->ad_id);
    int ad = ad_fmt ? bcf_get_allelic_depth(ad_fmt, args->gt0_arr, args->gt1_arr, args->ad0_arr, args->ad1_arr, args->nsmpl) : 0;
    bcf_fmt_t *baf_fmt = bcf_get_fmt_id(rec, args->baf_id);
    int baf = baf_fmt && baf_fmt->n == 1 && baf_fmt->type == BCF_BT_FLOAT ? 1 : 0;

    float ret[4];
    int ac_het = 0, ac_sex[] = {0, 0, 0, 0}, ac_het_sex[] = {0, 0}, ac_het_phase[] = {0, 0}, fmt_bal[] = {0,0}, fmt_bal_phase[] = {0, 0}, ad_het[] = {0, 0};

    for (int i=0; i<args->nsmpl; i++)
    {
        float curr_baf = NAN;

        // if genotype is missing, skip
        if ( args->gt0_arr[i] == bcf_int16_missing || args->gt0_arr[i] == bcf_int16_missing ) continue;

        int idx_fmt_sign = ( fmt_sign && args->fmt_sign_arr[i] != bcf_int8_missing && args->fmt_sign_arr[i] != 0 ) ? (1-args->fmt_sign_arr[i])/2 : -1;
        if ( idx_fmt_sign >= 0 ) fmt_bal[idx_fmt_sign]++;

        if ( args->sex && ( args->sex[i] == 1 || args->sex[i] == 2 ) )
        {
            if ( args->gt0_arr[i] == 0 && args->gt1_arr[i] == 0 ) ac_sex[args->sex[i]-1]++;
            else if ( args->gt0_arr[i] > 0 && args->gt1_arr[i] > 0 ) ac_sex[2+args->sex[i]-1]++;
        }

        // if genotype is not heterozygous, skip
        if ( args->gt0_arr[i] == args->gt1_arr[i] || ( args->gt0_arr[i] != 0 && args->gt1_arr[i] != 0 ) ) continue;

        int idx_gt_phase = gt_phase && (args->gt_phase_arr[i] == -1 || args->gt_phase_arr[i] == 1) ? (1-args->gt_phase_arr[i])/2 : -1;
        ac_het++;
        if ( args->sex && ( args->sex[i] == 1 || args->sex[i] == 2 ) ) ac_het_sex[args->sex[i]-1]++;
        if ( idx_gt_phase >= 0 ) ac_het_phase[idx_gt_phase]++;

        int idx_fmt_phase = ( idx_gt_phase >= 0 && idx_fmt_sign >= 0 ) ? (1-args->fmt_sign_arr[i]*args->gt_phase_arr[i])/2 : -1;
        if ( idx_fmt_phase >= 0 ) fmt_bal_phase[idx_fmt_phase]++;

        if ( ad )
        {
            int ref_cnt = args->ad0_arr[i];
            int alt_cnt = args->ad1_arr[i];
            ad_het[0] += ref_cnt;
            ad_het[1] += alt_cnt;
            curr_baf = ((float)alt_cnt + 0.5f) / ((float)ref_cnt + (float)alt_cnt + 1.0f);
        }
        if ( baf ) curr_baf = ((float *)baf_fmt->p)[i];
        if ( idx_gt_phase >= 0 && !isnan(curr_baf) )
        {
            args->baf_arr[idx_gt_phase][ac_het_phase[idx_gt_phase] - 1] = curr_baf;
        }
    }

    bcf_update_info_int32(args->out_hdr, rec, "AC_Het", &ac_het, 1);
    if ( args->sex )
    {
        bcf_update_info_int32(args->out_hdr, rec, "AC_Het_Sex", &ac_het_sex, 2);
        double left, right, fisher;
        ret[0] = 0.0f - (float)log10(kt_fisher_exact(ac_sex[0], ac_sex[1], ac_sex[2], ac_sex[3], &left, &right, &fisher));
        bcf_update_info_float(args->out_hdr, rec, "AC_Sex_Test", &ret, 1);
    }
    if ( args->phase )
    {
        bcf_update_info_int32(args->out_hdr, rec, "AC_Het_Phase", &ac_het_phase, 2);
        ret[0] = 0.0f - (float)log10(binom_exact(ac_het_phase[0], ac_het_phase[0] + ac_het_phase[1]));
        bcf_update_info_float(args->out_hdr, rec, "AC_Het_Phase_Test", &ret, 1);
    }
    if ( args->format )
    {
        bcf_update_info_int32(args->out_hdr, rec, "Bal", &fmt_bal, 2);
        ret[0] = 0.0f - (float)log10(binom_exact(fmt_bal[0], fmt_bal[0] + fmt_bal[1]));
        bcf_update_info_float(args->out_hdr, rec, "Bal_Test", &ret, 1);
        if ( args->phase )
        {
            bcf_update_info_int32(args->out_hdr, rec, "Bal_Phase", &fmt_bal_phase, 2);
            ret[0] = 0.0f - (float)log10(binom_exact(fmt_bal_phase[0], fmt_bal_phase[0] + fmt_bal_phase[1]));
            bcf_update_info_float(args->out_hdr, rec, "Bal_Phase_Test", &ret, 1);
        }
    }
    if ( args->ad )
    {
        bcf_update_info_int32(args->out_hdr, rec, "AD_Het", &ad_het, 2);
        ret[0] = 0.0f - (float)log10(binom_exact(ad_het[0], ad_het[0] + ad_het[1]));
        bcf_update_info_float(args->out_hdr, rec, "AD_Het_Test", &ret, 1);
    }
    if ( args->phase && ac_het_phase[0] && ac_het_phase[1] )
    {
        ret[0] = get_median( args->baf_arr[0], ac_het_phase[0], NULL );
        ret[1] = get_median( args->baf_arr[1], ac_het_phase[1], NULL );
        ret[2] = 0.0f - (float)log10(welch_t_test( args->baf_arr[0], args->baf_arr[1], ac_het_phase[0], ac_het_phase[1] ));
        ret[3] = 0.0f - (float)log10(mann_whitney_u( args->baf_arr[0], args->baf_arr[1], ac_het_phase[0], ac_het_phase[1] ));
        bcf_update_info_float(args->out_hdr, rec, "BAF_Phase_Test", &ret, 4);
    }

ret:
    // remove all samples if sites_only was selected
    if ( bcf_hdr_nsamples(args->out_hdr) == 0 ) bcf_subset(args->out_hdr, rec, 0, NULL);
    return rec;
}

void destroy(void)
{
    binom_exact(-1, -1);
    free(args->sex);
    free(args->gt_phase_arr);
    free(args->fmt_sign_arr);
    free(args->gt0_arr);
    free(args->gt1_arr);
    free(args->ad0_arr);
    free(args->ad1_arr);
    free(args->baf_arr[0]);
    free(args->baf_arr[1]);
    free(args);
}
