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

// kv macro from https://github.com/attractivechaos/klib/blob/master/kvec.h by Attractive Chaos <attractor@live.co.uk>
#define kvec_t(type) struct { size_t n, m; type *a; }
#define kv_push(type, v, x) do { \
    if ((v).n == (v).m) { \
        (v).m = (v).m? (v).m<<1 : 2; \
        (v).a = (type*)realloc((v).a, sizeof(type) * (v).m); \
    } \
    (v).a[(v).n++] = (x); \
} while (0)

#define SIGN(x) (((x) > 0) - ((x) < 0))

static inline float sqf(float x) { return x*x; }

typedef struct
{
    int ase, ad, tx, ws;
    int ngt_arr, nad_arr, nbdev_phase_arr;
    int *sex;
    int32_t *gt_arr, *ad_arr, *bdev_phase_arr;
    faidx_t *fai;
    bcf_hdr_t *in_hdr, *out_hdr;
}
args_t;

args_t *args;

const char *about(void)
{
    return "MOsaic CHromosomal Alterations tools\n";
}

const char *usage(void)
{
    return
"\n"
"About:   tools for the MOsaic CHromosomal Alterations pipeline.\n"
"\n"
"Usage:   bcftools +mochatools [General Options] -- [Plugin Options]\n"
"\n"
"General options:\n"
"   run \"bcftools plugin\" for a list of common options\n"
"\n"
"Plugin options:\n"
"   -b, --binom-ase               performs binomial test for asymmetry of B Allele Frequency (Bdev_Phase)\n"
"   -a, --ad-het                  performs binomial test for reference / alternate allelic depth (AD)\n"
"   -x  --sex <file>              file including information about sex of sample\n"
"   -t, --tx-bias                 perform transmission bias tests (requires absolute phasing)\n"
"   -f, --fasta-ref <file>        reference sequence to compute GC and CpG content\n"
"   -w, --window-size <int>       Window size in bp used to compute the GC and CpG content [200]\n"
"   -s, --samples [^]<list>       comma separated list of samples to include (or exclude with \"^\" prefix)\n"
"   -S, --samples-file [^]<file>  file of samples to include (or exclude with \"^\" prefix)\n"
"       --force-samples           only warn about unknown subset samples\n"
"   -G, --drop-genotypes          drop individual genotype information (after running statistical tests)\n"
"\n";
}

static void parse_sex(args_t *args, char *fname)
{
    htsFile *fp = hts_open(fname, "r");
    if ( !fp ) error("Could not read: %s\n", fname);

    kstring_t str = {0,0,0};
    if ( hts_getline(fp, KS_SEP_LINE, &str) <= 0 ) error("Empty file: %s\n", fname);

    args->sex = (int *)calloc(bcf_hdr_nsamples(args->in_hdr), sizeof(int));

    int moff = 0, *off = NULL;
    char *tmp = NULL;
    do
    {
        int ncols = ksplit_core(str.s, 0, &moff, &off);
        if ( ncols<2 ) error("Could not parse the sex file: %s\n", str.s);

        int sample = bcf_hdr_id2int(args->in_hdr, BCF_DT_SAMPLE, &str.s[off[0]]);
        if ( sample>=0 ) args->sex[sample] = (int)strtol(&str.s[off[1]], &tmp, 0);
    } while ( hts_getline(fp, KS_SEP_LINE, &str)>=0 );

    free(str.s);
    free(off);
    hts_close(fp);
}

int init(int argc, char **argv, bcf_hdr_t *in, bcf_hdr_t *out)
{
    args = (args_t *)calloc(1, sizeof(args_t));
    args->ws = 200;
    args->in_hdr = in;
    args->out_hdr = out;
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
        {"binom-ase", no_argument, NULL, 'b'},
        {"ad-het", no_argument, NULL, 'a'},
        {"sex", required_argument, NULL, 'x'},
        {"tx-bias", no_argument, NULL, 't'},
        {"fasta-ref", required_argument, NULL, 'f'},
        {"window-size", required_argument, NULL, 'w'},
        {"samples", required_argument, NULL, 's'},
        {"samples-file", required_argument, NULL, 'S'},
        {"force-samples", no_argument, NULL, 1},
        {"drop-genotypes", no_argument,NULL, 'G'},
        {NULL,0,NULL,0}
    };

    while ((c = getopt_long(argc, argv, "h?bax:tf:w:s:S:G",loptions,NULL)) >= 0)
    {
        switch (c)
        {
            case 'b': args->ase = 1; break;
            case 'a': args->ad = 1; break;
            case 'x': sex_fname = optarg; break;
            case 't': args->tx = 1; break;
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

    // this ugly workaround is required to make sure we can set samples on both headers even when sample_is_file is true
    if ( sample_names )
    {
        int nsmpl;
        char **smpl = hts_readlist(sample_names[0]=='^' ? sample_names + 1 : sample_names, sample_is_file, &nsmpl);
        kstring_t tmp = {0, 0, NULL};
        if (sample_names[0]=='^') ksprintf(&tmp, "^%s", smpl[0]);
        else ksprintf(&tmp, "%s", smpl[0]);
        for (int i=1; i<nsmpl; i++) ksprintf(&tmp, ",%s", smpl[i]);
        int ret = bcf_hdr_set_samples(args->in_hdr, tmp.s, 0);
        if ( ret<0 ) error("Error parsing the sample list\n");
        else if ( ret>0 )
        {
            if ( force_samples )
                fprintf(stderr, "Warn: subset called for sample that does not exist in header: \"%s\"... skipping\n", smpl[ret-1]);
            else
                error("Error: subset called for sample that does not exist in header: \"%s\". Use \"--force-samples\" to ignore this error.\n", smpl[ret-1]);
        }
        if ( bcf_hdr_nsamples(args->in_hdr) == 0 )
            error("Error: subsetting has removed all samples\n");
        bcf_hdr_set_samples(args->out_hdr, tmp.s, 0);
        free(tmp.s);
        for (int i=0; i<nsmpl; i++) free(smpl[i]);
        free(smpl);
    }

    if ( sex_fname ) parse_sex(args, sex_fname);

    if ( ref_fname )
    {
        args->fai = fai_load(ref_fname);
        if ( !args->fai ) error("Failed to load the fai index: %s\n", ref_fname);
        bcf_hdr_append(args->out_hdr, "##INFO=<ID=GC,Number=1,Type=Float,Description=\"GC ratio content around the variant\">");
        bcf_hdr_append(args->out_hdr, "##INFO=<ID=CpG,Number=1,Type=Float,Description=\"CpG ratio content around the variant\">");
    }

    // check for what fields are present
    int gt = bcf_hdr_id2int(args->in_hdr,BCF_DT_ID,"GT");
    int ad = bcf_hdr_id2int(args->in_hdr,BCF_DT_ID,"AD");
    int baf = bcf_hdr_id2int(args->in_hdr,BCF_DT_ID,"BAF");
    int bdev_phase = bcf_hdr_id2int(args->in_hdr,BCF_DT_ID,"Bdev_Phase");

    if ( args->ase && bdev_phase < 0 ) error("Error: Bdev_Phase format field is not present, cannot perform --binom-ase analysis\n");
    if ( args->ad && ( gt < 0 || ad < 0 ) ) error("Error: Either GT or AD format fields are not present, cannot perform --ad-het analysis\n");
    if ( args->tx && ( gt < 0 || (ad < 0 && baf < 0 && bdev_phase < 0 ) ) ) error("Error: Either GT or AD/BAF/Bdev_Phase format fields are not present, cannot perform --tx-bias analysis\n");

    if ( args->ase )
    {
        bcf_hdr_append(args->out_hdr, "##INFO=<ID=ASE,Number=2,Type=Integer,Description=\"Reference alternate allelic shift counts\">");
        bcf_hdr_append(args->out_hdr, "##INFO=<ID=ASE_Test,Number=1,Type=Float,Description=\"Reference alternate allelic shift binomial test -log10(P)\">");
        if ( args->tx )
        {
            bcf_hdr_append(args->out_hdr, "##INFO=<ID=ASE_Tx,Number=2,Type=Integer,Description=\"Paternal maternal allelic shift counts\">");
            bcf_hdr_append(args->out_hdr, "##INFO=<ID=ASE_Tx_Test,Number=1,Type=Float,Description=\"Paternal maternal allelic shift binomial test -log10(P)\">");
        }
    }

    bcf_hdr_append(args->out_hdr, "##INFO=<ID=AC_Het,Number=1,Type=Integer,Description=\"Number of heterozygous genotypes\">");
    if ( args->sex )
    {
        bcf_hdr_append(args->out_hdr, "##INFO=<ID=AC_Het_Sex,Number=2,Type=Integer,Description=\"Number of heterozygous genotypes by sex\">");
        bcf_hdr_append(args->out_hdr, "##INFO=<ID=AC_Sex_Test,Number=1,Type=Float,Description=\"Fisher's exact test for alternate alleles and sex\">");
    }

    if ( args->ad && ad >= 0 )
    {
        bcf_hdr_append(args->out_hdr, "##INFO=<ID=AD_Het,Number=2,Type=Integer,Description=\"Allelic depths for the reference and alternate alleles across heterozygous genotypes\">");
        bcf_hdr_append(args->out_hdr, "##INFO=<ID=AD_Het_Test,Number=1,Type=Float,Description=\"Binomial test for reference and alternate allelic depth across heterozygous genotypes -log10(P)\">");
    }
    if ( args->tx )
    {
        bcf_hdr_append(args->out_hdr, "##INFO=<ID=AC_Het_Tx,Number=2,Type=Integer,Description=\"Number of heterozygous genotypes by transmission type\">");
        bcf_hdr_append(args->out_hdr, "##INFO=<ID=AC_Het_Tx_Test,Number=1,Type=Float,Description=\"Binomial test for allelic transmission bias across heterozygous genotypes -log10(P)\">");
        if ( ad >= 0 || baf >= 0 )
            bcf_hdr_append(args->out_hdr, "##INFO=<ID=BAF_Tx_Test,Number=4,Type=Float,Description=\"Welch's t-test and Mann-Whitney U test for allelic transmission ratios across heterozygous genotypes\">");
    }

    if ( sites_only ) bcf_hdr_set_samples(args->out_hdr, NULL, 0);

    return 0;
}

// Petr Danecek's implementation in bcftools/mcall.c
double binom_dist(int N, double p, int k);

// Giulio Genovese's implementation in bcftools/vcfmocha.c
static int sample_mean_var(const float *x, int n, float *mu, float *s2)
{
    if (n<2) return -1;
    *mu = 0;
    *s2 = 0;
    int j = 0;
    for (int i=0; i<n; i++)
    {
        if ( !isnan(x[i]) )
        {
            *mu += x[i];
            *s2 += sqf(x[i]);
            j++;
        }
    }
    if ( j <= 1 ) return -1;
    *mu /= (float)j;
    *s2 -= sqf(*mu)*(float)j;
    *s2 /= (float)(j-1);
    return 0;
}

static float welch_t_test(float *a, float *b, int na, int nb)
{
    float mua, mub, sa2, sb2, t, v;
    if (na < 2 || nb < 2) return HUGE_VAL;
    sample_mean_var(a, na, &mua, &sa2);
    sample_mean_var(b, nb, &mub, &sb2);
    t = ( mua - mub ) / sqrtf( sa2 / na + sb2 / nb );
    v = ( sa2 / na + sb2 / nb );
    v *= v;
    v /= sqf(sa2) / na / na / (na - 1) + sqf(sb2) / nb / nb / (nb - 1);
    return kf_betai(v/2.0f, 0.5, v/(v+sqf(t)));
}

// Petr Danecek's and James Bonfield's implementation in bcftools/bam2bcf.c
double mann_whitney_1947_cdf(int n, int m, int U);

static int cmpfunc(const void * a, const void * b) {
   return ( *(float*)a > *(float*)b ) - ( *(float*)a < *(float*)b );
}

// it currently does not handle nans
// adapted from Petr Danecek's implementation in calc_mwu_bias_cdf() in bcftools/bam2bcf.c
static float mann_whitney_u(float *a, float *b, int na, int nb)
{
    qsort (a, na, sizeof(float), cmpfunc);
    qsort (b, nb, sizeof(float), cmpfunc);

    int i = 0, j = 0, ca, cb;
    float U = 0, ties = 0;
    while ( i<na || j<nb )
    {
        float curr = (j==nb || (i<na && a[i]<b[j])) ? a[i] : b[j];
        for (ca=0; i<na && a[i]==curr; i++) ca++;
        for (cb=0; j<nb && b[j]==curr; j++) cb++;
        U += ca * (j - cb*0.5);
        if ( ca && cb )
        {
            float tie = ca + cb;
            ties += (sqf(tie)-1)*tie;
        }
    }
    if ( !na || !nb ) return HUGE_VAL;

    float U_min = ((float)na * nb) - U;
    if ( U < U_min ) U_min = U;

    if ( na==1 ) return 2.0f * (floorf(U_min)+1.0f) / (float)(nb+1);
    if ( nb==1 ) return 2.0f * (floorf(U_min)+1.0f) / (float)(na+1);

    // Normal approximation, very good for na>=8 && nb>=8 and reasonable if na<8 or nb<8
    if ( na>=8 || nb>=8 )
    {
        float mean = ((float)na*nb)*0.5f;
        // Correction for ties:
        float N = na+nb;
        float var2 = (sqf(N)-1)*N-ties;
        if ( var2==0 ) return 1.0f;
        var2 *= ((float)na*nb)/N/(N-1)/12.0f;
        // No correction for ties:
        // float var2 = ((float)na*nb)*(na+nb+1)/12.0f;
        float z = (U_min - mean)/sqrtf(2.0f*var2); // z is N(0,1)
        // return 2.0 - kf_erfc(z);  // which is 1 + erf(z)
        return (float)kf_erfc(-z); // which is 1 - erf(-z)
    }

    // Exact calculation
    float pval = 2.0f * (float)mann_whitney_1947_cdf(na,nb,U_min);
    return pval > 1.0f ? 1.0f : pval;
}

// Giulio Genovese's implementation in bcftools/vcfmocha.c
float get_median(const float *v, int n, const int *imap);

bcf1_t *process(bcf1_t *rec)
{
    // compute GC and CpG content for each site
    if (args->fai)
    {
        int fa_len;
        int at_cnt = 0, cg_cnt = 0, cpg_cnt = 0;
        const char *ref = rec->d.allele[0];
        char *fa = faidx_fetch_seq(args->fai, bcf_seqname(args->in_hdr,rec), rec->pos-args->ws, rec->pos+strlen(ref)-1+args->ws, &fa_len);
        if ( !fa ) error("fai_fetch_seq failed at %s:%d\n", bcf_hdr_id2name(args->in_hdr,rec->rid), rec->pos+1);
        for (int i=0; i<fa_len; i++)
        {
            if ( (int)fa[i]>96 ) fa[i] -= 32;
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

    // extract genotypes
    int nsmpl = bcf_hdr_nsamples(args->in_hdr);
    if ( nsmpl == 0 ) goto end;
    int ngt = bcf_get_genotypes(args->in_hdr, rec, &args->gt_arr, &args->ngt_arr);
    if ( ngt <= 0 ) goto end;
    int max_ploidy = ngt / nsmpl;
    if( rec->n_allele <= 1 || max_ploidy != 2 ) goto end;

    float ret[4];
    kvec_t(float) kv_baf[2] = { {0, 0, NULL}, {0, 0, NULL} };
    int ac_het = 0, ac_sex[] = {0, 0, 0, 0}, ac_het_sex[] = {0, 0}, ac_het_tx[] = {0, 0}, ase[] = {0,0}, ase_tx[] = {0, 0}, ad_het[] = {0, 0};

    int nbdev_phase = args->ase ? bcf_get_format_int32(args->in_hdr, rec, "Bdev_Phase", &args->bdev_phase_arr, &args->nbdev_phase_arr) : 0;
    int nad = args->ad ? bcf_get_format_int32(args->in_hdr, rec, "AD", &args->ad_arr, &args->nad_arr) : 0;

    for (int i=0; i<nsmpl; i++)
    {
        float baf = NAN;
        int32_t *ptr = args->gt_arr + i * max_ploidy;
        // if genotype is missing, skip
        if ( bcf_gt_is_missing(ptr[0]) || bcf_gt_is_missing(ptr[1]) || ptr[1]==bcf_int32_vector_end ) continue;
        if ( args->sex && ( args->sex[i] == 1 || args->sex[i] == 2 ) )
        {
            if ( bcf_gt_allele(ptr[0]) == 0 && bcf_gt_allele(ptr[1]) == 0 ) ac_sex[args->sex[i]-1]++;
            else if ( bcf_gt_allele(ptr[0]) > 0 && bcf_gt_allele(ptr[1]) > 0 ) ac_sex[2+args->sex[i]-1]++;
        }
        // if genotype is not heterozygous, skip
        if ( bcf_gt_allele(ptr[0]) == bcf_gt_allele(ptr[1]) ||
           ( bcf_gt_allele(ptr[0]) != 0 && bcf_gt_allele(ptr[1]) != 0 ) ) continue;
        int8_t gt_phase = (int8_t)SIGN((ptr[0]>>1) - (ptr[1]>>1)) * (int8_t)bcf_gt_is_phased(ptr[1]);
        ac_het++;
        if ( args->sex && ( args->sex[i] == 1 || args->sex[i] == 2 ) ) ac_het_sex[args->sex[i]-1]++;
        if ( args->tx && gt_phase ) ac_het_tx[(1-gt_phase)/2]++;

        if ( nbdev_phase )
        {
            int8_t bdev_phase = SIGN(args->bdev_phase_arr[i]);
            if ( bdev_phase )
            {
                ase[(1-bdev_phase)/2]++;
                if ( args->tx && gt_phase ) ase_tx[(1-bdev_phase*gt_phase)/2]++;
            }
        }

        if ( nad )
        {
            int nalleles = nad / nsmpl;
            int ref_cnt = args->ad_arr[i*nalleles];
            int alt_cnt = 0;
            for (int j=1; j<nalleles; j++)
                alt_cnt += args->ad_arr[i*nalleles+j];
            ad_het[0] += ref_cnt;
            ad_het[1] += alt_cnt;
            baf = ((float)alt_cnt + 0.5f) / ((float)ref_cnt + (float)alt_cnt + 1.0f);
        }
        bcf_fmt_t *baf_fmt = bcf_get_fmt(args->in_hdr, rec, "BAF");
        if ( baf_fmt ) baf = ((float*)(baf_fmt->p + baf_fmt->size * i))[0];
        if ( args->tx && gt_phase && !isnan(baf) ) kv_push(float, kv_baf[(1-gt_phase)/2], baf);
    }

    bcf_update_info_int32(args->out_hdr, rec, "AC_Het", &ac_het, 1);
    if ( args->sex )
    {
        bcf_update_info_int32(args->out_hdr, rec, "AC_Het_Sex", &ac_het_sex, 2);
        double left, right, fisher;
        ret[0] = 0.0f - log10f(kt_fisher_exact(ac_sex[0], ac_sex[1], ac_sex[2], ac_sex[3], &left, &right, &fisher));
        bcf_update_info_float(args->out_hdr, rec, "AC_Sex_Test", &ret, 1);
    }
    if ( args->tx )
    {
        bcf_update_info_int32(args->out_hdr, rec, "AC_Het_Tx", &ac_het_tx, 2);
        ret[0] = 0.0f - log10f(binom_dist((ac_het_tx[0] + ac_het_tx[1]), 0.5, ac_het_tx[0]));
        bcf_update_info_float(args->out_hdr, rec, "AC_Het_Tx_Test", &ret, 1);
    }
    if ( args->ase )
    {
        bcf_update_info_int32(args->out_hdr, rec, "ASE", &ase, 2);
        ret[0] = 0.0f - log10f(binom_dist(ase[0] + ase[1], 0.5, ase[0]));
        bcf_update_info_float(args->out_hdr, rec, "ASE_Test", &ret, 1);
        if ( args->tx )
        {
            bcf_update_info_int32(args->out_hdr, rec, "ASE_Tx", &ase_tx, 2);
            ret[0] = 0.0f - log10f(binom_dist(ase_tx[0] + ase_tx[1], 0.5, ase_tx[0]));
            bcf_update_info_float(args->out_hdr, rec, "ASE_Tx_Test", &ret, 1);
        }
    }
    if ( args->ad )
    {
        bcf_update_info_int32(args->out_hdr, rec, "AD_Het", &ad_het, 2);
        ret[0] = 0.0f - log10f(binom_dist(ad_het[0] + ad_het[1], 0.5, ad_het[0]));
        bcf_update_info_float(args->out_hdr, rec, "AD_Het_Test", &ret, 1);
    }
    if ( args->tx && ( kv_baf[0].a || kv_baf[0].a ) )
    {
        ret[0] = get_median( kv_baf[0].a, kv_baf[0].n, NULL );
        ret[1] = get_median( kv_baf[1].a, kv_baf[1].n, NULL );
        ret[2] = 0.0f - log10f(welch_t_test(kv_baf[0].a, kv_baf[1].a, kv_baf[0].n, kv_baf[1].n));
        ret[3] = 0.0f - log10f(mann_whitney_u(kv_baf[0].a, kv_baf[1].a, kv_baf[0].n, kv_baf[1].n));
        bcf_update_info_float(args->out_hdr, rec, "BAF_Tx_Test", &ret, 4);
    }
    free(kv_baf[0].a);
    free(kv_baf[1].a);

end:
    // remove all samples if sites_only was selected
    if ( bcf_hdr_nsamples(args->out_hdr) == 0 ) bcf_subset(args->out_hdr, rec, 0, NULL);
    return rec;
}

void destroy(void)
{
    free(args->sex);
    free(args->gt_arr);
    free(args->ad_arr);
    free(args->bdev_phase_arr);
    free(args);
}
