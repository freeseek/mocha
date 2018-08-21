/* The MIT License

   Copyright (c) 2015-2018 Giulio Genovese

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
#include <ctype.h>
#include <math.h>
#include <htslib/kseq.h>
#include <htslib/vcf.h>
#include <htslib/synced_bcf_reader.h>
#include <htslib/ksort.h>
#include "regidx.h"
#include "kmin.h"
#include "bcftools.h"

// kv macro from https://github.com/attractivechaos/klib/blob/master/kvec.h by Attractive Chaos <attractor@live.co.uk>
#define kvec_t(type) struct { size_t n, m; type *a; }
#define kv_push(type, v, x) do { \
    if ((v).n == (v).m) { \
        (v).m = (v).m? (v).m<<1 : 2; \
        (v).a = (type *)realloc((v).a, sizeof(type) * (v).m); \
    } \
    (v).a[(v).n++] = (x); \
} while (0)

typedef kvec_t(int8_t) kv_int8_t;
typedef kvec_t(int16_t) kv_int16_t;
typedef kvec_t(int) kv_int;
typedef kvec_t(float) kv_float;

#define SIGN(x) (((x) > 0) - ((x) < 0))

#define MOCHA_VERSION "2018-08-21"

#define FLT_INCLUDE  (1<<0)
#define FLT_EXCLUDE  (1<<1)
#define WGS_DATA     (1<<2)
#define NO_LOG       (1<<3)
#define NO_ANNOT     (1<<4)

#define LRR 0
#define BAF 1
#define AD0 0
#define AD1 1

#define LRR_BAF   0
#define BAF_PHASE 1

#define SEX_UNK 0
#define SEX_MAL 1
#define SEX_FEM 2

#define MOCHA_UNK 0
#define MOCHA_DEL 1
#define MOCHA_DUP 2
#define MOCHA_UPD 3
#define MOCHA_CNP_DEL 4
#define MOCHA_CNP_DUP 5
#define MOCHA_CNP_CNV 6

#define MOCHA_NAN 0
#define MOCHA_ARM 1
#define MOCHA_TEL 2

#define MAX_ORDER 5

// this macro from ksort.h defines the function
// void ks_introsort_int(size_t n, int a[]);
KSORT_INIT_GENERIC(int)

// this macro from ksort.h defines the function
// float ks_ksmall_float(size_t n, float arr[], size_t kk);
KSORT_INIT_GENERIC(float)

static inline float sqf(float x) { return x*x; }
static inline double sq(double x) { return x*x; }


// default values for the model
const float xy_prob_default = 1e-09f;
const float err_prob_default = 1e-04f;
const float flip_prob_default = 1e-02f;
const float telomere_prob_default = 1e-02f;
const float centromere_prob_default = 1e-04f;
static const char *cnf_default = "1.0,1.5,3.0,4.0";
static const char *bdev_default = "6.0,8.0,10.0,15.0,20.0,30.0,50.0,80.0,100.0,150.0,200.0";
const float lrr_bias_default = 0.2f;
const int min_dist_default = 400;
const int median_baf_adjust_default = 5;
const int order_lrr_gc_default = 2;

/****************************************
 * CONVERT FLAOT TO INT16 AND VICEVERSA *
 ****************************************/

#define INT16_SCALE 1000 // most BAF values from Illumina and UKBB are scaled to 1000
#define INT16_NAN (int16_t)0x8000

inline int16_t float_to_int16(float in)
{
    return isnan(in) ? INT16_NAN : (int16_t)roundf(INT16_SCALE * in);
}

inline float int16_to_float(int16_t in)
{
    return in == INT16_NAN ? NAN : ((float)in) / INT16_SCALE;
}

static float *realloc_int16_to_float(int16_t *in, int n)
{
    float *out = (float *)malloc(n * sizeof(float));
    if ( !out ) return NULL;
    for (int i=0; i<n; i++) out[i] = int16_to_float(in[i]);
    free(in);
    return out;
}

/******************************************
 * LRR AND COVERAGE POLYNOMIAL REGRESSION *
 ******************************************/

// the following alternative code snippets were considered to perform GC regression:
// https://stackoverflow.com/questions/5083465/fast-efficient-least-squares-fit-algorithm-in-c
// https://github.com/natedomin/polyfit/blob/master/polyfit.c
// https://software.intel.com/sites/products/documentation/doclib/mkl_sa/11/mkl_lapack_examples/sgels_ex.c.htm
// this function needs to use doubles internally when dealing with WGS data
float *polyfit(const float *lrr, const float *gc, int n, const int *imap, int order)
{
    int m = order + 1;
    if (n < m || order > MAX_ORDER) return NULL;
    double B[MAX_ORDER+1] = {0.0};
    double P[((MAX_ORDER + 1) * 2)+1] = {0.0};
    double A[(MAX_ORDER + 1)*2*(MAX_ORDER + 1)] = {0.0};

    // identify the column vector
    for (int i=0; i<n; i++)
    {
        float x = imap ? gc[ imap[i] ]: gc[i];
        float y = lrr[i];
        if ( isnan(x) || isnan(y) ) continue;
        float powx = 1.0f;

        for (int j=0; j<m; j++)
        {
            B[j] += (double)(y * powx);
            powx *= x;
        }
    }

    // initialize the PowX array
    P[0] = (float)n;

    // compute the sum of the powers of X
    for (int i=0; i<n; i++)
    {
        float x = imap ? gc[ imap[i] ]: gc[i];
        if ( isnan(x) ) continue;
        float powx = x;

        for (int j=1; j<((2 * m) + 1); j++)
        {
            P[j] += (double)powx;
            powx *= x;
        }
    }

    // initialize the reduction matrix
    for (int i=0; i<m; i++)
    {
        for (int j=0; j<m; j++)
        {
            A[(i * (2 * m)) + j] = P[i+j];
        }

        A[(i*(2 * m)) + (i + m)] = 1.0;
    }

    // move the identity matrix portion of the redux matrix to the
    // left side (find the inverse of the left side of the redux matrix)
    for (int i=0; i<m; i++)
    {
        double x = A[(i * (2 * m)) + i];
        if (x != 0)
        {
            for (int k=0; k<(2 * m); k++)
            {
                A[(i * (2 * m)) + k] /= x;
            }

            for (int j=0; j<m; j++)
            {
                if (i != j)
                {
                    double y = A[(j * (2 * m)) + i];
                    for (int k=0; k<(2 * m); k++)
                    {
                        A[(j * (2 * m)) + k] -=
                            y * A[(i * (2 * m)) + k];
                    }
                }
            }
        }
        else
        {
            // cannot work with singular matrices
            return NULL;
        }
    }

    // calculate coefficients
    float *coeffs = (float *)malloc(m * sizeof(float));
    if ( !coeffs ) return NULL;
    for (int i=0; i<m; i++)
    {
        for (int j=0; j<m; j++)
        {
            double x = 0.0;
            for (int k=0; k<m; k++)
            {
                x += (A[(i * (2 * m)) + (k + m)] * B[k]);
            }
            coeffs[i] = (float)x;
        }
    }

    return coeffs;
}

static float *ad_to_lrr(int16_t *ad0, int16_t *ad1, int n)
{
    // this function keeps a list of logarithms of integers to minimize log calls
    static float *int_logf = NULL;
    static int nint_logf = 0, mint_logf = 0;
    if (ad0 == NULL && ad1 == NULL && n == 0) { free(int_logf); return NULL; }

    float *lrr = (float *)malloc(n * sizeof(n));
    for(int i=0; i<n; i++)
    {
        int cov = (int)(ad0[i]==INT16_NAN ? 0 : ad0[i]) + (int)(ad1[i]==INT16_NAN ? 0 : ad1[i]);
        if (cov==0)
        {
            lrr[i] = NAN;
        }
        else
        {
            if (cov > nint_logf)
            {
                hts_expand(float, cov, mint_logf, int_logf);
                for (int j=nint_logf; j<cov; j++) int_logf[j] = logf(j+1);
                nint_logf = cov;

            }
            lrr[i] = int_logf[cov-1];
        }
    }
    return lrr;
}

static float *ad_to_baf(int16_t *ad0, int16_t *ad1, int n)
{
    float *baf = (float *)malloc(n * sizeof(n));
    for(int i=0; i<n; i++)
        baf[i] = (ad0[i]==INT16_NAN || ad1[i]==INT16_NAN) ? NAN : (float)ad1[i] / (float)(ad0[i] + ad1[i]);
    return baf;
}

static void adjust_lrr(float *lrr, const float *gc, int n, const int *imap, const float *coeffs, int order)
{
    for (int i=0; i<n; i++)
    {
        float x = imap ? gc[ imap[i] ] : gc[i];
        float powx = 1.0f;
        for (int j=0; j<=order; j++)
        {
            lrr[i] -= coeffs[j] * powx;
            powx *= x;
        }
    }
}

// computes total sum of squares
// this function needs to use doubles internally when dealing with WGS data
static float get_tss(const float *v, int n)
{
    double mean = 0.0;
    int j = 0;
    for (int i = 0; i < n; i++)
    {
        if (!isnan(v[i]))
        {
            mean += (double)v[i];
            j++;
        }
    }
    if ( j <= 1 ) return NAN;
    mean /= (double)j;

    double tss = 0.0;
    for (int i = 0; i < n; i++)
    {
        if (!isnan(v[i])) tss += sq((double)v[i] - mean);
    }
    return (float)tss;
}

/*********************************
 * HMM AND OPTIMIZATION METHODS  *
 *********************************/

// compute Viterbi path from probabilities
static int8_t *retrace_viterbi(int T, int N, const float *prb, const int8_t *ptr)
{
    int i, t;
    int8_t *path = (int8_t *)malloc(T * sizeof(int8_t));

    // initialize last path state
    path[T-1] = 0;
    for (i=1; i<N; i++)
        if (prb[(int)path[T-1]] < prb[i])
            path[T-1] = (int8_t)i;

    // compute best path by tracing back the Markov chain
    for (t=T-1; t>0; t--)
        path[t-1] = ptr[(t-1) * N + (int)path[t]];

    return path;
}

// rescale Viterbi probabilities to avoid underflow issues
static void rescale_prb(float *prb, int N)
{
    int i;
    float sum;
    for (sum=0, i=0; i<N; sum += prb[i++]);
    for (i=0; i<N; prb[i++] /= sum);
}

// compute the Viterbi path from BAF
// T is the length of the hidden Markov model
// m is the number of possible BAF deviations
// TODO find better terminology for all variables and make sure they match the variables in the manuscript
static int8_t *get_viterbi(const float *emis, int T, int m, float xy_prob, float flip_prob, float telomere_prob, float centromere_prob, int last_p, int first_q)
{
    int t, i, j, changeidx;
    float sqrt_xy_prob = sqrtf(xy_prob);
    float sqrt_centromere_prob = sqrtf(centromere_prob);

    // determine the number of hidden states based on whether phase information is used
    int N = 1 + m + ( isnan(flip_prob) ? 0 : m);

    // allocate memory necessary for running the algorithm
    float *prb = (float *)malloc(N * sizeof(float));
    float *newprb = (float *)malloc(N * sizeof(float));
    int8_t *ptr = (int8_t *)malloc(N * (T-1) * sizeof(int8_t));
    int8_t *path;

    // initialize and rescale the first state
    prb[0] = emis[0];
    for (i=1; i<N; i++) prb[i] = xy_prob / (last_p == 0 ? sqrt_centromere_prob : telomere_prob) * emis[i];
    rescale_prb(prb, N);

    // TODO do I want to include this probability switch as a parameter?
    float switch_alternate = xy_prob * sqrt_xy_prob; // switching across alternate states incurs a penalty

    // compute best probabilities at each position
    for (t=1; t<T; t++)
    {
        // this causes a penalty for mosaic chromosomal calls across the centromeres
        float exit_prob = t>last_p ? xy_prob * sqrt_centromere_prob : xy_prob / sqrt_centromere_prob;
        float enter_prob = t<first_q ? xy_prob * sqrt_centromere_prob : xy_prob / sqrt_centromere_prob;

        for (i=0; i<N; i++)
        {
            newprb[i] = prb[i];
            ptr[(t-1) * N + i] = (int8_t)i;
        }

        // compute whether a state switch should be considered for null state
        for (i=1; i<N; i++)
        {
             if (newprb[0] < prb[i] * exit_prob)
             {
                 newprb[0] = prb[i] * exit_prob;
                 ptr[(t-1) * N] = ptr[(t-1) * N + i];
             }
             if (newprb[i] < prb[0] * enter_prob)
             {
                 newprb[i] = prb[0] * enter_prob;
                 ptr[(t-1) * N + i] = ptr[(t-1) * N];
             }
        }

        // compute whether a state switch should be considered for each other state
        // it will run twice if and only if phasing is used
        for (j=0; j==0 || (!isnan(flip_prob) && j==m); j+=m)
        {
            float changeprb = prb[0] * enter_prob; changeidx = 0;
            for (i=0; i<m; i++)
            {
                if (changeprb < prb[1+j+i] * switch_alternate)
                {
                    changeprb = prb[1+j+i] * switch_alternate;
                    changeidx = 1+j+i;
                }
            }
            for (i=0; i<m; i++)
            {
                if (newprb[1+j+i] < changeprb)
                {
                    newprb[1+j+i] = changeprb;
                    ptr[(t-1) * N + 1+j+i] = ptr[(t-1) * N + changeidx];
                }
            }
        }

        // compute whether a phase flip should be considered for non-null states
        if (!isnan(flip_prob))
        {
            for (i=0; i<m; i++)
            {
                if (newprb[1+i] < newprb[1+m+i] * flip_prob)
                {
                    newprb[1+i] = newprb[1+m+i] * flip_prob;
                    ptr[(t-1) * N + 1+i] = ptr[(t-1) * N + 1+m+i];
                }
            if (newprb[1+m+i] < newprb[1+i] * flip_prob)
                {
                    newprb[1+m+i] = newprb[1+i] * flip_prob;
                    ptr[(t-1) * N + 1+m+i] = ptr[(t-1) * N + 1+i];
                }
            }
        }

        // update and rescale the current state
        newprb[0] *= emis[t*N];
        for (i=0; i<m; i++)
        {
            newprb[1+i] *= emis[t*N + 1+i];
            if (!isnan(flip_prob)) newprb[1+m+i] *= emis[t*N + 1+m+i];
        }
        for (i=0; i<N; i++) prb[i] = newprb[i];
        rescale_prb(prb, N);
    }

    // add closing cost to the last state
    for (i=1; i<N; i++) prb[i] *= xy_prob / (first_q == T ? sqrt_centromere_prob : telomere_prob);
    rescale_prb(prb, N);

    path = retrace_viterbi(T, N, prb, ptr);

    // free memory
    free(prb);
    free(newprb);
    free(ptr);

    // symmetrize the path
    if (!isnan(flip_prob))
        for (i=0; i<T; i++)
            if (path[i]>m)
                path[i] = (int8_t)m - path[i];

    return path;
}

/*********************************
 * LRR AND BAF LIKELIHOODS       *
 *********************************/

static inline float normal_pdf(float x, float m, float s, float w)
{
    static const float log_inv_sqrt_2pi = -0.918938533204672669541f;
    float a = (x - m) / s;
    return expf( ( log_inv_sqrt_2pi - 0.5f * a * a ) * w );
}

// lrr_bias is used in a different way from what done by Petr Danecek in bcftools/vcfcnv.c
// TODO lrr_sd_bias is just a constant here, it could be removed
static inline float lkl_lrr_baf(float lrr, float baf, float ldev, float bdev, float lrr_sd, float baf_sd, float lrr_bias, float lrr_sd_bias, float err_prob)
{
    float ret = isnan( lrr ) ? 1.0f : normal_pdf( lrr, ldev, lrr_sd, lrr_bias ) / lrr_sd_bias;
    if ( !isnan(baf) ) ret *= ( normal_pdf( baf - 0.5f,  bdev, baf_sd, 1.0f ) +
                                normal_pdf( baf - 0.5f, -bdev, baf_sd, 1.0f ) ) * 0.5f / baf_sd;
    return ret < err_prob ? err_prob : ret;
}

static inline float lkl_baf_phase(float baf, int8_t phase, float bdev, float baf_sd, float err_prob)
{
    if ( isnan( baf ) ) return 1.0f;
    float ret;
    if ( phase == 0 ) ret = ( normal_pdf( baf - 0.5f,  bdev, baf_sd, 1.0f ) +
                              normal_pdf( baf - 0.5f, -bdev, baf_sd, 1.0f ) ) * 0.5f / baf_sd;
    else ret = normal_pdf( ( baf - 0.5f ) * (float)SIGN( phase ), bdev, baf_sd, 1.0f ) / baf_sd;
    return ret < err_prob ? err_prob : ret;
}

// precomupute emission probabilities
static float *get_lrr_baf_emis(const float *lrr, const float *baf, int T, const int *imap, float err_prob, float lrr_bias, float lrr_hap2dip, float lrr_sd, float baf_sd, const float *cnf, int m)
{
    float lrr_sd_bias = powf( lrr_sd, lrr_bias );
    float *ldev = (float *)malloc(m * sizeof(float));
    for (int i=0; i<m; i++) ldev[i] = ( logf(cnf[i]) / (float)M_LN2 - 1.0f ) * lrr_hap2dip;
    int N = 1 + m;
    float *emis = (float *)malloc(N * T * sizeof(float));
    for (int t=0; t<T; t++)
    {
        float x = imap ? lrr[ imap[t] ] : lrr[t];
        float y = imap ? baf[ imap[t] ] : baf[t];
        emis[t*N] = lkl_lrr_baf( x, y, 0.0f, 0.0f, lrr_sd, baf_sd, lrr_bias, lrr_sd_bias, err_prob );
        for (int i=0; i<m; i++)
        {
            float bdev = fabsf( 0.5f - 1 / cnf[i] );
            emis[t*N + 1+i] = lkl_lrr_baf( x, y, ldev[i], bdev, lrr_sd, baf_sd, lrr_bias, lrr_sd_bias, err_prob );
        }
    }
    free(ldev);
    return emis;
}

// precomupute emission probabilities
static float *get_baf_phase_emis(const float *baf, const int8_t *gt_phase, int T, const int *imap, float err_prob, float flip_prob, float baf_sd, const float *bdev, int m)
{
    int N = 1 + m + ( isnan(flip_prob) ? 0 : m);
    float *emis = (float *)malloc(N * T * sizeof(float));
    for (int t=0; t<T; t++)
    {
        float x = imap ? baf[ imap[t] ] : baf[t];
        int8_t p = imap ? gt_phase[ imap[t] ] : gt_phase[t];
        emis[t*N] = lkl_baf_phase( x, (int8_t)1, 0.0f, baf_sd, err_prob );
        for (int i=0; i<m; i++)
        {
            if( isnan(flip_prob) )
                emis[t*N + 1+i] = lkl_baf_phase( x, (int8_t)0, bdev[i], baf_sd, err_prob );
            if( !isnan(flip_prob) )
            {
                emis[t*N + 1+i  ] = lkl_baf_phase( x, p, bdev[i], baf_sd, err_prob );
                if ( p == 0 ) emis[t*N + 1+m+i] = emis[t*N + 1+i];
                else emis[t*N + 1+m+i] = lkl_baf_phase( x, p, -bdev[i], baf_sd, err_prob );
            }
        }
    }
    return emis;
}

// TODO see if the sumlog.c code can be replaced by MIT license code
#include "sumlog.c" // http://homepages.inf.ed.ac.uk/imurray2/code/ by Iain Murray <i.murray@ed.ac.uk>

// return the log10 likelihood for a segment
static double log10_lkl_lrr_baf(const float *lrr,
                                const float *baf,
                                int n,
                                const int *imap,
                                float err_prob,
                                float lrr_bias,
                                float lrr_hap2dip,
                                float lrr_sd,
                                float baf_sd,
                                double cnf)
{
    if ( n==0 || cnf < 0.0 || cnf > 4.0 ) return -INFINITY; // kmin_brent does not handle NAN
    float lrr_sd_bias = powf( lrr_sd, lrr_bias );
    float ldev = ( logf((float)cnf) / (float)M_LN2 - 1.0f ) * lrr_hap2dip;
    float bdev = fabsf( 0.5f - 1.0f / (float)cnf );
    float *v = (float *)calloc(n, sizeof(float));
    for (int i=0; i<n; i++)
    {
        float x = imap ? lrr[ imap[i] ] : lrr[i];
        float y = imap ? baf[ imap[i] ] : baf[i];
        v[i] = lkl_lrr_baf( x, y, ldev, bdev, lrr_sd, baf_sd, lrr_bias, lrr_sd_bias, err_prob );
    }
    float ret = 0.0f;
    sumlogf(v, n, &ret);
    free(v);
    return (double)ret / M_LN10;
}

// return the log10 likelihood for a segment
static double log10_lkl_baf_phase(const float *baf,
                                  const int8_t *gt_phase,
                                  int n,
                                  const int *imap,
                                  const int8_t *bdev_phase,
                                  float err_prob,
                                  float flip_prob,
                                  float baf_sd,
                                  double bdev)
{
    if ( n==0 || bdev < -0.5 || bdev > 0.5 ) return -INFINITY; // kmin_brent does not handle NAN
    float *v = (float *)calloc(n, sizeof(float));
    for (int i=0; i<n; i++)
    {
        float x = imap ? baf[ imap[i] ] : baf[i];
        int8_t p = imap ? gt_phase[ imap[i] ] : gt_phase[i];
        if ( bdev_phase ) p *= (int8_t)SIGN( bdev_phase[i] ); // notice bdev_phase has no imap
        v[i] = lkl_baf_phase( x, isnan(flip_prob) ? 0 : p, (float)bdev, baf_sd, err_prob );
    }
    float ret = 0.0f;
    sumlogf(v, n, &ret);
    free(v);
    return (double)ret / M_LN10;
}

// TODO need a better title here
static float compare_models(const float *baf,
                            const int8_t *gt_phase,
                            int n,
                            const int *imap,
                            float xy_prob,
                            float err_prob,
                            float flip_prob,
                            float telomere_prob,
                            float baf_sd,
                            const float *bdev,
                            int m)
{
    if ( n == 0 ) return NAN;
    float *emis = get_baf_phase_emis(baf, gt_phase, n, imap, err_prob, flip_prob, baf_sd, bdev, m);
    int8_t *path = get_viterbi(emis, n, m, xy_prob, flip_prob, telomere_prob, 1.0f, 0, 0); // TODO can't I pass these values instead of 0 0?
    free(emis);
    int nflips = 0;
    for (int i=1; i<n; i++) if ( path[i-1] && path[i] && path[i-1] != path[i] ) nflips++;
    double f(double x, void *data) { return -log10_lkl_baf_phase(baf, gt_phase, n, imap, path, err_prob, flip_prob, baf_sd, x); }
    double x, fx = kmin_brent(f, -0.5, 0.5, NULL, KMIN_EPS, &x);
    float lod = (float)(f(0.0f, NULL) - fx) + (float)nflips * log10f(flip_prob);
    free(path);
    return lod;
}

/*********************************
 * WGS AD LIKELIHOODS            *
 *********************************/

static inline float lkl_lrr_ad(float lrr,
                               int16_t ad0,
                               int16_t ad1,
                               float ldev,
                               float lrr_sd,
                               float lrr_bias,
                               float lrr_sd_bias,
                               float *log_gamma,
                               float *lod_gamma_alpha,
                               float *lod_gamma_beta,
                               float *lod_gamma_alpha_beta,
                               float err_prob)
{
    float ret = isnan( lrr ) ? 1.0f : normal_pdf( lrr, ldev, lrr_sd, lrr_bias ) / lrr_sd_bias;
    if ( ad0!=INT16_NAN && ad1!=INT16_NAN )
        ret *= ( expf( lod_gamma_alpha[ad0] + lod_gamma_beta[ad1] - lod_gamma_alpha_beta[ad0 + ad1] + log_gamma[ad0 + ad1] - log_gamma[ad0] - log_gamma[ad1] ) +
                 expf( lod_gamma_alpha[ad1] + lod_gamma_beta[ad0] - lod_gamma_alpha_beta[ad0 + ad1] + log_gamma[ad0 + ad1] - log_gamma[ad0] - log_gamma[ad1] ) ) * 0.5f;
    return ret < err_prob ? err_prob : ret;
}

static inline float lkl_ad_phase(int16_t ad0,
                                 int16_t ad1,
                                 int8_t phase,
                                 float *log_gamma,
                                 float *lod_gamma_alpha,
                                 float *lod_gamma_beta,
                                 float *lod_gamma_alpha_beta,
                                 float err_prob)
{
    if ( ad0==INT16_NAN || ad1==INT16_NAN ) return 1.0f;
    float ret;
    if (phase == 0)
        ret = ( expf( lod_gamma_alpha[ad0] + lod_gamma_beta[ad1] - lod_gamma_alpha_beta[ad0 + ad1] + log_gamma[ad0 + ad1] - log_gamma[ad0] - log_gamma[ad1] ) +
                expf( lod_gamma_alpha[ad1] + lod_gamma_beta[ad0] - lod_gamma_alpha_beta[ad0 + ad1] + log_gamma[ad0 + ad1] - log_gamma[ad0] - log_gamma[ad1] ) ) * 0.5f;
    else if (phase>0) ret = expf( lod_gamma_alpha[ad0] + lod_gamma_beta[ad1] - lod_gamma_alpha_beta[ad0 + ad1] + log_gamma[ad0 + ad1] - log_gamma[ad0] - log_gamma[ad1]);
    else ret = expf( lod_gamma_alpha[ad1] + lod_gamma_beta[ad0] - lod_gamma_alpha_beta[ad0 + ad1] + log_gamma[ad0 + ad1] - log_gamma[ad0] - log_gamma[ad1]);
    return ret < err_prob ? err_prob : ret;
}

// precompute table of gammas of integer values such that
// log_gamma[n] = log(n!)
// see https://en.wikipedia.org/wiki/Gamma_function
static float *precompute_log_gammas(const int16_t *ad0, const int16_t *ad1, int n, const int *imap)
{
    static float *log_gamma = NULL;
    static int nlog_gamma = 1, mlog_gamma = 0;
    if (ad0 == NULL && ad1 == NULL && n == 0 && imap == NULL) { free(log_gamma); return NULL; }

    int max = 0;
    for (int i=0; i<n; i++)
    {
        int a = imap ? ad0[ imap[i] ] : ad0[i];
        int b = imap ? ad1[ imap[i] ] : ad1[i];
        if ( a!=INT16_NAN && b!=INT16_NAN )
        {
            if (max < a + b) max = a + b;
        }
    }

    if (max + 1 > nlog_gamma)
    {
        hts_expand(float, max + 1, mlog_gamma, log_gamma);
        log_gamma[0] = 0.0f;
        for (int i=nlog_gamma; i<=max; i++) log_gamma[i] = log_gamma[i-1] + logf(i);
        nlog_gamma = max + 1;
    }

    return log_gamma;
}

// precompute table of values for the beta binomial log likelihoods such that
// lod_gamma_alpha[n] = log( \gamma(n+\alpha) / \gamma(\alpha) )
// lod_gamma_beta[n] = log( \gamma(n+\beta) / \gamma(\beta) )
// lod_gamma_alpha_beta[n] = log( \gamma(n+\alpha+\beta) / \gamma(\alpha+\beta) )
// see https://en.wikipedia.org/wiki/Beta-binomial_distribution#As_a_compound_distribution
// TODO change these to double as you need the extra precision
static void precompute_lod_gammas(const int16_t *ad0,
                                  const int16_t *ad1,
                                  int n,
                                  const int *imap,
                                  float alpha,
                                  float beta,
                                  float **ptr_lod_gamma_alpha,
                                  float **ptr_lod_gamma_beta,
                                  float **ptr_lod_gamma_alpha_beta)
{
    static float *lod_gamma_alpha = NULL;
    static float *lod_gamma_beta = NULL;
    static float *lod_gamma_alpha_beta = NULL;
    static float curr_alpha = NAN, curr_beta = NAN;
    static int nlod_gamma_alpha = 1, mlod_gamma_alpha = 0;
    static int nlod_gamma_beta = 1, mlod_gamma_beta = 0;
    static int nlod_gamma_alpha_beta = 1, mlod_gamma_alpha_beta = 0;
    if (ad0 == NULL && ad1 == NULL && n == 0 && imap == NULL)
    {
        free(lod_gamma_alpha);
        free(lod_gamma_beta);
        free(lod_gamma_alpha_beta);
        return;
    }

    if (curr_alpha != alpha || curr_beta != beta)
    {
        nlod_gamma_alpha = 1;
        nlod_gamma_beta = 1;
        nlod_gamma_alpha_beta = 1;
        curr_alpha = alpha;
        curr_beta = beta;
    }

    int max1 = 0, max2 = 0;
    for (int i=0; i<n; i++)
    {
        int a = imap ? ad0[ imap[i] ] : ad0[i];
        int b = imap ? ad1[ imap[i] ] : ad1[i];
        if ( a!=INT16_NAN && b!=INT16_NAN )
        {
            if (a > max1) max1 = a;
            if (b > max1) max1 = b;
            if (a + b > max2) max2 = a + b;
        }
    }

    if (max1 + 1 > nlod_gamma_alpha)
    {
        hts_expand(float, max1 + 1, mlod_gamma_alpha, lod_gamma_alpha);
        lod_gamma_alpha[0] = 0.0f;
        float tmp = alpha;
        for (int i=nlod_gamma_alpha; i<=max1; i++) lod_gamma_alpha[i] = lod_gamma_alpha[i-1] + logf(tmp++);
        nlod_gamma_alpha = max1 + 1;
    }

    if (max1 + 1 > nlod_gamma_beta)
    {
        hts_expand(float, max1 + 1, mlod_gamma_beta, lod_gamma_beta);
        lod_gamma_beta[0] = 0.0f;
        float tmp = beta;
        for (int i=nlod_gamma_beta; i<=max1; i++) lod_gamma_beta[i] = lod_gamma_beta[i-1] + logf(tmp++);
        nlod_gamma_beta = max1 + 1;
    }

    if (max2 + 1 > nlod_gamma_alpha_beta)
    {
        hts_expand(float, max2 + 1, mlod_gamma_alpha_beta, lod_gamma_alpha_beta);
        lod_gamma_alpha_beta[0] = 0.0f;
        float tmp = alpha + beta;
        for (int i=nlod_gamma_alpha_beta; i<=max2; i++) lod_gamma_alpha_beta[i] = lod_gamma_alpha_beta[i-1] + logf(tmp++);
        nlod_gamma_alpha_beta = max2 + 1;
    }

    *ptr_lod_gamma_alpha = lod_gamma_alpha;
    *ptr_lod_gamma_beta = lod_gamma_beta;
    *ptr_lod_gamma_alpha_beta = lod_gamma_alpha_beta;
}

static float *get_lrr_ad_emis(const float *lrr,
                              const int16_t *ad0,
                              const int16_t *ad1,
                              int T,
                              const int *imap,
                              float err_prob,
                              float lrr_bias,
                              float lrr_hap2dip,
                              float lrr_sd,
                              float ad_rho,
                              const float *cnf,
                              int m)
{
    float *log_gamma = precompute_log_gammas(ad0, ad1, T, imap);
    float *lod_gamma_alpha, *lod_gamma_beta, *lod_gamma_alpha_beta;
    float lrr_sd_bias = powf( lrr_sd, lrr_bias );
    float *ldev = (float *)malloc((1+m) * sizeof(float));
    ldev[0] = 0.0f;
    for (int i=1; i<1+m; i++) ldev[i] = ( logf(cnf[i-1]) / (float)M_LN2 - 1.0f ) * lrr_hap2dip;
    float *emis = (float *)malloc((1+m) * T * sizeof(float));
    for (int i=0; i<1+m; i++)
    {
        float alpha = ( 0.5f + ( i==0 ? 0.0f : fabsf( 0.5f - 1.0f / cnf[i-1] ) ) ) * ( 1.0f - ad_rho ) / ad_rho;
        float beta = ( 0.5f - ( i==0 ? 0.0f : fabsf( 0.5f - 1.0f / cnf[i-1] ) ) ) * ( 1.0f - ad_rho ) / ad_rho;
        precompute_lod_gammas(ad0, ad1, T, imap, alpha, beta, &lod_gamma_alpha, &lod_gamma_beta, &lod_gamma_alpha_beta);
        for (int t=0; t<T; t++)
        {
            float x = imap ? lrr[ imap[t] ] : lrr[t];
            int16_t a = imap ? ad0[ imap[t] ] : ad0[t];
            int16_t b = imap ? ad1[ imap[t] ] : ad1[t];
            emis[t*(1+m) + i] = lkl_lrr_ad( x, a, b, ldev[i], lrr_sd, lrr_bias, lrr_sd_bias, log_gamma, lod_gamma_alpha, lod_gamma_beta, lod_gamma_alpha_beta, err_prob );
        }
    }
    free(ldev);
    return emis;
}

static float *get_ad_phase_emis(const int16_t *ad0,
                                const int16_t *ad1,
                                const int8_t *gt_phase,
                                int T,
                                const int *imap,
                                float err_prob,
                                float flip_prob,
                                float ad_rho,
                                const float *bdev,
                                int m)
{
    float *log_gamma = precompute_log_gammas(ad0, ad1, T, imap);
    float *lod_gamma_alpha, *lod_gamma_beta, *lod_gamma_alpha_beta;
    float *emis = (float *)malloc((1+2*m) * T * sizeof(float));
    for (int i=0; i<1+m; i++)
    {
        float alpha = ( 0.5f + ( i==0 ? 0.0f : bdev[i-1] ) ) * ( 1.0f - ad_rho ) / ad_rho;
        float beta = ( 0.5f - ( i==0 ? 0.0f : bdev[i-1] ) ) * ( 1.0f - ad_rho ) / ad_rho;
        precompute_lod_gammas(ad0, ad1, T, imap, alpha, beta, &lod_gamma_alpha, &lod_gamma_beta, &lod_gamma_alpha_beta);
        for (int t=0; t<T; t++)
        {
            int16_t a = imap ? ad0[ imap[t] ] : ad0[t];
            int16_t b = imap ? ad1[ imap[t] ] : ad1[t];
            int8_t p = imap ? gt_phase[ imap[t] ] : gt_phase[t];
            emis[t*(1+2*m) + i] = lkl_ad_phase( a, b, p, log_gamma, lod_gamma_alpha, lod_gamma_beta, lod_gamma_alpha_beta, err_prob );
            if (i>0) emis[t*(1+2*m) + m + i] = lkl_ad_phase( a, b, p, log_gamma, lod_gamma_beta, lod_gamma_alpha, lod_gamma_alpha_beta, err_prob );
// if (i>0 && p==0) { fprintf(stderr, "%.4f %.4f\n", emis[t*(1+2*m) + i], emis[t*(1+2*m) + m + i]); exit(-1); }
        }
    }
    return emis;
}

// return the log10 likelihood for a segment
static double log10_lkl_lrr_ad(const float *lrr,
                               const int16_t *ad0,
                               const int16_t *ad1,
                               int n,
                               const int *imap,
                               float err_prob,
                               float lrr_bias,
                               float lrr_hap2dip,
                               float lrr_sd,
                               float ad_rho,
                               double cnf)
{
    if ( n==0 || cnf < 0.0 || cnf > 4.0 ) return -INFINITY; // kmin_brent does not handle NAN
    float *log_gamma = precompute_log_gammas(ad0, ad1, n, imap);
    float *lod_gamma_alpha, *lod_gamma_beta, *lod_gamma_alpha_beta;
    float lrr_sd_bias = powf( lrr_sd, lrr_bias );
    float ldev = ( logf((float)cnf) / (float)M_LN2 - 1.0f ) * lrr_hap2dip;
    float alpha = ( 0.5f + fabsf( 0.5f - 1.0f / cnf ) ) * ( 1.0f - ad_rho ) / ad_rho;
    float beta = ( 0.5f - fabsf( 0.5f - 1.0f / cnf ) ) * ( 1.0f - ad_rho ) / ad_rho;
    precompute_lod_gammas(ad0, ad1, n, imap, alpha, beta, &lod_gamma_alpha, &lod_gamma_beta, &lod_gamma_alpha_beta);
    float *v = (float *)calloc(n, sizeof(float));
    for (int i=0; i<n; i++)
    {
        float x = imap ? lrr[ imap[i] ] : lrr[i];
        int16_t a = imap ? ad0[ imap[i] ] : ad0[i];
        int16_t b = imap ? ad1[ imap[i] ] : ad1[i];
        v[i] = lkl_lrr_ad( x, a, b, ldev, lrr_sd, lrr_bias, lrr_sd_bias, log_gamma, lod_gamma_alpha, lod_gamma_beta, lod_gamma_alpha_beta, err_prob );
    }
    float ret = 0.0f;
    sumlogf(v, n, &ret);
    free(v);
    return (double)ret / M_LN10;
}

// return the log10 likelihood for a segment
static double log10_lkl_ad_phase(const int16_t *ad0,
                                 const int16_t *ad1,
                                 const int8_t *gt_phase,
                                 int n,
                                 const int *imap,
                                 const int8_t *bdev_phase,
                                 float err_prob,
                                 float flip_prob,
                                 float ad_rho,
                                 double bdev)
{
    if ( n==0 || bdev < -0.5 || bdev > 0.5 ) return -INFINITY; // kmin_brent does not handle NAN
    float *log_gamma = precompute_log_gammas(ad0, ad1, n, imap);
    float *lod_gamma_alpha, *lod_gamma_beta, *lod_gamma_alpha_beta;
    float alpha = ( 0.5f + bdev ) * ( 1.0f - ad_rho ) / ad_rho;
    float beta = ( 0.5f - bdev ) * ( 1.0f - ad_rho ) / ad_rho;
    precompute_lod_gammas(ad0, ad1, n, imap, alpha, beta, &lod_gamma_alpha, &lod_gamma_beta, &lod_gamma_alpha_beta);
    float *v = (float *)calloc(n, sizeof(float));
    for (int i=0; i<n; i++)
    {
        int16_t a = imap ? ad0[ imap[i] ] : ad0[i];
        int16_t b = imap ? ad1[ imap[i] ] : ad1[i];
        int8_t p = imap ? gt_phase[ imap[i] ] : gt_phase[i];
        if ( bdev_phase ) p *= (int8_t)SIGN( bdev_phase[i] ); // notice bdev_phase has no imap
        v[i] = lkl_ad_phase( a, b, p, log_gamma, lod_gamma_alpha, lod_gamma_beta, lod_gamma_alpha_beta, err_prob );
    }
    float ret = 0.0f;
    sumlogf(v, n, &ret);
    free(v);
    return (double)ret / M_LN10;
}

// TODO need a better title here
static float compare_wgs_models(const int16_t *ad0,
                            const int16_t *ad1,
                            const int8_t *gt_phase,
                            int n,
                            const int *imap,
                            float xy_prob,
                            float err_prob,
                            float flip_prob,
                            float telomere_prob,
                            float ad_rho,
                            const float *bdev,
                            int m)
{
    if ( n == 0 ) return NAN;
    float *emis = get_ad_phase_emis(ad0, ad1, gt_phase, n, imap, err_prob, flip_prob, ad_rho, bdev, m);
    int8_t *path = get_viterbi(emis, n, m, xy_prob, flip_prob, telomere_prob, 1.0f, 0, 0); // TODO can't I pass these values instead of 0 0?
    free(emis);
    int nflips = 0;
    for (int i=1; i<n; i++) if ( path[i-1] && path[i] && path[i-1] != path[i] ) nflips++;
    double f(double x, void *data) { return -log10_lkl_ad_phase(ad0, ad1, gt_phase, n, imap, path, err_prob, flip_prob, ad_rho, x); }
    double x, fx = kmin_brent(f, -0.5, 0.5, NULL, KMIN_EPS, &x);
    float lod = (float)(f(0.0f, NULL) - fx) + (float)nflips * log10f(flip_prob);
    free(path);
    return lod;
}

static double log10_lkl_beta_binomial(const int16_t *ad0,
                                      const int16_t *ad1,
                                      int n,
                                      double ad_rho)
{
    if ( n==0 || ad_rho <= 0.0 || ad_rho >= 1.0 ) return -INFINITY;
    float alpha = 0.5f * ( 1.0f - (float)ad_rho ) / (float)ad_rho;
    float *log_gamma = precompute_log_gammas(ad0, ad1, n, NULL);
    float *lod_gamma_alpha, *lod_gamma_beta, *lod_gamma_alpha_beta;
    precompute_lod_gammas(ad0, ad1, n, NULL, alpha, alpha, &lod_gamma_alpha, &lod_gamma_beta, &lod_gamma_alpha_beta);
    float ret = 0.0f;
    for (int i=0; i<n; i++)
    {
        int16_t a = ad0[i];
        int16_t b = ad1[i];
        if ( a!=INT16_NAN && b!=INT16_NAN )
            ret += lod_gamma_alpha[a] + lod_gamma_beta[b] - lod_gamma_alpha_beta[a+b] +
                   log_gamma[a+b] - log_gamma[a] - log_gamma[b];
    }
    return (double)ret / M_LN10;
}

/*********************************
 * BASIC STATISTICS FUNCTIONS    *
 *********************************/

// iterator of non-NaN values
static inline float next_not_nan(const float *v, const int *imap, int n, int *i)
{
    float x = NAN;
    while (*i < n)
    {
         x = imap ? v[ imap[*i] ] : v[*i];
         if ( !isnan(x) ) break;
         (*i)++;
    }
    return x;
}

// compute BAF phase concordance for a float array with iterator
static float get_baf_conc(const float *baf, const int8_t *gt_phase, int n, const int *imap)
{
    int i, a = 0, b = 0;
    float prev = NAN, next = NAN;
    for (i=0; i<n; i++)
    {
        prev = ( ( imap ? baf[ imap[i] ] : baf[i] ) - 0.5f ) * ( imap ? gt_phase[ imap[i] ] : gt_phase[i] );
        if (!isnan(prev) && prev != 0.0f) break;
    }
    if ( i == n ) return NAN;

    for (i++; i<n; i++)
    {
        next = ( ( imap ? baf[ imap[i] ] : baf[i] ) - 0.5f ) * ( imap ? gt_phase[ imap[i] ] : gt_phase[i] );
        if (!isnan(next) && next != 0.0f)
        {
            if (prev * next > 0.0f) a++;
            else b++;
            prev = next;
        }
    }
    return (float)a/(float)(a+b);
}

static float get_sample_sd(const float *v, int n, const int *imap);

// compute (adjusted) LRR autocorrelation for a float array with iterator
static float get_lrr_auto(const float *lrr, int n, const int *imap)
{
    float value;
    double mean = 0.0;
    int i = 0, j = 0;
    for (value = next_not_nan(lrr, imap, n, &i); i < n; i++, value = next_not_nan(lrr, imap, n, &i))
    {
        mean += (double)value;
        j++;
    }
    if ( j <= 1 ) return NAN;
    mean /= (double)j;

    double var = 0.0;
    i = 0;
    for (value = next_not_nan(lrr, imap, n, &i); i < n; i++, value = next_not_nan(lrr, imap, n, &i))
    {
        var += sq((double)value - mean);
    }

    double auto_corr = 0.0;
    i = 0;
    double prev = (double)next_not_nan(lrr, imap, n, &i) - mean, next;
    for (i++, value = next_not_nan(lrr, imap, n, &i); i < n ; i++, value = next_not_nan(lrr, imap, n, &i))
    {
        next = (double)value - mean;
        auto_corr += prev * next;
        prev = next;
    }
    auto_corr /= var;
    return auto_corr;
}

// compute the n50 of a vector
static int get_n50(const int *v, int n, const int *imap)
{
    if ( n <= 1 ) return -1;
    int i;
    int sum, sum2;
    int *w = (int *)malloc((n-1) * sizeof(int));

    for (i = 0, sum = 0; i < n-1; i++)
    {
        w[i] = imap ? v[imap[i+1]] - v[imap[i]] : v[i+1] - v[i];
        sum += w[i];
    }
    sum /= 2;

    ks_introsort_int((size_t)n-1, w);

    for (i = 0, sum2 = 0; sum2 < sum && i < n-1; i++) sum2 += w[i];
    int n50 = w[i-1];
    free(w);
    return n50;
}

// compute sample standard deviation of a float array (with iterator)
// sqrt ( ( \sum x^2 - (\sum x)^2 / N ) / ( N - 1 ) )
// TODO this is not okay as the above way to compute the sd is error prone with floats!!!
static float get_sample_sd(const float *v, int n, const int *imap)
{
    // float s = 0.0, s2 = 0.0;
    double mean = 0.0;
    int j = 0;
    for (int i = 0; i < n; i++)
    {
        double tmp = (double)(imap ? v[imap[i]] : v[i]);
        if ( !isnan(tmp) )
        {
            mean += tmp;
            j++;
        }
    }
    if ( j <= 1 ) return NAN;
    mean /= (double)j;

    double s2 = 0.0;
    for (int i = 0; i < n; i++)
    {
        double tmp = (double)(imap ? v[imap[i]] : v[i]);
        if ( !isnan(tmp) ) s2 += sq(tmp - mean);
    }
    s2 /= (double)(j-1);

    return (float)sqrt(s2);
}

// compute standard error of mean a float array (with iterator)
// sqrt ( ( \sum x^2 - (\sum x)^2 / N ) / ( N - 1 ) / N )
static float get_se_mean(const float *v, int n, const int *imap)
{
    int j = 0;
    for (int i = 0; i < n; i++)
    {
        float tmp = imap ? v[imap[i]] : v[i];
        if ( !isnan(tmp) ) j++;
    }
    if ( j <= 1 ) return NAN;

    return get_sample_sd(v, n, imap) / sqrtf( j );
}

// compute the median of a vector using the ksort library (with iterator)
float get_median(const float *v, int n, const int *imap)
{
    if ( n == 0 ) return NAN;
    float tmp, *w = (float *)malloc(n * sizeof(float));
    int j = 0;
    for (int i = 0; i < n; i++)
    {
        tmp = imap ? v[imap[i]] : v[i];
        if (!isnan(tmp)) w[j++] = tmp;
    }
    if ( j == 0 ) { free(w); return NAN; }
    float ret = ks_ksmall_float((size_t)j, w, (size_t)j/2);
    if (j%2==0) ret = (ret + w[j/2-1]) * 0.5f;
    free(w);
    return ret;
}

/*********************************
 * GENOME REFERENCE CONTIG RULES *
 *********************************/

// structure defining regions of interest in the genome
typedef struct
{
    int *length;
    int *cen_beg;
    int *cen_end;
    int x_rid;
    int x_nonpar_beg;
    int x_nonpar_end;
    int x_xtr_beg;
    int x_xtr_end;
    int y_rid;
    int y_nonpar_beg;
    int y_nonpar_end;
    int y_xtr_beg;
    int y_xtr_end;
    int mt_rid;
} genome_param_t;

// inspired by Petr Danecek's implementation of bcftools/vcfcall.c and bcftools/plugins/mendelian.c
typedef struct
{
    const char *alias, *about, *rules;
}
rules_predef_t;

// the following definitions were derived from these files:
// http://hgdownload.cse.ucsc.edu/goldenPath/hg18/database/chr{{1..22},{X,Y}}_gap.txt.gz
// http://hgdownload.cse.ucsc.edu/goldenPath/hg18/database/genomicSuperDups.txt.gz
// http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/gap.txt.gz
// http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/genomicSuperDups.txt.gz
// http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/centromeres.txt.gz
// http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/genomicSuperDups.txt.gz
static rules_predef_t rules_predefs[] =
{
    { .alias = "NCBI36",
      .about = "Human Genome reference assembly NCBI36 / hg18, both chr naming conventions",
      .rules =
            "   1:121236957-123476957     centromere\n"
            "   2:91689898-94689898       centromere\n"
            "   3:90587544-93487544       centromere\n"
            "   4:49354874-52354874       centromere\n"
            "   5:46441398-49441398       centromere\n"
            "   6:58938125-61938125       centromere\n"
            "   7:58058273-61058273       centromere\n"
            "   8:43958052-46958052       centromere\n"
            "   9:47107499-50107499       centromere\n"
            "   10:39244941-41624941      centromere\n"
            "   11:51450781-54450781      centromere\n"
            "   12:34747961-36142961      centromere\n"
            "   13:16000000-17868000      centromere\n"
            "   14:15070000-18070000      centromere\n"
            "   15:15260000-18260000      centromere\n"
            "   16:35143302-36943302      centromere\n"
            "   17:22187133-22287133      centromere\n"
            "   18:15400898-16764896      centromere\n"
            "   19:26923622-29923622      centromere\n"
            "   20:26267569-28033230      centromere\n"
            "   21:10260000-13260000      centromere\n"
            "   22:11330000-14330000      centromere\n"
            "   X:58598737-61598737       centromere\n"
            "   Y:11253954-12308578       centromere\n"
            "   X:2709520-154583483       X_nonpar\n"
            "   X:88343457-92262165       X_xtr\n"
            "   Y:2709520-57442674        Y_nonpar\n"
            "   Y:2977958-6676600         Y_xtr\n"
            "   MT:1-16571                mitochondria\n"
            "\n"
            "   chr1:121236957-123476957  centromere\n"
            "   chr2:91689898-94689898    centromere\n"
            "   chr3:90587544-93487544    centromere\n"
            "   chr4:49354874-52354874    centromere\n"
            "   chr5:46441398-49441398    centromere\n"
            "   chr6:58938125-61938125    centromere\n"
            "   chr7:58058273-61058273    centromere\n"
            "   chr8:43958052-46958052    centromere\n"
            "   chr9:47107499-50107499    centromere\n"
            "   chr10:39244941-41624941   centromere\n"
            "   chr11:51450781-54450781   centromere\n"
            "   chr12:34747961-36142961   centromere\n"
            "   chr13:16000000-17868000   centromere\n"
            "   chr14:15070000-18070000   centromere\n"
            "   chr15:15260000-18260000   centromere\n"
            "   chr16:35143302-36943302   centromere\n"
            "   chr17:22187133-22287133   centromere\n"
            "   chr18:15400898-16764896   centromere\n"
            "   chr19:26923622-29923622   centromere\n"
            "   chr20:26267569-28033230   centromere\n"
            "   chr21:10260000-13260000   centromere\n"
            "   chr22:11330000-14330000   centromere\n"
            "   chrX:58598737-61598737    centromere\n"
            "   chrY:11253954-12308578    centromere\n"
            "   chrX:2709520-154583483    X_nonpar\n"
            "   chrX:88343457-92262165    X_xtr\n"
            "   chrY:2709520-57442674     Y_nonpar\n"
            "   chrY:2977958-6676600      Y_xtr\n"
            "   chrM:1-16571              mitochondria\n"
    },
    { .alias = "GRCh37",
      .about = "Human Genome reference assembly GRCh37 / hg19, both chr naming conventions",
      .rules =
            "   1:121535434-124535434     centromere\n"
            "   2:92326171-95326171       centromere\n"
            "   3:90504854-93504854       centromere\n"
            "   4:49660117-52660117       centromere\n"
            "   5:46405641-49405641       centromere\n"
            "   6:58830166-61830166       centromere\n"
            "   7:58054331-61054331       centromere\n"
            "   8:43838887-46838887       centromere\n"
            "   9:47367679-50367679       centromere\n"
            "   10:39254935-42254935      centromere\n"
            "   11:51644205-54644205      centromere\n"
            "   12:34856694-37856694      centromere\n"
            "   13:16000000-19000000      centromere\n"
            "   14:16000000-19000000      centromere\n"
            "   15:17000000-20000000      centromere\n"
            "   16:35335801-38335801      centromere\n"
            "   17:22263006-25263006      centromere\n"
            "   18:15460898-18460898      centromere\n"
            "   19:24681782-27681782      centromere\n"
            "   20:26369569-29369569      centromere\n"
            "   21:11288129-14288129      centromere\n"
            "   22:13000000-16000000      centromere\n"
            "   X:58632012-61632012       centromere\n"
            "   Y:10104553-13104553       centromere\n"
            "   X:2699520-154930289       X_nonpar\n"
            "   X:88456801-92375509       X_xtr\n"
            "   Y:2649520-59033286        Y_nonpar\n"
            "   Y:2917958-6616600         Y_xtr\n"
            "   MT:1-16569                mitochondria\n"
            "\n"
            "   chr1:121535434-124535434  centromere\n"
            "   chr2:92326171-95326171    centromere\n"
            "   chr3:90504854-93504854    centromere\n"
            "   chr4:49660117-52660117    centromere\n"
            "   chr5:46405641-49405641    centromere\n"
            "   chr6:58830166-61830166    centromere\n"
            "   chr7:58054331-61054331    centromere\n"
            "   chr8:43838887-46838887    centromere\n"
            "   chr9:47367679-50367679    centromere\n"
            "   chr10:39254935-42254935   centromere\n"
            "   chr11:51644205-54644205   centromere\n"
            "   chr12:34856694-37856694   centromere\n"
            "   chr13:16000000-19000000   centromere\n"
            "   chr14:16000000-19000000   centromere\n"
            "   chr15:17000000-20000000   centromere\n"
            "   chr16:35335801-38335801   centromere\n"
            "   chr17:22263006-25263006   centromere\n"
            "   chr18:15460898-18460898   centromere\n"
            "   chr19:24681782-27681782   centromere\n"
            "   chr20:26369569-29369569   centromere\n"
            "   chr21:11288129-14288129   centromere\n"
            "   chr22:13000000-16000000   centromere\n"
            "   chrX:58632012-61632012    centromere\n"
            "   chrY:10104553-13104553    centromere\n"
            "   chrX:2699520-154930289    X_nonpar\n"
            "   chrX:88456801-92375509    X_xtr\n"
            "   chrY:2649520-59033286     Y_nonpar\n"
            "   chrY:2917958-6616600      Y_xtr\n"
            "   chrM:1-16571              mitochondria\n"
    },
    { .alias = "GRCh38",
      .about = "Human Genome reference assembly GRCh38 / hg38, both chr naming conventions",
      .rules =
            "   1:122026459-124932724     centromere\n"
            "   2:92188145-94090557       centromere\n"
            "   3:90772458-93655574       centromere\n"
            "   4:49712061-51743951       centromere\n"
            "   5:46485900-50059807       centromere\n"
            "   6:58553888-59829934       centromere\n"
            "   7:58169653-61528020       centromere\n"
            "   8:44033744-45877265       centromere\n"
            "   9:43389635-45518558       centromere\n"
            "   10:39686682-41593521      centromere\n"
            "   11:51078348-54425074      centromere\n"
            "   12:34769407-37185252      centromere\n"
            "   13:16000000-18051248      centromere\n"
            "   14:16000000-18173523      centromere\n"
            "   15:17083673-19725254      centromere\n"
            "   16:36311158-38265669      centromere\n"
            "   17:22813679-26616164      centromere\n"
            "   18:15460899-20861206      centromere\n"
            "   19:24498980-27190874      centromere\n"
            "   20:26436232-30038348      centromere\n"
            "   21:10864560-12915808      centromere\n"
            "   22:12954788-15054318      centromere\n"
            "   X:58605579-62412542       centromere\n"
            "   Y:10316944-10544039       centromere\n"
            "   X:2781479-155700628       X_nonpar\n"
            "   X:89201802-93120510       X_xtr\n"
            "   Y:2781479-56887139        Y_nonpar\n"
            "   Y:3049917-6748559         Y_xtr\n"
            "   MT:1-16569                mitochondria\n"
            "\n"
            "   chr1:122026459-124932724  centromere\n"
            "   chr2:92188145-94090557    centromere\n"
            "   chr3:90772458-93655574    centromere\n"
            "   chr4:49712061-51743951    centromere\n"
            "   chr5:46485900-50059807    centromere\n"
            "   chr6:58553888-59829934    centromere\n"
            "   chr7:58169653-61528020    centromere\n"
            "   chr8:44033744-45877265    centromere\n"
            "   chr9:43389635-45518558    centromere\n"
            "   chr10:39686682-41593521   centromere\n"
            "   chr11:51078348-54425074   centromere\n"
            "   chr12:34769407-37185252   centromere\n"
            "   chr13:16000000-18051248   centromere\n"
            "   chr14:16000000-18173523   centromere\n"
            "   chr15:17083673-19725254   centromere\n"
            "   chr16:36311158-38265669   centromere\n"
            "   chr17:22813679-26616164   centromere\n"
            "   chr18:15460899-20861206   centromere\n"
            "   chr19:24498980-27190874   centromere\n"
            "   chr20:26436232-30038348   centromere\n"
            "   chr21:10864560-12915808   centromere\n"
            "   chr22:12954788-15054318   centromere\n"
            "   chrX:58605579-62412542    centromere\n"
            "   chrY:10316944-10544039    centromere\n"
            "   chrX:2781479-155700628    X_nonpar\n"
            "   chrX:89201802-93120510    X_xtr\n"
            "   chrY:2781479-56887139     Y_nonpar\n"
            "   chrY:3049917-6748559      Y_xtr\n"
            "   chrM:1-16569              mitochondria\n"
    },
    { .alias = NULL,
      .about = NULL,
      .rules = NULL,
    }
};

// adapted from Petr Danecek's implementation of parse_rules() in bcftools/plugins/mendelian.c
static void push_rule(genome_param_t *genome_param, char *line, const bcf_hdr_t *hdr)
{
    // eat any leading spaces
    char *ss = (char *)line;
    while ( *ss && isspace(*ss) ) ss++;
    if ( !*ss ) return; // skip empty lines

    // chromosome name, beg, end
    char *tmp, *se = ss;
    while ( *se && !isspace(*se) && *se!=':' ) se++;
    if ( *se != ':' ) error("Could not parse the region: %s\n", line);
    *se = '\0'; // terminates the chromosome string
    int rid = bcf_hdr_name2id(hdr, ss);
    if (rid < 0) return;
    *se = ':'; // restores the separator
    ss = ++se;
    while ( *se && !isspace(*se) && *se!='-' ) se++;
    if ( *se != '-' ) error("Could not parse the region: %s\n", line);
    int beg = (int)strtol(ss, &tmp, 10);
    if ( tmp==ss ) error("Could not parse the region: %s\n", line);
    ss = ++se;
    int end = (int)strtol(ss, &tmp, 10);
    if ( tmp==ss || beg > end ) error("Could not parse the region: %s\n", line);

    // skip region
    while ( *ss && !isspace(*ss) ) ss++;
    while ( *ss && isspace(*ss) ) ss++;

    if ( strncmp(ss, "centromere", 10) == 0 )
    {
        if ( genome_param->cen_beg[rid] != 0 || genome_param->cen_end[rid] != 0 )
            error("Second centromere rule %s\n", line);
        genome_param->cen_beg[rid] = beg;
        genome_param->cen_end[rid] = end;
    }
    else if ( strncmp(ss, "X_nonpar", 8) == 0 )
    {
        if ( ( genome_param->x_xtr_beg != 0 || genome_param->x_xtr_end != 0 ) && genome_param->x_rid != rid )
            error("Chromosome X XTR and nonPAR regions declared on different contigs: %s\n", line);
        genome_param->x_rid = rid;
        if ( genome_param->x_nonpar_beg != 0 || genome_param->x_nonpar_end != 0 )
            error("Second chromosome X nonPAR rule: %s\n", line);
        genome_param->x_nonpar_beg = beg;
        genome_param->x_nonpar_end = end;
    }
    else if ( strncmp(ss, "X_xtr", 5) == 0 )
    {
        if ( ( genome_param->x_nonpar_beg != 0 || genome_param->x_nonpar_end != 0 ) && genome_param->x_rid != rid )
            error("Chromosome X nonPAR and XTR regions declared on different contigs: %s\n", line);
        genome_param->x_rid = rid;
        if ( genome_param->x_xtr_beg != 0 || genome_param->x_xtr_end != 0 )
            error("Second chromosome X XTR rule: %s\n", line);
        genome_param->x_xtr_beg = beg;
        genome_param->x_xtr_end = end;
    }
    else if ( strncmp(ss, "Y_nonpar", 8) == 0 )
    {
        if ( ( genome_param->y_xtr_beg != 0 || genome_param->y_xtr_end != 0 ) && genome_param->y_rid != rid )
            error("Chromosome Y XTR and nonPAR regions declared on different contigs: %s\n", line);
        genome_param->y_rid = rid;
        if ( genome_param->y_nonpar_beg != 0 || genome_param->y_nonpar_end != 0 )
            error("Second chromosome Y nonPAR rule: %s\n", line);
        genome_param->y_nonpar_beg = beg;
        genome_param->y_nonpar_end = end;
    }
    else if ( strncmp(ss, "Y_xtr", 5) == 0 )
    {
        if ( ( genome_param->y_nonpar_beg != 0 || genome_param->y_nonpar_end != 0 ) && genome_param->y_rid != rid )
            error("Chromosome Y nonPAR and XTR regions declared on different contigs: %s\n", line);
        genome_param->y_rid = rid;
        if ( genome_param->y_xtr_beg != 0 || genome_param->y_xtr_end != 0 )
            error("Second chromosome Y XTR rule: %s\n", line);
        genome_param->y_xtr_beg = beg;
        genome_param->y_xtr_end = end;
    }
    else if ( strncmp(ss, "mitochondria", 12) == 0 )
    {
        if ( genome_param->mt_rid != -1 )
            error("Second mitochondria rule %s\n", line);
        genome_param->mt_rid = rid;
    }
}

// split file into lines
// adapted from Petr Danecek's implementation of regidx_init_string() in bcftools/regidx.c
static void genome_init_file(genome_param_t *genome_param, const char *fname, const bcf_hdr_t *hdr)
{
    if ( !fname ) return;
    kstring_t tmp = {0, 0, NULL};
    htsFile *fp = hts_open(fname, "r");
    if ( !fp ) error("Failed to open %s: %s\n", fname, strerror(errno));
    while ( hts_getline(fp, KS_SEP_LINE, &tmp) > 0 )
        push_rule(genome_param, tmp.s, hdr);
    free(tmp.s);
    hts_close(fp);
}

// split string into lines
// adapted from Petr Danecek's implementation of regidx_init_string() in bcftools/regidx.c
static void genome_init_string(genome_param_t *genome_param, const char *str, const bcf_hdr_t *hdr)
{
    if ( !str ) return;
    kstring_t tmp = {0, 0, NULL};
    const char *ss = str;
    while ( *ss )
    {
        while ( *ss && isspace(*ss) ) ss++;
        const char *se = ss;
        while ( *se && *se!='\r' && *se!='\n' ) se++; // equivalent to KS_SEP_LINE
        tmp.l = 0;
        kputsn(ss, (int)(se-ss), &tmp);
        push_rule(genome_param, tmp.s, hdr);
        while ( *se && isspace(*se) ) se++;
        ss = se;
    }
    free(tmp.s);
}

// adapted from Petr Danecek's implementation of init_rules() in bcftools/regidx.c
static void genome_init_alias(genome_param_t *genome_param, char *alias, const bcf_hdr_t *hdr)
{
    const rules_predef_t *rules = rules_predefs;

    int detailed = 0, len = strlen(alias);
    if ( alias[len-1]=='?' ) { detailed = 1; alias[len-1] = '\0'; }

    while ( rules->alias && strcasecmp(alias, rules->alias) ) rules++;

    if ( !rules->alias )
    {
        fprintf(stderr, "\nPRE-DEFINED REFERENCE GENOME RULES\n");
        fprintf(stderr, "\b * Columns are: CHROM:BEG-END centromere/nonpar\n");
        fprintf(stderr, " * Coordinates are 1-based inclusive.\n");
        for ( rules = rules_predefs; rules->alias; rules++ )
        {
            fprintf(stderr, "\n%s\n   .. %s\n\n", rules->alias, rules->about);
            if ( detailed ) fprintf(stderr, "%s\n", rules->rules);
        }
        if ( !detailed )
        {
            fprintf(stderr, "\nRun as --rules <assembly> (e.g. --rules GRCh37).\n");
            fprintf(stderr, "To see the detailed rules definition, append a question mark (e.g. --rules GRCh37?).\n");
            fprintf(stderr, "\n");
        }
        exit(1);
    }
    else if ( detailed )
    {
        fprintf(stderr, "\n%s\n   .. %s\n\n", rules->alias, rules->about);
        fprintf(stderr, "%s", rules->rules);
        fprintf(stderr, "\n");
        exit(1);
    }
    return genome_init_string(genome_param, rules->rules, hdr);
}

int cnp_parse(const char *line, char **chr_beg, char **chr_end, uint32_t *beg, uint32_t *end, void *payload, void *usr)
{
    // Use the standard parser for CHROM,FROM,TO
    int ret = regidx_parse_bed(line, chr_beg, chr_end, beg, end, NULL, NULL);
    if ( ret!=0 ) return ret;

    // Skip the fields that were parsed above
    char *ss = (char*) line;
    while ( *ss && isspace(*ss) ) ss++;
    for (int i=0; i<3; i++)
    {
        while ( *ss && !isspace(*ss) ) ss++;
        if ( !*ss ) return -2;  // wrong number of fields
        while ( *ss && isspace(*ss) ) ss++;
    }
    if ( !*ss ) return -2;

    // Parse the payload
    int *dat = (int*)payload;
    if ( strncmp(ss, "DEL", 3) == 0 ) *dat = MOCHA_CNP_DEL;
    else if ( strncmp(ss, "DUP", 3) == 0 ) *dat = MOCHA_CNP_DUP;
    else if ( strncmp(ss, "CNV", 3) == 0 ) *dat = MOCHA_CNP_CNV;
    else *dat = MOCHA_UNK;
    return 0;
}

/*********************************
 * SAMPLE STRUCTURES AND METHODS *
 *********************************/

typedef struct
{
    float xy_prob;
    float err_prob;
    float flip_prob;
    float centromere_prob;
    float telomere_prob;
    float *cnf, *bdev;
    int cnf_n, bdev_n;
    int min_dist;
    float lrr_hap2dip;
    float lrr_auto2sex;
    float lrr_bias;
    int median_baf_adjust;
    int order_lrr_gc;
} model_param_t;

typedef struct
{
    int sample_idx;
    int sex;
    int rid;
    int beg_pos;
    int end_pos;
    int length;
    int8_t p_arm;
    int8_t q_arm;
    int nsites;
    int nhets;
    int n50_hets;
    float cnf;
    float bdev;
    float bdev_se;
    float ldev;
    float ldev_se;
    float lod_lrr_baf;
    float lod_baf_phase;
    int nflips;
    float baf_conc;
    int8_t type;
    float cf;
} mocha_t;

typedef kvec_t(mocha_t) kv_mocha_t;

typedef struct
{
    int rid;
    int n;
    int *pos;
    float *gc;
} contig_param_t;

// this structure will contain at most 9*size bytes as data
// will be freed when bdev and bdev_phase get populated
typedef struct
{
    int size;
    int *vcf_imap;
    int16_t *data[2];
    int8_t *gt_phase;
    int16_t *ldev;
    int16_t *bdev;
    int8_t *bdev_phase;
} contig_t;

typedef struct
{
    int idx;
    int sex;
    float baf_sd;
    float lrr_median;
    float lrr_sd;
    float adjlrr_sd;
    float baf_conc;
    float lrr_auto;
    int nhets;
    int x_nonpar_nhets;
    float x_nonpar_lrr_median;
    float y_nonpar_lrr_median;
    float mt_lrr_median;
    contig_t contig;

    kv_float kv_lrr_median;
    kv_float kv_lrr_sd;
    kv_float kv_baf_sd;
    kv_float kv_baf_conc;
    kv_float kv_lrr_auto;

    // LRR polynomial regression parameters
    float coeffs[MAX_ORDER+1];
    kv_float kv_coeffs;
    float rel_ess;
    kv_float kv_rel_ess;
} sample_t;

static void mocha_print_ucsc(const mocha_t *mocha, int n, FILE *restrict stream, const bcf_hdr_t *hdr)
{
    if (stream == NULL) return;
    const char *name[4];
    name[MOCHA_UNK] = "mCA_undetermined";
    name[MOCHA_DEL] = "mCA_loss";
    name[MOCHA_DUP] = "mCA_gain";
    name[MOCHA_UPD] = "mCA_neutral";
    const char *desc[4];
    desc[MOCHA_UNK] = "Undetermined";
    desc[MOCHA_DEL] = "Deletions";
    desc[MOCHA_DUP] = "Duplications";
    desc[MOCHA_UPD] = "CNN-LOHs";
    const char *color[4];
    color[MOCHA_UNK] = "127,127,127";
    color[MOCHA_DEL] = "255,0,0";
    color[MOCHA_DUP] = "0,0,255";
    color[MOCHA_UPD] = "0,255,0";
    for (int i=0; i<4; i++)
    {
        fprintf(stream, "track name=%s description=\"%s\" visibility=4 priority=1 itemRgb=\"On\"\n", name[i], desc[i]);
        for (int j=0; j<n; j++)
        {
            if (i == mocha[j].type)
            {
                const char *sample_name = bcf_hdr_int2id(hdr, BCF_DT_SAMPLE, mocha[j].sample_idx);
                const char *seq_name = bcf_hdr_id2name(hdr, mocha[j].rid);
                if ( strncmp(seq_name, "chr", 3) == 0 ) fprintf(stream, "%s\t%d\t%d\t%s\t0\t.\t%d\t%d\t%s\n", seq_name,
                    mocha[j].beg_pos, mocha[j].end_pos, sample_name, mocha[j].beg_pos, mocha[j].end_pos, color[i]);
                else fprintf(stream, "chr%s\t%d\t%d\t%s\t0\t.\t%d\t%d\t%s\n", seq_name, mocha[j].beg_pos,
                    mocha[j].end_pos, sample_name, mocha[j].beg_pos, mocha[j].end_pos, color[i]);
            }
        }
    }
    if (stream != stdout && stream != stderr) fclose(stream);
}

static void mocha_print(const mocha_t *mocha, int n, FILE *restrict stream, const bcf_hdr_t *hdr, int flags)
{
    if (stream == NULL) return;
    char sex[3];
    sex[SEX_UNK] = 'U';
    sex[SEX_MAL] = 'M';
    sex[SEX_FEM] = 'F';
    const char *type[6];
    type[MOCHA_UNK] = "Undetermined";
    type[MOCHA_DEL] = "Deletion";
    type[MOCHA_DUP] = "Duplication";
    type[MOCHA_UPD] = "CNN-LOH";
    type[MOCHA_CNP_DEL] = "CNP Deletion";
    type[MOCHA_CNP_DUP] = "CNP Duplication";
    char arm_type[3];
    arm_type[MOCHA_NAN] = 'N';
    arm_type[MOCHA_ARM] = 'Y';
    arm_type[MOCHA_TEL] = 'T';
    if ( flags & WGS_DATA )
    {
        fprintf(stream, "SAMPLE\tSEX\tCHROM\tBEG\tEND\tLENGTH\tP\tQ\tNSITES\tHETS\tN50_HETS\tCNF\tBDEV\tBDEV_SE\tREL_COV\tREL_COV_SE\tLOD_LRR_BAF\tLOD_BAF_PHASE\tFLIPS\tBAF_CONC\tTYPE\tCF\n");
        for (int i=0; i<n; i++)
        {
            const char *sample_name = bcf_hdr_int2id(hdr, BCF_DT_SAMPLE, mocha->sample_idx);
            const char *seq_name = bcf_hdr_id2name(hdr, mocha->rid);
            fprintf(stream, "%s\t%c\t%s\t%d\t%d\t%d\t%c\t%c\t%d\t%d\t%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.2f\t%d\t%.4f\t%s\t%.4f\n",
                sample_name, sex[mocha->sex], seq_name, mocha->beg_pos, mocha->end_pos, mocha->length,
                arm_type[mocha->p_arm], arm_type[mocha->q_arm], mocha->nsites, mocha->nhets, mocha->n50_hets, mocha->cnf,
                mocha->bdev, mocha->bdev_se, 2.0f * expf(mocha->ldev), 2.0f * expf(mocha->ldev) * mocha->ldev_se,
                mocha->lod_lrr_baf, mocha->lod_baf_phase, mocha->nflips, mocha->baf_conc, type[mocha->type], mocha->cf);
            mocha++;
        }
    }
    else
    {
        fprintf(stream, "SAMPLE\tSEX\tCHROM\tBEG\tEND\tLENGTH\tP\tQ\tNSITES\tHETS\tN50_HETS\tCNF\tBDEV\tBDEV_SE\tLDEV\tLDEV_SE\tLOD_LRR_BAF\tLOD_BAF_PHASE\tFLIPS\tBAF_CONC\tTYPE\tCF\n");
        for (int i=0; i<n; i++)
        {
            const char *sample_name = bcf_hdr_int2id(hdr, BCF_DT_SAMPLE, mocha->sample_idx);
            const char *seq_name = bcf_hdr_id2name(hdr, mocha->rid);
            fprintf(stream, "%s\t%c\t%s\t%d\t%d\t%d\t%c\t%c\t%d\t%d\t%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.2f\t%d\t%.4f\t%s\t%.4f\n",
                sample_name, sex[mocha->sex], seq_name, mocha->beg_pos, mocha->end_pos, mocha->length,
                arm_type[mocha->p_arm], arm_type[mocha->q_arm], mocha->nsites, mocha->nhets, mocha->n50_hets, mocha->cnf,
                mocha->bdev, mocha->bdev_se, mocha->ldev, mocha->ldev_se, mocha->lod_lrr_baf, mocha->lod_baf_phase,
                mocha->nflips, mocha->baf_conc, type[mocha->type], mocha->cf);
            mocha++;
        }
    }
    if (stream != stdout && stream != stderr) fclose(stream);
}

// this function checks whether the edge of a CNP is indeed diploid
// by looking for evidence for the diploid model over the CNP model (duplication or deletion)
static int cnp_edge_is_not_cn2(const float *lrr, const float *baf, int n, int a, int b, float xy_prob, float err_prob, float lrr_bias, float lrr_hap2dip, float lrr_sd, float baf_sd, double cnf)
{
    float lrr_sd_bias = powf( lrr_sd, lrr_bias );
    float ldev = ( logf((float)cnf) / (float)M_LN2 - 1.0f ) * lrr_hap2dip;
    float bdev = fabsf( 0.5f - 1.0f / (float)cnf );

    // test left edge
    float lkl = 1.0f;
    for (int i=a-1; i>=0; i--)
    {
        lkl *= lkl_lrr_baf( lrr[i], baf[i], ldev, bdev, lrr_sd, baf_sd, lrr_bias, lrr_sd_bias, err_prob );
        lkl /= lkl_lrr_baf( lrr[i], baf[i], 0, 0, lrr_sd, baf_sd, lrr_bias, lrr_sd_bias, err_prob );
        if ( lkl > 1/xy_prob ) return -1;
        if ( lkl < xy_prob ) break;
    }

    // test right edge
    lkl = 1.0f;
    for (int i=b+1; i<n; i++)
    {
        lkl *= lkl_lrr_baf( lrr[i], baf[i], ldev, bdev, lrr_sd, baf_sd, lrr_bias, lrr_sd_bias, err_prob );
        lkl /= lkl_lrr_baf( lrr[i], baf[i], 0, 0, lrr_sd, baf_sd, lrr_bias, lrr_sd_bias, err_prob );
        if ( lkl > 1/xy_prob ) return -1;
        if ( lkl < xy_prob ) break;
    }

    return 0;
}

// this function returns two values (a, b) such that:
// (i) a <= b (ii) pos[a-1] < beg (iii) beg <= pos[a] < end (iv) beg <= pos[b] < end (v) pos[b+1] >= end
static int get_cnp_edges(const int *pos, int n, int beg, int end, int *a, int *b)
{
    if ( pos[ 0 ] >= end || pos[ n-1 ] < beg ) return -1;

    int i = 0, j = n-1, k;
    while( j - i > 1 )
    {
        k = ( i + j ) / 2;
        if ( pos[ k ] < beg ) i = k; else j = k;
    }
    if( pos[ j ] >= end ) return -1;

    *a = j;

    i = j;
    j = n-1;
    while( j - i > 1 )
    {
        k = ( i + j ) / 2;
        if ( pos[ k ] < end ) i = k; else j = k;
    }

    *b = i;
    return 0;
}

// classify mosaic chromosomal alteration type based on LRR and BAF
static int8_t mocha_type(float ldev, float ldev_se, float bdev, float bdev_se, float lrr_hap2dip, int8_t p_arm, int8_t q_arm)
{
    float z2_upd = sqf( ldev / ldev_se ) * 0.5f;
    if (p_arm != MOCHA_TEL && q_arm != MOCHA_TEL) z2_upd += 2.0f; // penalty for not ending in a telomere
    if (ldev > 0)
    {
        if ( z2_upd > 5.0 ) return MOCHA_DUP;
        if ( !isnan(bdev) )
        {
            float z2_dup = sqf( (ldev - 2.0f * M_LOG2E * lrr_hap2dip * bdev ) / ldev_se ) * 0.5f;
            if (z2_upd > z2_dup + 3.0f) return MOCHA_DUP;
            if (z2_dup > z2_upd + 3.0f) return MOCHA_UPD;
        }
    }
    else
    {
        if ( z2_upd > 5.0 ) return MOCHA_DEL;
        if ( !isnan(bdev) )
        {
            float z2_del = sqf( (ldev + 2.0f * M_LOG2E * lrr_hap2dip * bdev) / ldev_se ) * 0.5f;
            if (z2_upd > z2_del + 3.0f) return MOCHA_DEL;
            if (z2_del > z2_upd + 3.0f) return MOCHA_UPD;
        }
    }
    return MOCHA_UNK;
}

// best estimate for cell fraction using the following formula (if BDEV is available):
// BDEV = | 1 / 2 - 1 / CNF |
// CNF = 2 / ( 1 + 2 x BDEV ) for deletions
// CNF = 2 / ( 1 - 2 x BDEV ) for duplications
static float mocha_cellfraction(float ldev, float ldev_se, float bdev, float bdev_se, int8_t type, float lrr_hap2dip)
{
    if (isnan(bdev))
    {
        switch (type)
        {
            case MOCHA_DEL:
                return -ldev / lrr_hap2dip;
            case MOCHA_DUP:
                return ldev / lrr_hap2dip;
            default:
                return NAN;
        }
    }
    else
    {
        switch (type)
        {
            case MOCHA_UNK:
                return 4.0f * bdev; // here it assumes it is either a deletion or a duplication
            case MOCHA_DEL:
                return 4.0f * bdev / (1.0f + 2.0f * bdev);
            case MOCHA_DUP:
                return 4.0f * bdev / (1.0f - 2.0f * bdev);
            case MOCHA_UPD:
                return 2.0f * bdev;
            default:
                return NAN;
        }
    }
}

static void get_mocha_stats(const int *pos, const float *lrr, const float *baf, const int8_t *gt_phase, int n, int a, int b, int cen_beg, int cen_end, int length, mocha_t *mocha)
{
    mocha->nsites = b + 1 - a;

    if ( a == 0 )
        if ( pos[ a ] < cen_beg )
            mocha->beg_pos = 0;
        else
            mocha->beg_pos = cen_end;
    else mocha->beg_pos = pos[ a ];

    if ( b == n-1 )
        if ( pos[ b ] >= cen_end )
            mocha->end_pos = length;
        else
            mocha->end_pos = cen_beg;
    else mocha->end_pos = pos[ b ];

    mocha->length = mocha->end_pos - mocha->beg_pos;

    if ( mocha->beg_pos == 0 ) mocha->p_arm = MOCHA_TEL;
    else if ( mocha->beg_pos < cen_beg ) mocha->p_arm = MOCHA_ARM;
    else mocha->p_arm = MOCHA_NAN;

    if ( mocha->end_pos == length ) mocha->q_arm = MOCHA_TEL;
    else if ( mocha->end_pos > cen_end ) mocha->q_arm = MOCHA_ARM;
    else mocha->q_arm = MOCHA_NAN;

    mocha->ldev_se = get_se_mean( lrr + a, b + 1 - a, NULL );
    mocha->nhets = 0;
    for ( int i=a; i<=b; i++ ) if ( !isnan( baf[i] ) ) mocha->nhets++;
    mocha->baf_conc = get_baf_conc( baf + a, gt_phase + a, b + 1 - a, NULL );
    mocha->n50_hets = get_n50( pos + a, b + 1 - a, NULL );
    mocha->nflips = -1;
    mocha->bdev = NAN;
    mocha->bdev_se = NAN;
    mocha->lod_baf_phase = NAN;
}

// return segments called by the HMM or state with consecutive call
static int get_path_segs(const int8_t *path, int n, int except, int **beg, int **end, int *nseg)
{
    int beg_m = 0, end_m = 0, a = 0, b = 0;
    *beg = NULL;
    *end = NULL;
    *nseg = 0;
    for (b=0; b<n; b++)
    {
        // check whether it is the end of a segment
        if ( b != n-1 )
        {
            int curr = abs(path[b]);
            int next = abs(path[b+1]);
            if ( curr == next ) continue;

            // if two consecutive segments have consecutive HMM states
            if ( abs( curr - next ) == 1 )
            {
                // if the consecutive HMM states are non-zero and don't correspond to consecutive deletions and duplications
                if ( curr && next && ( ( curr != except && curr != except+1 ) || ( next != except && next != except+1 ) ) )
                {
                    free(*beg);
                    free(*end);
                    return curr < next ? curr : next;
                }
            }
        }

        if ( path[b] )
        {
            (*nseg)++;
            hts_expand(int, *nseg, beg_m, *beg);
            (*beg)[(*nseg)-1] = a;
            hts_expand(int, *nseg, end_m, *end);
            (*end)[(*nseg)-1] = b;
        }
        a = b + 1;
    }
    return 0;
}

// process one contig for one sample
static void sample_contig_run(sample_t *self, kv_mocha_t *kv_mocha, const contig_param_t *contig_param, const genome_param_t *genome_param, const model_param_t *model_param, regitr_t *cnp_itr, int flags)
{
    // do nothing if chromosome Y or MT are being tested
    if ( contig_param->rid == genome_param->y_rid || contig_param->rid == genome_param->mt_rid )
    {
        free(self->contig.data[BAF]);
        free(self->contig.data[LRR]);
        free(self->contig.gt_phase);
        self->contig.ldev = NULL;
        self->contig.bdev = NULL;
        self->contig.bdev_phase = NULL;
        return;
    }

    mocha_t mocha;
    mocha.sample_idx = self->idx;
    mocha.sex = self->sex;
    mocha.rid = contig_param->rid;

    int cen_beg = genome_param->cen_beg[contig_param->rid];
    int cen_end = genome_param->cen_end[contig_param->rid];
    int length = genome_param->length[contig_param->rid];

    // declutter code by copying these values onto the stack
    int size = self->contig.size;
    int8_t *gt_phase = self->contig.gt_phase;
    float *lrr, *baf;
    int16_t *ad0, *ad1;
    if ( flags & WGS_DATA )
    {
        ad0 = self->contig.data[AD0];
        ad1 = self->contig.data[AD1];
        lrr = ad_to_lrr(ad0, ad1, self->contig.size);
        baf = ad_to_baf(ad0, ad1, self->contig.size);
    }
    else
    {
        ad0 = NULL;
        ad1 = NULL;
        lrr = realloc_int16_to_float(self->contig.data[LRR], self->contig.size);
        baf = realloc_int16_to_float(self->contig.data[BAF], self->contig.size);
    }

    if ( model_param->order_lrr_gc >= 0 )
        adjust_lrr(lrr, contig_param->gc, self->contig.size, self->contig.vcf_imap, self->coeffs, model_param->order_lrr_gc);
    int16_t *ldev = (int16_t *)calloc(self->contig.size, sizeof(int16_t));
    int16_t *bdev = (int16_t *)calloc(self->contig.size, sizeof(int16_t));
    int8_t *bdev_phase = (int8_t *)calloc(self->contig.size, sizeof(int8_t));
    int *pos = (int *)malloc(self->contig.size * sizeof(int));
    for (int i=0; i<self->contig.size; i++) pos[i] = contig_param->pos[ self->contig.vcf_imap [i] ];

    // TODO do I need special normalization for the X nonPAR region
    if ( contig_param->rid == genome_param->x_rid )
    {
        for (int i=0; i<self->contig.size; i++)
        {
            if ( pos[i] > genome_param->x_nonpar_beg && pos[i] < genome_param->x_nonpar_end )
            {
                lrr[i] = ( self->sex == SEX_MAL ) ? NAN : lrr[i] - model_param->lrr_auto2sex;
                baf[i] = ( self->sex == SEX_MAL ) ? NAN : baf[i];
            }
        }
    }

    if ( cnp_itr )
    {
        while ( regitr_overlap( cnp_itr ) )
        {
            int a, b;
            if ( get_cnp_edges( pos, size, cnp_itr->beg, cnp_itr->end, &a, &b) == 0 )
            {
                float cn2_lod = log10_lkl_lrr_baf(lrr + a, baf + a, b + 1 - a, NULL, model_param->err_prob, model_param->lrr_bias, model_param->lrr_hap2dip, self->adjlrr_sd, self->baf_sd, 2.0f);
                int cnp_type = regitr_payload(cnp_itr, int);
                mocha.type = MOCHA_UNK;
                mocha.ldev = get_median( lrr + a, b + 1 - a, NULL );
                if ( mocha.ldev > 0 && ( cnp_type == MOCHA_CNP_DUP || cnp_type == MOCHA_CNP_CNV ) )
                {
                    mocha.lod_lrr_baf = log10_lkl_lrr_baf(lrr + a, baf + a, b + 1 - a, NULL, model_param->err_prob, model_param->lrr_bias, model_param->lrr_hap2dip, self->adjlrr_sd, self->baf_sd, 3.0f) - cn2_lod;
                    if ( mocha.lod_lrr_baf > -log10f( model_param->xy_prob ) )
                    {
                        mocha.cnf = 3.0f;
                        mocha.type = MOCHA_CNP_DUP;
                        mocha.cf = NAN;
                    }
                }
                else if ( mocha.ldev <= 0 && ( cnp_type == MOCHA_CNP_DEL || cnp_type == MOCHA_CNP_CNV ) )
                {
                    mocha.lod_lrr_baf = log10_lkl_lrr_baf(lrr + a, baf + a, b + 1 - a, NULL, model_param->err_prob, model_param->lrr_bias, model_param->lrr_hap2dip, self->adjlrr_sd, self->baf_sd, 1.0f) - cn2_lod;
                    if ( mocha.lod_lrr_baf > -log10f( model_param->xy_prob ) )
                    {
                        mocha.cnf = 1.0f;
                        mocha.type = MOCHA_CNP_DEL;
                        mocha.cf = NAN;
                    }
                }
                if ( mocha.type == MOCHA_CNP_DUP || mocha.type == MOCHA_CNP_DEL )
                {
                    if ( cnp_edge_is_not_cn2(lrr, baf, size, a, b, model_param->xy_prob, model_param->err_prob, model_param->lrr_bias, model_param->lrr_hap2dip, self->adjlrr_sd, self->baf_sd, mocha.cnf) ) continue;
                    get_mocha_stats( pos, lrr, baf, gt_phase, size, a, b, cen_beg, cen_end, length, &mocha);
                    // compute bdev, if possible
                    if ( mocha.nhets > 0 )
                    {
                        double f(double x, void *data) { return -log10_lkl_baf_phase(baf + a, gt_phase + a, b + 1 - a, NULL, NULL, model_param->err_prob, NAN, self->baf_sd, x); }
                        double x;
                        kmin_brent(f, 0.0, 0.5, NULL, KMIN_EPS, &x);
                        mocha.bdev = fabsf((float)x);
                    }
                    else mocha.bdev = NAN;
                    kv_push(mocha_t, *kv_mocha, mocha);
                    for (int j=a; j<=b; j++)
                    {
                        // TODO add other stuff here, like setting ldev and bdev
                        lrr[j] = NAN; // do not use the data again
                        baf[j] = NAN; // do not use the data again
                    }
                }
            }
        }
    }

    kv_int imap = {0, 0, NULL};
    for (int model=0; model<2; model++)
    {
        // select data to use from the contig, depending on which model is being used
        int last_p = 0, first_q = 0;
        imap.n = 0;
        for (int i=0; i<size; i++)
            if ( ( model == LRR_BAF && !isnan( lrr[i] ) ) || ( model == BAF_PHASE && !isnan( baf[i] ) ) )
            {
                if ( pos[ i ] < cen_beg ) last_p++;
                if ( pos[ i ] < cen_end ) first_q++;
                kv_push(int, imap, i);
            }
        if ( imap.n == 0 ) continue;

        // compute emission probabilities and Viterbi path according to model
        int except = 0;
        kv_float hs = {0, 0, NULL};
        if ( model == LRR_BAF ) for (int i=0; i<model_param->cnf_n; i++) { kv_push(float, hs, model_param->cnf[i]); if ( model_param->cnf[i] < 2.0f ) except++; }
        else if ( model == BAF_PHASE ) for (int i=0; i<model_param->bdev_n; i++) kv_push(float, hs, model_param->bdev[i]);
        int8_t *path;
        int ret, *beg, *end, nseg;
        do
        {
            float *emis;
            if ( flags & WGS_DATA )
            {
                emis = model==LRR_BAF ? get_lrr_ad_emis( lrr, ad0, ad1, imap.n, imap.a, model_param->err_prob, model_param->lrr_bias, model_param->lrr_hap2dip, self->adjlrr_sd, self->baf_sd, hs.a, hs.n )
                                      : get_ad_phase_emis( ad0, ad1, gt_phase, imap.n, imap.a, model_param->err_prob, model_param->flip_prob, self->baf_sd, hs.a, hs.n );
            }
            else
            {
                emis = model==LRR_BAF ? get_lrr_baf_emis( lrr, baf, imap.n, imap.a, model_param->err_prob, model_param->lrr_bias, model_param->lrr_hap2dip, self->adjlrr_sd, self->baf_sd, hs.a, hs.n )
                                      : get_baf_phase_emis( baf, gt_phase, imap.n, imap.a, model_param->err_prob, model_param->flip_prob, self->baf_sd, hs.a, hs.n );
            }
            path = model==LRR_BAF ? get_viterbi(emis, imap.n, hs.n, model_param->xy_prob, NAN, model_param->telomere_prob, model_param->centromere_prob, last_p, first_q)
                                  : get_viterbi(emis, imap.n, hs.n, model_param->xy_prob, model_param->flip_prob, model_param->telomere_prob, model_param->centromere_prob, last_p, first_q);
            free(emis);

            ret = get_path_segs(path, imap.n, except, &beg, &end, &nseg);

            if ( ret ) // two consecutive hidden states were used, hinting that testing of a middle state might be necessary
            {
                free(path);
                kv_push(float, hs, (hs.a[ret-1] + hs.a[ret]) * 0.5f);
                ks_introsort_float(hs.n, hs.a);
            }
        }
        while ( ret );
        free(hs.a);

        // loop through all the segments called by the Viterbi algorithm
        for (int i=0; i<nseg; i++)
        {
            // compute edges of the call
            int a = imap.a[ beg[i] ];
            if ( beg[i] == 0 ) while ( a>0 && ldev[a-1]==0 && bdev[a-1]==0 ) a--; // extend call towards p telomere
            int b = imap.a[ end[i] ];
            if ( end[i] == imap.n-1 ) while ( b<size-1 && ldev[b+1]==0 && bdev[b+1]==0 ) b++; // extend call towards q telomere
            mocha.ldev = get_median( lrr + a, b + 1 - a, NULL );
            get_mocha_stats( pos, lrr, baf, gt_phase, size, a, b, cen_beg, cen_end, length, &mocha);

            double f(double x, void *data)
            {
                if ( flags & WGS_DATA ) return -log10_lkl_lrr_ad(lrr + a, ad0 + a, ad1 + a, mocha.nsites, NULL, model_param->err_prob, model_param->lrr_bias, model_param->lrr_hap2dip, self->adjlrr_sd, self->baf_sd, x);
                else return -log10_lkl_lrr_baf(lrr + a, baf + a, mocha.nsites, NULL, model_param->err_prob, model_param->lrr_bias, model_param->lrr_hap2dip, self->adjlrr_sd, self->baf_sd, x);
            }
            double x, fx = kmin_brent(f, ( mocha.ldev > 0 ) ? 2.0 : 0, ( mocha.ldev > 0 ) ? 4.0 : 2.0, NULL, KMIN_EPS, &x);
            mocha.cnf = x;
            mocha.lod_lrr_baf = (float)(f(2.0, NULL) - fx);

            if ( model == LRR_BAF )
            {
                // here you need to check whether the call would have been better with the phased model
                kv_int hets_imap = {0, 0, NULL};
                for (int j=beg[i]; j<=end[i]; j++)
                    if ( !isnan( baf[ imap.a[j] ] ) )
                        kv_push(int, hets_imap, imap.a[j]);
                // TODO here it needs to pass information about the centromeres
                if ( flags & WGS_DATA )
                {
                    mocha.lod_baf_phase = compare_wgs_models(ad0, ad1, gt_phase, hets_imap.n, hets_imap.a, model_param->xy_prob,
                        model_param->err_prob, model_param->flip_prob, model_param->telomere_prob, self->baf_sd, model_param->bdev, model_param->bdev_n);
                }
                else
                {
                    mocha.lod_baf_phase = compare_models(baf, gt_phase, hets_imap.n, hets_imap.a, model_param->xy_prob,
                        model_param->err_prob, model_param->flip_prob, model_param->telomere_prob, self->baf_sd, model_param->bdev, model_param->bdev_n);
                }
                if (mocha.lod_baf_phase > mocha.lod_lrr_baf) { free(hets_imap.a); continue; }

                // compute bdev, if possible
                if ( hets_imap.n > 0 )
                {
                    double f(double x, void *data)
                    {
                        if ( flags & WGS_DATA ) return -log10_lkl_ad_phase(ad0, ad1, gt_phase, hets_imap.n, hets_imap.a, NULL, model_param->err_prob, NAN, self->baf_sd, x);
                        else return -log10_lkl_baf_phase(baf, gt_phase, hets_imap.n, hets_imap.a, NULL, model_param->err_prob, NAN, self->baf_sd, x);
                    }
                    fx = kmin_brent(f, 0.0, 0.5, NULL, KMIN_EPS, &x); // since f is symmetric around 0, it is important not to use an interval symmetric around 0
                    mocha.bdev = fabsf((float)x);
                }
                else mocha.bdev = NAN;
                mocha.bdev_se = NAN;
                for (int j=0; j<hets_imap.n; j++) bdev_phase[ hets_imap.a[j] ] = (int8_t)SIGN( baf[ hets_imap.a[j] ] - 0.5f );
                free(hets_imap.a);
            }
            else
            {
                // penalizes the LOD by the number of phase flips
                mocha.nflips = 0;
                for (int j=beg[i]; j<end[i]; j++) if ( path[j] != path[j+1] ) mocha.nflips++;

                double f(double x, void *data)
                {
                    if ( flags & WGS_DATA ) return -log10_lkl_ad_phase(ad0, ad1, gt_phase, mocha.nhets, imap.a + beg[i], path + beg[i], model_param->err_prob, model_param->flip_prob, self->baf_sd, x);
                    else return -log10_lkl_baf_phase(baf, gt_phase, mocha.nhets, imap.a + beg[i], path + beg[i], model_param->err_prob, model_param->flip_prob, self->baf_sd, x);
                }
                double x, fx = kmin_brent(f, -0.5f, 0.5f, NULL, KMIN_EPS, &x);
                mocha.bdev = fabsf((float)x);
                mocha.lod_baf_phase = (float)(f(0.0, NULL) - fx) + (float)mocha.nflips * log10f(model_param->flip_prob);

                float *pbaf = (float *)malloc(mocha.nhets * sizeof(float));
                for (int j=0; j<mocha.nhets; j++) pbaf[ j ] = ( baf[ imap.a[beg[i]+j] ] - 0.5f ) * (float)SIGN( path[beg[i]+j] );
                mocha.bdev_se = get_se_mean( pbaf, mocha.nhets, NULL );
                free(pbaf);
                for (int j=beg[i]; j<=end[i]; j++) bdev_phase[ imap.a[j] ] = (int8_t)SIGN( path[j] ) * gt_phase[ imap.a[j] ];
            }

            mocha.type = mocha_type(mocha.ldev, mocha.ldev_se, mocha.bdev, mocha.bdev_se, model_param->lrr_hap2dip, mocha.p_arm, mocha.q_arm);
            mocha.cf = mocha_cellfraction(mocha.ldev, mocha.ldev_se, mocha.bdev, mocha.bdev_se, mocha.type, model_param->lrr_hap2dip);
            kv_push(mocha_t, *kv_mocha, mocha);

            // update information that will be stored in the output VCF and make remaining sites NAN
            for (int j=a; j<=b; j++)
            {
                if ( ldev[j] == 0 ) ldev[j] = float_to_int16( mocha.ldev );
                if ( bdev[j] == 0 ) bdev[j] = float_to_int16( mocha.bdev );
                lrr[j] = NAN; // do not use the data again
                baf[j] = NAN; // do not use the data again
            }
        }
        free(path);
        free(beg);
        free(end);
    }
    free(imap.a);

    free(pos);
    free(ad0);
    free(ad1);
    free(lrr);
    free(baf);
    free(self->contig.gt_phase);
    self->contig.ldev = ldev;
    self->contig.bdev = bdev;
    self->contig.bdev_phase = bdev_phase;
}

// computes the medoid contig for LRR regression
// TODO weight the coefficients appropriately
static int get_medoid(const float *coeffs, int n, int order)
{
    n /= order + 1;
    int medoid_idx = -1;
    float prev = INFINITY;
    for (int i=0; i<n; i++)
    {
        float next = 0.0f;
        for (int j=0; j<n; j++)
        {
            if (i != j)
            {
                for (int k=0; k<=order; k++)
                {
                    next += fabsf( coeffs[i * (order+1) + k] - coeffs[j * (order+1) + k] );
                }
            }
        }
        if (next < prev)
        {
            prev = next;
            medoid_idx = i;
        }
    }
    return medoid_idx;
}

// groups numbers in two separate distributions
static float get_cutoff(const float *v, int n)
{
    if ( n <= 1 ) return NAN;
    float *w = (float *)malloc(n * sizeof(float));
    int j = 0;
    for (int i=0; i<n; i++)
    {
        if ( !isnan(v[i]) ) w[j++] = v[i];
    }
    if ( j <= 1 ) { free(w); return NAN; }
    ks_introsort_float((size_t)j, w);

    // run k-means clustering EM
    int k = j/2;
    while (k > 0 && w[k-1] - (w[(k-1)/2] + w[k/2])*0.5f > (w[(j+k-1)/2] + w[(j+k)/2])*0.5f - w[k-1]) k--;
    while (k < j && w[k  ] - (w[(k-1)/2] + w[k/2])*0.5f < (w[(j+k-1)/2] + w[(j+k)/2])*0.5f - w[k  ]) k++;
    float cutoff = ( k>0 && k<j ) ? ( w[k-1] + w[k] ) * 0.5f : NAN;
    free(w);
    return cutoff;
}

// this function computes the median of contig stats and then clears the contig stats
static float sample_stats(sample_t *self, int n, model_param_t *model_param, int flags)
{
    int j = 0;
    float *tmp = (float *)malloc(n * sizeof(float));

    for (int i=0; i<n; i++)
    {
        self[i].lrr_median = get_median( self[i].kv_lrr_median.a, self[i].kv_lrr_median.n, NULL );
        self[i].lrr_sd = get_median( self[i].kv_lrr_sd.a, self[i].kv_lrr_sd.n, NULL );
        self[i].baf_sd = get_median( self[i].kv_baf_sd.a, self[i].kv_baf_sd.n, NULL );
        self[i].baf_conc = get_median( self[i].kv_baf_conc.a, self[i].kv_baf_conc.n, NULL );
        self[i].lrr_auto = get_median( self[i].kv_lrr_auto.a, self[i].kv_lrr_auto.n, NULL );
        free(self[i].kv_lrr_median.a);
        free(self[i].kv_lrr_sd.a);
        free(self[i].kv_baf_sd.a);
        free(self[i].kv_baf_conc.a);
        free(self[i].kv_lrr_auto.a);

        self[i].adjlrr_sd = self[i].lrr_sd;
        if ( model_param->order_lrr_gc == 0 )
        {
            self[i].coeffs[0] = get_median( self[i].kv_coeffs.a, self[i].kv_coeffs.n, NULL );
        }
        else if ( model_param->order_lrr_gc > 0 && self[i].kv_coeffs.n > 0 )
        {
            int medoid_idx = get_medoid( self[i].kv_coeffs.a, self[i].kv_coeffs.n, model_param->order_lrr_gc );
            for (int j=0; j<=model_param->order_lrr_gc; j++)
            {
                self[i].coeffs[j] = self[i].kv_coeffs.a[medoid_idx * (model_param->order_lrr_gc+1) + j];
            }
            free(self[i].kv_coeffs.a);
            self[i].rel_ess = get_median( self[i].kv_rel_ess.a, self[i].kv_rel_ess.n, NULL );
            free(self[i].kv_rel_ess.a);
            self[i].adjlrr_sd *= sqrtf( 1.0f - self[i].rel_ess ); // not perfect, but good enough(?)
        }
    }

    j = 0;
    for (int i=0; i<n; i++)
        if( !isnan(self[i].x_nonpar_lrr_median) )
            tmp[j++] = self[i].x_nonpar_lrr_median - self[i].lrr_median;
    float cutoff = get_cutoff( tmp, j );
    for (int i=0; i<n; i++)
    {
        if( self[i].x_nonpar_lrr_median - self[i].lrr_median < cutoff ) self[i].sex = SEX_MAL;
        else if( self[i].x_nonpar_lrr_median - self[i].lrr_median > cutoff ) self[i].sex = SEX_FEM;
    }

    if ( flags & WGS_DATA )
    {
        model_param->lrr_hap2dip = (float)M_LN2;
        model_param->lrr_auto2sex = 0.0f;
    }

    if ( isnan(model_param->lrr_hap2dip) || isnan(model_param->lrr_auto2sex) )
    {
        j = 0;
        for (int i=0; i<n; i++) if ( self[i].sex == SEX_MAL ) tmp[j++] = self[i].x_nonpar_lrr_median - self[i].lrr_median;
        float lrr_males = get_median( tmp, j, NULL );
        j = 0;
        for (int i=0; i<n; i++) if ( self[i].sex == SEX_FEM ) tmp[j++] = self[i].x_nonpar_lrr_median - self[i].lrr_median;
        float lrr_females = get_median( tmp, j, NULL );
        if ( isnan(model_param->lrr_hap2dip) ) model_param->lrr_hap2dip = lrr_females - lrr_males;
        if ( isnan(model_param->lrr_auto2sex) ) model_param->lrr_auto2sex = lrr_females;
    }
    free(tmp);
    return cutoff;
}

// this function computes several contig stats and then clears the contig data from the sample
static void sample_contig_stats(sample_t *self, const contig_param_t *contig_param, const genome_param_t *genome_param, const model_param_t *model_param, int flags)
{
    if (self->contig.size == 0) return;

    float *lrr, *baf;
    if ( flags & WGS_DATA )
    {
        lrr = ad_to_lrr(self->contig.data[AD0], self->contig.data[AD1], self->contig.size);
        baf = ad_to_baf(self->contig.data[AD0], self->contig.data[AD1], self->contig.size);
    }
    else
    {
        lrr = realloc_int16_to_float(self->contig.data[LRR], self->contig.size);
        baf = realloc_int16_to_float(self->contig.data[BAF], self->contig.size);
    }

    if ( contig_param->rid == genome_param->x_rid )
    {
        kv_int imap = {0, 0, NULL};
        for (int i=0; i<self->contig.size; i++)
        {
            if ( !isnan( baf[i] ) ) self->nhets++;
            int pos = contig_param->pos[ self->contig.vcf_imap[i] ];
            if ( pos > genome_param->x_nonpar_beg && pos < genome_param->x_nonpar_end &&
                  ( pos < genome_param->x_xtr_beg || pos > genome_param->x_xtr_end ) )
            {
                if ( !isnan( baf[i] ) ) self->x_nonpar_nhets++;
                kv_push(int, imap, i);
            }
        }
        self->x_nonpar_lrr_median = get_median( lrr, imap.n, imap.a );
        free(imap.a);
    }
    else if ( contig_param->rid == genome_param->y_rid )
    {
        kv_int imap = {0, 0, NULL};
        for (int i=0; i<self->contig.size; i++)
        {
            if ( !isnan( baf[i] ) ) self->nhets++;
            int pos = contig_param->pos[ self->contig.vcf_imap[i] ];
            if ( pos > genome_param->y_nonpar_beg && pos < genome_param->y_nonpar_end &&
                  ( pos < genome_param->y_xtr_beg || pos > genome_param->y_xtr_end ) )
                kv_push(int, imap, i);
        }
        self->y_nonpar_lrr_median = get_median( lrr, imap.n, imap.a );
        free(imap.a);
    }
    else if ( contig_param->rid == genome_param->mt_rid )
    {
        self->mt_lrr_median = get_median(lrr, self->contig.size, NULL );
    }
    else
    {
        if ( flags & WGS_DATA )
        {
            double f(double x, void *data) { return -log10_lkl_beta_binomial(self->contig.data[AD0], self->contig.data[AD1], self->contig.size, x); }
            double x; kmin_brent(f, 0.0, 0.5, NULL, KMIN_EPS, &x);
            kv_push(float, self->kv_baf_sd, (float)x );
        }
        else
        {
            kv_push(float, self->kv_baf_sd, get_sample_sd( baf, self->contig.size, NULL ) );
        }
        for (int i=0; i<self->contig.size; i++) if ( !isnan( baf[i] ) ) self->nhets++;
        kv_push(float, self->kv_lrr_median, get_median(lrr, self->contig.size, NULL ) );
        kv_push(float, self->kv_lrr_sd, get_sample_sd( lrr, self->contig.size, NULL ) );

        kv_push(float, self->kv_baf_conc, get_baf_conc( baf, self->contig.gt_phase, self->contig.size, NULL ) );
        if ( model_param->order_lrr_gc == 0 )
        {
            kv_push(float, self->kv_coeffs, get_median( lrr, self->contig.size, NULL ) );
        }
        // performs polynomial regression for LRR
        else if ( model_param->order_lrr_gc > 0 )
        {
            float tss = get_tss(lrr, self->contig.size);
            float *coeffs = polyfit(lrr, contig_param->gc, self->contig.size, self->contig.vcf_imap, model_param->order_lrr_gc);
            adjust_lrr(lrr, contig_param->gc, self->contig.size, self->contig.vcf_imap, coeffs, model_param->order_lrr_gc);
            coeffs[0] += get_median( lrr, self->contig.size, NULL ); // further adjusts by median
            for (int i=0; i<=model_param->order_lrr_gc; i++)
                kv_push(float, self->kv_coeffs, coeffs[i] );
            free(coeffs);
            float rss = get_tss(lrr, self->contig.size);
            kv_push(float, self->kv_rel_ess, 1.0f - rss / tss );
        }
        // compute autocorrelation after GC correction
        kv_push(float, self->kv_lrr_auto, get_lrr_auto( lrr, self->contig.size, NULL ) );
    }

    free(self->contig.vcf_imap);
    if ( flags & WGS_DATA )
    {
        free(self->contig.data[AD0]);
        free(self->contig.data[AD1]);
    }
    free(lrr);
    free(baf);
    free(self->contig.gt_phase);
}

static void sample_print(const sample_t *self, int n, FILE *restrict stream, const bcf_hdr_t *hdr, int flags)
{
    if (stream == NULL) return;
    char sex[3];
    sex[SEX_UNK] = 'U';
    sex[SEX_MAL] = 'M';
    sex[SEX_FEM] = 'F';
    if ( flags & WGS_DATA )
    {
        fprintf(stream, "SAMPLE\tCOV_MEDIAN\tCOV_SD\tBAF_CORR\tBAF_CONC\tNHETS\tX_NONPAR_NHETS\tX_NONPAR_COV_MEDIAN\tY_NONPAR_COV_MEDIAN\tMT_COV_MEDIAN\tSEX\tREL_ESS\tCOV_AUTO\n");
        for (int i=0; i<n; i++)
        {
            const char *sample_name = bcf_hdr_int2id(hdr, BCF_DT_SAMPLE, self[i].idx);
            fprintf(stream, "%s\t%.4f\t%.4f\t%.4f\t%.4f\t%d\t%d\t%.4f\t%.4f\t%.4f\t%c\t%.4f\t%.4f\n", sample_name, expf(self[i].lrr_median),
            expf(self[i].lrr_median)*self[i].lrr_sd, self[i].baf_sd, self[i].baf_conc, self[i].nhets, self[i].x_nonpar_nhets,
            expf(self[i].x_nonpar_lrr_median), expf(self[i].y_nonpar_lrr_median), expf(self[i].mt_lrr_median), sex[self[i].sex],
            self[i].rel_ess, self[i].lrr_auto);
        }
    }
    else
    {
        fprintf(stream, "SAMPLE\tLRR_MEDIAN\tLRR_SD\tBAF_SD\tBAF_CONC\tNHETS\tX_NONPAR_NHETS\tX_NONPAR_LRR_MEDIAN\tY_NONPAR_LRR_MEDIAN\tMT_LRR_MEDIAN\tSEX\tREL_ESS\tLRR_AUTO\n");
        for (int i=0; i<n; i++)
        {
            const char *sample_name = bcf_hdr_int2id(hdr, BCF_DT_SAMPLE, self[i].idx);
            fprintf(stream, "%s\t%.4f\t%.4f\t%.4f\t%.4f\t%d\t%d\t%.4f\t%.4f\t%.4f\t%c\t%.4f\t%.4f\n", sample_name, self[i].lrr_median,
            self[i].lrr_sd, self[i].baf_sd, self[i].baf_conc, self[i].nhets, self[i].x_nonpar_nhets,
            self[i].x_nonpar_lrr_median, self[i].y_nonpar_lrr_median, self[i].mt_lrr_median, sex[self[i].sex],
            self[i].rel_ess, self[i].lrr_auto);
        }
    }
    if (stream != stdout && stream != stderr) fclose(stream);
}

/*********************************
 * VCF READ AND WRITE METHODS    *
 *********************************/

// moves synced bcf reader to first useful line for reader0
static int bcf_sr_next_line_reader0(bcf_srs_t *sr)
{
    int nret = bcf_sr_next_line(sr);
    while ( nret > 0 && !bcf_sr_has_line(sr, 0) ) nret = bcf_sr_next_line(sr);
    return nret;
}

// write one contig
static int put_contig(bcf_srs_t *sr, int rid, const sample_t *sample, const contig_param_t contig_param, int flags, htsFile *out_fh, bcf_hdr_t *out_hdr)
{
    bcf_hdr_t *hdr = bcf_sr_get_header(sr, 0);
    bcf_sr_seek(sr, bcf_hdr_id2name( hdr, rid ), 0);
    int nsmpl = bcf_hdr_nsamples(out_hdr);

    int *synced_iter = (int *)calloc(nsmpl, sizeof(int));
    float *bdev = (float *)calloc(nsmpl, sizeof(float));
    float *ldev = (float *)calloc(nsmpl, sizeof(float));
    int *bdev_phase = (int *)calloc(nsmpl, sizeof(int));

    int i;
    for (i=0; bcf_sr_next_line_reader0(sr); i++)
    {
        bcf1_t *line = bcf_sr_get_line(sr, 0);
        if (rid != line->rid) break;

        for (int j=0; j<nsmpl; j++)
        {
            while ( synced_iter[j] < sample[j].contig.size-1 && sample[j].contig.vcf_imap[ synced_iter[j] ] < i ) synced_iter[j]++;
            if ( sample[j].contig.vcf_imap[ synced_iter[j] ] == i )
            {
                if ( sample[j].contig.ldev ) ldev[ sample[j].idx ] = int16_to_float( sample[j].contig.ldev[ synced_iter[j] ] );
                if ( sample[j].contig.bdev ) bdev[ sample[j].idx ] = int16_to_float( sample[j].contig.bdev[ synced_iter[j] ] );
                if ( sample[j].contig.bdev_phase ) bdev_phase[ sample[j].idx ] = sample[j].contig.bdev_phase[ synced_iter[j] ];
            }
            else
            {
                // if no match variant found, match the end of the contig or keep conservative
                if ( i==0 && sample[j].contig.bdev )
                    bdev[ sample[j].idx ] = int16_to_float( sample[j].contig.bdev[0] );
                if ( i==0 && sample[j].contig.ldev )
                    ldev[ sample[j].idx ] = int16_to_float( sample[j].contig.ldev[0] );
                if ( sample[j].contig.bdev && int16_to_float( sample[j].contig.bdev[ synced_iter[j] ] ) == 0.0f )
                    bdev[ sample[j].idx ] = 0.0f;
                if ( sample[j].contig.ldev && int16_to_float( sample[j].contig.ldev[ synced_iter[j] ] ) == 0.0f )
                    ldev[ sample[j].idx ] = 0.0f;
                if ( sample[j].contig.bdev_phase ) bdev_phase[ sample[j].idx ] = 0;
            }
        }
        if ( !(flags & NO_ANNOT) )
        {
            bcf_update_format_float(out_hdr, line, "Ldev", ldev, (int)nsmpl);
            bcf_update_format_float(out_hdr, line, "Bdev", bdev, (int)nsmpl);
        }
        bcf_update_format_int32(out_hdr, line, "Bdev_Phase", bdev_phase, (int)nsmpl);

        bcf_write(out_fh, out_hdr, line);
    }

    free(synced_iter);
    free(ldev);
    free(bdev);
    free(bdev_phase);

    return i;
}

// write header
static bcf_hdr_t *print_hdr(htsFile *out_fh, bcf_hdr_t *hdr, int argc, char *argv[], int record_cmd_line, int flags)
{
    bcf_hdr_t *out_hdr = bcf_hdr_dup(hdr);
    if ( !(flags & NO_ANNOT) )
    {
        if ( bcf_hdr_id2int(out_hdr, BCF_DT_ID, "Ldev") < 0 )
            bcf_hdr_append(out_hdr, "##FORMAT=<ID=Ldev,Number=1,Type=Float,Description=\"LRR deviation due to chromosomal alteration\">");
        if ( bcf_hdr_id2int(out_hdr, BCF_DT_ID, "Bdev") < 0 )
            bcf_hdr_append(out_hdr, "##FORMAT=<ID=Bdev,Number=1,Type=Float,Description=\"BAF deviation due to chromosomal alteration\">");
    }
    if ( bcf_hdr_id2int(out_hdr, BCF_DT_ID, "Bdev_Phase") < 0 )
        bcf_hdr_append(out_hdr, "##FORMAT=<ID=Bdev_Phase,Number=1,Type=Integer,Description=\"BAF deviation phase, if available\">");
    if (record_cmd_line) bcf_hdr_append_version(out_hdr, argc, argv, "bcftools_mocha");
    bcf_hdr_write(out_fh, out_hdr);
    return out_hdr;
}

// read one contig
static int get_contig(bcf_srs_t *sr, int rid, sample_t *sample, contig_param_t *contig_param, int flags, float median_baf_adjust, int min_dist)
{
    bcf_hdr_t *hdr = bcf_sr_get_header(sr, 0);
    bcf_sr_seek(sr, bcf_hdr_id2name( hdr, rid ), 0);

    bcf_fmt_t *baf_fmt = NULL, *lrr_fmt = NULL;
    bcf_info_t *info;
    int nsmpl = bcf_hdr_nsamples(hdr);

    kv_int kv_pos = {0, 0, NULL};
    kv_float kv_gc = {0, 0, NULL};

    int *size = (int *)calloc(nsmpl, sizeof(int));
    int *next_het = (int *)calloc(nsmpl, sizeof(int));
    int *next_hom = (int *)calloc(nsmpl, sizeof(int));
    kv_int *kv_vcf_imap = (kv_int *)calloc(nsmpl, sizeof(kv_int));
    kv_int16_t *kv_data[2];
    for (int i=0; i<2; i++) kv_data[i] = (kv_int16_t *)calloc(nsmpl, sizeof(kv_int16_t));
    kv_int8_t *kv_gt_phase = (kv_int8_t *)calloc(nsmpl, sizeof(kv_int8_t));

    int i;
    int nad_arr = 0, *ad_arr = NULL, nad = 0, nalleles = 0;
    int ngt_arr = 0, *gt_arr = NULL, ngt = 0, max_ploidy = 0;
    for (i=0; bcf_sr_next_line_reader0(sr); i++)
    {
        bcf1_t *line = bcf_sr_get_line(sr, 0);
        if (rid != line->rid) break;
        int pos = line->pos + 1;
        kv_push(int, kv_pos, pos);
        if ( ( info = bcf_get_info( hdr, line, "GC" ) ) ) kv_push(float, kv_gc, info->v1.f );
        else kv_push(float, kv_gc, NAN);

        // if failing inclusion/exclusion requirement, skip line
        if ( (flags & FLT_EXCLUDE) && bcf_sr_get_line(sr, 1) ) continue;
        if ( (flags & FLT_INCLUDE) && !bcf_sr_get_line(sr, 1) ) continue;

        // if there are no genotypes, skip line
        ngt = bcf_get_genotypes(hdr, line, &gt_arr, &ngt_arr);
        if ( ngt <= 0 ) continue;
        max_ploidy = ngt / nsmpl;

        // if neither AD nor LRR and BAF formats are present, skip line
        if ( flags & WGS_DATA )
        {
            nad = bcf_get_format_int32(hdr, line, "AD", &ad_arr, &nad_arr);
            if ( nad <= 0 ) continue;
            nalleles = nad / nsmpl;
        }
        else
        {
            if ( !(lrr_fmt = bcf_get_fmt(hdr, line, "LRR")) || !(baf_fmt = bcf_get_fmt(hdr, line, "BAF")) ) continue;

            // make nan all BAF values for non heterozygous SNPs (we do not use those BAFs)
            int nbaf = 0;
            for (int j=0; j<nsmpl; j++)
            {
                int *ptr = gt_arr + max_ploidy * sample[j].idx;
                if ( bcf_gt_is_missing(ptr[0]) || bcf_gt_is_missing(ptr[1]) || ptr[1]==bcf_int32_vector_end || (ptr[0]>>1) == (ptr[1]>>1) )
                    ((float *)(baf_fmt->p + baf_fmt->size * sample[j].idx))[0] = NAN;
                if ( !isnan(((float *)(baf_fmt->p + baf_fmt->size * sample[j].idx))[0]) ) nbaf++;
            }

            if ( median_baf_adjust >= 0 && nbaf >= median_baf_adjust )
            {
                float baf_median = get_median( (float *)baf_fmt->p, nsmpl, NULL ) - 0.5f;
                for (int j=0; j<nsmpl; j++)
                    ((float *)(baf_fmt->p + baf_fmt->size * sample[j].idx))[0] -= baf_median;
            }
        }

        for (int j=0; j<nsmpl; j++)
        {
            // if genotype is not diploid, skip line
            int *ptr = gt_arr + max_ploidy * sample[j].idx;
            if ( bcf_gt_is_missing(ptr[0]) || bcf_gt_is_missing(ptr[1]) || ptr[1]==bcf_int32_vector_end ) continue;
            int gt0 = bcf_gt_allele(ptr[0]), gt1 = bcf_gt_allele(ptr[1]);
            int8_t gt_phase = (int8_t)( SIGN( gt0 - gt1 ) * bcf_gt_is_phased(ptr[1]) );

            if ( flags & WGS_DATA )
            {
                assert(gt0 >= 0 && gt0 < nalleles && gt1 >= 0 && gt1 < nalleles);
                int *ad_ptr = ad_arr + nalleles * sample[j].idx;
                // missing allelic depth information
                if ( ad_ptr[gt0] == 0 && ad_ptr[gt1] == 0) continue;
                // site too close to last het site or hom site too close to last site
                if ( pos < next_het[j] || ( (gt0 == gt1) && (pos < next_hom[j]) ) ) continue;
                if ( pos < next_hom[j] )
                {
                    // substitute the last hom site with the current het site
                    kv_vcf_imap[j].n--;
                    kv_gt_phase[j].n--;
                    kv_data[AD0][j].n--;
                    kv_data[AD1][j].n--;
                }
                else
                {
                    size[j]++;
                }

                if ( gt0 != gt1 ) next_het[j] = pos + min_dist;
                next_hom[j] = pos + min_dist;

                kv_push(int, kv_vcf_imap[j], i);
                kv_push(int8_t, kv_gt_phase[j], gt_phase );
                int16_t ad0 = (int16_t)ad_ptr[gt0];
                int16_t ad1 = (gt0 == gt1) ? (int16_t)INT16_NAN : (int16_t)ad_ptr[gt1];
                kv_push(int16_t, kv_data[AD0][j], ad0);
                kv_push(int16_t, kv_data[AD1][j], ad1);
            }
            else
            {
                size[j]++;
                kv_push(int, kv_vcf_imap[j], i);
                kv_push(int8_t, kv_gt_phase[j], gt_phase );
                float lrr = ((float *)(lrr_fmt->p + lrr_fmt->size * sample[j].idx))[0];
                float baf = ((float *)(baf_fmt->p + baf_fmt->size * sample[j].idx))[0];
                kv_push(int16_t, kv_data[LRR][j], float_to_int16(lrr));
                kv_push(int16_t, kv_data[BAF][j], float_to_int16(baf));
            }
        }
    }
    free(ad_arr);
    free(gt_arr);
    contig_param->n = kv_pos.n;

    if ( i > 0 )
    {
        contig_param->pos = (int *)realloc(kv_pos.a, kv_pos.n * sizeof(int));
        contig_param->gc = (float *)realloc(kv_gc.a, kv_gc.n * sizeof(float));
        for (int j=0; j<nsmpl; j++)
        {
            sample[j].contig.size = size[j];
            sample[j].contig.vcf_imap = (int *)realloc(kv_vcf_imap[j].a, size[j] * sizeof(int));
            sample[j].contig.gt_phase = (int8_t *)realloc(kv_gt_phase[j].a, size[j] * sizeof(int8_t));
            for (int i=0; i<2; i++) sample[j].contig.data[i] = (int16_t *)realloc(kv_data[i][j].a, size[j] * sizeof(int16_t));
        }
    }

    free(size);
    free(next_het);
    free(next_hom);
    free(kv_vcf_imap);
    free(kv_data[LRR]);
    free(kv_data[BAF]);
    free(kv_gt_phase);

    return i;
}

/*********************************
 * MAIN PART OF THE COMMAND      *
 *********************************/

// TODO add examples at the end of the usage section
static void usage()
{
    fprintf(stderr, "\n");
    fprintf(stderr, "About:   MOsaic CHromosomal Alterations caller, requires phased genotypes (GT)\n");
    fprintf(stderr, "         and either B-allele frequency (BAF) and Log R Ratio intensity (LRR)\n");
    fprintf(stderr, "         or allelic depth coverage (AD). (version %s)\n", MOCHA_VERSION);
    fprintf(stderr, "Usage:   bcftools mocha [OPTIONS] <in.vcf>\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Required options:\n");
    fprintf(stderr, "    -r, --rules <assembly>[?]         predefined genome reference rules, 'list' to print available settings, append '?' for details\n");
    fprintf(stderr, "    -R, --rules-file <file>           genome reference rules, space/tab-delimited CHROM:FROM-TO,TYPE\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "General Options:\n");
    fprintf(stderr, "    -s, --samples [^]<list>           comma separated list of samples to include (or exclude with \"^\" prefix)\n");
    fprintf(stderr, "    -S, --samples-file [^]<file>      file of samples to include (or exclude with \"^\" prefix)\n");
    fprintf(stderr, "        --force-samples               only warn about unknown subset samples\n");
    fprintf(stderr, "    -v, --variants [^]<file>          tabix-indexed [compressed] VCF/BCF file containing variants\n");
    fprintf(stderr, "                                      to include (or exclude with \"^\" prefix) in the analysis\n");
    fprintf(stderr, "        --threads <int>               number of extra output compression threads [0]\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Output Options:\n");
    fprintf(stderr, "    -o, --output <file>               write output to a file [no output]\n");
    fprintf(stderr, "    -O, --output-type <b|u|z|v>       b: compressed BCF, u: uncompressed BCF, z: compressed VCF, v: uncompressed VCF [v]\n");
    fprintf(stderr, "        --no-version                  do not append version and command line to the header\n");
    fprintf(stderr, "    -a  --no-annotations              omit Ldev and Bdev FORMAT from output VCF (requires --output)\n");
    fprintf(stderr, "    -m, --mosaic-calls <file>         write mosaic chromosomal alterations to a file [standard output]\n");
    fprintf(stderr, "    -g, --genome-stats <file>         write sample genome-wide statistics to a file [no output]\n");
    fprintf(stderr, "    -u, --bed-ucsc <file>             write UCSC bed track to a file [no output]\n");
    fprintf(stderr, "    -l  --no-log                      suppress progress report on standard error\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "HMM Options:\n");
    fprintf(stderr, "    -x, --xy-prob <float>             transition probability [%.0e]\n", xy_prob_default);
    fprintf(stderr, "    -z, --err-prob <float>            uniform error probability [%.0e]\n", err_prob_default);
    fprintf(stderr, "    -f, --flip-prob <float>           phase flip probability [%.0e]\n", flip_prob_default);
    fprintf(stderr, "    -t, --telomere-advantage <float>  telomere advantage [%.0e]\n", telomere_prob_default);
    fprintf(stderr, "    -c, --centromere-penalty <float>  centromere penalty [%.0e]\n", centromere_prob_default);
    fprintf(stderr, "    -p  --cnp <file>                  list of regions to genotype in BED format\n");
    fprintf(stderr, "    -n, --cnf <list>                  comma separated list of copy number fractions for LRR+BAF model\n");
    fprintf(stderr, "                                      [%s]\n", cnf_default);
    fprintf(stderr, "    -b, --bdev <list>                 comma separated list of inverse BAF deviations for BAF+phase model\n");
    fprintf(stderr, "                                      [%s]\n", bdev_default);
    fprintf(stderr, "    -d, --min-dist <int>              minimum base pair distance between consecutive sites for WGS data [%d]\n", min_dist_default);
    fprintf(stderr, "        --LRR-hap2dip <float>         LRR difference between haploid and diploid [estimated from X nonPAR]\n");
    fprintf(stderr, "        --LRR-auto2sex <float>        LRR difference between autosomes and diploid sex chromosomes [estimated from X nonPAR]\n");
    fprintf(stderr, "        --LRR-weight <float>          relative contribution from LRR for LRR+BAF model [%g]\n", lrr_bias_default);
    fprintf(stderr, "        --median-BAF-adjust <int>     minimum number of heterozygous genotypes required to perform\n");
    fprintf(stderr, "                                      median BAF adjustment (-1 for no BAF adjustment) [%d]\n", median_baf_adjust_default);
    fprintf(stderr, "        --order-LRR-GC <int>          order of polynomial in local GC content to be used for polynomial\n");
    fprintf(stderr, "                                      regression of LRR (-1 for no LRR adjustment, %d maximum) [%d]\n", MAX_ORDER, order_lrr_gc_default);
    fprintf(stderr, "\n");
    fprintf(stderr, "Examples:\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Publications:\n");
    fprintf(stderr, "    Loh P., Genovese G., McCarroll S., Price A. et al. Insights about clonal expansions from 8,342 mosaic\n");
    fprintf(stderr, "    chromosomal alterations. Nature 559, 350–355 (2018). [PMID: 29995854] [DOI: 10.1038/s41586-018-0321-x]\n");
    fprintf(stderr, "\n");
    exit(1);
}

static float *readlistf(const char *str, int *n, float min, float max)
{
    char *tmp, **list = hts_readlist(str, 0, n);
    if ( *n >= 128 ) error("Cannot handle list of 128 or more parameters: %s\n", str);
    float *ret = (float *)malloc(*n * sizeof(float));
    for (int i=0; i<*n; i++)
    {
        ret[i] = strtof(list[i], &tmp);
        if ( *tmp ) error("Could not parse: %s\n", list[i]);
        if ( min!=max && (ret[i]<min || ret[i]>max) )
            error("Expected values from the interval [%f,%f], found %s\n", min, max, list[i]);
        free(list[i]);
    }
    free(list);
    return ret;
}

static FILE *get_file_handle(const char *str)
{
    FILE *ret;
    if ( strcmp(str, "-") == 0 )
        ret = stdout;
    else
    {
        ret = fopen(str, "w");
        if ( !ret ) error("Failed to open %s: %s\n", str, strerror(errno));
    }
    return ret;
}

int main_vcfmocha(int argc, char *argv[])
{
    // program options
    int rules_is_file = 0;
    int sample_is_file = 0;
    int force_samples = 0;
    int flags = 0;
    int output_type = FT_VCF;
    int n_threads = 0;
    int record_cmd_line = 1;
    char *sample_names = NULL;
    char *output_fname = NULL;
    char *mocha_fname = NULL;
    char *stats_fname = NULL;
    char *ucsc_fname = NULL;
    char *cnp_fname = NULL;
    regidx_t *cnp_idx = NULL;
    regitr_t *cnp_itr = NULL;
    char *filter_fname = NULL;
    char *rules = NULL;
    sample_t *sample = NULL;
    genome_param_t genome_param;
    memset(&genome_param, 0, sizeof(genome_param));
    genome_param.x_rid = -1;
    genome_param.y_rid = -1;
    genome_param.mt_rid = -1;
    contig_param_t contig_param;
    memset(&contig_param, 0, sizeof(contig_param));
    FILE *out_fm = stdout;
    FILE *out_fg = NULL;
    FILE *out_fu = NULL;
    bcf_srs_t *sr = NULL;
    bcf_hdr_t *hdr = NULL;
    bcf_hdr_t *out_hdr = NULL;
    htsFile *out_fh = NULL;
    kv_mocha_t kv_mocha = {0, 0, NULL};
    const char *cnf = cnf_default;
    const char *bdev = bdev_default;

    // model parameters
    model_param_t model_param;
    model_param.xy_prob = xy_prob_default;
    model_param.err_prob = err_prob_default;
    model_param.flip_prob = flip_prob_default;
    model_param.telomere_prob = telomere_prob_default;
    model_param.centromere_prob = centromere_prob_default;
    model_param.min_dist = min_dist_default;
    model_param.lrr_bias = lrr_bias_default;
    model_param.lrr_hap2dip = NAN;
    model_param.lrr_auto2sex = NAN;
    model_param.median_baf_adjust = median_baf_adjust_default;
    model_param.order_lrr_gc = order_lrr_gc_default;

    int c;
    char *tmp = NULL;

    static struct option loptions[] =
    {
        {"samples", required_argument, NULL, 's'},
        {"samples-file", required_argument, NULL, 'S'},
        {"force-samples", no_argument, NULL, 1},
        {"variants", required_argument, NULL,'v'},
        {"no-annotations", no_argument, NULL,'a'},
        {"mosaic-calls", required_argument, NULL, 'm'},
        {"genome-stats", required_argument, NULL, 'g'},
        {"bed-ucsc", required_argument, NULL, 'u'},
        {"rules", required_argument, NULL, 'r'},
        {"rules-file", required_argument, NULL, 'R'},
        {"xy-prob", required_argument, NULL, 'x'},
        {"err-prob", required_argument, NULL, 'z'},
        {"flip-prob", required_argument, NULL, 'f'},
        {"centromere-penalty", required_argument, NULL, 'c'},
        {"telomere-advantage", required_argument, NULL, 't'},
        {"cnf", required_argument, NULL, 'c'},
        {"bdev", required_argument, NULL, 'b'},
        {"min-dist", required_argument, NULL, 'd'},
        {"LRR-hap2dip", required_argument, NULL, 2},
        {"LRR-auto2sex", required_argument, NULL, 3},
        {"LRR-weight", required_argument, NULL, 4},
        {"cnp", required_argument, NULL, 'p'},
        {"median-BAF-adjust", required_argument, NULL, 5},
        {"order-LRR-GC", required_argument, NULL, 6},
        {"output", required_argument, NULL, 'o'},
        {"output-type", required_argument, NULL, 'O'},
        {"no-version", no_argument, NULL, 8},
        {"threads", required_argument, NULL, 9},
        {"no-log", no_argument, NULL, 'd'},
        {NULL, 0, NULL, 0}
    };
    while ((c = getopt_long(argc, argv, "h?r:R:s:S:v:o:O:am:g:u:lx:z:f:t:c:d:p:n:b:", loptions, NULL)) >= 0)
    {
        switch (c)
        {
            case 'o': output_fname = optarg; break;
            case 'O':
                switch (optarg[0]) {
                    case 'b': output_type = FT_BCF_GZ; break;
                    case 'u': output_type = FT_BCF; break;
                    case 'z': output_type = FT_VCF_GZ; break;
                    case 'v': output_type = FT_VCF; break;
                    default: error("The output type \"%s\" not recognised\n", optarg);
                };
                break;
            case 'v':
                if (optarg[0]=='^')
                {
                    filter_fname = optarg + 1;
                    flags |= FLT_EXCLUDE;
                }
                else
                {
                    filter_fname = optarg;
                    flags |= FLT_INCLUDE;
                }
                break;
            case 'a': flags |= NO_ANNOT; break;
            case 'x':
                model_param.xy_prob = strtof(optarg, &tmp);
                if ( *tmp ) error("Could not parse: --xy-prob %s\n", optarg);
                break;
            case 'z':
                model_param.err_prob = strtof(optarg, &tmp);
                if ( *tmp ) error("Could not parse: --err-prob %s\n", optarg);
                break;
            case 'f':
                model_param.flip_prob = strtof(optarg, &tmp);
                if ( *tmp ) error("Could not parse: --flip-prob %s\n", optarg);
                break;
            case 't':
                model_param.telomere_prob = strtof(optarg, &tmp);
                if ( *tmp ) error("Could not parse: --telomere-advantage %s\n", optarg);
                break;
            case 'c':
                model_param.centromere_prob = strtof(optarg, &tmp);
                if ( *tmp ) error("Could not parse: --centromere-penalty %s\n", optarg);
                break;
            case 'd':
                model_param.min_dist = (int)strtol(optarg, &tmp, 0);
                if ( *tmp ) error("Could not parse: --min-dist %s\n", optarg);
                break;
            case 'n': cnf = optarg; break;
            case 'b': bdev = optarg; break;
            case 'm': mocha_fname = optarg; break;
            case 'g': stats_fname = optarg; break;
            case 'u': ucsc_fname = optarg; break;
            case 'p': cnp_fname = optarg; break;
            case 's': sample_names = optarg; break;
            case 'S': sample_names = optarg; sample_is_file = 1; break;
            case 'r': rules = optarg; break;
            case 'R': rules = optarg; rules_is_file = 1; break;
            case  1 : force_samples = 1; break;
            case  8 : record_cmd_line = 0; break;
            case  9 :
                n_threads = (int)strtol(optarg, &tmp, 0);
                if ( *tmp ) error("Could not parse: --threads %s\n", optarg);
                break;
            case  2 :
                model_param.lrr_hap2dip = strtof(optarg, &tmp);
                if ( *tmp ) error("Could not parse: --LRR-hap2dip %s\n", optarg);
                break;
            case  3 :
                model_param.lrr_auto2sex = strtof(optarg, &tmp);
                if ( *tmp ) error("Could not parse: --LRR-auto2sex %s\n", optarg);
                break;
            case  4 :
                model_param.lrr_bias = strtof(optarg, &tmp);
                if ( *tmp ) error("Could not parse: --LRR-weight %s\n", optarg);
                break;
            case  5 :
                model_param.median_baf_adjust = (int)strtol(optarg, &tmp, 0);
                if ( *tmp ) error("Could not parse: --median-BAF-adjust %s\n", optarg);
                break;
            case  6 :
                model_param.order_lrr_gc = (int)strtol(optarg, &tmp, 0);
                if ( *tmp ) error("Could not parse: --order-LRR-GC %s\n", optarg);
                break;
            case 'l': flags |= NO_LOG; break;
            case 'h':
            case '?': usage(); break;
            default: error("Unknown argument: %s\n", optarg);
        }
    }

    if ( !rules )
    {
        fprintf(stderr, "Genome reference assembly was not specified with --rules or --rules-file\n");
        usage();
    }
    int len = strlen(rules);
    if ( !rules_is_file && ( strncmp(rules, "list", 4) == 0 || rules[len-1]=='?' ) ) genome_init_alias(NULL, rules, NULL);

    if ( output_fname == NULL && ( flags & NO_ANNOT ) )
    {
        fprintf(stderr, "Option --no-annotations requires option --output\n");
        usage();
    }

    if (model_param.order_lrr_gc > MAX_ORDER)
    {
        fprintf(stderr, "Polynomial order must not be greater than %d: --order-LRR-GC %d\n", MAX_ORDER, model_param.order_lrr_gc);
        usage();
    }

    // parse parameters defining hidden states
    model_param.cnf = readlistf(cnf, &model_param.cnf_n, 0.0f, 4.0f);
    ks_introsort_float((size_t)model_param.cnf_n, model_param.cnf);
    model_param.bdev = readlistf(bdev, &model_param.bdev_n, 2.0f, INFINITY);
    for (int i=0; i<model_param.bdev_n; i++) model_param.bdev[i] = 1.0f / model_param.bdev[i]; // compute inverses
    ks_introsort_float((size_t)model_param.bdev_n, model_param.bdev);

    // output tables with mosaic chromosomal alteration calls
    if ( mocha_fname ) out_fm = get_file_handle( mocha_fname );
    if ( stats_fname ) out_fg = get_file_handle( stats_fname );
    if ( ucsc_fname  ) out_fu = get_file_handle( ucsc_fname  );

    // read list of regions to genotype
    if ( cnp_fname )
    {
        cnp_idx = regidx_init(cnp_fname, cnp_parse, NULL, sizeof(int), NULL);
        if ( !cnp_idx ) error("Error: failed to initialize CNP regions: --cnp %s\n", cnp_fname);
        cnp_itr = regitr_init(cnp_idx);
    }

    // input VCF
    char *input_fname = NULL;
    if ( optind>=argc )
    {
        if ( !isatty(fileno((FILE *)stdin)) ) input_fname = "-";
    }
    else input_fname = argv[optind];
    if ( !input_fname ) usage();

    // create synced reader object
    sr = bcf_sr_init();
    bcf_sr_set_opt(sr, BCF_SR_REQUIRE_IDX);
    if ( bcf_sr_set_threads(sr, n_threads)<0 ) error("Failed to create threads\n");
    if ( !bcf_sr_add_reader(sr, input_fname) ) error("Failed to open %s: %s\n", input_fname, bcf_sr_strerror(sr->errnum));
    if ( filter_fname )
        if ( !bcf_sr_add_reader(sr, filter_fname) ) error("Failed to open %s: %s\n", filter_fname, bcf_sr_strerror(sr->errnum));

    // check whether the necessary information has been included in the VCF
    hdr = bcf_sr_get_header(sr, 0);
    if ( bcf_hdr_nsamples(hdr) == 0 )
        error("Error: input VCF file has no samples\n");
    if ( bcf_hdr_id2int( hdr, BCF_DT_ID, "GT" ) < 0 )
        error("Error: input VCF file has no GT format field\n");
    if ( !( bcf_hdr_id2int( hdr, BCF_DT_ID, "AD" ) < 0 ) )
        flags |= WGS_DATA;
    else if ( (bcf_hdr_id2int( hdr, BCF_DT_ID, "LRR" ) < 0 || bcf_hdr_id2int( hdr, BCF_DT_ID, "BAF" ) < 0) )
        error("Error: input VCF file must contain either the AD format field or the LRR and BAF format fields\n");
    if ( model_param.order_lrr_gc > 0 && ( bcf_hdr_id2int( hdr, BCF_DT_ID, "GC" ) < 0 ) )
        error("Error: input VCF has no GC info field\n");

    // initialize genome parameters
    genome_param.length = (int *)calloc(hdr->n[BCF_DT_CTG], sizeof(int));
    for (int rid=0; rid<hdr->n[BCF_DT_CTG]; rid++)
        genome_param.length[rid] = hdr->id[BCF_DT_CTG][rid].val->info[0];
    genome_param.cen_beg = (int *)calloc(hdr->n[BCF_DT_CTG], sizeof(int));
    genome_param.cen_end = (int *)calloc(hdr->n[BCF_DT_CTG], sizeof(int));
    if ( rules_is_file )
        genome_init_file(&genome_param, rules, hdr);
    else
        genome_init_alias(&genome_param, rules, hdr);
    if ( !(flags & NO_LOG) ) fprintf(stderr, "Using genome assembly rules from %s\n", rules);

    // subset VCF file
    if (sample_names)
    {
        int ret = bcf_hdr_set_samples(hdr, sample_names, sample_is_file);
        if ( ret<0 ) error("Error parsing the sample list\n");
        else if ( ret>0 )
        {
            if ( force_samples )
                fprintf(stderr, "Warn: sample #%d not found in the header... skipping\n", ret);
            else
                error("Error: sample #%d not found in the header. Use \"--force-samples\" to ignore this error\n", ret);
        }
        if ( bcf_hdr_nsamples(hdr) == 0 )
            error("Error: subsetting has removed all samples\n");
    }

    // output VCF
    if (output_fname)
    {
        out_fh = hts_open(output_fname, hts_bcf_wmode(output_type));
        if ( out_fh == NULL ) error("Can't write to \"%s\": %s\n", output_fname, strerror(errno));
        if ( n_threads ) hts_set_opt(out_fh, HTS_OPT_THREAD_POOL, sr->p);
        out_hdr = print_hdr(out_fh, hdr, argc, argv, record_cmd_line, flags);
    }

    int nsmpl = bcf_hdr_nsamples(hdr);
    if ( nsmpl == 0 ) error("Subsetting has removed all samples\n");
    sample = (sample_t *)calloc(nsmpl, sizeof(sample_t));
    for (int i=0; i<nsmpl; i++)
    {
        sample[i].idx = i; // TODO is this right?
        sample[i].x_nonpar_lrr_median = NAN;
        sample[i].y_nonpar_lrr_median = NAN;
        sample[i].mt_lrr_median = NAN;
    }

    for (int rid=0; rid < hdr->n[BCF_DT_CTG]; rid++)
    {
        contig_param.rid = rid;
        int nret = get_contig(sr, rid, sample, &contig_param, flags, model_param.median_baf_adjust, model_param.min_dist);
        if ( nret<=0 ) continue;
        if ( !(flags & NO_LOG) ) fprintf(stderr, "Read %d variants from contig %s\n", nret, bcf_hdr_id2name( hdr, rid ));
        if( genome_param.length[rid] < contig_param.pos[contig_param.n-1] )
            genome_param.length[rid] = contig_param.pos[contig_param.n-1];
        for (int j=0; j<nsmpl; j++) sample_contig_stats(sample + j, &contig_param, &genome_param, &model_param, flags);
        free(contig_param.pos);
        free(contig_param.gc);
    }

    if ( !(flags & NO_LOG) ) fprintf(stderr, "Compute and print sample statistics\n");
    float cutoff = sample_stats(sample, nsmpl, &model_param, flags);
    if ( !(flags & NO_LOG) ) fprintf(stderr, "Estimated cutoff: %.4f\n", cutoff);
    sample_print(sample, nsmpl, out_fg, hdr, flags);

    if ( isnan( model_param.lrr_hap2dip ) )
        error("Error: Unable to estimate LRR-hap2dip from X nonPAR. Make sure "
              "the regions are present in the VCF or specify the parameter\n");
    if ( isnan( model_param.lrr_auto2sex ) )
        error("Error: Unable to estimate LRR-auto2sex from autosomes and X nonPAR. Make sure "
              "the regions are present in the VCF or specify the parameter\n");

    if ( !(flags & NO_LOG) ) fprintf(stderr, "Estimated parameters: LRR-hap2dip=%.4f LRR-auto2sex=%.4f LRR-dip2trip=%.4f\n",
        model_param.lrr_hap2dip, model_param.lrr_auto2sex, model_param.lrr_hap2dip * (log2f(3.0f) - 1.0f));

    for (int rid=0; rid < hdr->n[BCF_DT_CTG]; rid++)
    {
        contig_param.rid = rid;
        int nret = get_contig(sr, rid, sample, &contig_param, flags, model_param.median_baf_adjust, model_param.min_dist);
        if ( nret<=0 ) continue;
        if ( !(flags & NO_LOG) ) fprintf(stderr, "Read %d variants from contig %s\n", nret, bcf_hdr_id2name( hdr, rid ));
        for (int j=0; j<nsmpl; j++)
        {
            if ( cnp_idx ) regidx_overlap(cnp_idx, bcf_hdr_id2name( hdr, rid ), 0, genome_param.length[rid], cnp_itr);
            sample_contig_run(sample + j, &kv_mocha, &contig_param, &genome_param, &model_param, cnp_itr, flags);
        }

        if (output_fname)
        {
            nret = put_contig(sr, rid, sample, contig_param, flags, out_fh, out_hdr);
            if ( !(flags & NO_LOG) ) fprintf(stderr, "Written %d variants for contig %s\n", nret, bcf_hdr_id2name( hdr, rid ));
        }
        free(contig_param.pos);
        free(contig_param.gc);
        for (int j=0; j<nsmpl; j++)
        {
            free(sample[j].contig.vcf_imap);
            free(sample[j].contig.ldev);
            free(sample[j].contig.bdev);
            free(sample[j].contig.bdev_phase);
        }
    }

    // clear parameter values
    free(model_param.cnf);
    free(model_param.bdev);
    free(genome_param.length);
    free(genome_param.cen_beg);
    free(genome_param.cen_end);

    // free precomputed tables
    ad_to_lrr(NULL, NULL, 0);
    precompute_log_gammas(NULL, NULL, 0, NULL);
    precompute_lod_gammas(NULL, NULL, 0, NULL, NAN, NAN, NULL, NULL, NULL);

    // write table with mosaic chromosomal alterations (and UCSC bed track)
    mocha_print(kv_mocha.a, kv_mocha.n, out_fm, hdr, flags);
    mocha_print_ucsc(kv_mocha.a, kv_mocha.n, out_fu, hdr);
    free(kv_mocha.a);

    // close output VCF
    if (output_fname)
    {
        bcf_hdr_destroy(out_hdr);
        hts_close(out_fh);
    }

    // clean up
    if ( cnp_idx ) regidx_destroy(cnp_idx);
    if ( cnp_itr ) regitr_destroy(cnp_itr);
    bcf_sr_destroy(sr);
    free(sample);
    return 0;
}
