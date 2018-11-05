/* The MIT License

   Copyright (C) 2015-2018 Giulio Genovese

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

/****************************************
 * CONSTANT DEFINITIONS                 *
 ****************************************/

#define SIGN(x) (((x) > 0) - ((x) < 0))

#define MOCHA_VERSION "2018-11-05"

#define FLT_INCLUDE      (1<<0)
#define FLT_EXCLUDE      (1<<1)
#define WGS_DATA         (1<<2)
#define NO_LOG           (1<<3)
#define NO_ANNOT         (1<<4)
#define USE_SHORT_ARMS   (1<<5)
#define USE_CENTROMERES  (1<<6)

#define LRR 0
#define BAF 1
#define AD0 0
#define AD1 1
#define LDEV 0
#define BDEV 1

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

/****************************************
 * DATA STRUCTURES                      *
 ****************************************/

// structure defining regions of interest in the genome
typedef struct
{
    int *length;
    int *cen_beg;
    int *cen_end;
    int *is_short_arm;
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
} genome_t;

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
    int flags;
    genome_t genome;
    regidx_t *cnp_idx;
    regitr_t *cnp_itr;

    int rid;
    int n;
    int *pos_arr;
    int m_pos;
    float *gc_arr;
    int m_gc;
} model_t;

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

typedef struct
{
    int n, m;
    mocha_t *a;
} mocha_table_t;

typedef struct
{
    float lrr_median;
    float lrr_sd;
    float lrr_auto;
    float baf_sd;
    float baf_conc;
    float baf_auto;
    float coeffs[MAX_ORDER+1];
    float rel_ess;
} stats_t;

typedef struct
{
    int idx;
    int sex;
    float adjlrr_sd;
    int nhets;
    int x_nonpar_nhets;
    float x_nonpar_lrr_median;
    float y_nonpar_lrr_median;
    float mt_lrr_median;
    stats_t stats;
    stats_t *stats_arr;
    int m_stats, n_stats;

    int n;
    int *vcf_imap_arr;
    int m_vcf_imap;
    int16_t *data_arr[2];
    int m_data[2];
    int8_t *phase_arr;
    int m_phase;
} sample_t;

/****************************************
 * INLINE FUNCTIONS AND CONSTANTS       *
 ****************************************/

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
static const char *short_arm_chrs_default = "13,14,15,21,22,chr13,chr14,chr15,chr21,chr22";
const float lrr_bias_default = 0.2f;
const int min_dist_default = 400;
const int median_baf_adjust_default = 5;
const int order_lrr_gc_default = 2;

/****************************************
 * CONVERT FLOAT TO INT16 AND VICEVERSA *
 ****************************************/

#define INT16_SCALE 1000 // most BAF values from Illumina and UKBB are scaled to 1000

inline int16_t float_to_int16(float in)
{
    return isnan(in) ? bcf_int16_missing : (int16_t)roundf(INT16_SCALE * in);
}

inline float int16_to_float(int16_t in)
{
    return in == bcf_int16_missing ? NAN : ((float)in) / INT16_SCALE;
}

/******************************************
 * LRR AND COVERAGE POLYNOMIAL REGRESSION *
 ******************************************/

// the following alternative code snippets were considered to perform GC regression:
// https://stackoverflow.com/questions/5083465/fast-efficient-least-squares-fit-algorithm-in-c
// https://github.com/natedomin/polyfit/blob/master/polyfit.c
// https://software.intel.com/sites/products/documentation/doclib/mkl_sa/11/mkl_lapack_examples/sgels_ex.c.htm
// this function needs to use doubles internally when dealing with WGS data
static int polyfit(const float *lrr,
                   const float *gc,
                   int n,
                   const int *imap,
                   int order,
                   float *coeffs)
{
    int m = order + 1;
    if (n < m || order > MAX_ORDER) return -1;
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
            return -1;
        }
    }

    // calculate coefficients
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

    return 0;
}

static void ad_to_lrr_baf(const int16_t *ad0,
                           const int16_t *ad1,
                           float *lrr,
                           float *baf,
                           int n)
{
    // this function keeps a list of logarithms of integers to minimize log calls
    static float *logf_arr = NULL;
    static int n_logf = 0, m_logf = 0;
    if (ad0 == NULL && ad1 == NULL && n == 0) { free(logf_arr); return; }

    for(int i=0; i<n; i++)
    {
        int cov = (int)(ad0[i]==bcf_int16_missing ? 0 : ad0[i]) + (int)(ad1[i]==bcf_int16_missing ? 0 : ad1[i]);
        if (cov==0)
        {
            lrr[i] = NAN;
            baf[i] = NAN;
        }
        else
        {
            if (cov > n_logf)
            {
                hts_expand(float, cov, m_logf, logf_arr);
                for (int j=n_logf; j<cov; j++) logf_arr[j] = logf(j+1);
                n_logf = cov;

            }
            lrr[i] = logf_arr[cov-1];
            baf[i] = (ad0[i]==bcf_int16_missing || ad1[i]==bcf_int16_missing) ? NAN : (float)ad1[i] / (float)cov;
        }
    }
}

static void adjust_lrr(float *lrr,
                       const float *gc,
                       int n,
                       const int *imap,
                       const float *coeffs,
                       int order)
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
static int8_t *retrace_viterbi(int T,
                               int N,
                               const float *prb,
                               const int8_t *ptr)
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
static void rescale_prb(float *prb,
                        int N)
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
static int8_t *get_viterbi(const float *emis,
                           int T,
                           int m,
                           float xy_prob,
                           float flip_prob,
                           float telomere_prob,
                           float centromere_prob,
                           int last_p,
                           int first_q)
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

static inline float normal_pdf(float x,
                               float m,
                               float s,
                               float w)
{
    static const float log_inv_sqrt_2pi = -0.918938533204672669541f;
    float a = (x - m) / s;
    return expf( ( log_inv_sqrt_2pi - 0.5f * a * a ) * w );
}

// lrr_bias is used in a different way from what done by Petr Danecek in bcftools/vcfcnv.c
// TODO lrr_sd_bias is just a constant here, it could be removed
static inline float lkl_lrr_baf(float lrr,
                                float baf,
                                float ldev,
                                float bdev,
                                float lrr_sd,
                                float baf_sd,
                                float lrr_bias,
                                float lrr_sd_bias,
                                float err_prob)
{
    float ret = isnan( lrr ) ? 1.0f : normal_pdf( lrr, ldev, lrr_sd, lrr_bias ) / lrr_sd_bias;
    if ( !isnan(baf) ) ret *= ( normal_pdf( baf - 0.5f,  bdev, baf_sd, 1.0f ) +
                                normal_pdf( baf - 0.5f, -bdev, baf_sd, 1.0f ) ) * 0.5f / baf_sd;
    return ret < err_prob ? err_prob : ret;
}

static inline float lkl_baf_phase(float baf,
                                  int8_t phase,
                                  float bdev,
                                  float baf_sd,
                                  float err_prob)
{
    if ( isnan( baf ) ) return 1.0f;
    float ret;
    if ( phase == 0 ) ret = ( normal_pdf( baf - 0.5f,  bdev, baf_sd, 1.0f ) +
                              normal_pdf( baf - 0.5f, -bdev, baf_sd, 1.0f ) ) * 0.5f / baf_sd;
    else ret = normal_pdf( ( baf - 0.5f ) * (float)SIGN( phase ), bdev, baf_sd, 1.0f ) / baf_sd;
    return ret < err_prob ? err_prob : ret;
}

// precomupute emission probabilities
static float *get_lrr_baf_emis(const float *lrr,
                               const float *baf,
                               int T,
                               const int *imap,
                               float err_prob,
                               float lrr_bias,
                               float lrr_hap2dip,
                               float lrr_sd,
                               float baf_sd,
                               const float *cnf,
                               int m)
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
static float *get_baf_phase_emis(const float *baf,
                                 const int8_t *gt_phase,
                                 int T,
                                 const int *imap,
                                 float err_prob,
                                 float flip_prob,
                                 float baf_sd,
                                 const float *bdev,
                                 int m)
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
    int8_t *path = get_viterbi(emis, n, m, xy_prob, flip_prob, telomere_prob, 1.0f, 0, 0); // TODO can I not pass these values instead of 0 0?
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
                               double *lod_gamma_alpha,
                               double *lod_gamma_beta,
                               double *lod_gamma_alpha_beta,
                               float err_prob)
{
    float ret = isnan( lrr ) ? 1.0f : normal_pdf( lrr, ldev, lrr_sd, lrr_bias ) / lrr_sd_bias;
    if ( ad0!=bcf_int16_missing && ad1!=bcf_int16_missing )
        ret *= ( expf( lod_gamma_alpha[ad0] + lod_gamma_beta[ad1] - lod_gamma_alpha_beta[ad0 + ad1] ) +
                 expf( lod_gamma_alpha[ad1] + lod_gamma_beta[ad0] - lod_gamma_alpha_beta[ad0 + ad1] ) ) * 0.5f;
    return ret < err_prob ? err_prob : ret;
}

static inline float lkl_ad_phase(int16_t ad0,
                                 int16_t ad1,
                                 int8_t phase,
                                 double *lod_gamma_alpha,
                                 double *lod_gamma_beta,
                                 double *lod_gamma_alpha_beta,
                                 float err_prob)
{
    if ( ad0==bcf_int16_missing || ad1==bcf_int16_missing ) return 1.0f;
    float ret;
    if (phase == 0)
        ret = ( expf( lod_gamma_alpha[ad0] + lod_gamma_beta[ad1] - lod_gamma_alpha_beta[ad0 + ad1] ) +
                expf( lod_gamma_alpha[ad1] + lod_gamma_beta[ad0] - lod_gamma_alpha_beta[ad0 + ad1] ) ) * 0.5f;
    else if (phase>0) ret = expf( lod_gamma_alpha[ad0] + lod_gamma_beta[ad1] - lod_gamma_alpha_beta[ad0 + ad1] );
    else ret = expf( lod_gamma_alpha[ad1] + lod_gamma_beta[ad0] - lod_gamma_alpha_beta[ad0 + ad1] );
    return ret < err_prob ? err_prob : ret;
}

// precompute table of values for the beta binomial log likelihoods such that
// f(n, x) = log( \gamma(n+x) / \gamma(x) / n! )
// lod_gamma_alpha[n] = f(n, \alpha)
// lod_gamma_beta[n] = f(n, \beta)
// lod_gamma_alpha_beta[n] = f(n, \alpha + \beta)
// see https://en.wikipedia.org/wiki/Beta-binomial_distribution#As_a_compound_distribution
static void precompute_lod_gammas(const int16_t *ad0,
                                  const int16_t *ad1,
                                  int n,
                                  const int *imap,
                                  float alpha,
                                  float beta,
                                  double **lod_gamma_alpha,
                                  double **lod_gamma_beta,
                                  double **lod_gamma_alpha_beta)
{
    static double *internal_lod_gamma_alpha = NULL;
    static double *internal_lod_gamma_beta = NULL;
    static double *internal_lod_gamma_alpha_beta = NULL;
    static float internal_alpha = NAN, internal_beta = NAN;
    static int mlod_gamma_alpha = 0;
    static int mlod_gamma_beta = 0;
    static int mlod_gamma_alpha_beta = 0;
    static int n1 = 1, n2 = 1;
    if (ad0 == NULL && ad1 == NULL && n == 0 && imap == NULL)
    {
        free(lod_gamma_alpha);
        free(lod_gamma_beta);
        free(lod_gamma_alpha_beta);
        return;
    }

    if (internal_alpha != alpha || internal_beta != beta)
    {
        n1 = 1;
        n2 = 1;
        internal_alpha = alpha;
        internal_beta = beta;
    }

    int max1 = 0, max2 = 0;
    for (int i=0; i<n; i++)
    {
        int a = imap ? ad0[ imap[i] ] : ad0[i];
        int b = imap ? ad1[ imap[i] ] : ad1[i];
        if ( a!=bcf_int16_missing && b!=bcf_int16_missing )
        {
            if (a > max1) max1 = a;
            if (b > max1) max1 = b;
            if (a + b > max2) max2 = a + b;
        }
    }

    if (max1 + 1 > n1)
    {
        hts_expand(double, max1 + 1, mlod_gamma_alpha, internal_lod_gamma_alpha);
        hts_expand(double, max1 + 1, mlod_gamma_beta, internal_lod_gamma_beta);
        internal_lod_gamma_alpha[0] = 0.0f;
        internal_lod_gamma_beta[0] = 0.0f;
        while (n1 <= max1)
        {
            internal_lod_gamma_alpha[n1] = internal_lod_gamma_alpha[n1-1] + logf((alpha + (float)n1 - 1.0f) / (float)n1);
            internal_lod_gamma_beta[n1] = internal_lod_gamma_beta[n1-1] + logf((beta + (float)n1 - 1.0f) / (float)n1);
            n1++;
        }
    }

    if (max2 + 1 > n2)
    {
        hts_expand(double, max2 + 1, mlod_gamma_alpha_beta, internal_lod_gamma_alpha_beta);
        internal_lod_gamma_alpha_beta[0] = 0.0f;
        while (n2 <= max2)
        {
            internal_lod_gamma_alpha_beta[n2] = internal_lod_gamma_alpha_beta[n2-1] + logf((alpha + beta + n2 - 1) / n2);
            n2++;
        }
    }

    *lod_gamma_alpha = internal_lod_gamma_alpha;
    *lod_gamma_beta = internal_lod_gamma_beta;
    *lod_gamma_alpha_beta = internal_lod_gamma_alpha_beta;
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
    double *lod_gamma_alpha, *lod_gamma_beta, *lod_gamma_alpha_beta;
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
            emis[t*(1+m) + i] = lkl_lrr_ad( x, a, b, ldev[i], lrr_sd, lrr_bias, lrr_sd_bias, lod_gamma_alpha, lod_gamma_beta, lod_gamma_alpha_beta, err_prob );
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
    double *lod_gamma_alpha, *lod_gamma_beta, *lod_gamma_alpha_beta;
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
            emis[t*(1+2*m) + i] = lkl_ad_phase( a, b, p, lod_gamma_alpha, lod_gamma_beta, lod_gamma_alpha_beta, err_prob );
            if (i>0) emis[t*(1+2*m) + m + i] = lkl_ad_phase( a, b, p, lod_gamma_beta, lod_gamma_alpha, lod_gamma_alpha_beta, err_prob );
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
    double *lod_gamma_alpha, *lod_gamma_beta, *lod_gamma_alpha_beta;
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
        v[i] = lkl_lrr_ad( x, a, b, ldev, lrr_sd, lrr_bias, lrr_sd_bias, lod_gamma_alpha, lod_gamma_beta, lod_gamma_alpha_beta, err_prob );
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
    double *lod_gamma_alpha, *lod_gamma_beta, *lod_gamma_alpha_beta;
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
        v[i] = lkl_ad_phase( a, b, p, lod_gamma_alpha, lod_gamma_beta, lod_gamma_alpha_beta, err_prob );
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
    int8_t *path = get_viterbi(emis, n, m, xy_prob, flip_prob, telomere_prob, 1.0f, 0, 0); // TODO can I not pass these values instead of 0 0?
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
    double *lod_gamma_alpha, *lod_gamma_beta, *lod_gamma_alpha_beta;
    precompute_lod_gammas(ad0, ad1, n, NULL, alpha, alpha, &lod_gamma_alpha, &lod_gamma_beta, &lod_gamma_alpha_beta);
    float ret = 0.0f;
    for (int i=0; i<n; i++)
    {
        int16_t a = ad0[i];
        int16_t b = ad1[i];
        if ( a!=bcf_int16_missing && b!=bcf_int16_missing )
            ret += lod_gamma_alpha[a] + lod_gamma_beta[b] - lod_gamma_alpha_beta[a+b];
    }
    return (double)ret / M_LN10;
}

/*********************************
 * BASIC STATISTICS FUNCTIONS    *
 *********************************/

// iterator of non-NaN values
static inline float next_not_nan(const float *v,
                                 const int *imap,
                                 int n,
                                 int *i)
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
static float get_baf_conc(const float *baf,
                          const int8_t *gt_phase,
                          int n,
                          const int *imap)
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

// compute phased BAF autocorrelation for a float array with iterator
static float get_baf_auto_corr(const float *baf,
                               const int8_t *gt_phase,
                               int n,
                               const int *imap)
{
    double var = 0.0, auto_corr = 0.0;
    float prev = NAN, next = NAN;
    for (int i=0; i<n; i++)
    {
        next = ( ( imap ? baf[ imap[i] ] : baf[i] ) - 0.5f ) * ( imap ? gt_phase[ imap[i] ] : gt_phase[i] );
        if (!isnan(next))
        {
            var += sq((double)next);
            if (!isnan(prev)) auto_corr += prev * next;
            prev = next;
        }
    }
    auto_corr /= var;
    return auto_corr;
}

static float get_sample_sd(const float *v, int n, const int *imap);

// compute (adjusted) LRR autocorrelation for a float array with iterator
static float get_lrr_auto_corr(const float *lrr,
                               int n,
                               const int *imap)
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
static int get_n50(const int *v,
                   int n,
                   const int *imap)
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
// TODO this is not completely okay as the above way to compute the sd is error prone with floats
static float get_sample_sd(const float *v,
                           int n,
                           const int *imap)
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
static float get_se_mean(const float *v,
                         int n,
                         const int *imap)
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
float get_median(const float *v,
                 int n,
                 const int *imap)
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
static void push_rule(genome_t *genome,
                      char *line,
                      const bcf_hdr_t *hdr)
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
        if ( genome->cen_beg[rid] != 0 || genome->cen_end[rid] != 0 )
            error("Second centromere rule %s\n", line);
        genome->cen_beg[rid] = beg;
        genome->cen_end[rid] = end;
    }
    else if ( strncmp(ss, "X_nonpar", 8) == 0 )
    {
        if ( ( genome->x_xtr_beg != 0 || genome->x_xtr_end != 0 ) && genome->x_rid != rid )
            error("Chromosome X XTR and nonPAR regions declared on different contigs: %s\n", line);
        genome->x_rid = rid;
        if ( genome->x_nonpar_beg != 0 || genome->x_nonpar_end != 0 )
            error("Second chromosome X nonPAR rule: %s\n", line);
        genome->x_nonpar_beg = beg;
        genome->x_nonpar_end = end;
    }
    else if ( strncmp(ss, "X_xtr", 5) == 0 )
    {
        if ( ( genome->x_nonpar_beg != 0 || genome->x_nonpar_end != 0 ) && genome->x_rid != rid )
            error("Chromosome X nonPAR and XTR regions declared on different contigs: %s\n", line);
        genome->x_rid = rid;
        if ( genome->x_xtr_beg != 0 || genome->x_xtr_end != 0 )
            error("Second chromosome X XTR rule: %s\n", line);
        genome->x_xtr_beg = beg;
        genome->x_xtr_end = end;
    }
    else if ( strncmp(ss, "Y_nonpar", 8) == 0 )
    {
        if ( ( genome->y_xtr_beg != 0 || genome->y_xtr_end != 0 ) && genome->y_rid != rid )
            error("Chromosome Y XTR and nonPAR regions declared on different contigs: %s\n", line);
        genome->y_rid = rid;
        if ( genome->y_nonpar_beg != 0 || genome->y_nonpar_end != 0 )
            error("Second chromosome Y nonPAR rule: %s\n", line);
        genome->y_nonpar_beg = beg;
        genome->y_nonpar_end = end;
    }
    else if ( strncmp(ss, "Y_xtr", 5) == 0 )
    {
        if ( ( genome->y_nonpar_beg != 0 || genome->y_nonpar_end != 0 ) && genome->y_rid != rid )
            error("Chromosome Y nonPAR and XTR regions declared on different contigs: %s\n", line);
        genome->y_rid = rid;
        if ( genome->y_xtr_beg != 0 || genome->y_xtr_end != 0 )
            error("Second chromosome Y XTR rule: %s\n", line);
        genome->y_xtr_beg = beg;
        genome->y_xtr_end = end;
    }
    else if ( strncmp(ss, "mitochondria", 12) == 0 )
    {
        if ( genome->mt_rid != -1 )
            error("Second mitochondria rule %s\n", line);
        genome->mt_rid = rid;
    }
}

// split file into lines
// adapted from Petr Danecek's implementation of regidx_init_string() in bcftools/regidx.c
static void genome_init_file(genome_t *genome,
                             const char *fname,
                             const bcf_hdr_t *hdr)
{
    if ( !fname ) return;
    kstring_t tmp = {0, 0, NULL};
    htsFile *fp = hts_open(fname, "r");
    if ( !fp ) error("Failed to open %s: %s\n", fname, strerror(errno));
    while ( hts_getline(fp, KS_SEP_LINE, &tmp) > 0 )
        push_rule(genome, tmp.s, hdr);
    free(tmp.s);
    hts_close(fp);
}

// split string into lines
// adapted from Petr Danecek's implementation of regidx_init_string() in bcftools/regidx.c
static void genome_init_string(genome_t *genome,
                               const char *str,
                               const bcf_hdr_t *hdr)
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
        push_rule(genome, tmp.s, hdr);
        while ( *se && isspace(*se) ) se++;
        ss = se;
    }
    free(tmp.s);
}

// adapted from Petr Danecek's implementation of init_rules() in bcftools/regidx.c
static void genome_init_alias(genome_t *genome,
                              char *alias,
                              const bcf_hdr_t *hdr)
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
    return genome_init_string(genome, rules->rules, hdr);
}

static int readlist_short_arms(genome_t *genome,
                               const bcf_hdr_t *hdr,
                               const char *str)
{
    int n;
    char **list = hts_readlist(str, 0, &n);
    if (!list) return 0;
    for (int i=0; i<n; i++)
    {
        int rid = bcf_hdr_name2id(hdr, list[i]);
        free(list[i]);
        if (rid < 0) continue;
        genome->is_short_arm[rid] = 1;
    }
    free(list);
    return 1;
}

static int cnp_parse(const char *line,
                     char **chr_beg,
                     char **chr_end,
                     uint32_t *beg,
                     uint32_t *end,
                     void *payload,
                     void *usr)
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
 * SAMPLE METHODS                *
 *********************************/

static void mocha_print_ucsc(const mocha_t *mocha,
                             int n,
                             FILE *restrict stream,
                             const bcf_hdr_t *hdr)
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
                if ( strncmp(seq_name, "chr", 3) == 0 ) fprintf(stream, "%s\t%d\t%d\t%s\t%d\t.\t%d\t%d\t%s\n", seq_name,
                    mocha[j].beg_pos, mocha[j].end_pos, sample_name, (int)(1e3*mocha[j].cf), mocha[j].beg_pos, mocha[j].end_pos, color[i]);
                else fprintf(stream, "chr%s\t%d\t%d\t%s\t%d\t.\t%d\t%d\t%s\n", seq_name, mocha[j].beg_pos,
                    mocha[j].end_pos, sample_name, (int)(1e3*mocha[j].cf), mocha[j].beg_pos, mocha[j].end_pos, color[i]);
            }
        }
    }
    if (stream != stdout && stream != stderr) fclose(stream);
}

static void mocha_print(const mocha_t *mocha,
                        int n,
                        FILE *restrict stream,
                        const bcf_hdr_t *hdr,
                        int flags,
                        char *rules)
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
        fprintf(stream, "SAMPLE\tSEX\tCHROM\tBEG_%s\tEND_%s\tLENGTH\tP_ARM\tQ_ARM\tNSITES\tNHETS\tN50_HETS\tBDEV\tBDEV_SE\tREL_COV\tREL_COV_SE\tLOD_LRR_BAF\tLOD_BAF_PHASE\tNFLIPS\tBAF_CONC\tTYPE\tCF\n", rules, rules);
        for (int i=0; i<n; i++)
        {
            const char *sample_name = bcf_hdr_int2id(hdr, BCF_DT_SAMPLE, mocha->sample_idx);
            const char *seq_name = bcf_hdr_id2name(hdr, mocha->rid);
            fprintf(stream, "%s\t%c\t%s\t%d\t%d\t%d\t%c\t%c\t%d\t%d\t%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.2f\t%d\t%.4f\t%s\t%.4f\n",
                sample_name, sex[mocha->sex], seq_name, mocha->beg_pos, mocha->end_pos, mocha->length,
                arm_type[mocha->p_arm], arm_type[mocha->q_arm], mocha->nsites, mocha->nhets, mocha->n50_hets,
                mocha->bdev, mocha->bdev_se, 2.0f * expf(mocha->ldev), 2.0f * expf(mocha->ldev) * mocha->ldev_se,
                mocha->lod_lrr_baf, mocha->lod_baf_phase, mocha->nflips, mocha->baf_conc, type[mocha->type], mocha->cf);
            mocha++;
        }
    }
    else
    {
        fprintf(stream, "SAMPLE\tSEX\tCHROM\tBEG_%s\tEND_%s\tLENGTH\tP_ARM\tQ_ARM\tNSITES\tNHETS\tN50_HETS\tBDEV\tBDEV_SE\tLDEV\tLDEV_SE\tLOD_LRR_BAF\tLOD_BAF_PHASE\tNFLIPS\tBAF_CONC\tTYPE\tCF\n", rules, rules);
        for (int i=0; i<n; i++)
        {
            const char *sample_name = bcf_hdr_int2id(hdr, BCF_DT_SAMPLE, mocha->sample_idx);
            const char *seq_name = bcf_hdr_id2name(hdr, mocha->rid);
            fprintf(stream, "%s\t%c\t%s\t%d\t%d\t%d\t%c\t%c\t%d\t%d\t%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.2f\t%d\t%.4f\t%s\t%.4f\n",
                sample_name, sex[mocha->sex], seq_name, mocha->beg_pos, mocha->end_pos, mocha->length,
                arm_type[mocha->p_arm], arm_type[mocha->q_arm], mocha->nsites, mocha->nhets, mocha->n50_hets,
                mocha->bdev, mocha->bdev_se, mocha->ldev, mocha->ldev_se, mocha->lod_lrr_baf, mocha->lod_baf_phase,
                mocha->nflips, mocha->baf_conc, type[mocha->type], mocha->cf);
            mocha++;
        }
    }
    if (stream != stdout && stream != stderr) fclose(stream);
}

// this function checks whether the edge of a CNP is indeed diploid
// by looking for evidence for the diploid model over the CNP model (duplication or deletion)
static int cnp_edge_is_not_cn2(const float *lrr,
                               const float *baf,
                               int n,
                               int a,
                               int b,
                               float xy_prob,
                               float err_prob,
                               float lrr_bias,
                               float lrr_hap2dip,
                               float lrr_sd,
                               float baf_sd,
                               int type)
{
    float lrr_sd_bias = powf( lrr_sd, lrr_bias );
    float ldev, bdev;
    switch (type)
    {
        case MOCHA_CNP_DUP:
            ldev = log2f(1.5f) * lrr_hap2dip;
            bdev = 1.0f / 6.0f;
            break;
        case MOCHA_CNP_DEL:
            ldev = -lrr_hap2dip;
            bdev = 0.5f;
           break;
        default:
            error("CNP type must be either duplication or deletion\n");
    }

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
static int get_cnp_edges(const int *pos,
                         int n,
                         int beg,
                         int end,
                         int *a,
                         int *b)
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
static int8_t mocha_type(float ldev,
                         float ldev_se,
                         float bdev,
                         float bdev_se,
                         float lrr_hap2dip,
                         int8_t p_arm,
                         int8_t q_arm)
{
    float z2_upd = sqf( ldev / ldev_se ) * 0.5f;
    if (p_arm != MOCHA_TEL && q_arm != MOCHA_TEL) z2_upd += 2.0f; // penalty for not ending in a telomere
    if (ldev > 0)
    {
        if ( z2_upd > 5.0 ) return MOCHA_DUP;
        if ( !(isnan(bdev) || isnan(bdev_se)) )
        {
            float z2_dup = sqf( (ldev - 2.0f * M_LOG2E * lrr_hap2dip * bdev ) / ldev_se ) * 0.5f;
            if (z2_upd > z2_dup + 3.0f) return MOCHA_DUP;
            if (z2_dup > z2_upd + 3.0f) return MOCHA_UPD;
        }
    }
    else
    {
        if ( z2_upd > 5.0 ) return MOCHA_DEL;
        if ( !(isnan(bdev) || isnan(bdev_se)) )
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
static float mocha_cell_fraction(float ldev,
                                 float ldev_se,
                                 float bdev,
                                 float bdev_se,
                                 int8_t type,
                                 float lrr_hap2dip)
{
    if (isnan(bdev) || isnan(bdev_se))
    {
        switch (type)
        {
            case MOCHA_DEL:
                return -2.0f * (expf( ldev / lrr_hap2dip * (float)M_LN2 ) - 1.0f);
            case MOCHA_DUP:
                return 2.0f * (expf( ldev / lrr_hap2dip * (float)M_LN2 ) - 1.0f);
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

static void get_mocha_stats(const int *pos,
                            const float *lrr,
                            const float *baf,
                            const int8_t *gt_phase,
                            int n,
                            int a,
                            int b,
                            int cen_beg,
                            int cen_end,
                            int length,
                            mocha_t *mocha)
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
static int get_path_segs(const int8_t *path,
                         int n,
                         int except,
                         int **beg,
                         int **end,
                         int *nseg)
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
                // if the consecutive HMM states are non-zero and do not correspond to consecutive deletions and duplications
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
static void sample_run(sample_t *self,
                       mocha_table_t *mocha_table,
                       const model_t *model)
{
    // do nothing if chromosome Y or MT are being tested
    if ( model->rid == model->genome.y_rid || model->rid == model->genome.mt_rid ) return;

    mocha_t mocha;
    mocha.sample_idx = self->idx;
    mocha.sex = self->sex;
    mocha.rid = model->rid;

    int cen_beg = model->genome.cen_beg[model->rid];
    int cen_end = model->genome.cen_end[model->rid];
    int length = model->genome.length[model->rid];

    // declutter code by copying these values onto the stack
    int n = self->n;
    int8_t *gt_phase = self->phase_arr;
    int16_t *ad0 = self->data_arr[AD0];
    int16_t *ad1 = self->data_arr[AD1];
    float *lrr = (float *)malloc(n * sizeof(float));
    float *baf = (float *)malloc(n * sizeof(float));
    if ( model->flags & WGS_DATA )
    {
        ad_to_lrr_baf(ad0, ad1, lrr, baf, n);
    }
    else
    {
        for (int i=0; i<n; i++)
        {
            lrr[i] = int16_to_float(self->data_arr[LRR][i]);
            baf[i] = int16_to_float(self->data_arr[BAF][i]);
        }
    }

    if ( model->order_lrr_gc >= 0 )
        adjust_lrr(lrr, model->gc_arr, n, self->vcf_imap_arr, self->stats.coeffs, model->order_lrr_gc);

    int8_t *bdev_phase = (int8_t *)calloc(n, sizeof(int8_t));
    int *pos = (int *)malloc(n * sizeof(int));
    for (int i=0; i<n; i++) pos[i] = model->pos_arr[ self->vcf_imap_arr[i] ];
    int *imap_arr = (int *)malloc(n * sizeof(int));
    int *hets_imap_arr = (int *)malloc(n * sizeof(int));
    float *pbaf_arr = (float *)malloc(n * sizeof(float));

    int16_t *ldev = (int16_t *)calloc(n, sizeof(int16_t));
    int16_t *bdev = (int16_t *)calloc(n, sizeof(int16_t));

    // TODO do I need special normalization for the X nonPAR region
    if ( model->rid == model->genome.x_rid )
    {
        for (int i=0; i<n; i++)
        {
            if ( pos[i] > model->genome.x_nonpar_beg && pos[i] < model->genome.x_nonpar_end )
            {
                lrr[i] = ( self->sex == SEX_MAL ) ? NAN : lrr[i] - model->lrr_auto2sex;
                baf[i] = ( self->sex == SEX_MAL ) ? NAN : baf[i];
            }
        }
    }

    if ( model->cnp_itr )
    {
        while ( regitr_overlap( model->cnp_itr ) )
        {
            int a, b;
            if ( get_cnp_edges( pos, n, model->cnp_itr->beg, model->cnp_itr->end, &a, &b) == 0 )
            {
                float cn2_lod = log10_lkl_lrr_baf(lrr + a, baf + a, b + 1 - a, NULL, model->err_prob, model->lrr_bias, model->lrr_hap2dip, self->adjlrr_sd, self->stats.baf_sd, 2.0f);
                int cnp_type = regitr_payload(model->cnp_itr, int);
                mocha.type = MOCHA_UNK;
                mocha.ldev = get_median( lrr + a, b + 1 - a, NULL );
                if ( mocha.ldev > 0 && ( cnp_type == MOCHA_CNP_DUP || cnp_type == MOCHA_CNP_CNV ) )
                {
                    mocha.lod_lrr_baf = log10_lkl_lrr_baf(lrr + a, baf + a, b + 1 - a, NULL, model->err_prob, model->lrr_bias, model->lrr_hap2dip, self->adjlrr_sd, self->stats.baf_sd, 3.0f) - cn2_lod;
                    if ( mocha.lod_lrr_baf > -log10f( model->xy_prob ) )
                    {
                        mocha.type = MOCHA_CNP_DUP;
                        mocha.cf = NAN;
                    }
                }
                else if ( mocha.ldev <= 0 && ( cnp_type == MOCHA_CNP_DEL || cnp_type == MOCHA_CNP_CNV ) )
                {
                    mocha.lod_lrr_baf = log10_lkl_lrr_baf(lrr + a, baf + a, b + 1 - a, NULL, model->err_prob, model->lrr_bias, model->lrr_hap2dip, self->adjlrr_sd, self->stats.baf_sd, 1.0f) - cn2_lod;
                    if ( mocha.lod_lrr_baf > -log10f( model->xy_prob ) )
                    {
                        mocha.type = MOCHA_CNP_DEL;
                        mocha.cf = NAN;
                    }
                }
                if ( mocha.type == MOCHA_CNP_DUP || mocha.type == MOCHA_CNP_DEL )
                {
                    if ( cnp_edge_is_not_cn2(lrr, baf, n, a, b, model->xy_prob, model->err_prob, model->lrr_bias, model->lrr_hap2dip, self->adjlrr_sd, self->stats.baf_sd, mocha.type) ) continue;
                    get_mocha_stats( pos, lrr, baf, gt_phase, n, a, b, cen_beg, cen_end, length, &mocha);
                    // compute bdev, if possible
                    if ( mocha.nhets > 0 )
                    {
                        double f(double x, void *data) { return -log10_lkl_baf_phase(baf + a, gt_phase + a, b + 1 - a, NULL, NULL, model->err_prob, NAN, self->stats.baf_sd, x); }
                        double x;
                        kmin_brent(f, 0.0, 0.5, NULL, KMIN_EPS, &x);
                        mocha.bdev = fabsf((float)x);
                    }
                    else mocha.bdev = NAN;
                    mocha_table->n++;
                    hts_expand(mocha_t, mocha_table->n, mocha_table->m, mocha_table->a);
                    mocha_table->a[mocha_table->n - 1] = mocha;
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

    float *hs_arr = NULL;
    int n_hs = 0, m_hs = 0;
    for (int hmm_model=0; hmm_model<2; hmm_model++)
    {
        // select data to use from the contig, depending on which HMM model is being used
        int last_p = 0, first_q = 0;
        int n_imap = 0;
        for (int i=0; i<n; i++)
            if ( ( hmm_model == LRR_BAF && !isnan( lrr[i] ) ) || ( hmm_model == BAF_PHASE && !isnan( baf[i] ) ) )
            {
                if ( pos[ i ] < cen_beg ) last_p++;
                if ( pos[ i ] < cen_end ) first_q++;
                n_imap++;
                imap_arr[n_imap - 1] = i;
            }
        if ( n_imap == 0 ) continue;

        // compute emission probabilities and Viterbi path according to HMM model
        int except = 0;
        if ( hmm_model == LRR_BAF )
        {
            n_hs = model->cnf_n;
            hts_expand(float, n_hs, m_hs, hs_arr);
            for (int i=0; i<model->cnf_n; i++)
            {
                hs_arr[i] = model->cnf[i];
                if ( model->cnf[i] < 2.0f ) except++;
            }
        }
        else if ( hmm_model == BAF_PHASE )
        {
            n_hs = model->bdev_n;
            hts_expand(float, n_hs, m_hs, hs_arr);
            for (int i=0; i<model->bdev_n; i++)
                hs_arr[i] = model->bdev[i];
        }
        int8_t *path;
        int ret, *beg, *end, nseg;
        do
        {
            float *emis;
            if ( model->flags & WGS_DATA )
            {
                emis = hmm_model==LRR_BAF ? get_lrr_ad_emis( lrr, ad0, ad1, n_imap, imap_arr, model->err_prob, model->lrr_bias, model->lrr_hap2dip, self->adjlrr_sd, self->stats.baf_sd, hs_arr, n_hs )
                                          : get_ad_phase_emis( ad0, ad1, gt_phase, n_imap, imap_arr, model->err_prob, model->flip_prob, self->stats.baf_sd, hs_arr, n_hs );
            }
            else
            {
                emis = hmm_model==LRR_BAF ? get_lrr_baf_emis( lrr, baf, n_imap, imap_arr, model->err_prob, model->lrr_bias, model->lrr_hap2dip, self->adjlrr_sd, self->stats.baf_sd, hs_arr, n_hs )
                                          : get_baf_phase_emis( baf, gt_phase, n_imap, imap_arr, model->err_prob, model->flip_prob, self->stats.baf_sd, hs_arr, n_hs );
            }
            path = hmm_model==LRR_BAF ? get_viterbi(emis, n_imap, n_hs, model->xy_prob, NAN, model->telomere_prob, model->centromere_prob, last_p, first_q)
                                      : get_viterbi(emis, n_imap, n_hs, model->xy_prob, model->flip_prob, model->telomere_prob, model->centromere_prob, last_p, first_q);
            free(emis);

            ret = get_path_segs(path, n_imap, except, &beg, &end, &nseg);

            if ( ret ) // two consecutive hidden states were used, hinting that testing of a middle state might be necessary
            {
                free(path);
                n_hs++;
                hts_expand(float, n_hs, m_hs, hs_arr);
                hs_arr[n_hs - 1] = (hs_arr[ret-1] + hs_arr[ret]) * 0.5f;
                ks_introsort_float(n_hs, hs_arr);
            }
        }
        while ( ret );

        // loop through all the segments called by the Viterbi algorithm
        for (int i=0; i<nseg; i++)
        {
            // compute edges of the call
            int a = imap_arr[ beg[i] ];
            if ( beg[i] == 0 ) while ( a>0 && ldev[a-1]==0 && bdev[a-1]==0 ) a--; // extend call towards p telomere
            int b = imap_arr[ end[i] ];
            if ( end[i] == n_imap-1 ) while ( b<n-1 && ldev[b+1]==0 && bdev[b+1]==0 ) b++; // extend call towards q telomere
            mocha.ldev = get_median( lrr + a, b + 1 - a, NULL );
            get_mocha_stats( pos, lrr, baf, gt_phase, n, a, b, cen_beg, cen_end, length, &mocha);

            double f(double x, void *data)
            {
                if ( model->flags & WGS_DATA ) return -log10_lkl_lrr_ad(lrr + a, ad0 + a, ad1 + a, mocha.nsites, NULL, model->err_prob, model->lrr_bias, model->lrr_hap2dip, self->adjlrr_sd, self->stats.baf_sd, x);
                else return -log10_lkl_lrr_baf(lrr + a, baf + a, mocha.nsites, NULL, model->err_prob, model->lrr_bias, model->lrr_hap2dip, self->adjlrr_sd, self->stats.baf_sd, x);
            }
            double x, fx = kmin_brent(f, ( mocha.ldev > 0 ) ? 2.0 : 0, ( mocha.ldev > 0 ) ? 4.0 : 2.0, NULL, KMIN_EPS, &x);
            mocha.lod_lrr_baf = (float)(f(2.0, NULL) - fx);

            if ( hmm_model == LRR_BAF )
            {
                // here you need to check whether the call would have been better with the phased HMM model
                int n_hets_imap = 0;
                for (int j=beg[i]; j<=end[i]; j++)
                    if ( !isnan( baf[ imap_arr[j] ] ) )
                    {
                        n_hets_imap++;
                        hets_imap_arr[n_hets_imap - 1] = imap_arr[j];
                    }
                // TODO here it needs to pass information about the centromeres
                if ( model->flags & WGS_DATA )
                {
                    mocha.lod_baf_phase = compare_wgs_models(ad0, ad1, gt_phase, n_hets_imap, hets_imap_arr, model->xy_prob,
                        model->err_prob, model->flip_prob, model->telomere_prob, self->stats.baf_sd, model->bdev, model->bdev_n);
                }
                else
                {
                    mocha.lod_baf_phase = compare_models(baf, gt_phase, n_hets_imap, hets_imap_arr, model->xy_prob,
                        model->err_prob, model->flip_prob, model->telomere_prob, self->stats.baf_sd, model->bdev, model->bdev_n);
                }
                if (mocha.lod_baf_phase > mocha.lod_lrr_baf) continue;

                // compute bdev, if possible
                if ( n_hets_imap > 0 )
                {
                    double f(double x, void *data)
                    {
                        if ( model->flags & WGS_DATA ) return -log10_lkl_ad_phase(ad0, ad1, gt_phase, n_hets_imap, hets_imap_arr, NULL, model->err_prob, NAN, self->stats.baf_sd, x);
                        else return -log10_lkl_baf_phase(baf, gt_phase, n_hets_imap, hets_imap_arr, NULL, model->err_prob, NAN, self->stats.baf_sd, x);
                    }
                    fx = kmin_brent(f, 0.0, 0.5, NULL, KMIN_EPS, &x); // since f is symmetric around 0, it is important not to use an interval symmetric around 0
                    mocha.bdev = fabsf((float)x);
                }
                else mocha.bdev = NAN;
                mocha.bdev_se = NAN;
                for (int j=0; j<n_hets_imap; j++) bdev_phase[ hets_imap_arr[j] ] = (int8_t)SIGN( baf[ hets_imap_arr[j] ] - 0.5f );
            }
            else
            {
                // penalizes the LOD by the number of phase flips
                mocha.nflips = 0;
                for (int j=beg[i]; j<end[i]; j++) if ( path[j] != path[j+1] ) mocha.nflips++;

                double f(double x, void *data)
                {
                    if ( model->flags & WGS_DATA ) return -log10_lkl_ad_phase(ad0, ad1, gt_phase, mocha.nhets, imap_arr + beg[i], path + beg[i], model->err_prob, model->flip_prob, self->stats.baf_sd, x);
                    else return -log10_lkl_baf_phase(baf, gt_phase, mocha.nhets, imap_arr + beg[i], path + beg[i], model->err_prob, model->flip_prob, self->stats.baf_sd, x);
                }
                double x, fx = kmin_brent(f, -0.5f, 0.5f, NULL, KMIN_EPS, &x);
                mocha.bdev = fabsf((float)x);
                mocha.lod_baf_phase = (float)(f(0.0, NULL) - fx) + (float)mocha.nflips * log10f(model->flip_prob);

                for (int j=0; j<mocha.nhets; j++) pbaf_arr[ j ] = ( baf[ imap_arr[beg[i]+j] ] - 0.5f ) * (float)SIGN( path[beg[i]+j] );
                mocha.bdev_se = get_se_mean( pbaf_arr, mocha.nhets, NULL );
                for (int j=beg[i]; j<=end[i]; j++) bdev_phase[ imap_arr[j] ] = (int8_t)SIGN( path[j] ) * gt_phase[ imap_arr[j] ];
            }

            mocha.type = mocha_type(mocha.ldev, mocha.ldev_se, mocha.bdev, mocha.bdev_se, model->lrr_hap2dip, mocha.p_arm, mocha.q_arm);
            mocha.cf = mocha_cell_fraction(mocha.ldev, mocha.ldev_se, mocha.bdev, mocha.bdev_se, mocha.type, model->lrr_hap2dip);
            mocha_table->n++;
            hts_expand(mocha_t, mocha_table->n, mocha_table->m, mocha_table->a);
            mocha_table->a[mocha_table->n - 1] = mocha;

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

    // clean up
    free(hs_arr);
    free(imap_arr);
    free(hets_imap_arr);
    free(pbaf_arr);
    free(pos);
    free(lrr);
    free(baf);
    memcpy(self->data_arr[LDEV], ldev, n * sizeof(int16_t));
    memcpy(self->data_arr[BDEV], bdev, n * sizeof(int16_t));
    memcpy(self->phase_arr, bdev_phase, n * sizeof(int8_t));
    free(ldev);
    free(bdev);
    free(bdev_phase);
}

// computes the medoid contig for LRR regression
// TODO weight the coefficients appropriately
static int get_medoid(const float *coeffs,
                      int n,
                      int order)
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

// this function computes the median of contig stats
static float sample_summary(sample_t *self,
                            int n,
                            model_t *model)
{
    float *tmp_arr = (float *)malloc(n * sizeof(float));
    int m_tmp = n;

    for (int i=0; i<n; i++)
    {
        hts_expand(float, self[i].n_stats * (model->order_lrr_gc + 1), m_tmp, tmp_arr);
        for (int j=0; j<self[i].n_stats; j++) tmp_arr[j] = self[i].stats_arr[j].lrr_median;
        self[i].stats.lrr_median = get_median( tmp_arr, self[i].n_stats, NULL );
        for (int j=0; j<self[i].n_stats; j++) tmp_arr[j] = self[i].stats_arr[j].lrr_sd;
        self[i].stats.lrr_sd = get_median( tmp_arr, self[i].n_stats, NULL );
        for (int j=0; j<self[i].n_stats; j++) tmp_arr[j] = self[i].stats_arr[j].lrr_auto;
        self[i].stats.lrr_auto = get_median( tmp_arr, self[i].n_stats, NULL );
        for (int j=0; j<self[i].n_stats; j++) tmp_arr[j] = self[i].stats_arr[j].baf_sd;
        self[i].stats.baf_sd = get_median( tmp_arr, self[i].n_stats, NULL );
        for (int j=0; j<self[i].n_stats; j++) tmp_arr[j] = self[i].stats_arr[j].baf_conc;
        self[i].stats.baf_conc = get_median( tmp_arr, self[i].n_stats, NULL );
        for (int j=0; j<self[i].n_stats; j++) tmp_arr[j] = self[i].stats_arr[j].baf_auto;
        self[i].stats.baf_auto = get_median( tmp_arr, self[i].n_stats, NULL );

        self[i].adjlrr_sd = self[i].stats.lrr_sd;
        if ( model->order_lrr_gc == 0 )
        {
            for (int j=0; j<self[i].n_stats; j++) tmp_arr[j] = self[i].stats_arr[j].coeffs[0];
            self[i].stats.coeffs[0] = get_median( tmp_arr, self[i].n_stats, NULL );
        }
        else if ( model->order_lrr_gc > 0 && self[i].n_stats > 0 )
        {
            for (int j=0; j<self[i].n_stats; j++)
                for (int k=0; k<=model->order_lrr_gc; k++)
                    tmp_arr[j * (model->order_lrr_gc + 1) + k] = self[i].stats_arr[j].coeffs[k];
            int medoid_idx = get_medoid( tmp_arr, self[i].n_stats, model->order_lrr_gc );
            for (int k=0; k<=model->order_lrr_gc; k++)
                self[i].stats.coeffs[k] = tmp_arr[medoid_idx * (model->order_lrr_gc + 1) + k];
            for (int j=0; j<self[i].n_stats; j++) tmp_arr[j] = self[i].stats_arr[j].rel_ess;
            self[i].stats.rel_ess = get_median( tmp_arr, self[i].n_stats, NULL );
            self[i].adjlrr_sd *= sqrtf( 1.0f - self[i].stats.rel_ess ); // not perfect, but good enough(?)
        }
        free(self[i].stats_arr);
    }

    int j = 0;
    for (int i=0; i<n; i++)
        if( !isnan(self[i].x_nonpar_lrr_median) )
            tmp_arr[j++] = self[i].x_nonpar_lrr_median - self[i].stats.lrr_median;
    float cutoff = get_cutoff( tmp_arr, j );
    for (int i=0; i<n; i++)
    {
        if( self[i].x_nonpar_lrr_median - self[i].stats.lrr_median < cutoff ) self[i].sex = SEX_MAL;
        else if( self[i].x_nonpar_lrr_median - self[i].stats.lrr_median > cutoff ) self[i].sex = SEX_FEM;
    }

    if ( model->flags & WGS_DATA )
    {
        model->lrr_hap2dip = (float)M_LN2;
        model->lrr_auto2sex = 0.0f;
    }

    if ( isnan(model->lrr_hap2dip) || isnan(model->lrr_auto2sex) )
    {
        j = 0;
        for (int i=0; i<n; i++) if ( self[i].sex == SEX_MAL ) tmp_arr[j++] = self[i].x_nonpar_lrr_median - self[i].stats.lrr_median;
        float lrr_males = get_median( tmp_arr, j, NULL );
        j = 0;
        for (int i=0; i<n; i++) if ( self[i].sex == SEX_FEM ) tmp_arr[j++] = self[i].x_nonpar_lrr_median - self[i].stats.lrr_median;
        float lrr_females = get_median( tmp_arr, j, NULL );
        if ( isnan(model->lrr_hap2dip) ) model->lrr_hap2dip = lrr_females - lrr_males;
        if ( isnan(model->lrr_auto2sex) ) model->lrr_auto2sex = lrr_females;
    }
    free(tmp_arr);
    return cutoff;
}

// this function computes several contig stats and then clears the contig data from the sample
static void sample_stats(sample_t *self, const model_t *model)
{
    int n = self->n;
    if (n == 0) return;

    int16_t *ad0 = self->data_arr[AD0];
    int16_t *ad1 = self->data_arr[AD1];
    float *lrr = (float *)malloc(n * sizeof(float));
    float *baf = (float *)malloc(n * sizeof(float));
    if ( model->flags & WGS_DATA )
    {
        ad_to_lrr_baf(ad0, ad1, lrr, baf, n);
    }
    else
    {
        for (int i=0; i<n; i++)
        {
            lrr[i] = int16_to_float(self->data_arr[LRR][i]);
            baf[i] = int16_to_float(self->data_arr[BAF][i]);
        }
    }
    int *imap_arr = (int *)malloc(n * sizeof(int));

    if ( model->rid == model->genome.x_rid )
    {
        int n_imap = 0;
        for (int i=0; i<n; i++)
        {
            if ( !isnan( baf[i] ) ) self->nhets++;
            int pos = model->pos_arr[ self->vcf_imap_arr[i] ];
            if ( pos > model->genome.x_nonpar_beg && pos < model->genome.x_nonpar_end &&
                  ( pos < model->genome.x_xtr_beg || pos > model->genome.x_xtr_end ) )
            {
                if ( !isnan( baf[i] ) ) self->x_nonpar_nhets++;
                n_imap++;
                imap_arr[n_imap - 1] = i;
            }
        }
        self->x_nonpar_lrr_median = get_median( lrr, n_imap, imap_arr );
    }
    else if ( model->rid == model->genome.y_rid )
    {
        int n_imap = 0;
        for (int i=0; i<n; i++)
        {
            if ( !isnan( baf[i] ) ) self->nhets++;
            int pos = model->pos_arr[ self->vcf_imap_arr[i] ];
            if ( pos > model->genome.y_nonpar_beg && pos < model->genome.y_nonpar_end &&
                  ( pos < model->genome.y_xtr_beg || pos > model->genome.y_xtr_end ) )
            {
                n_imap++;
                imap_arr[n_imap - 1] = i;
            }
        }
        self->y_nonpar_lrr_median = get_median( lrr, n_imap, imap_arr );
    }
    else if ( model->rid == model->genome.mt_rid )
    {
        self->mt_lrr_median = get_median(lrr, n, NULL );
    }
    else
    {
        // expand arrays if necessary
        self->n_stats++;
        hts_expand(stats_t, self->n_stats, self->m_stats, self->stats_arr);

        if ( model->flags & WGS_DATA )
        {
            double f(double x, void *data) { return -log10_lkl_beta_binomial(ad0, ad1, n, x); }
            double x; kmin_brent(f, 0.0, 0.5, NULL, KMIN_EPS, &x);
            self->stats_arr[self->n_stats - 1].baf_sd = (float)x;
        }
        else
        {
            self->stats_arr[self->n_stats - 1].baf_sd = get_sample_sd( baf, n, NULL );
        }
        for (int i=0; i<n; i++) if ( !isnan( baf[i] ) ) self->nhets++;
        self->stats_arr[self->n_stats - 1].lrr_median = get_median(lrr, n, NULL );
        self->stats_arr[self->n_stats - 1].lrr_sd = get_sample_sd( lrr, n, NULL );

        self->stats_arr[self->n_stats - 1].baf_conc = get_baf_conc( baf, self->phase_arr, n, NULL );
        self->stats_arr[self->n_stats - 1].baf_auto = get_baf_auto_corr( baf, self->phase_arr, n, NULL );
        if ( model->order_lrr_gc == 0 )
        {
            self->stats_arr[self->n_stats - 1].coeffs[0] = get_median( lrr, n, NULL );
        }
        // performs polynomial regression for LRR
        else if ( model->order_lrr_gc > 0 )
        {
            float tss = get_tss(lrr, n);
            polyfit(lrr, model->gc_arr, n, self->vcf_imap_arr, model->order_lrr_gc, self->stats_arr[self->n_stats - 1].coeffs);
            adjust_lrr(lrr, model->gc_arr, n, self->vcf_imap_arr, self->stats_arr[self->n_stats - 1].coeffs, model->order_lrr_gc);
            self->stats_arr[self->n_stats - 1].coeffs[0] += get_median( lrr, n, NULL ); // further adjusts by median
            float rss = get_tss(lrr, n);
            self->stats_arr[self->n_stats - 1].rel_ess =  1.0f - rss / tss;
        }
        // compute autocorrelation after GC correction
        self->stats_arr[self->n_stats - 1].lrr_auto = get_lrr_auto_corr( lrr, n, NULL );
    }

    free(lrr);
    free(baf);
    free(imap_arr);
}

static void sample_print(const sample_t *self,
                         int n,
                         FILE *restrict stream,
                         const bcf_hdr_t *hdr,
                         int flags)
{
    if (stream == NULL) return;
    char sex[3];
    sex[SEX_UNK] = 'U';
    sex[SEX_MAL] = 'M';
    sex[SEX_FEM] = 'F';
    if ( flags & WGS_DATA )
    {
        fprintf(stream, "SAMPLE\tCOV_MEDIAN\tCOV_SD\tCOV_AUTO\tBAF_CORR\tBAF_CONC\tBAF_AUTO\tNHETS\tX_NONPAR_NHETS\tX_NONPAR_COV_MEDIAN\tY_NONPAR_COV_MEDIAN\tMT_COV_MEDIAN\tSEX\tREL_ESS\n");
        for (int i=0; i<n; i++)
        {
            const char *sample_name = bcf_hdr_int2id(hdr, BCF_DT_SAMPLE, self[i].idx);
            fprintf(stream, "%s\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%d\t%d\t%.4f\t%.4f\t%.4f\t%c\t%.4f\n", sample_name,
            expf(self[i].stats.lrr_median), expf(self[i].stats.lrr_median)*self[i].stats.lrr_sd, self[i].stats.lrr_auto, self[i].stats.baf_sd,
            self[i].stats.baf_conc, self[i].stats.baf_auto, self[i].nhets, self[i].x_nonpar_nhets, expf(self[i].x_nonpar_lrr_median),
            expf(self[i].y_nonpar_lrr_median), expf(self[i].mt_lrr_median), sex[self[i].sex], self[i].stats.rel_ess);
        }
    }
    else
    {
        fprintf(stream, "SAMPLE\tLRR_MEDIAN\tLRR_SD\tLRR_AUTO\tBAF_SD\tBAF_CONC\tBAF_AUTO\tNHETS\tX_NONPAR_NHETS\tX_NONPAR_LRR_MEDIAN\tY_NONPAR_LRR_MEDIAN\tMT_LRR_MEDIAN\tSEX\tREL_ESS\n");
        for (int i=0; i<n; i++)
        {
            const char *sample_name = bcf_hdr_int2id(hdr, BCF_DT_SAMPLE, self[i].idx);
            fprintf(stream, "%s\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%d\t%d\t%.4f\t%.4f\t%.4f\t%c\t%.4f\n", sample_name,
            self[i].stats.lrr_median, self[i].stats.lrr_sd, self[i].stats.lrr_auto, self[i].stats.baf_sd, self[i].stats.baf_conc, self[i].stats.baf_auto,
            self[i].nhets, self[i].x_nonpar_nhets, self[i].x_nonpar_lrr_median, self[i].y_nonpar_lrr_median,
            self[i].mt_lrr_median, sex[self[i].sex], self[i].stats.rel_ess);
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
static int put_contig(bcf_srs_t *sr,
                      const sample_t *sample,
                      const model_t *model,
                      htsFile *out_fh,
                      bcf_hdr_t *out_hdr)
{
    int rid = model->rid;
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
            while ( synced_iter[j] < sample[j].n - 1 && sample[j].vcf_imap_arr[ synced_iter[j] ] < i ) synced_iter[j]++;
            if ( sample[j].vcf_imap_arr[ synced_iter[j] ] == i )
            {
                if ( sample[j].data_arr[LDEV] ) ldev[ sample[j].idx ] = int16_to_float( sample[j].data_arr[LDEV][ synced_iter[j] ] );
                if ( sample[j].data_arr[BDEV] ) bdev[ sample[j].idx ] = int16_to_float( sample[j].data_arr[BDEV][ synced_iter[j] ] );
                if ( sample[j].phase_arr ) bdev_phase[ sample[j].idx ] = sample[j].phase_arr[ synced_iter[j] ];
            }
            else
            {
                // if no match variant found, match the end of the contig or keep conservative
                if ( i==0 && sample[j].data_arr[BDEV] )
                    bdev[ sample[j].idx ] = int16_to_float( sample[j].data_arr[BDEV][0] );
                if ( i==0 && sample[j].data_arr[LDEV] )
                    ldev[ sample[j].idx ] = int16_to_float( sample[j].data_arr[LDEV][0] );
                if ( sample[j].data_arr[BDEV] && int16_to_float( sample[j].data_arr[BDEV][ synced_iter[j] ] ) == 0.0f )
                    bdev[ sample[j].idx ] = 0.0f;
                if ( sample[j].data_arr[LDEV] && int16_to_float( sample[j].data_arr[LDEV][ synced_iter[j] ] ) == 0.0f )
                    ldev[ sample[j].idx ] = 0.0f;
                if ( sample[j].phase_arr ) bdev_phase[ sample[j].idx ] = 0;
            }
        }
        if ( !(model->flags & NO_ANNOT) )
        {
            bcf_update_format_float(out_hdr, line, "Ldev", ldev, (int)nsmpl);
            bcf_update_format_float(out_hdr, line, "Bdev", bdev, (int)nsmpl);
        }
        bcf_update_format_int32(out_hdr, line, "Bdev_Phase", bdev_phase, (int)nsmpl);

        if ( bcf_write(out_fh, out_hdr, line) < 0 ) error("Unable to write to output VCF file\n");
    }

    free(synced_iter);
    free(ldev);
    free(bdev);
    free(bdev_phase);

    return i;
}

// write header
static bcf_hdr_t *print_hdr(htsFile *out_fh,
                            bcf_hdr_t *hdr,
                            int argc,
                            char *argv[],
                            int record_cmd_line,
                            int flags)
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
    if ( bcf_hdr_write(out_fh, out_hdr) < 0 ) error("Unable to write to output VCF file\n");
    return out_hdr;
}

// retrieve phase information from BCF record
// bcf_int8_missing if phase does not apply
// 0 if phase is not available
// 1 if higher number allele received from the mother
// -1 if higher number allele received from the father
// assumes little endian architecture
int bcf_get_genotype_phase(bcf_fmt_t *fmt,
                           int8_t *gt_phase_arr,
                           int nsmpl)
{
    // bcf_fmt_t *fmt = bcf_get_fmt_id(line, id);
    if ( !fmt || fmt->n != 2 ) return 0;

    #define BRANCH(type_t, bcf_type_vector_end) { \
        type_t *p = (type_t *)fmt->p; \
        for (int i=0; i<nsmpl; i++, p+=2) \
        { \
            if ( p[0]==bcf_type_vector_end || bcf_gt_is_missing(p[0]) || \
                 p[1]==bcf_type_vector_end || bcf_gt_is_missing(p[1]) ) \
            { \
                gt_phase_arr[i] = bcf_int8_missing; \
            } \
            else \
            { \
                type_t gt0 = bcf_gt_allele(p[0]) > 0; \
                type_t gt1 = bcf_gt_allele(p[1]) > 0; \
                if ( gt0 == gt1 ) gt_phase_arr[i] = bcf_int8_missing; \
                else if ( !bcf_gt_is_phased(p[1]) ) gt_phase_arr[i] = 0; \
                else if ( gt1 > gt0 ) gt_phase_arr[i] = 1; \
                else gt_phase_arr[i] = -1; \
            } \
        } \
    }
    switch (fmt->type) {
        case BCF_BT_INT8:  BRANCH(int8_t, bcf_int8_vector_end); break;
        case BCF_BT_INT16: BRANCH(int16_t, bcf_int16_vector_end); break;
        case BCF_BT_INT32: BRANCH(int32_t, bcf_int32_vector_end); break;
        default: error("Unexpected type %d", fmt->type);
    }
    #undef BRANCH

    return 1;
}

// retrive genotype alleles information from BCF record
// assumes little endian architecture
int bcf_get_genotype_alleles(const bcf_fmt_t *fmt,
                             int16_t *gt0_arr,
                             int16_t *gt1_arr,
                             int nsmpl)
{
    if ( !fmt || fmt->n != 2 ) return 0;

    // temporarily store genotype alleles in AD array
    #define BRANCH(type_t, bcf_type_vector_end) { \
        type_t *p = (type_t *)fmt->p; \
        for (int i=0; i<nsmpl; i++, p+=2) \
        { \
            if ( p[0]==bcf_type_vector_end || bcf_gt_is_missing(p[0]) || \
                 p[1]==bcf_type_vector_end || bcf_gt_is_missing(p[1]) ) \
            { \
                gt0_arr[i] = bcf_int16_missing; \
                gt1_arr[i] = bcf_int16_missing; \
            } \
            else \
            { \
                gt0_arr[i] = (int16_t)bcf_gt_allele(p[0]); \
                gt1_arr[i] = (int16_t)bcf_gt_allele(p[1]); \
            } \
        } \
    }
    switch (fmt->type) {
        case BCF_BT_INT8:  BRANCH(int8_t, bcf_int8_vector_end); break;
        case BCF_BT_INT16: BRANCH(int16_t, bcf_int16_vector_end); break;
        case BCF_BT_INT32: BRANCH(int32_t, bcf_int32_vector_end); break;
        default: error("Unexpected type %d", fmt->type);
    }
    #undef BRANCH

    return 1;
}

// retrive allelic depth information from BCF record
// assumes little endian architecture
int bcf_get_allelic_depth(const bcf_fmt_t *fmt,
                          const int16_t *gt0_arr,
                          const int16_t *gt1_arr,
                          int16_t *ad0_arr,
                          int16_t *ad1_arr,
                          int nsmpl)
{
    if ( !fmt ) return 0;
    int nalleles = fmt->n;

    #define BRANCH(type_t, bcf_type_vector_end, bcf_type_missing) { \
        type_t *p = (type_t *)fmt->p; \
        for (int i=0; i<nsmpl; i++, p+=nalleles) \
        { \
            if ( ( gt0_arr[i] != bcf_int16_missing && (int)gt0_arr[i] >= nalleles ) || \
                 ( gt1_arr[i] != bcf_int16_missing && (int)gt1_arr[i] >= nalleles ) ) \
                error("Error: found VCF record with GT alleles %d and %d and %d number of alleles\n", gt0_arr[i], gt1_arr[i], nalleles); \
            if ( gt0_arr[i] == bcf_int16_missing || gt1_arr[i] == bcf_int16_missing || \
                 p[gt0_arr[i]]==bcf_type_vector_end || p[gt0_arr[i]]==bcf_type_missing || \
                 p[gt1_arr[i]]==bcf_type_vector_end || p[gt1_arr[i]]==bcf_type_missing ) \
            { \
                ad0_arr[i] = bcf_int16_missing; \
                ad1_arr[i] = bcf_int16_missing; \
            } \
            else \
            { \
                type_t gt0 = gt0_arr[i] > 0; \
                type_t gt1 = gt1_arr[i] > 0; \
                if ( gt0 == gt1 ) \
                { \
                    ad0_arr[i] = (int16_t)(gt0_arr[i] == gt1_arr[i] ? p[gt0_arr[i]] : p[gt0_arr[i]] + p[gt1_arr[i]]); \
                    ad1_arr[i] = bcf_int16_missing; \
                } \
                else \
                { \
                    ad0_arr[i] = (int16_t)p[gt0_arr[i]]; \
                    ad1_arr[i] = (int16_t)p[gt1_arr[i]]; \
                } \
            } \
        } \
    }
    switch (fmt->type) {
        case BCF_BT_INT8:  BRANCH(int8_t, bcf_int8_vector_end, bcf_int8_missing); break;
        case BCF_BT_INT16: BRANCH(int16_t, bcf_int16_vector_end, bcf_int16_missing); break;
        case BCF_BT_INT32: BRANCH(int32_t, bcf_int32_vector_end, bcf_int32_missing); break;
        default: error("Unexpected type %d", fmt->type);
    }
    #undef BRANCH

    return 1;
}

// read one contig
static int get_contig(bcf_srs_t *sr,
                      sample_t *sample,
                      model_t *model)
{
    int rid = model->rid;
    bcf_hdr_t *hdr = bcf_sr_get_header(sr, 0);
    bcf_sr_seek(sr, bcf_hdr_id2name( hdr, rid ), 0);

    bcf_fmt_t *baf_fmt = NULL, *lrr_fmt = NULL;
    bcf_info_t *info;
    int nsmpl = bcf_hdr_nsamples(hdr);

    int i;
    int gt_id = bcf_hdr_id2int(hdr, BCF_DT_ID, "GT");
    int baf_id = bcf_hdr_id2int(hdr, BCF_DT_ID, "BAF");
    int lrr_id = bcf_hdr_id2int(hdr, BCF_DT_ID, "LRR");
    int ad_id = bcf_hdr_id2int(hdr, BCF_DT_ID, "AD");
    int gc_id = bcf_hdr_id2int(hdr, BCF_DT_ID, "GC");

    // maybe this should be assert instead
    if ( gt_id < 0 ) error("Error: input VCF file has no GT format field\n");
    if ( (model->flags & WGS_DATA) && ad_id < 0 ) error("Error: input VCF file has no AD format field\n");
    if ( !(model->flags & WGS_DATA) && ( baf_id < 0 || lrr_id <0 ) ) error("Error: input VCF file has no BAF or LRR format field\n");

    int8_t *phase_arr = (int8_t *)malloc(nsmpl * sizeof(int8_t));
    int16_t *gt0 = (int16_t *)malloc(nsmpl * sizeof(int16_t));
    int16_t *gt1 = (int16_t *)malloc(nsmpl * sizeof(int16_t));
    int16_t *ad0 = (int16_t *)malloc(nsmpl * sizeof(int16_t));
    int16_t *ad1 = (int16_t *)malloc(nsmpl * sizeof(int16_t));
    int *last_het_pos = (int *)calloc(nsmpl, sizeof(int));
    int *last_pos = (int *)calloc(nsmpl, sizeof(int));

    model->n = 0;
    for (int j=0; j<nsmpl; j++) sample[j].n = 0;

    for (i=0; bcf_sr_next_line_reader0(sr); i++)
    {
        bcf1_t *line = bcf_sr_get_line(sr, 0);
        if ( rid != line->rid ) break;
        int pos = line->pos + 1;

        hts_expand(int, i+1, model->m_pos, model->pos_arr);
        model->pos_arr[i] = pos;
        hts_expand(float, i+1, model->m_gc, model->gc_arr);
        if ( gc_id>=0 && ( info = bcf_get_info_id( line, gc_id ) ) ) model->gc_arr[i] = info->v1.f;
        else model->gc_arr[i] = NAN;

        // if failing inclusion/exclusion requirement, skip line
        if ( ( model->flags & FLT_EXCLUDE ) && bcf_sr_get_line(sr, 1) ) continue;
        if ( ( model->flags & FLT_INCLUDE ) && !bcf_sr_get_line(sr, 1) ) continue;

        // if site falls in short arm or centromere regions skip line
        if ( !( model->flags & USE_SHORT_ARMS ) && model->genome.is_short_arm[rid] && pos < model->genome.cen_beg[rid] ) continue;
        if ( !( model->flags & USE_CENTROMERES ) && pos > model->genome.cen_beg[rid] && pos < model->genome.cen_end[rid] ) continue;

        // if there are no genotypes, skip line
        bcf_fmt_t *gt_fmt = bcf_get_fmt_id(line, gt_id);
        if ( !bcf_get_genotype_phase(gt_fmt, phase_arr, nsmpl) ) continue;

        // if neither AD nor LRR and BAF formats are present, skip line
        if ( model->flags & WGS_DATA )
        {
            if ( !bcf_get_genotype_alleles(gt_fmt, gt0, gt1, nsmpl) ) continue;
            if ( !bcf_get_allelic_depth(bcf_get_fmt_id(line, ad_id), gt0, gt1, ad0, ad1, nsmpl) ) continue;
        }
        else
        {
            if ( !(lrr_fmt = bcf_get_fmt_id(line, lrr_id)) || !(baf_fmt = bcf_get_fmt_id(line, baf_id)) ) continue;

            // make nan all BAF values for non heterozygous SNPs (we do not use those BAFs)
            int nbaf = 0;
            for (int j=0; j<nsmpl; j++)
            {
                if ( phase_arr[j] == bcf_int8_missing ) ((float *)(baf_fmt->p + baf_fmt->size * sample[j].idx))[0] = NAN;
                if ( !isnan(((float *)(baf_fmt->p + baf_fmt->size * sample[j].idx))[0]) ) nbaf++;
            }

            if ( model->median_baf_adjust >= 0 && nbaf >= model->median_baf_adjust )
            {
                float baf_median = get_median( (float *)baf_fmt->p, nsmpl, NULL ) - 0.5f;
                for (int j=0; j<nsmpl; j++)
                    ((float *)(baf_fmt->p + baf_fmt->size * sample[j].idx))[0] -= baf_median;
            }
        }

        // read line in memory
        model->n++;
        for (int j=0; j<nsmpl; j++)
        {
            if ( model->flags & WGS_DATA )
            {
                // site too close to last het site or hom site too close to last site
                if ( ( pos < last_het_pos[j] + model->min_dist ) ||
                   ( ( phase_arr[j] == bcf_int8_missing ) && ( pos < last_pos[j] + model->min_dist ) ) ) continue;

                // substitute the last hom site with the current het site
                if ( pos < last_pos[j] + model->min_dist ) sample[j].n--;

                if ( phase_arr[j] != bcf_int8_missing ) last_het_pos[j] = pos;
                last_pos[j] = pos;

                sample[j].n++;
                hts_expand(int, sample[j].n, sample[j].m_vcf_imap, sample[j].vcf_imap_arr);
                hts_expand(int8_t, sample[j].n, sample[j].m_phase, sample[j].phase_arr);
                hts_expand(int16_t, sample[j].n, sample[j].m_data[AD0], sample[j].data_arr[AD0]);
                hts_expand(int16_t, sample[j].n, sample[j].m_data[AD1], sample[j].data_arr[AD1]);
                sample[j].vcf_imap_arr[sample[j].n - 1] = i;
                sample[j].phase_arr[sample[j].n - 1] = phase_arr[j];
                sample[j].data_arr[AD0][sample[j].n - 1] = ad0[j];
                sample[j].data_arr[AD1][sample[j].n - 1] = ad1[j];

            }
            else
            {
                sample[j].n++;
                hts_expand(int, sample[j].n, sample[j].m_vcf_imap, sample[j].vcf_imap_arr);
                hts_expand(int8_t, sample[j].n, sample[j].m_phase, sample[j].phase_arr);
                hts_expand(int16_t, sample[j].n, sample[j].m_data[LRR], sample[j].data_arr[LRR]);
                hts_expand(int16_t, sample[j].n, sample[j].m_data[BAF], sample[j].data_arr[BAF]);
                sample[j].vcf_imap_arr[sample[j].n - 1] = i;
                sample[j].phase_arr[sample[j].n - 1] = phase_arr[j];
                float lrr = ((float *)(lrr_fmt->p + lrr_fmt->size * sample[j].idx))[0];
                float baf = ((float *)(baf_fmt->p + baf_fmt->size * sample[j].idx))[0];
                sample[j].data_arr[LRR][sample[j].n - 1] = float_to_int16(lrr);
                sample[j].data_arr[BAF][sample[j].n - 1] = float_to_int16(baf);

            }
        }
    }
    free(phase_arr);
    free(gt0);
    free(gt1);
    free(ad0);
    free(ad1);
    free(last_het_pos);
    free(last_pos);

    return i;
}

/*********************************
 * MAIN PART OF THE COMMAND      *
 *********************************/

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
    fprintf(stderr, "    -u, --ucsc-bed <file>             write UCSC bed track to a file [no output]\n");
    fprintf(stderr, "    -l  --no-log                      suppress progress report on standard error\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "HMM Options:\n");
    fprintf(stderr, "    -x, --xy-prob <float>             transition probability [%.0e]\n", xy_prob_default);
    fprintf(stderr, "    -z, --err-prob <float>            uniform error probability [%.0e]\n", err_prob_default);
    fprintf(stderr, "    -f, --flip-prob <float>           phase flip probability [%.0e]\n", flip_prob_default);
    fprintf(stderr, "    -t, --telomere-advantage <float>  telomere advantage [%.0e]\n", telomere_prob_default);
    fprintf(stderr, "    -c, --centromere-penalty <float>  centromere penalty [%.0e]\n", centromere_prob_default);
    fprintf(stderr, "        --short_arm_chrs <list>       list of chromosomes with short arms [%s]\n", short_arm_chrs_default);
    fprintf(stderr, "        --use_short_arms              use variants in short arms\n");
    fprintf(stderr, "        --use_centromeres             use variants in centromeres\n");
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
    fprintf(stderr, "    bcftools mocha -r GRCh37 input.bcf -v ^exclude.bcf -g stats.tsv -m mocha.tsv -p cnp.grch37.bed\n");
    fprintf(stderr, "    bcftools mocha -r GRCh38 input.bcf -Ob -o output.bcf -g stats.tsv -m mocha.tsv -n 1.0 --LRR-weight 0.5\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Publications:\n");
    fprintf(stderr, "    Loh P., Genovese G., McCarroll S., Price A. et al. Insights about clonal expansions\n");
    fprintf(stderr, "    from 8,342 mosaic chromosomal alterations. Nature 559, 350–355 (2018)\n");
    fprintf(stderr, "    [PubMed: http://www.ncbi.nlm.nih.gov/pubmed/29995854]\n");
    fprintf(stderr, "    [DOI: http://doi.org/10.1038/s41586-018-0321-x]\n");
    fprintf(stderr, "\n");
    exit(1);
}

static float *readlistf(const char *str,
                        int *n,
                        float min,
                        float max)
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
    int output_type = FT_VCF;
    int n_threads = 0;
    int record_cmd_line = 1;
    char *sample_names = NULL;
    char *output_fname = NULL;
    char *mocha_fname = NULL;
    char *stats_fname = NULL;
    char *ucsc_fname = NULL;
    char *cnp_fname = NULL;
    char *filter_fname = NULL;
    char *rules = NULL;
    sample_t *sample = NULL;
    FILE *out_fm = stdout;
    FILE *out_fg = NULL;
    FILE *out_fu = NULL;
    bcf_srs_t *sr = NULL;
    bcf_hdr_t *hdr = NULL;
    bcf_hdr_t *out_hdr = NULL;
    htsFile *out_fh = NULL;
    mocha_table_t mocha_table = {0, 0, NULL};
    const char *short_arm_chrs = short_arm_chrs_default;
    const char *cnf = cnf_default;
    const char *bdev = bdev_default;

    // model parameters
    model_t model;
    memset(&model, 0, sizeof(model_t));
    model.xy_prob = xy_prob_default;
    model.err_prob = err_prob_default;
    model.flip_prob = flip_prob_default;
    model.telomere_prob = telomere_prob_default;
    model.centromere_prob = centromere_prob_default;
    model.min_dist = min_dist_default;
    model.lrr_bias = lrr_bias_default;
    model.lrr_hap2dip = NAN;
    model.lrr_auto2sex = NAN;
    model.median_baf_adjust = median_baf_adjust_default;
    model.order_lrr_gc = order_lrr_gc_default;
    model.genome.x_rid = -1;
    model.genome.y_rid = -1;
    model.genome.mt_rid = -1;

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
        {"ucsc-bed", required_argument, NULL, 'u'},
        {"rules", required_argument, NULL, 'r'},
        {"rules-file", required_argument, NULL, 'R'},
        {"xy-prob", required_argument, NULL, 'x'},
        {"err-prob", required_argument, NULL, 'z'},
        {"flip-prob", required_argument, NULL, 'f'},
        {"telomere-advantage", required_argument, NULL, 't'},
        {"centromere-penalty", required_argument, NULL, 'c'},
        {"short_arm_chrs", required_argument, NULL, 10},
        {"use_short_arms", no_argument, NULL, 11},
        {"use_centromeres", no_argument, NULL, 12},
        {"cnf", required_argument, NULL, 'n'},
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
                    model.flags |= FLT_EXCLUDE;
                }
                else
                {
                    filter_fname = optarg;
                    model.flags |= FLT_INCLUDE;
                }
                break;
            case 'a': model.flags |= NO_ANNOT; break;
            case 'x':
                model.xy_prob = strtof(optarg, &tmp);
                if ( *tmp ) error("Could not parse: --xy-prob %s\n", optarg);
                break;
            case 'z':
                model.err_prob = strtof(optarg, &tmp);
                if ( *tmp ) error("Could not parse: --err-prob %s\n", optarg);
                break;
            case 'f':
                model.flip_prob = strtof(optarg, &tmp);
                if ( *tmp ) error("Could not parse: --flip-prob %s\n", optarg);
                break;
            case 't':
                model.telomere_prob = strtof(optarg, &tmp);
                if ( *tmp ) error("Could not parse: --telomere-advantage %s\n", optarg);
                break;
            case 'c':
                model.centromere_prob = strtof(optarg, &tmp);
                if ( *tmp ) error("Could not parse: --centromere-penalty %s\n", optarg);
                break;
            case 'd':
                model.min_dist = (int)strtol(optarg, &tmp, 0);
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
                model.lrr_hap2dip = strtof(optarg, &tmp);
                if ( *tmp ) error("Could not parse: --LRR-hap2dip %s\n", optarg);
                break;
            case  3 :
                model.lrr_auto2sex = strtof(optarg, &tmp);
                if ( *tmp ) error("Could not parse: --LRR-auto2sex %s\n", optarg);
                break;
            case  4 :
                model.lrr_bias = strtof(optarg, &tmp);
                if ( *tmp ) error("Could not parse: --LRR-weight %s\n", optarg);
                break;
            case  5 :
                model.median_baf_adjust = (int)strtol(optarg, &tmp, 0);
                if ( *tmp ) error("Could not parse: --median-BAF-adjust %s\n", optarg);
                break;
            case  6 :
                model.order_lrr_gc = (int)strtol(optarg, &tmp, 0);
                if ( *tmp ) error("Could not parse: --order-LRR-GC %s\n", optarg);
                break;
            case 10 : short_arm_chrs = optarg; break;
            case 11 : model.flags |= USE_SHORT_ARMS; break;
            case 12 : model.flags |= USE_CENTROMERES; break;
            case 'l': model.flags |= NO_LOG; break;
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

    if ( output_fname == NULL && ( model.flags & NO_ANNOT ) )
    {
        fprintf(stderr, "Option --no-annotations requires option --output\n");
        usage();
    }

    if (model.order_lrr_gc > MAX_ORDER)
    {
        fprintf(stderr, "Polynomial order must not be greater than %d: --order-LRR-GC %d\n", MAX_ORDER, model.order_lrr_gc);
        usage();
    }

    // parse parameters defining hidden states
    model.cnf = readlistf(cnf, &model.cnf_n, 0.0f, 4.0f);
    ks_introsort_float((size_t)model.cnf_n, model.cnf);
    model.bdev = readlistf(bdev, &model.bdev_n, 2.0f, INFINITY);
    for (int i=0; i<model.bdev_n; i++) model.bdev[i] = 1.0f / model.bdev[i]; // compute inverses
    ks_introsort_float((size_t)model.bdev_n, model.bdev);

    // output tables with mosaic chromosomal alteration calls
    if ( mocha_fname ) out_fm = get_file_handle( mocha_fname );
    if ( stats_fname ) out_fg = get_file_handle( stats_fname );
    if ( ucsc_fname ) out_fu = get_file_handle( ucsc_fname );

    // read list of regions to genotype
    if ( cnp_fname )
    {
        model.cnp_idx = regidx_init(cnp_fname, cnp_parse, NULL, sizeof(int), NULL);
        if ( !model.cnp_idx ) error("Error: failed to initialize CNP regions: --cnp %s\n", cnp_fname);
        model.cnp_itr = regitr_init(model.cnp_idx);
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
        model.flags |= WGS_DATA;
    else if ( (bcf_hdr_id2int( hdr, BCF_DT_ID, "LRR" ) < 0 || bcf_hdr_id2int( hdr, BCF_DT_ID, "BAF" ) < 0) )
        error("Error: input VCF file must contain either the AD format field or the LRR and BAF format fields\n");
    if ( model.order_lrr_gc > 0 && ( bcf_hdr_id2int( hdr, BCF_DT_ID, "GC" ) < 0 ) )
        error("Error: input VCF has no GC info field\n");

    // initialize genome parameters
    model.genome.length = (int *)calloc(hdr->n[BCF_DT_CTG], sizeof(int));
    for (int rid=0; rid<hdr->n[BCF_DT_CTG]; rid++)
        model.genome.length[rid] = hdr->id[BCF_DT_CTG][rid].val->info[0];
    model.genome.cen_beg = (int *)calloc(hdr->n[BCF_DT_CTG], sizeof(int));
    model.genome.cen_end = (int *)calloc(hdr->n[BCF_DT_CTG], sizeof(int));
    model.genome.is_short_arm = (int *)calloc(hdr->n[BCF_DT_CTG], sizeof(int));
    if ( rules_is_file )
        genome_init_file(&model.genome, rules, hdr);
    else
        genome_init_alias(&model.genome, rules, hdr);
    if ( !(model.flags & NO_LOG) ) fprintf(stderr, "Using genome assembly rules from %s\n", rules);
    readlist_short_arms(&model.genome, hdr, short_arm_chrs);

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
        if ( out_fh == NULL ) error("Cannot write to \"%s\": %s\n", output_fname, strerror(errno));
        if ( n_threads ) hts_set_opt(out_fh, HTS_OPT_THREAD_POOL, sr->p);
        out_hdr = print_hdr(out_fh, hdr, argc, argv, record_cmd_line, model.flags);
    }

    int nsmpl = bcf_hdr_nsamples(hdr);
    if ( nsmpl == 0 ) error("Subsetting has removed all samples\n");
    if ( !(model.flags & NO_LOG) ) fprintf(stderr, "Loading %d samples from the VCF file\n", nsmpl);
    if ( nsmpl < 30 && !(model.flags & WGS_DATA) && (model.median_baf_adjust >= 0) )
        fprintf(stderr, "Warning: median BAF adjustment with only %d samples is not recommended\n", nsmpl);

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
        model.rid = rid;
        int nret = get_contig(sr, sample, &model);
        if ( nret<=0 ) continue;
        if ( !(model.flags & NO_LOG) ) fprintf(stderr, "Read %d variants from contig %s\n", nret, bcf_hdr_id2name( hdr, rid ));
        if( model.genome.length[rid] < model.pos_arr[model.n-1] )
            model.genome.length[rid] = model.pos_arr[model.n-1];
        for (int j=0; j<nsmpl; j++) sample_stats(sample + j, &model);
    }

    if ( !(model.flags & NO_LOG) ) fprintf(stderr, "Compute and print sample statistics\n");
    float cutoff = sample_summary(sample, nsmpl, &model);
    if ( !(model.flags & NO_LOG) ) fprintf(stderr, "Estimated X nonPAR LRR cutoff for sex determination: %.4f\n", cutoff);
    int cnt[3] = {0, 0, 0}; for (int i=0; i<nsmpl; i++) cnt[sample[i].sex]++;
    if ( !(model.flags & NO_LOG) ) fprintf(stderr, "Estimated %d samples of unknown sex, %d males and %d females\n", cnt[SEX_UNK], cnt[SEX_MAL], cnt[SEX_FEM]);
    sample_print(sample, nsmpl, out_fg, hdr, model.flags);

    if ( isnan( model.lrr_hap2dip ) )
        error("Error: Unable to estimate LRR-hap2dip from X nonPAR. Make sure "
              "the regions are present in the VCF or specify the parameter\n");
    if ( isnan( model.lrr_auto2sex ) )
        error("Error: Unable to estimate LRR-auto2sex from autosomes and X nonPAR. Make sure "
              "the regions are present in the VCF or specify the parameter\n");

    if ( !(model.flags & NO_LOG) ) fprintf(stderr, "Estimated parameters: LRR-hap2dip=%.4f LRR-auto2sex=%.4f LRR-dip2trip=%.4f\n",
        model.lrr_hap2dip, model.lrr_auto2sex, model.lrr_hap2dip * (log2f(3.0f) - 1.0f));

    for (int rid=0; rid < hdr->n[BCF_DT_CTG]; rid++)
    {
        model.rid = rid;
        int nret = get_contig(sr, sample, &model);
        if ( nret<=0 ) continue;
        if ( !(model.flags & NO_LOG) ) fprintf(stderr, "Read %d variants from contig %s\n", nret, bcf_hdr_id2name( hdr, rid ));
        for (int j=0; j<nsmpl; j++)
        {
            if ( model.cnp_idx ) regidx_overlap(model.cnp_idx, bcf_hdr_id2name( hdr, rid ), 0, model.genome.length[rid], model.cnp_itr);
            sample_run(sample + j, &mocha_table, &model);
        }

        if (output_fname)
        {
            nret = put_contig(sr, sample, &model, out_fh, out_hdr);
            if ( !(model.flags & NO_LOG) ) fprintf(stderr, "Written %d variants for contig %s\n", nret, bcf_hdr_id2name( hdr, rid ));
        }
    }

    // clear sample data
    for (int j=0; j<nsmpl; j++)
    {
        free(sample[j].vcf_imap_arr);
        free(sample[j].data_arr[BDEV]);
        free(sample[j].data_arr[LDEV]);
        free(sample[j].phase_arr);
    }

    // clear model data
    free(model.pos_arr);
    free(model.gc_arr);
    free(model.cnf);
    free(model.bdev);
    free(model.genome.length);
    free(model.genome.cen_beg);
    free(model.genome.cen_end);
    free(model.genome.is_short_arm);

    // free precomputed tables
    ad_to_lrr_baf(NULL, NULL, NULL, NULL, 0);
    precompute_lod_gammas(NULL, NULL, 0, NULL, NAN, NAN, NULL, NULL, NULL);

    // write table with mosaic chromosomal alterations (and UCSC bed track)
    mocha_print(mocha_table.a, mocha_table.n, out_fm, hdr, model.flags, rules);
    mocha_print_ucsc(mocha_table.a, mocha_table.n, out_fu, hdr);
    free(mocha_table.a);

    // close output VCF
    if (output_fname)
    {
        bcf_hdr_destroy(out_hdr);
        hts_close(out_fh);
    }

    // clean up
    if ( model.cnp_idx ) regidx_destroy(model.cnp_idx);
    if ( model.cnp_itr ) regitr_destroy(model.cnp_itr);
    bcf_sr_destroy(sr);
    free(sample);
    return 0;
}
