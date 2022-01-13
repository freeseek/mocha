/* The MIT License

   Copyright (C) 2015-2022 Giulio Genovese

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
#include "mocha.h"
#include "genome_rules.h"
#include "beta_binom.h"
#include "bcftools.h"
#include "filter.h"
#include "tsv2vcf.h"

#define MOCHA_VERSION "2022-01-12"

/****************************************
 * CONSTANT DEFINITIONS                 *
 ****************************************/

#define SIGN(x) (((x) > 0) - ((x) < 0))

#define BDEV_LRR_BAF_DFLT "-2.0,-4.0,-6.0,10.0,6.0,4.0"
#define BDEV_BAF_PHASE_DFLT "6.0,8.0,10.0,15.0,20.0,30.0,50.0,80.0,130.0,210.0,340.0,550.0"
#define MIN_DST_DFLT "400"
#define ADJ_BAF_LRR_DFLT "5"
#define REGRESS_BAF_LRR_DFLT "15"
#define LRR_GC_ORDER_DFLT "2"
#define MAX_ORDER 5
#define XY_MAJOR_PL_DFLT "65.0"
#define XY_MINOR_PL_DFLT "35.0"
#define AUTO_TEL_PL_DFLT "20.0"
#define CHRX_TEL_PL_DFLT "8.0"
#define CHRY_TEL_PL_DFLT "6.0"
#define ERR_PL_DFLT "15.0"
#define FLIP_PL_DFLT "20.0"
#define SHORT_ARM_CHRS_DFLT "13,14,15,21,22,chr13,chr14,chr15,chr21,chr22"
#define LRR_BIAS_DFLT "0.2"
// https://www.illumina.com/documents/products/technotes/technote_cnv_algorithms.pdf
#define LRR_HAP2DIP_DFLT "0.45"

#define FLT_INCLUDE (1 << 0)
#define FLT_EXCLUDE (1 << 1)
#define WGS_DATA (1 << 2)
#define NO_LOG (1 << 3)
#define NO_ANNOT (1 << 4)
#define USE_SHORT_ARMS (1 << 5)
#define USE_CENTROMERES (1 << 6)
#define USE_MALES_XTR (1 << 7)
#define USE_MALES_PAR2 (1 << 8)
#define USE_NO_RULES_CHRS (1 << 9)

#define LRR 0
#define BAF 1
#define AD0 0
#define AD1 1
#define LDEV 0
#define BDEV 1

#define LRR_BAF 0
#define BAF_PHASE 1

#define MOCHA_UNDET 0
#define MOCHA_LOSS 1
#define MOCHA_GAIN 2
#define MOCHA_CNLOH 3
#define MOCHA_CNP_LOSS 4
#define MOCHA_CNP_GAIN 5
#define MOCHA_CNP_CNV 6

#define MOCHA_NOT 0
#define MOCHA_CEN 1
#define MOCHA_ARM 2
#define MOCHA_TEL 3

#define GT_NC 0
#define GT_AA 1
#define GT_AB 2
#define GT_BB 3

/****************************************
 * DATA STRUCTURES                      *
 ****************************************/

typedef struct {
    int pos;
    int allele_a;
    int allele_b;
    float adjust[3][3]; // shift
} locus_t;

typedef struct {
    int allele_a_id, allele_b_id, gc_id, gt_id, ad_id, baf_id, lrr_id;
    float xy_major_log_prb;
    float xy_minor_log_prb;
    float auto_tel_log_prb;
    float chrX_tel_log_prb;
    float chrY_tel_log_prb;
    float err_log_prb;
    float flip_log_prb;
    float *bdev_lrr_baf, *bdev_baf_phase;
    int bdev_lrr_baf_n, bdev_baf_phase_n;
    int min_dst;
    float lrr_cutoff;
    float lrr_hap2dip;
    float lrr_bias;
    int adj_baf_lrr;
    int regress_baf_lrr;
    int lrr_gc_order;
    int flags;
    int filter_logic;
    filter_t *filter;
    genome_rules_t *genome_rules;
    regidx_t *cnp_idx;
    regitr_t *cnp_itr;
    regidx_t *mhc_idx;
    regidx_t *kir_idx;

    int rid;
    int n;
    locus_t *locus_arr;
    int m_locus;
    float *gc_arr;
    int m_gc;
    int n_flipped;
} model_t;

typedef struct {
    int sample_idx;
    int computed_gender;
    int rid;
    int beg_pos;
    int end_pos;
    int length;
    int8_t p_arm;
    int8_t q_arm;
    int n_sites;
    int n_hets;
    int n50_hets;
    float bdev;
    float bdev_se;
    float ldev;
    float ldev_se;
    float lod_lrr_baf;
    float lod_baf_phase;
    int n_flips;
    float baf_conc;
    float lod_baf_conc;
    int8_t type;
    float cf;
} mocha_t;

typedef struct {
    int n, m;
    mocha_t *a;
} mocha_table_t;

typedef struct {
    float call_rate;
    float lrr_median;
    float lrr_sd;
    float lrr_auto;
    float dispersion; // either rho(AD0, AD1) for WGS model or sd(BAF)
    float baf_conc;
    float baf_auto;
    float lrr_gc_rel_ess;
    float coeffs[MAX_ORDER + 1];
} stats_t;

typedef struct {
    int idx;
    int computed_gender;
    float adjlrr_sd;
    int n_sites;
    int n_missing_gts;
    int n_hets;
    int x_nonpar_n_hets;
    int par1_n_hets;
    int xtr_n_hets;
    int par2_n_hets;
    float x_nonpar_dispersion; // either rho(AD0, AD1) for WGS model or sd(BAF)
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

static inline float sqf(float x) { return x * x; }
static inline double sq(double x) { return x * x; }
// the x == y is necessary in case x == -INFINITY
static inline float log_mean_expf(float x, float y) {
    return x == y ? x : (x > y ? x + logf(1.0f + expf(y - x)) : y + logf(1.0f + expf(x - y))) - (float)M_LN2;
}

static beta_binom_t *beta_binom_null, *beta_binom_alt;
static inline float beta_binom_log_lkl(const beta_binom_t *self, int16_t ad0, int16_t ad1) {
    return ad0 == bcf_int16_missing || ad1 == bcf_int16_missing ? 0.0f : beta_binom_log_unsafe(self, ad0, ad1);
}

/****************************************
 * CONVERT FLOAT TO INT16 AND VICEVERSA *
 ****************************************/

#define INT16_SCALE 1000 // BAF values from Illumina are scaled to 1000

static inline int16_t float_to_int16(float in) {
    return isnan(in) ? bcf_int16_missing : (int16_t)roundf(INT16_SCALE * in);
}

static inline float int16_to_float(int16_t in) { return in == bcf_int16_missing ? NAN : ((float)in) / INT16_SCALE; }

/******************************************
 * LRR AND COVERAGE POLYNOMIAL REGRESSION *
 ******************************************/

// https://stackoverflow.com/questions/5083465/fast-efficient-least-squares-fit-algorithm-in-c
// https://github.com/natedomin/polyfit/blob/master/polyfit.c
// https://software.intel.com/sites/products/documentation/doclib/mkl_sa/11/mkl_lapack_examples/sgels_ex.c.htm
// this function needs to use doubles internally when dealing with WGS data
static int polyfit(const float *lrr, const float *gc, int n, const int *imap, int order, float *coeffs) {
    int m = order + 1;
    if (n < m || order > MAX_ORDER) return -1;
    double B[MAX_ORDER + 1] = {0.0};
    double P[((MAX_ORDER + 1) * 2) + 1] = {0.0};
    double A[(MAX_ORDER + 1) * 2 * (MAX_ORDER + 1)] = {0.0};

    // identify the column vector
    for (int i = 0; i < n; i++) {
        float x = imap ? gc[imap[i]] : gc[i];
        float y = lrr[i];
        if (isnan(x) || isnan(y)) continue;
        float powx = 1.0f;

        for (int j = 0; j < m; j++) {
            B[j] += (double)(y * powx);
            powx *= x;
        }
    }

    // initialize the PowX array
    P[0] = (float)n;

    // compute the sum of the powers of X
    for (int i = 0; i < n; i++) {
        float x = imap ? gc[imap[i]] : gc[i];
        if (isnan(x)) continue;
        float powx = x;

        for (int j = 1; j < ((2 * m) + 1); j++) {
            P[j] += (double)powx;
            powx *= x;
        }
    }

    // initialize the reduction matrix
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < m; j++) {
            A[(i * (2 * m)) + j] = P[i + j];
        }

        A[(i * (2 * m)) + (i + m)] = 1.0;
    }

    // move the identity matrix portion of the redux matrix to the
    // left side (find the inverse of the left side of the redux matrix)
    for (int i = 0; i < m; i++) {
        double x = A[(i * (2 * m)) + i];
        if (x != 0) {
            for (int k = 0; k < (2 * m); k++) {
                A[(i * (2 * m)) + k] /= x;
            }

            for (int j = 0; j < m; j++) {
                if (i != j) {
                    double y = A[(j * (2 * m)) + i];
                    for (int k = 0; k < (2 * m); k++) {
                        A[(j * (2 * m)) + k] -= y * A[(i * (2 * m)) + k];
                    }
                }
            }
        } else {
            // cannot work with singular matrices
            return -1;
        }
    }

    // calculate coefficients
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < m; j++) {
            double x = 0.0;
            for (int k = 0; k < m; k++) {
                x += (A[(i * (2 * m)) + (k + m)] * B[k]);
            }
            coeffs[i] = (float)x;
        }
    }

    return 0;
}

static void ad_to_lrr_baf(const int16_t *ad0, const int16_t *ad1, float *lrr, float *baf, int n) {
    // this function keeps a list of logarithms of integers to minimize log calls
    static float *logf_arr = NULL;
    static int n_logf = 0, m_logf = 0;
    if (ad0 == NULL && ad1 == NULL && n == 0) {
        free(logf_arr);
        return;
    }

    for (int i = 0; i < n; i++) {
        if (ad0[i] == bcf_int16_missing && ad1[i] == bcf_int16_missing) {
            lrr[i] = NAN;
            baf[i] = NAN;
            continue;
        }
        int cov = (int)(ad0[i] == bcf_int16_missing ? 0 : ad0[i]) + (int)(ad1[i] == bcf_int16_missing ? 0 : ad1[i]);
        if (cov == 0) {
            lrr[i] = 0;
            baf[i] = NAN;
        } else {
            if (cov > n_logf) {
                hts_expand(float, cov, m_logf, logf_arr);
                for (int j = n_logf; j < cov; j++) logf_arr[j] = logf(j + 1);
                n_logf = cov;
            }
            lrr[i] = logf_arr[cov - 1];
            baf[i] = (ad0[i] == bcf_int16_missing || ad1[i] == bcf_int16_missing) ? NAN : (float)ad1[i] / (float)cov;
        }
    }
}

static void adjust_lrr(float *lrr, const float *gc, int n, const int *imap, const float *coeffs, int order) {
    for (int i = 0; i < n; i++) {
        float x = imap ? gc[imap[i]] : gc[i];
        float powx = 1.0f;
        for (int j = 0; j <= order; j++) {
            lrr[i] -= coeffs[j] * powx;
            powx *= x;
        }
    }
}

// computes total sum of squares
// this function needs to use doubles internally when dealing with WGS data
static float get_tss(const float *v, int n) {
    double mean = 0.0;
    int j = 0;
    for (int i = 0; i < n; i++) {
        if (!isnan(v[i])) {
            mean += (double)v[i];
            j++;
        }
    }
    if (j <= 1) return NAN;
    mean /= (double)j;

    double tss = 0.0;
    for (int i = 0; i < n; i++) {
        if (!isnan(v[i])) tss += sq((double)v[i] - mean);
    }
    return (float)tss;
}

/*********************************
 * HMM AND OPTIMIZATION METHODS  *
 *********************************/

// compute Viterbi path from probabilities
static int8_t *retrace_viterbi(int T, int N, const float *log_prb, const int8_t *ptr) {
    int i, t;
    int8_t *path = (int8_t *)malloc(T * sizeof(int8_t));

    // initialize last path state
    path[T - 1] = 0;
    for (i = 1; i < N; i++)
        if (log_prb[(int)path[T - 1]] < log_prb[i]) path[T - 1] = (int8_t)i;

    // compute best path by tracing back the Markov chain
    for (t = T - 1; t > 0; t--) path[t - 1] = ptr[(t - 1) * N + (int)path[t]];

    return path;
}

// rescale Viterbi log probabilities to avoid underflow issues
static inline void rescale_log_prb(float *log_prb, int n) {
    float max = -INFINITY;
    for (int i = 0; i < n; i++) max = max > log_prb[i] ? max : log_prb[i];
    for (int i = 0; i < n; i++) log_prb[i] -= max;
}

// compute the Viterbi path from BAF
// n is the length of the hidden Markov model
// m is the number of possible BAF deviations
static int8_t *log_viterbi_run(const float *emis_log_lkl, int T, int m, float xy_major_log_prb, float xy_minor_log_prb,
                               float tel_log_prb, float flip_log_prb, int last_p, int first_q) {
    int t, i, j, changeidx;

    // determine the number of hidden states based on whether phase information is used
    int N = 1 + m + (isnan(flip_log_prb) ? 0 : m);

    // allocate memory necessary for running the algorithm
    float *log_prb = (float *)malloc(N * sizeof(float));
    float *new_log_prb = (float *)malloc(N * sizeof(float));
    int8_t *ptr = (int8_t *)malloc(N * (T - 1) * sizeof(int8_t));
    int8_t *path;

    // initialize and rescale the first state
    log_prb[0] = emis_log_lkl[0];
    for (i = 1; i < N; i++) log_prb[i] = tel_log_prb + emis_log_lkl[i];
    rescale_log_prb(log_prb, N);

    // compute best probabilities at each position
    for (t = 1; t < T; t++) {
        // this causes a penalty for mosaic chromosomal calls across the centromeres
        float exit_log_prb = t > last_p ? xy_major_log_prb : xy_minor_log_prb;
        float enter_log_prb = t < first_q ? xy_major_log_prb : xy_minor_log_prb;

        for (i = 0; i < N; i++) {
            new_log_prb[i] = log_prb[i];
            ptr[(t - 1) * N + i] = (int8_t)i;
        }

        // compute whether a state switch should be considered for null state
        for (i = 1; i < N; i++) {
            if (new_log_prb[0] < log_prb[i] + exit_log_prb) {
                new_log_prb[0] = log_prb[i] + exit_log_prb;
                ptr[(t - 1) * N] = ptr[(t - 1) * N + i];
            }
            if (new_log_prb[i] < log_prb[0] + enter_log_prb) {
                new_log_prb[i] = log_prb[0] + enter_log_prb;
                ptr[(t - 1) * N + i] = ptr[(t - 1) * N];
            }
        }

        // compute whether a state switch should be considered for each other state
        // it will run twice if and only if phasing is used
        for (j = 0; j == 0 || (!isnan(flip_log_prb) && j == m); j += m) {
            float change_log_prb = log_prb[0] + enter_log_prb;
            changeidx = 0;
            for (i = 0; i < m; i++) {
                if (change_log_prb < log_prb[1 + j + i] + (xy_major_log_prb + xy_minor_log_prb) * 0.75f) {
                    change_log_prb = log_prb[1 + j + i] + (xy_major_log_prb + xy_minor_log_prb) * 0.75f;
                    changeidx = 1 + j + i;
                }
            }
            for (i = 0; i < m; i++) {
                if (new_log_prb[1 + j + i] < change_log_prb) {
                    new_log_prb[1 + j + i] = change_log_prb;
                    ptr[(t - 1) * N + 1 + j + i] = ptr[(t - 1) * N + changeidx];
                }
            }
        }

        // compute whether a phase flip should be considered for non-null states
        if (!isnan(flip_log_prb)) {
            for (i = 0; i < m; i++) {
                if (new_log_prb[1 + i] < new_log_prb[1 + m + i] + flip_log_prb) {
                    new_log_prb[1 + i] = new_log_prb[1 + m + i] + flip_log_prb;
                    ptr[(t - 1) * N + 1 + i] = ptr[(t - 1) * N + 1 + m + i];
                }
                if (new_log_prb[1 + m + i] < new_log_prb[1 + i] + flip_log_prb) {
                    new_log_prb[1 + m + i] = new_log_prb[1 + i] + flip_log_prb;
                    ptr[(t - 1) * N + 1 + m + i] = ptr[(t - 1) * N + 1 + i];
                }
            }
        }

        // update and rescale the current state
        new_log_prb[0] += emis_log_lkl[t * N];
        for (i = 0; i < m; i++) {
            new_log_prb[1 + i] += emis_log_lkl[t * N + 1 + i];
            if (!isnan(flip_log_prb)) new_log_prb[1 + m + i] += emis_log_lkl[t * N + 1 + m + i];
        }
        for (i = 0; i < N; i++) log_prb[i] = new_log_prb[i];
        rescale_log_prb(log_prb, N);
    }

    // add closing cost to the last state
    for (i = 1; i < N; i++) log_prb[i] += tel_log_prb;
    rescale_log_prb(log_prb, N);

    path = retrace_viterbi(T, N, log_prb, ptr);

    // free memory
    free(log_prb);
    free(new_log_prb);
    free(ptr);

    // symmetrize the path
    if (!isnan(flip_log_prb))
        for (i = 0; i < T; i++)
            if (path[i] > m) path[i] = (int8_t)m - path[i];

    return path;
}

/*********************************
 * LRR AND BAF LIKELIHOODS       *
 *********************************/

// rescale emission probabilities to avoid problems with outliers
// TODO find a better name for this function
static void rescale_emis_log_lkl(float *log_prb, int n, float err_log_prb) {
    float min_thr = -INFINITY;
    for (int i = 0; i < n; i++)
        if (min_thr < log_prb[i]) min_thr = log_prb[i];
    min_thr += err_log_prb;
    for (int i = 0; i < n; i++)
        if (log_prb[i] < min_thr) log_prb[i] = min_thr;
}

static inline float norm_log_lkl(float x, float m, float s, float w) {
    return isnan(x) ? 0.0f : -0.5f * sqf((x - m) / s) * w;
}

// lrr_bias is used in a different way from what done by Petr Danecek in bcftools/vcfcnv.c
static inline float lrr_baf_log_lkl(float lrr, float baf, float ldev, float bdev, float lrr_sd, float baf_sd,
                                    float lrr_bias) {
    return norm_log_lkl(lrr, ldev, lrr_sd, lrr_bias)
           + log_mean_expf(norm_log_lkl(baf - 0.5f, bdev, baf_sd, 1.0f), norm_log_lkl(baf - 0.5f, -bdev, baf_sd, 1.0f));
}

static inline float baf_phase_log_lkl(float baf, int8_t phase, float bdev, float baf_sd) {
    return phase == 0 ? log_mean_expf(norm_log_lkl(baf - 0.5f, bdev, baf_sd, 1.0f),
                                      norm_log_lkl(baf - 0.5f, -bdev, baf_sd, 1.0f))
                      : norm_log_lkl(baf - 0.5f, (float)SIGN(phase) * bdev, baf_sd, 1.0f);
}

// precomupute emission probabilities
static float *lrr_baf_emis_log_lkl(const float *lrr, const float *baf, int T, const int *imap, float err_log_prb,
                                   float lrr_bias, float lrr_hap2dip, float lrr_sd, float baf_sd,
                                   const float *bdev_lrr_baf_arr, int m) {
    float *ldev = (float *)malloc(m * sizeof(float));
    for (int i = 0; i < m; i++) ldev[i] = -logf(1.0f - 2.0f * bdev_lrr_baf_arr[i]) / (float)M_LN2 * lrr_hap2dip;
    int N = 1 + 2 * m;
    float *emis_log_lkl = (float *)malloc(N * T * sizeof(float));
    for (int t = 0; t < T; t++) {
        float x = imap ? lrr[imap[t]] : lrr[t];
        float y = imap ? baf[imap[t]] : baf[t];
        emis_log_lkl[t * N] = lrr_baf_log_lkl(x, y, 0.0f, 0.0f, lrr_sd, baf_sd, lrr_bias);
        for (int i = 0; i < m; i++) {
            emis_log_lkl[t * N + 1 + i] = lrr_baf_log_lkl(x, y, ldev[i], bdev_lrr_baf_arr[i], lrr_sd, baf_sd, lrr_bias);
        }
        // add states to distinguish LRR waves from true mosaic gains/losses
        for (int i = 0; i < m; i++) {
            if (bdev_lrr_baf_arr[i] < -1.0f / 6.0f || bdev_lrr_baf_arr[i] >= 1.0f / 6.0f)
                emis_log_lkl[t * N + 1 + m + i] = emis_log_lkl[t * N] + err_log_prb;
            else
                emis_log_lkl[t * N + 1 + m + i] =
                    lrr_baf_log_lkl(x, y, bdev_lrr_baf_arr[i], 0.0f, lrr_sd, baf_sd, lrr_bias);
        }
        rescale_emis_log_lkl(&emis_log_lkl[t * N], N, err_log_prb);
    }
    free(ldev);
    return emis_log_lkl;
}

// precomupute emission probabilities
static float *baf_phase_emis_log_lkl(const float *baf, const int8_t *gt_phase, int T, const int *imap,
                                     float err_log_prb, float baf_sd, const float *bdev, int m) {
    int N = 1 + 2 * m;
    float *emis_log_lkl = (float *)malloc(N * T * sizeof(float));
    for (int t = 0; t < T; t++) {
        float x = imap ? baf[imap[t]] : baf[t];
        int8_t p = imap ? gt_phase[imap[t]] : gt_phase[t];
        emis_log_lkl[t * N] = baf_phase_log_lkl(x, (int8_t)1, 0.0f, baf_sd);
        for (int i = 0; i < m; i++) {
            emis_log_lkl[t * N + 1 + i] = baf_phase_log_lkl(x, p, bdev[i], baf_sd);
            if (p == 0)
                emis_log_lkl[t * N + 1 + m + i] = emis_log_lkl[t * N + 1 + i];
            else
                emis_log_lkl[t * N + 1 + m + i] = baf_phase_log_lkl(x, p, -bdev[i], baf_sd);
        }
        rescale_emis_log_lkl(&emis_log_lkl[t * N], N, err_log_prb);
    }
    return emis_log_lkl;
}

static int cnp_edge_is_not_cn2_lrr_baf(const float *lrr, const float *baf, int n, int a, int b, float xy_log_prb,
                                       float err_log_prb, float lrr_bias, float lrr_hap2dip, float lrr_sd, float baf_sd,
                                       float ldev, float bdev) {
    // test left edge
    float sum_log_lkl = 0.0f;
    for (int i = a - 1; i >= 0; i--) {
        float log_lkl = lrr_baf_log_lkl(lrr[i], baf[i], ldev, bdev, lrr_sd, baf_sd, lrr_bias)
                        - lrr_baf_log_lkl(lrr[i], baf[i], 0.0f, 0.0f, lrr_sd, baf_sd, lrr_bias);
        if (!isnan(err_log_prb)) {
            if (log_lkl < err_log_prb)
                log_lkl = err_log_prb;
            else if (log_lkl > -err_log_prb)
                log_lkl = -err_log_prb;
        }
        sum_log_lkl += log_lkl;
        if (sum_log_lkl > -xy_log_prb) return -1;
        if (sum_log_lkl < xy_log_prb) break;
    }

    // test right edge
    sum_log_lkl = 0.0f;
    for (int i = b + 1; i < n; i++) {
        float log_lkl = lrr_baf_log_lkl(lrr[i], baf[i], ldev, bdev, lrr_sd, baf_sd, lrr_bias)
                        - lrr_baf_log_lkl(lrr[i], baf[i], 0, 0, lrr_sd, baf_sd, lrr_bias);
        if (!isnan(err_log_prb)) {
            if (log_lkl < err_log_prb)
                log_lkl = err_log_prb;
            else if (log_lkl > -err_log_prb)
                log_lkl = -err_log_prb;
        }
        sum_log_lkl += log_lkl;
        if (sum_log_lkl > -xy_log_prb) return -1;
        if (sum_log_lkl < xy_log_prb) break;
    }

    return 0;
}

// return the LOD likelihood for a segment
typedef struct {
    const float *lrr_arr;
    const float *baf_arr;
    int n;
    const int *imap;
    float err_log_prb;
    float lrr_bias;
    float lrr_hap2dip;
    float lrr_sd;
    float baf_sd;
} minus_lrr_baf_lod_t;

static double minus_lrr_baf_lod(double bdev_lrr_baf, void *ap) {
    minus_lrr_baf_lod_t *data = (minus_lrr_baf_lod_t *)ap;
    if (data->n == 0 || bdev_lrr_baf < -0.5 || bdev_lrr_baf > 0.25) return INFINITY; // kmin_brent does not handle NAN
    float ldev = -logf(1.0f - 2.0f * (float)bdev_lrr_baf) / (float)M_LN2 * data->lrr_hap2dip;
    float ret = 0.0f;
    for (int i = 0; i < data->n; i++) {
        float lrr = data->imap ? data->lrr_arr[data->imap[i]] : data->lrr_arr[i];
        float baf = data->imap ? data->baf_arr[data->imap[i]] : data->baf_arr[i];
        float log_lkl = lrr_baf_log_lkl(lrr, baf, ldev, (float)bdev_lrr_baf, data->lrr_sd, data->baf_sd, data->lrr_bias)
                        - lrr_baf_log_lkl(lrr, baf, 0.0f, 0.0f, data->lrr_sd, data->baf_sd, data->lrr_bias);
        if (!isnan(data->err_log_prb)) {
            if (log_lkl < data->err_log_prb)
                log_lkl = data->err_log_prb;
            else if (log_lkl > -data->err_log_prb)
                log_lkl = -data->err_log_prb;
        }
        ret += log_lkl;
    }
    return -(double)ret * M_LOG10E;
}

// return the LOD likelihood for a segment
typedef struct {
    const float *baf_arr;
    const int8_t *gt_phase;
    int n;
    const int *imap;
    const int8_t *as;
    float err_log_prb;
    float baf_sd;
} minus_baf_phase_lod_t;

static double minus_baf_phase_lod(double bdev, void *ap) {
    minus_baf_phase_lod_t *data = (minus_baf_phase_lod_t *)ap;
    if (data->n == 0 || bdev < 0.0 || bdev > 0.5) return INFINITY; // kmin_brent does not handle NAN

    float ret = 0.0f;
    for (int i = 0; i < data->n; i++) {
        float baf = data->imap ? data->baf_arr[data->imap[i]] : data->baf_arr[i];
        int8_t p = data->imap ? data->gt_phase[data->imap[i]] : data->gt_phase[i];
        if (data->as) p *= (int8_t)SIGN(data->as[i]); // notice as has no imap
        float log_lkl =
            baf_phase_log_lkl(baf, p, (float)bdev, data->baf_sd) - baf_phase_log_lkl(baf, 0, 0.0f, data->baf_sd);
        if (!isnan(data->err_log_prb)) {
            if (log_lkl < data->err_log_prb)
                log_lkl = data->err_log_prb;
            else if (log_lkl > -data->err_log_prb)
                log_lkl = -data->err_log_prb;
        }
        ret += log_lkl;
    }
    return -(double)ret * M_LOG10E;
}

// TODO find a better title for this function
static float compare_models(const float *baf, const int8_t *gt_phase, int n, const int *imap, float xy_log_prb,
                            float err_log_prb, float flip_log_prb, float tel_log_prb, float baf_sd, const float *bdev,
                            int m) {
    if (n == 0) return NAN;
    float *emis_log_lkl = baf_phase_emis_log_lkl(baf, gt_phase, n, imap, err_log_prb, baf_sd, bdev, m);
    int8_t *path = log_viterbi_run(emis_log_lkl, n, m, xy_log_prb, xy_log_prb, tel_log_prb, flip_log_prb, 0, 0);
    free(emis_log_lkl);
    int n_flips = 0;
    for (int i = 1; i < n; i++)
        if (path[i - 1] && path[i] && path[i - 1] != path[i]) n_flips++;
    minus_baf_phase_lod_t data = {baf, gt_phase, n, imap, path, err_log_prb, baf_sd};
    double x, fx = kmin_brent(minus_baf_phase_lod, 0.1, 0.2, (void *)&data, KMIN_EPS, &x);
    free(path);
    return -(float)fx + (float)n_flips * flip_log_prb * (float)M_LOG10E;
}

typedef struct {
    const float *baf_arr;
    int n;
    const int *imap;
    float baf_sd;
} minus_baf_log_lkl_t;

static double minus_baf_log_lkl(double bdev, void *ap) {
    minus_baf_log_lkl_t *data = (minus_baf_log_lkl_t *)ap;
    if (data->n == 0 || bdev < 0.0 || bdev > 0.5) return INFINITY; // kmin_brent does not handle NAN

    double ret = 0.0;
    for (int i = 0; i < data->n; i++) {
        float baf = data->imap ? data->baf_arr[data->imap[i]] : data->baf_arr[i];
        if (isnan(baf)) continue;
        float log_lkl = log_mean_expf(norm_log_lkl(baf - 0.5f, (float)bdev, data->baf_sd, 1.0f),
                                      norm_log_lkl(baf - 0.5f, -(float)bdev, data->baf_sd, 1.0f));
        ret += (double)log_lkl;
    }
    return -ret * M_LOG10E;
}

static float get_sample_mean(const float *v, int n, const int *imap) {
    float mean = 0.0;
    int j = 0;
    for (int i = 0; i < n; i++) {
        float tmp = imap ? v[imap[i]] : v[i];
        if (!isnan(tmp)) {
            mean += tmp;
            j++;
        }
    }
    if (j <= 1) return NAN;
    return mean /= (float)j;
}

static float get_baf_bdev(const float *baf_arr, int n, const int *imap, float baf_sd) {
    double bdev = 0.0;
    int j = 0;
    for (int i = 0; i < n; i++) {
        float baf = imap ? baf_arr[imap[i]] : baf_arr[i];
        if (isnan(baf)) continue;
        bdev += fabs((double)baf - 0.5);
        j++;
    }
    if (j == 0) return NAN;
    bdev /= j;
    // simple method to compute bdev should work well for germline duplications
    if ((float)bdev > 2.0f * baf_sd) return (float)bdev;
    minus_baf_log_lkl_t data = {baf_arr, n, imap, baf_sd};
    kmin_brent(minus_baf_log_lkl, 0.1, 0.2, (void *)&data, KMIN_EPS, &bdev);
    return (float)bdev < 1e-4 ? (float)NAN : (float)bdev;
}

static void get_max_sum(const int16_t *ad0, const int16_t *ad1, int n, const int *imap, int *n1, int *n2) {
    *n1 = 0;
    *n2 = 0;
    for (int i = 0; i < n; i++) {
        int a = imap ? ad0[imap[i]] : ad0[i];
        int b = imap ? ad1[imap[i]] : ad1[i];
        if (a != bcf_int16_missing && b != bcf_int16_missing) {
            if (a > *n1) *n1 = a;
            if (b > *n1) *n1 = b;
            if (a + b > *n2) *n2 = a + b;
        }
    }
}

typedef struct {
    const int16_t *ad0_arr;
    const int16_t *ad1_arr;
    int n;
    const int *imap;
    float ad_rho;
} minus_ad_log_lkl_t;

static double minus_ad_log_lkl(double bdev, void *ap) {
    minus_ad_log_lkl_t *data = (minus_ad_log_lkl_t *)ap;
    if (data->n == 0 || bdev < 0.0 || bdev > 0.5) return INFINITY; // kmin_brent does not handle NAN

    int n1, n2;
    get_max_sum(data->ad0_arr, data->ad1_arr, data->n, data->imap, &n1, &n2);
    beta_binom_update(beta_binom_alt, 0.5f + (float)bdev, data->ad_rho, n1, n2);

    double ret = 0.0;
    for (int i = 0; i < data->n; i++) {
        int16_t ad0 = data->imap ? data->ad0_arr[data->imap[i]] : data->ad0_arr[i];
        int16_t ad1 = data->imap ? data->ad1_arr[data->imap[i]] : data->ad1_arr[i];
        float log_lkl =
            log_mean_expf(beta_binom_log_lkl(beta_binom_alt, ad0, ad1), beta_binom_log_lkl(beta_binom_alt, ad1, ad0));
        ret += (double)log_lkl;
    }
    return -(double)ret * M_LOG10E;
}

static float get_ad_bdev(const int16_t *ad0_arr, const int16_t *ad1_arr, int n, const int *imap, float ad_rho) {
    double bdev = 0.0;
    minus_ad_log_lkl_t data = {ad0_arr, ad1_arr, n, imap, ad_rho};
    kmin_brent(minus_ad_log_lkl, 0.1, 0.2, (void *)&data, KMIN_EPS, &bdev);
    return (float)bdev < 1e-4 ? (float)NAN : (float)bdev;
}

/*********************************
 * WGS AD LIKELIHOODS            *
 *********************************/

static inline float lrr_ad_log_lkl(float lrr, int16_t ad0, int16_t ad1, float ldev, float lrr_sd, float lrr_bias,
                                   const beta_binom_t *beta_binom) {
    return norm_log_lkl(lrr, ldev, lrr_sd, lrr_bias)
           + log_mean_expf(beta_binom_log_lkl(beta_binom, ad0, ad1), beta_binom_log_lkl(beta_binom, ad1, ad0));
}

static inline float ad_phase_log_lkl(int16_t ad0, int16_t ad1, int8_t phase, const beta_binom_t *beta_binom) {
    return phase == 0
               ? log_mean_expf(beta_binom_log_lkl(beta_binom, ad0, ad1), beta_binom_log_lkl(beta_binom, ad1, ad0))
               : (phase > 0 ? beta_binom_log_lkl(beta_binom, ad0, ad1) : beta_binom_log_lkl(beta_binom, ad1, ad0));
}

static float *lrr_ad_emis_log_lkl(const float *lrr, const int16_t *ad0, const int16_t *ad1, int T, const int *imap,
                                  float err_log_prb, float lrr_bias, float lrr_hap2dip, float lrr_sd, float ad_rho,
                                  const float *bdev_lrr_baf_arr, int m) {
    int N = 1 + 2 * m;
    int n1, n2;
    get_max_sum(ad0, ad1, T, NULL, &n1, &n2);
    float *emis_log_lkl = (float *)malloc(N * T * sizeof(float));
    for (int i = 0; i < 1 + m; i++) {
        float ldev = i == 0 ? 0.0f : -logf(1.0f - 2.0f * bdev_lrr_baf_arr[i - 1]) / (float)M_LN2 * lrr_hap2dip;
        float bdev = i == 0 ? 0.0f : fabsf(bdev_lrr_baf_arr[i - 1]);
        beta_binom_t *beta_binom = i == 0 ? beta_binom_null : beta_binom_alt;
        beta_binom_update(beta_binom, 0.5f + bdev, ad_rho, n1, n2);

        for (int t = 0; t < T; t++) {
            float x = imap ? lrr[imap[t]] : lrr[t];
            int16_t a = imap ? ad0[imap[t]] : ad0[t];
            int16_t b = imap ? ad1[imap[t]] : ad1[t];
            emis_log_lkl[t * N + i] = lrr_ad_log_lkl(x, a, b, ldev, lrr_sd, lrr_bias, beta_binom);
        }
    }
    // generate states that should attract LRR waves with no BAF signal
    for (int i = 0; i < m; i++) {
        // do not make extreme states compete with LRR waves
        if (bdev_lrr_baf_arr[i] < -1.0f / 6.0f || bdev_lrr_baf_arr[i] >= 1.0f / 6.0f) {
            for (int t = 0; t < T; t++) emis_log_lkl[t * N + 1 + m + i] = emis_log_lkl[t * N] + err_log_prb;
        } else {
            float ldev = -logf(1.0f - 2.0f * bdev_lrr_baf_arr[i]) / (float)M_LN2 * lrr_hap2dip;
            for (int t = 0; t < T; t++) {
                float x = imap ? lrr[imap[t]] : lrr[t];
                int16_t a = imap ? ad0[imap[t]] : ad0[t];
                int16_t b = imap ? ad1[imap[t]] : ad1[t];
                emis_log_lkl[t * N + 1 + m + i] = lrr_ad_log_lkl(x, a, b, ldev, lrr_sd, lrr_bias, beta_binom_null);
            }
        }
    }
    for (int t = 0; t < T; t++) rescale_emis_log_lkl(&emis_log_lkl[t * N], N, err_log_prb);
    return emis_log_lkl;
}

static float *ad_phase_emis_log_lkl(const int16_t *ad0, const int16_t *ad1, const int8_t *gt_phase, int T,
                                    const int *imap, float err_log_prb, float ad_rho, const float *bdev_arr, int m) {
    int N = 1 + 2 * m;
    float *emis_log_lkl = (float *)malloc(N * T * sizeof(float));
    for (int i = 0; i < 1 + m; i++) {
        float bdev = i == 0 ? 0.0f : bdev_arr[i - 1];
        // TODO this function should come out of the loop
        int n1, n2;
        get_max_sum(ad0, ad1, T, imap, &n1, &n2);
        beta_binom_t *beta_binom = i == 0 ? beta_binom_null : beta_binom_alt;
        beta_binom_update(beta_binom, 0.5f + bdev, ad_rho, n1, n2);

        for (int t = 0; t < T; t++) {
            int16_t a = imap ? ad0[imap[t]] : ad0[t];
            int16_t b = imap ? ad1[imap[t]] : ad1[t];
            int8_t p = imap ? gt_phase[imap[t]] : gt_phase[t];
            emis_log_lkl[t * N + i] = ad_phase_log_lkl(a, b, p, beta_binom);
            if (i > 0) emis_log_lkl[t * N + m + i] = ad_phase_log_lkl(b, a, p, beta_binom);
        }
    }
    for (int t = 0; t < T; t++) rescale_emis_log_lkl(&emis_log_lkl[t * N], N, err_log_prb);
    return emis_log_lkl;
}

static int cnp_edge_is_not_cn2_lrr_ad(const float *lrr, int16_t *ad0, int16_t *ad1, int n, int a, int b,
                                      float xy_log_prb, float err_log_prb, float lrr_bias, float lrr_hap2dip,
                                      float lrr_sd, float ad_rho, float ldev, float bdev) {
    int n1, n2;
    get_max_sum(ad0, ad1, n, NULL, &n1, &n2);
    beta_binom_update(beta_binom_null, 0.5f, ad_rho, n1, n2);
    beta_binom_update(beta_binom_alt, 0.5f + bdev, ad_rho, n1, n2);

    // test left edge
    float sum_log_lkl = 0.0f;
    for (int i = a - 1; i >= 0; i--) {
        float log_lkl = lrr_ad_log_lkl(lrr[i], ad0[i], ad1[i], ldev, lrr_sd, lrr_bias, beta_binom_alt)
                        - lrr_ad_log_lkl(lrr[i], ad0[i], ad1[i], 0.0f, lrr_sd, lrr_bias, beta_binom_null);
        if (!isnan(err_log_prb)) {
            if (log_lkl < err_log_prb)
                log_lkl = err_log_prb;
            else if (log_lkl > -err_log_prb)
                log_lkl = -err_log_prb;
        }
        sum_log_lkl += log_lkl;
        if (sum_log_lkl > -xy_log_prb) return -1;
        if (sum_log_lkl < xy_log_prb) break;
    }

    // test right edge
    sum_log_lkl = 0.0f;
    for (int i = b + 1; i < n; i++) {
        float log_lkl = lrr_ad_log_lkl(lrr[i], ad0[i], ad1[i], ldev, lrr_sd, lrr_bias, beta_binom_alt)
                        - lrr_ad_log_lkl(lrr[i], ad0[i], ad1[i], 0.0f, lrr_sd, lrr_bias, beta_binom_null);
        if (!isnan(err_log_prb)) {
            if (log_lkl < err_log_prb)
                log_lkl = err_log_prb;
            else if (log_lkl > -err_log_prb)
                log_lkl = -err_log_prb;
        }
        sum_log_lkl += log_lkl;
        if (sum_log_lkl > -xy_log_prb) return -1;
        if (sum_log_lkl < xy_log_prb) break;
    }

    return 0;
}

// return the LOD likelihood for a segment
typedef struct {
    const float *lrr_arr;
    const int16_t *ad0_arr;
    const int16_t *ad1_arr;
    int n;
    const int *imap;
    float err_log_prb;
    float lrr_bias;
    float lrr_hap2dip;
    float lrr_sd;
    float ad_rho;
} minus_lrr_ad_lod_t;

static double minus_lrr_ad_lod(double bdev_lrr_baf, void *ap) {
    minus_lrr_ad_lod_t *data = (minus_lrr_ad_lod_t *)ap;
    if (data->n == 0 || bdev_lrr_baf < -0.5 || bdev_lrr_baf > 0.25) return INFINITY; // kmin_brent does not handle NAN

    float ldev = -logf(1.0f - 2.0f * (float)bdev_lrr_baf) / (float)M_LN2 * data->lrr_hap2dip;
    int n1, n2;
    get_max_sum(data->ad0_arr, data->ad1_arr, data->n, data->imap, &n1, &n2);
    beta_binom_update(beta_binom_null, 0.5f, data->ad_rho, n1, n2);
    beta_binom_update(beta_binom_alt, 0.5f + (float)bdev_lrr_baf, data->ad_rho, n1, n2);
    float ret = 0.0f;
    for (int i = 0; i < data->n; i++) {
        float lrr = data->imap ? data->lrr_arr[data->imap[i]] : data->lrr_arr[i];
        int16_t ad0 = data->imap ? data->ad0_arr[data->imap[i]] : data->ad0_arr[i];
        int16_t ad1 = data->imap ? data->ad1_arr[data->imap[i]] : data->ad1_arr[i];
        float log_lkl = lrr_ad_log_lkl(lrr, ad0, ad1, ldev, data->lrr_sd, data->lrr_bias, beta_binom_alt)
                        - lrr_ad_log_lkl(lrr, ad0, ad1, 0.0f, data->lrr_sd, data->lrr_bias, beta_binom_null);
        if (!isnan(data->err_log_prb)) {
            if (log_lkl < data->err_log_prb)
                log_lkl = data->err_log_prb;
            else if (log_lkl > -data->err_log_prb)
                log_lkl = -data->err_log_prb;
        }
        ret += log_lkl;
    }
    return -(double)ret * M_LOG10E;
}

// return the LOD likelihood for a segment
typedef struct {
    const int16_t *ad0_arr;
    const int16_t *ad1_arr;
    const int8_t *gt_phase;
    int n;
    const int *imap;
    const int8_t *as;
    float err_log_prb;
    float ad_rho;
} minus_ad_phase_lod_t;

static double minus_ad_phase_lod(double bdev, void *ap) {
    minus_ad_phase_lod_t *data = (minus_ad_phase_lod_t *)ap;
    if (data->n == 0 || bdev < 0.0 || bdev > 0.5) return INFINITY; // kmin_brent does not handle NAN

    int n1, n2;
    get_max_sum(data->ad0_arr, data->ad1_arr, data->n, data->imap, &n1, &n2);
    beta_binom_update(beta_binom_null, 0.5f, data->ad_rho, n1, n2);
    beta_binom_update(beta_binom_alt, 0.5f + (float)bdev, data->ad_rho, n1, n2);
    float ret = 0.0f;
    for (int i = 0; i < data->n; i++) {
        int16_t ad0 = data->imap ? data->ad0_arr[data->imap[i]] : data->ad0_arr[i];
        int16_t ad1 = data->imap ? data->ad1_arr[data->imap[i]] : data->ad1_arr[i];
        int8_t p = data->imap ? data->gt_phase[data->imap[i]] : data->gt_phase[i];
        if (data->as) p *= (int8_t)SIGN(data->as[i]); // notice as has no imap
        float log_lkl = ad_phase_log_lkl(ad0, ad1, p, beta_binom_alt) - ad_phase_log_lkl(ad0, ad1, 0, beta_binom_null);
        if (!isnan(data->err_log_prb)) {
            if (log_lkl < data->err_log_prb)
                log_lkl = data->err_log_prb;
            else if (log_lkl > -data->err_log_prb)
                log_lkl = -data->err_log_prb;
        }
        ret += log_lkl;
    }
    return -(double)ret * M_LOG10E;
}

// TODO find a better title for this function
static float compare_wgs_models(const int16_t *ad0, const int16_t *ad1, const int8_t *gt_phase, int n, const int *imap,
                                float xy_log_prb, float err_log_prb, float flip_log_prb, float tel_log_prb,
                                float ad_rho, const float *bdev, int m) {
    if (n == 0) return NAN;
    float *emis_log_lkl = ad_phase_emis_log_lkl(ad0, ad1, gt_phase, n, imap, err_log_prb, ad_rho, bdev, m);
    int8_t *path = log_viterbi_run(emis_log_lkl, n, m, xy_log_prb, xy_log_prb, tel_log_prb, flip_log_prb, 0,
                                   0); // TODO can I not pass these values instead of 0 0?
    free(emis_log_lkl);
    int n_flips = 0;
    for (int i = 1; i < n; i++)
        if (path[i - 1] && path[i] && path[i - 1] != path[i]) n_flips++;
    minus_ad_phase_lod_t data = {ad0, ad1, gt_phase, n, imap, path, err_log_prb, ad_rho};
    double x, fx = kmin_brent(minus_ad_phase_lod, 0.1, 0.2, (void *)&data, KMIN_EPS, &x);
    free(path);
    return -(float)fx + (float)n_flips * flip_log_prb * (float)M_LOG10E;
}

// TODO change this or integrate with ad_lod
typedef struct {
    const int16_t *ad0_arr;
    const int16_t *ad1_arr;
    int n;
    const int *imap;
} minus_lod_lkl_beta_binomial_t;

static double minus_lod_lkl_beta_binomial(double ad_rho, void *ap) {
    minus_lod_lkl_beta_binomial_t *data = (minus_lod_lkl_beta_binomial_t *)ap;
    if (data->n == 0 || ad_rho <= 0.0 || ad_rho >= 1.0) return INFINITY;
    float ret = 0.0f;
    int n1, n2;
    get_max_sum(data->ad0_arr, data->ad1_arr, data->n, data->imap, &n1, &n2);
    beta_binom_update(beta_binom_null, 0.5f, ad_rho, n1, n2);
    for (int i = 0; i < data->n; i++) {
        int16_t ad0 = data->imap ? data->ad0_arr[data->imap[i]] : data->ad0_arr[i];
        int16_t ad1 = data->imap ? data->ad1_arr[data->imap[i]] : data->ad1_arr[i];
        ret += beta_binom_log_lkl(beta_binom_null, ad0, ad1);
    }
    return -(double)ret * M_LOG10E;
}

/*********************************
 * BASIC STATISTICS FUNCTIONS    *
 *********************************/

// iterator of non-NaN values
static inline float next_not_missing(const float *v, const int *imap, int n, int *i) {
    float x = NAN;
    while (*i < n) {
        x = imap ? v[imap[*i]] : v[*i];
        if (!isnan(x)) break;
        (*i)++;
    }
    return x;
}

// compute BAF phase concordance for a float array with iterator
static void get_baf_conc(const float *baf, const int8_t *gt_phase, int n, int *conc, int *disc) {
    int i;
    float prev = NAN, next = NAN;
    *conc = 0, *disc = 0;
    for (i = 0; i < n; i++) {
        prev = (baf[i] - 0.5f) * (float)gt_phase[i];
        if (!isnan(prev) && prev != 0.0f) break;
    }
    if (i == n) return;

    for (i++; i < n; i++) {
        next = (baf[i] - 0.5f) * (float)gt_phase[i];
        if (!isnan(next) && next != 0.0f) {
            if (prev * next > 0.0f)
                (*conc)++;
            else if (prev * next < 0.0f)
                (*disc)++;
            prev = next;
        }
    }
}

// compute phased BAF autocorrelation for a float array with iterator
static float get_baf_auto_corr(const float *baf, const int8_t *gt_phase, int n) {
    double var = 0.0, auto_corr = 0.0;
    float prev = NAN, next = NAN;
    for (int i = 0; i < n; i++) {
        next = (baf[i] - 0.5f) * (float)gt_phase[i];
        if (!isnan(next)) {
            var += sq((double)next);
            if (!isnan(prev)) auto_corr += prev * next;
            prev = next;
        }
    }
    auto_corr /= var;
    return auto_corr;
}

// compute sample standard deviation of a float array (with iterator)
// sqrt ( ( \sum x^2 - (\sum x)^2 / N ) / ( N - 1 ) )
static float get_sample_sd(const float *v, int n, const int *imap) {
    float mean = get_sample_mean(v, n, imap);
    float s2 = 0.0;
    int j = 0;
    for (int i = 0; i < n; i++) {
        float tmp = imap ? v[imap[i]] : v[i];
        if (!isnan(tmp)) {
            s2 += sq(tmp - mean);
            j++;
        }
    }
    s2 /= (double)(j - 1);
    return (float)sqrt(s2);
}

// compute standard error of mean a float array (with iterator)
// sqrt ( ( \sum x^2 - (\sum x)^2 / N ) / ( N - 1 ) / N )
static float get_se_mean(const float *v, int n, const int *imap) {
    int j = 0;
    for (int i = 0; i < n; i++) {
        float tmp = imap ? v[imap[i]] : v[i];
        if (!isnan(tmp)) j++;
    }
    if (j <= 1) return NAN;

    return get_sample_sd(v, n, imap) / sqrtf(j);
}

// compute (adjusted) LRR autocorrelation for a float array with iterator
static float get_lrr_auto_corr(const float *lrr, int n, const int *imap) {
    float value;
    double mean = 0.0;
    int i = 0, j = 0;
    for (value = next_not_missing(lrr, imap, n, &i); i < n; i++, value = next_not_missing(lrr, imap, n, &i)) {
        mean += (double)value;
        j++;
    }
    if (j <= 1) return NAN;
    mean /= (double)j;

    double var = 0.0;
    i = 0;
    for (value = next_not_missing(lrr, imap, n, &i); i < n; i++, value = next_not_missing(lrr, imap, n, &i)) {
        var += sq((double)value - mean);
    }

    double auto_corr = 0.0;
    i = 0;
    double prev = (double)next_not_missing(lrr, imap, n, &i) - mean, next;
    for (i++, value = next_not_missing(lrr, imap, n, &i); i < n; i++, value = next_not_missing(lrr, imap, n, &i)) {
        next = (double)value - mean;
        auto_corr += prev * next;
        prev = next;
    }
    auto_corr /= var;
    return auto_corr;
}

// compute the n50 of the call when split at the heterozygous sites
static int get_n50_hets(const int *v, const float *baf, int n, int beg_pos, int end_pos) {
    int *w = (int *)malloc((n + 1) * sizeof(int));
    int i, j, sum, prev = beg_pos;
    for (i = 0, j = 0; i < n; i++) {
        if (!isnan(baf[i])) {
            w[j++] = v[i] - prev;
            prev = v[i];
        }
    }
    w[j++] = end_pos - prev;

    ks_introsort_int((size_t)j, w);

    for (i = 0, sum = 0; sum << 1 < end_pos - beg_pos && i < j; i++) sum += w[i];
    int n50 = w[i - 1];
    free(w);
    return n50;
}

/*********************************
 * SAMPLE METHODS                *
 *********************************/

static void mocha_print_ucsc(FILE *restrict stream, const mocha_t *mocha, int n, const bcf_hdr_t *hdr) {
    if (stream == NULL) return;
    const char *name[4];
    name[MOCHA_UNDET] = "mCA_undetermined";
    name[MOCHA_LOSS] = "mCA_loss";
    name[MOCHA_GAIN] = "mCA_gain";
    name[MOCHA_CNLOH] = "mCA_neutral";
    const char *desc[4];
    desc[MOCHA_UNDET] = "Undetermined";
    desc[MOCHA_LOSS] = "Losses";
    desc[MOCHA_GAIN] = "Gains";
    desc[MOCHA_CNLOH] = "CN-LOHs";
    kstring_t tmp = {0, 0, NULL};
    for (int i = 0; i < 4; i++) {
        fprintf(stream, "track name=%s description=\"%s\" visibility=4 priority=1 itemRgb=\"On\"\n", name[i], desc[i]);
        for (int j = 0; j < n; j++) {
            uint8_t red = 127, green = 127, blue = 127;
            switch (i) {
            case MOCHA_LOSS:
                red -= (uint8_t)(127.0f * sqf(mocha[j].cf));
                green -= (uint8_t)(127.0f * sqf(mocha[j].cf));
                blue += (uint8_t)(128.0f * sqf(mocha[j].cf));
                break;
            case MOCHA_GAIN:
                red += (uint8_t)(128.0f * sqf(mocha[j].cf));
                green -= (uint8_t)(127.0f * sqf(mocha[j].cf));
                blue -= (uint8_t)(127.0f * sqf(mocha[j].cf));
                break;
            case MOCHA_CNLOH:
                red += (uint8_t)(128.0f * mocha[j].cf);
                green += (uint8_t)(38.0f * mocha[j].cf);
                blue -= (uint8_t)(127.0f * mocha[j].cf);
                break;
            }
            if (i == mocha[j].type) {
                const char *sample_name = bcf_hdr_int2id(hdr, BCF_DT_SAMPLE, mocha[j].sample_idx);
                tmp.l = 0;
                kputs(sample_name, &tmp);
                for (char *p = tmp.s; (p = strchr(p, ' ')); ++p) *p = '_';
                const char *seq_name = bcf_hdr_id2name(hdr, mocha[j].rid);
                if (strncmp(seq_name, "chr", 3) == 0)
                    fprintf(stream, "%s\t%d\t%d\t%s\t%d\t.\t%d\t%d\t%d,%d,%d\n", seq_name, mocha[j].beg_pos,
                            mocha[j].end_pos, tmp.s, isnan(mocha[j].cf) ? 0 : (int)(1e3 * mocha[j].cf),
                            mocha[j].beg_pos, mocha[j].end_pos, red, green, blue);
                else
                    fprintf(stream, "chr%s\t%d\t%d\t%s\t%d\t.\t%d\t%d\t%d,%d,%d\n", seq_name, mocha[j].beg_pos,
                            mocha[j].end_pos, tmp.s, isnan(mocha[j].cf) ? 0 : (int)(1e3 * mocha[j].cf),
                            mocha[j].beg_pos, mocha[j].end_pos, red, green, blue);
            }
        }
    }
    free(tmp.s);
    if (stream != stdout && stream != stderr) fclose(stream);
}

static void mocha_print_calls(FILE *restrict stream, const mocha_t *mocha, int n, const bcf_hdr_t *hdr, int flags,
                              char *genome, float lrr_hap2dip) {
    if (stream == NULL) return;
    char gender[4];
    gender[GENDER_UNKNOWN] = 'U';
    gender[GENDER_MALE] = 'M';
    gender[GENDER_FEMALE] = 'F';
    gender[GENDER_KLINEFELTER] = 'K';
    const char *type[6];
    type[MOCHA_UNDET] = "Undetermined";
    type[MOCHA_LOSS] = "Loss";
    type[MOCHA_GAIN] = "Gain";
    type[MOCHA_CNLOH] = "CN-LOH";
    type[MOCHA_CNP_LOSS] = "CNP_Loss";
    type[MOCHA_CNP_GAIN] = "CNP_Gain";
    char arm_type[4];
    arm_type[MOCHA_NOT] = 'N';
    arm_type[MOCHA_CEN] = 'C';
    arm_type[MOCHA_ARM] = 'Y';
    arm_type[MOCHA_TEL] = 'T';
    fputs("sample_id", stream);
    fputs("\tcomputed_gender", stream);
    fputs("\tchrom", stream);
    fprintf(stream, "\tbeg_%s", genome);
    fprintf(stream, "\tend_%s", genome);
    fputs("\tlength", stream);
    fputs("\tp_arm", stream);
    fputs("\tq_arm", stream);
    fputs("\tn_sites", stream);
    fputs("\tn_hets", stream);
    fputs("\tn50_hets", stream);
    fputs("\tbdev", stream);
    fputs("\tbdev_se", stream);
    fputs("\trel_cov", stream);
    fputs("\trel_cov_se", stream);
    fputs("\tlod_lrr_baf", stream);
    fputs("\tlod_baf_phase", stream);
    fputs("\tn_flips", stream);
    fputs("\tbaf_conc", stream);
    fputs("\tlod_baf_conc", stream);
    fputs("\ttype", stream);
    fputs("\tcf", stream);
    fputc('\n', stream);
    for (int i = 0; i < n; i++) {
        fputs(bcf_hdr_int2id(hdr, BCF_DT_SAMPLE, mocha->sample_idx), stream);
        fprintf(stream, "\t%c", gender[mocha->computed_gender]);
        fprintf(stream, "\t%s", bcf_hdr_id2name(hdr, mocha->rid));
        fprintf(stream, "\t%d", mocha->beg_pos);
        fprintf(stream, "\t%d", mocha->end_pos);
        fprintf(stream, "\t%d", mocha->length);
        fprintf(stream, "\t%c", arm_type[mocha->p_arm]);
        fprintf(stream, "\t%c", arm_type[mocha->q_arm]);
        fprintf(stream, "\t%d", mocha->n_sites);
        fprintf(stream, "\t%d", mocha->n_hets);
        fprintf(stream, "\t%d", mocha->n50_hets);
        fprintf(stream, "\t%.4f", mocha->bdev);
        fprintf(stream, "\t%.4f", mocha->bdev_se);
        fprintf(stream, "\t%.4f", 2.0f * expf(mocha->ldev / lrr_hap2dip * (float)M_LN2));
        fprintf(stream, "\t%.4f",
                2.0f * expf(mocha->ldev / lrr_hap2dip * (float)M_LN2) * mocha->ldev_se / lrr_hap2dip * (float)M_LN2);
        fprintf(stream, "\t%.2f", mocha->lod_lrr_baf);
        fprintf(stream, "\t%.2f", mocha->lod_baf_phase);
        fprintf(stream, "\t%d", mocha->n_flips);
        fprintf(stream, "\t%.4f", mocha->baf_conc);
        fprintf(stream, "\t%.2f", mocha->lod_baf_conc);
        fprintf(stream, "\t%s", type[mocha->type]);
        fprintf(stream, "\t%.4f", mocha->cf);
        fputc('\n', stream);
        mocha++;
    }
    if (stream != stdout && stream != stderr) fclose(stream);
}

// this function returns two values (a, b) such that:
// (i) a <= b (ii) pos[a-1] < beg (iii) beg <= pos[a] < end (iv) beg <= pos[b] < end (v)
// pos[b+1] >= end
static int get_cnp_edges(const int *pos, int n, int beg, int end, int *a, int *b) {
    if (pos[0] >= end || pos[n - 1] < beg) return -1;

    int i = 0, j = n - 1, k;
    while (j - i > 1) {
        k = (i + j) / 2;
        if (pos[k] < beg)
            i = k;
        else
            j = k;
    }
    if (pos[j] >= end) return -1;

    *a = j;

    i = j;
    j = n - 1;
    while (j - i > 1) {
        k = (i + j) / 2;
        if (pos[k] < end)
            i = k;
        else
            j = k;
    }

    *b = i;
    return 0;
}

// classify mosaic chromosomal alteration type based on LRR and BAF
// LDEV = -log2( 1 - 2 x BDEV ) * LRR-hap2dip for gains
// LDEV = -log2( 1 + 2 x BDEV ) * LRR-hap2dip for losses
static int8_t mocha_type(float ldev, float ldev_se, float bdev, float bdev_se, int n_hets, float lrr_hap2dip,
                         int8_t p_arm, int8_t q_arm, int is_short_arm) {
    // a LOD score can be computed from a chi-squared statistic by dividing by 2ln(10) ~ 4.6
    float z2_cnloh = sqf(ldev / ldev_se);

    // equivalent of a 2 LOD score bonus for ending in one but not two telomeres (unless it is a short arm autosome)
    if (is_short_arm) {
        if (q_arm != MOCHA_TEL)
            z2_cnloh += 4.0f * M_LN10; // CN-LOH less likely
        else if (p_arm != MOCHA_CEN)
            z2_cnloh -= 4.0f * M_LN10; // CN-LOH more likely
    } else {
        if ((p_arm == MOCHA_TEL) == (q_arm == MOCHA_TEL))
            z2_cnloh += 4.0f * M_LN10; // CN-LOH less likely
        else
            z2_cnloh -= 4.0f * M_LN10; // CN-LOH more likely
    }

    // if one model has 4 LOD scores point more than the other model, select the better model
    if (ldev > 0) {
        if (n_hets < 5 || isnan(bdev)) {
            if (!isnan(bdev_se) && z2_cnloh < 8.0 * M_LN10)
                return MOCHA_CNLOH;
            else
                return MOCHA_GAIN;
        } else {
            float expected_ldev =
                -logf(1.0f - 2.0f * bdev > 2.0f / 3.0f ? 1.0f - 2.0f * bdev : 2.0f / 3.0f) * M_LOG2E * lrr_hap2dip;
            float z2_gain = sqf((ldev - expected_ldev) / ldev_se);
            if (z2_cnloh > z2_gain + 8.0f * M_LN10) return MOCHA_GAIN;
            if (z2_gain > z2_cnloh + 8.0f * M_LN10) return MOCHA_CNLOH;
        }
    } else {
        if (n_hets < 5 || isnan(bdev)) {
            if (!isnan(bdev_se) && z2_cnloh < 8.0 * M_LN10)
                return MOCHA_CNLOH;
            else
                return MOCHA_LOSS;
        } else {
            float expected_ldev = -logf(1.0f + 2.0f * bdev) * M_LOG2E * lrr_hap2dip;
            float z2_loss = sqf((ldev - expected_ldev) / ldev_se);
            if (z2_cnloh > z2_loss + 8.0f * M_LN10) return MOCHA_LOSS;
            if (z2_loss > z2_cnloh + 8.0f * M_LN10) return MOCHA_CNLOH;
        }
    }
    return MOCHA_UNDET;
}

// best estimate for cell fraction using the following formula (if BDEV is available):
// BDEV = | 1 / 2 - 1 / CNF |
// CNF = 2 / ( 1 + 2 x BDEV ) for losses
// CNF = 2 / ( 1 - 2 x BDEV ) for gains
// CNF = 2 x 2^( LDEV / LRR-hap2dip )
// LDEV = - LRR-hap2dip / ln(2) * ln( 1 + 2 x BDEV ) for losses
// LDEV = - LRR-hap2dip / ln(2) * ln( 1 - 2 x BDEV ) for gains
static float mocha_cell_fraction(float ldev, float ldev_se, float bdev, int n_hets, int8_t type, float lrr_hap2dip) {
    if (n_hets < 5 || isnan(bdev)) {
        switch (type) {
        case MOCHA_LOSS:
            return ldev < -lrr_hap2dip ? 1.0f : -2.0f * (expf(ldev / lrr_hap2dip * (float)M_LN2) - 1.0f);
        case MOCHA_GAIN:
            return ldev * (float)M_LN2 > lrr_hap2dip * logf(1.5f)
                       ? 1.0f
                       : 2.0f * (expf(ldev / lrr_hap2dip * (float)M_LN2) - 1.0f);
        default:
            return NAN;
        }
    } else {
        switch (type) {
        case MOCHA_LOSS:
            return 4.0f * bdev / (1.0f + 2.0f * bdev);
        case MOCHA_GAIN:
            return bdev > 1.0f / 6.0f ? 1.0f : 4.0f * bdev / (1.0f - 2.0f * bdev);
        case MOCHA_CNLOH:
            return 2.0f * bdev;
        case MOCHA_UNDET:
            return bdev < 0.05f ? 4.0f * bdev : NAN; // here it assumes it is either a loss or a gain
        default:
            return NAN;
        }
    }
}

static void get_mocha_stats(const int *pos, const float *lrr, const float *baf, const int8_t *gt_phase, int n, int a,
                            int b, int cen_beg, int cen_end, int length, float baf_conc, int is_short_arm,
                            mocha_t *mocha) {
    mocha->n_sites = b + 1 - a;

    if (a == 0) {
        mocha->p_arm = is_short_arm ? MOCHA_CEN : MOCHA_TEL;
        if (pos[a] < cen_beg)
            mocha->beg_pos = 0;
        else
            mocha->beg_pos = cen_end;
    } else {
        mocha->beg_pos = pos[a];
        if (mocha->beg_pos < cen_beg)
            mocha->p_arm = MOCHA_ARM;
        else if (mocha->beg_pos <= cen_end)
            mocha->p_arm = MOCHA_CEN;
        else
            mocha->p_arm = MOCHA_NOT;
    }

    if (b == n - 1) {
        mocha->q_arm = MOCHA_TEL;
        if (pos[b] >= cen_end)
            mocha->end_pos = length;
        else
            mocha->end_pos = cen_beg;
    } else {
        mocha->end_pos = pos[b];
        if (mocha->end_pos > cen_end)
            mocha->q_arm = MOCHA_ARM;
        else if (mocha->end_pos >= cen_beg)
            mocha->q_arm = MOCHA_CEN;
        else
            mocha->q_arm = MOCHA_NOT;
    }

    mocha->length = mocha->end_pos - mocha->beg_pos;

    mocha->ldev_se = get_se_mean(lrr + a, b + 1 - a, NULL);
    mocha->n_hets = 0;
    for (int i = a; i <= b; i++)
        if (!isnan(baf[i])) mocha->n_hets++;
    int conc, disc;
    get_baf_conc(baf + a, gt_phase + a, b + 1 - a, &conc, &disc);
    mocha->baf_conc = conc + disc > 0 ? (float)conc / (float)(conc + disc) : NAN;
    mocha->lod_baf_conc = ((mocha->baf_conc > 0 ? (float)conc * logf(mocha->baf_conc / baf_conc) : 0)
                           + (mocha->baf_conc < 1 ? (float)disc * logf((1 - mocha->baf_conc) / (1 - baf_conc)) : 0))
                          * (float)M_LOG10E;
    mocha->n50_hets = get_n50_hets(pos + a, baf + a, b + 1 - a, mocha->beg_pos, mocha->end_pos);
    mocha->n_flips = -1;
    mocha->bdev = NAN;
    mocha->bdev_se = NAN;
    mocha->lod_baf_phase = NAN;
}

// return segments called by the HMM or a suggestion of what state should be added to the HMM
static float get_path_segs(const int8_t *path, const float *hs_arr, int n, int hmm_model, int middle, int **beg,
                           int *m_beg, int **end, int *m_end, int *nseg) {
    int a = 0, b = 0;
    *nseg = 0;
    for (b = 0; b < n; b++) {
        // check whether it is the end of a segment
        if (b != n - 1) {
            int x = abs(path[b]);
            int y = abs(path[b + 1]);
            if (x == y) continue;

            // swap the two elements
            if (x > y) {
                x = abs(path[b + 1]);
                y = abs(path[b]);
            }

            if (y - x == 1) {
                if (hmm_model == LRR_BAF && x > 0 && x != middle)
                    return (hs_arr[x - 1] + hs_arr[y - 1]) * 0.5f;
                else if (hmm_model == BAF_PHASE)
                    return x == 0 ? hs_arr[0] * 0.5f : (hs_arr[x - 1] + hs_arr[y - 1]) * 0.5f;
            }
        }

        if (path[b]) {
            (*nseg)++;
            hts_expand(int, *nseg, *m_beg, *beg);
            (*beg)[(*nseg) - 1] = a;
            hts_expand(int, *nseg, *m_end, *end);
            (*end)[(*nseg) - 1] = b;
        }
        a = b + 1;
    }
    return 0;
}

// process one contig for one sample
static void sample_run(sample_t *self, mocha_table_t *mocha_table, const model_t *model) {
    // do nothing if chromosome Y or MT are being tested
    if (model->rid == model->genome_rules->y_rid || model->rid == model->genome_rules->mt_rid) {
        memset(self->data_arr[LDEV], 0, self->n * sizeof(int16_t));
        memset(self->data_arr[BDEV], 0, self->n * sizeof(int16_t));
        memset(self->phase_arr, 0, self->n * sizeof(int8_t));
        return;
    }

    mocha_t mocha;
    mocha.sample_idx = self->idx;
    mocha.computed_gender = self->computed_gender;
    mocha.rid = model->rid;

    int cen_beg = model->genome_rules->cen_beg[model->rid];
    int cen_end = model->genome_rules->cen_end[model->rid];
    int length = model->genome_rules->length[model->rid];
    if (length == 0) length = model->locus_arr[model->n - 1].pos;
    // incentive to extend to the telomere
    float tel_log_prb = model->rid == model->genome_rules->x_rid
                            ? (self->computed_gender == GENDER_MALE ? model->chrY_tel_log_prb : model->chrX_tel_log_prb)
                            : model->auto_tel_log_prb;

    // declutter code by copying these values onto the stack
    int n = self->n;
    int8_t *gt_phase = self->phase_arr;
    int16_t *ad0 = self->data_arr[AD0];
    int16_t *ad1 = self->data_arr[AD1];
    float *lrr = (float *)malloc(n * sizeof(float));
    float *baf = (float *)malloc(n * sizeof(float));
    if (model->flags & WGS_DATA) {
        ad_to_lrr_baf(ad0, ad1, lrr, baf, n);
    } else {
        for (int i = 0; i < n; i++) {
            lrr[i] = int16_to_float(self->data_arr[LRR][i]);
            baf[i] = int16_to_float(self->data_arr[BAF][i]);
        }
    }

    if (model->lrr_gc_order > 0 && n > model->lrr_gc_order)
        adjust_lrr(lrr, model->gc_arr, n, self->vcf_imap_arr, self->stats.coeffs, model->lrr_gc_order);
    else if (model->lrr_gc_order != -1)
        adjust_lrr(lrr, model->gc_arr, n, self->vcf_imap_arr, self->stats.coeffs, 0);

    int8_t *as = (int8_t *)calloc(n, sizeof(int8_t));
    int *pos = (int *)malloc(n * sizeof(int));
    for (int i = 0; i < n; i++) pos[i] = model->locus_arr[self->vcf_imap_arr[i]].pos;
    int *imap_arr = (int *)malloc(n * sizeof(int));
    int *hets_imap_arr = (int *)malloc(n * sizeof(int));
    float *pbaf_arr = (float *)malloc(n * sizeof(float));

    int16_t *ldev = (int16_t *)calloc(n, sizeof(int16_t));
    int16_t *bdev = (int16_t *)calloc(n, sizeof(int16_t));

    if (model->rid == model->genome_rules->x_rid) {
        if (self->computed_gender == GENDER_MALE) {
            for (int i = 0; i < n; i++) {
                if (pos[i] > model->genome_rules->x_nonpar_beg && pos[i] < model->genome_rules->x_xtr_beg)
                    lrr[i] = baf[i] = NAN;
                else if (pos[i] > model->genome_rules->x_xtr_end && pos[i] < model->genome_rules->x_nonpar_end)
                    lrr[i] = baf[i] = NAN;
                else if (!(model->flags & USE_MALES_XTR) && pos[i] >= model->genome_rules->x_xtr_beg
                         && pos[i] <= model->genome_rules->x_xtr_end)
                    lrr[i] = baf[i] = NAN;
                else if (!(model->flags & USE_MALES_PAR2) && pos[i] >= model->genome_rules->x_nonpar_end)
                    lrr[i] = baf[i] = NAN;
            }
        } else if (self->computed_gender == GENDER_KLINEFELTER) { // only analyze diploid region for XXY individuals
            for (int i = 0; i < n; i++) {
                if (pos[i] <= model->genome_rules->x_nonpar_beg)
                    lrr[i] = baf[i] = NAN;
                else if (pos[i] >= model->genome_rules->x_xtr_beg && pos[i] <= model->genome_rules->x_xtr_end)
                    lrr[i] = baf[i] = NAN;
                else if (pos[i] >= model->genome_rules->x_nonpar_end)
                    lrr[i] = baf[i] = NAN;
            }
        }
    }

    if (model->cnp_itr) {
        while (regitr_overlap(model->cnp_itr)) {
            int a, b;
            if (get_cnp_edges(pos, n, model->cnp_itr->beg, model->cnp_itr->end, &a, &b) == 0) {
                int cnp_type = regitr_payload(model->cnp_itr, int);
                float exp_ldev = NAN;
                float exp_bdev = NAN;
                mocha.type = MOCHA_UNDET;
                mocha.ldev = get_median(lrr + a, b + 1 - a, NULL);
                if (mocha.ldev > 0 && (cnp_type == MOCHA_CNP_GAIN || cnp_type == MOCHA_CNP_CNV)) {
                    if (model->flags & WGS_DATA) {
                        minus_lrr_ad_lod_t data = {lrr + a,
                                                   ad0 + a,
                                                   ad1 + a,
                                                   b + 1 - a,
                                                   NULL,
                                                   model->err_log_prb,
                                                   model->lrr_bias,
                                                   model->lrr_hap2dip,
                                                   self->adjlrr_sd,
                                                   self->stats.dispersion};
                        mocha.lod_lrr_baf = -minus_lrr_ad_lod(1.0f / 6.0f, (void *)&data);
                    } else {
                        minus_lrr_baf_lod_t data = {lrr + a,
                                                    baf + a,
                                                    b + 1 - a,
                                                    NULL,
                                                    model->err_log_prb,
                                                    model->lrr_bias,
                                                    model->lrr_hap2dip,
                                                    self->adjlrr_sd,
                                                    self->stats.dispersion};
                        mocha.lod_lrr_baf = -minus_lrr_baf_lod(1.0f / 6.0f, (void *)&data);
                    }
                    if (mocha.lod_lrr_baf
                        > -(model->xy_major_log_prb + model->xy_minor_log_prb) / 2.0f * (float)M_LOG10E) {
                        mocha.type = MOCHA_CNP_GAIN;
                        mocha.cf = NAN;
                        exp_ldev = log2f(1.5f) * model->lrr_hap2dip;
                        exp_bdev = 1.0f / 6.0f;
                    }
                } else if (mocha.ldev <= 0 && (cnp_type == MOCHA_CNP_LOSS || cnp_type == MOCHA_CNP_CNV)) {
                    if (model->flags & WGS_DATA) {
                        minus_lrr_ad_lod_t data = {lrr + a,
                                                   ad0 + a,
                                                   ad1 + a,
                                                   b + 1 - a,
                                                   NULL,
                                                   model->err_log_prb,
                                                   model->lrr_bias,
                                                   model->lrr_hap2dip,
                                                   self->adjlrr_sd,
                                                   self->stats.dispersion};
                        mocha.lod_lrr_baf = -minus_lrr_ad_lod(-0.5f, (void *)&data);
                    } else {
                        minus_lrr_baf_lod_t data = {lrr + a,
                                                    baf + a,
                                                    b + 1 - a,
                                                    NULL,
                                                    model->err_log_prb,
                                                    model->lrr_bias,
                                                    model->lrr_hap2dip,
                                                    self->adjlrr_sd,
                                                    self->stats.dispersion};
                        mocha.lod_lrr_baf = -minus_lrr_baf_lod(-0.5f, (void *)&data);
                    }
                    if (mocha.lod_lrr_baf
                        > -(model->xy_major_log_prb + model->xy_minor_log_prb) / 2.0f * (float)M_LOG10E) {
                        mocha.type = MOCHA_CNP_LOSS;
                        mocha.cf = NAN;
                        exp_ldev = -model->lrr_hap2dip;
                        exp_bdev = 0.5f;
                    }
                }
                if (mocha.type == MOCHA_CNP_GAIN || mocha.type == MOCHA_CNP_LOSS) {
                    if (model->flags & WGS_DATA) {
                        if (cnp_edge_is_not_cn2_lrr_ad(lrr, ad0, ad1, n, a, b,
                                                       (model->xy_major_log_prb + model->xy_minor_log_prb) / 2.0f,
                                                       model->err_log_prb, model->lrr_bias, model->lrr_hap2dip,
                                                       self->adjlrr_sd, self->stats.dispersion, exp_ldev, exp_bdev))
                            continue;
                    } else {
                        if (cnp_edge_is_not_cn2_lrr_baf(lrr, baf, n, a, b,
                                                        (model->xy_major_log_prb + model->xy_minor_log_prb) / 2.0f,
                                                        model->err_log_prb, model->lrr_bias, model->lrr_hap2dip,
                                                        self->adjlrr_sd, self->stats.dispersion, exp_ldev, exp_bdev))
                            continue;
                    }
                    get_mocha_stats(pos, lrr, baf, gt_phase, n, a, b, cen_beg, cen_end, length, self->stats.baf_conc,
                                    model->genome_rules->is_short_arm[model->rid], &mocha);
                    // compute bdev, if possible
                    if (mocha.n_hets > 0) {
                        mocha.bdev = model->flags & WGS_DATA
                                         ? get_ad_bdev(ad0 + a, ad1 + a, b + 1 - a, NULL, self->stats.dispersion)
                                         : get_baf_bdev(baf + a, b + 1 - a, NULL, self->stats.dispersion);
                    } else {
                        mocha.bdev = NAN;
                    }
                    mocha_table->n++;
                    hts_expand(mocha_t, mocha_table->n, mocha_table->m, mocha_table->a);
                    mocha_table->a[mocha_table->n - 1] = mocha;
                    for (int j = a; j <= b; j++) {
                        // TODO add other stuff here, like setting ldev
                        // and bdev
                        lrr[j] = NAN; // do not use the data again
                        baf[j] = NAN; // do not use the data again
                    }
                }
            }
        }
    }

    float *hs_arr = NULL;
    int n_hs = 0, m_hs = 0;
    for (int hmm_model = 0; hmm_model < 2; hmm_model++) {
        // select data to use from the contig, depending on which HMM model is being
        // used
        int last_p = 0, first_q = 0;
        int n_imap = 0;
        for (int i = 0; i < n; i++)
            if ((hmm_model == LRR_BAF && !isnan(lrr[i])) || (hmm_model == BAF_PHASE && !isnan(baf[i]))) {
                if (pos[i] < cen_beg) last_p++;
                if (pos[i] < cen_end) first_q++;
                imap_arr[n_imap] = i;
                n_imap++;
            }
        if (n_imap == 0) continue;

        // compute emission probabilities and Viterbi path according to HMM model
        int middle = 0;
        // TODO eliminate this redundancy
        if (hmm_model == LRR_BAF) {
            n_hs = model->bdev_lrr_baf_n;
            hts_expand(float, n_hs, m_hs, hs_arr);
            for (int i = 0; i < model->bdev_lrr_baf_n; i++) {
                hs_arr[i] = model->bdev_lrr_baf[i];
                if (model->bdev_lrr_baf[i] < 0.0f) middle++;
            }
        } else if (hmm_model == BAF_PHASE) {
            n_hs = model->bdev_baf_phase_n;
            hts_expand(float, n_hs, m_hs, hs_arr);
            for (int i = 0; i < model->bdev_baf_phase_n; i++) hs_arr[i] = model->bdev_baf_phase[i];
        }
        int8_t *path;
        float ret;
        int *beg = NULL, m_beg = 0, *end = NULL, m_end = 0, nseg;
        do {
            if (n_hs + (hmm_model == LRR_BAF ? n_hs : 0) > 50) error("Too many states being tested for the HMM\n");

            float *emis_log_lkl;
            if (model->flags & WGS_DATA) {
                emis_log_lkl =
                    hmm_model == LRR_BAF
                        ? lrr_ad_emis_log_lkl(lrr, ad0, ad1, n_imap, imap_arr, model->err_log_prb, model->lrr_bias,
                                              model->lrr_hap2dip, self->adjlrr_sd, self->stats.dispersion, hs_arr, n_hs)
                        : ad_phase_emis_log_lkl(ad0, ad1, gt_phase, n_imap, imap_arr, model->err_log_prb,
                                                self->stats.dispersion, hs_arr, n_hs);
            } else {
                emis_log_lkl = hmm_model == LRR_BAF
                                   ? lrr_baf_emis_log_lkl(lrr, baf, n_imap, imap_arr, model->err_log_prb,
                                                          model->lrr_bias, model->lrr_hap2dip, self->adjlrr_sd,
                                                          self->stats.dispersion, hs_arr, n_hs)
                                   : baf_phase_emis_log_lkl(baf, gt_phase, n_imap, imap_arr, model->err_log_prb,
                                                            self->stats.dispersion, hs_arr, n_hs);
            }
            path = log_viterbi_run(emis_log_lkl, n_imap, n_hs + (hmm_model == LRR_BAF ? n_hs : 0),
                                   model->xy_major_log_prb, model->xy_minor_log_prb, tel_log_prb,
                                   hmm_model == LRR_BAF ? NAN : model->flip_log_prb, last_p, first_q);
            free(emis_log_lkl);

            if (hmm_model == LRR_BAF)
                for (int i = 0; i < n_imap; i++)
                    if (path[i] > n_hs) path[i] = 0;
            ret = get_path_segs(path, hs_arr, n_imap, hmm_model, middle, &beg, &m_beg, &end, &m_end, &nseg);

            if (ret) // two consecutive hidden states were used, hinting that
                     // testing of a middle state might be necessary
            {
                free(path);
                n_hs++;
                if (middle && ret < 0.0f) middle++;
                hts_expand(float, n_hs, m_hs, hs_arr);
                hs_arr[n_hs - 1] = ret;
                ks_introsort_float(n_hs, hs_arr);
            }
        } while (ret);

        // loop through all the segments called by the Viterbi algorithm
        for (int i = 0; i < nseg; i++) {
            // compute edges of the call
            int a = imap_arr[beg[i]];
            if (beg[i] == 0)
                while (a > 0 && ldev[a - 1] == 0 && bdev[a - 1] == 0) a--; // extend call towards p telomere
            int b = imap_arr[end[i]];
            if (end[i] == n_imap - 1)
                while (b < n - 1 && ldev[b + 1] == 0 && bdev[b + 1] == 0) b++; // extend call towards q telomere

            // if an autosomal call spans the whole chromosome and seems an isochromosome event, split it in two
            if (model->rid != model->genome_rules->x_rid && a == 0 && b == n - 1 && last_p > 0 && first_q < n_imap) {
                float p_arm_ldev = get_median(lrr + 0, imap_arr[last_p - 1] + 1, NULL);
                float q_arm_ldev = get_median(lrr + imap_arr[first_q], n - imap_arr[first_q], NULL);
                float left_arm_ldev_se = get_se_mean(lrr + 0, imap_arr[last_p - 1] + 1, NULL);
                float right_arm_ldev_se = get_se_mean(lrr + imap_arr[first_q], n - imap_arr[first_q], NULL);
                float z2 = sqf(p_arm_ldev - q_arm_ldev) / (sqf(left_arm_ldev_se) + sqf(right_arm_ldev_se));
                if (p_arm_ldev * q_arm_ldev < 0.0f && z2 > 8.0f * M_LN10) {
                    end[i] = last_p - 1;
                    b = imap_arr[end[i]];
                    nseg++;
                    hts_expand(int, nseg, m_beg, beg);
                    hts_expand(int, nseg, m_end, end);
                    beg[nseg - 1] = first_q;
                    end[nseg - 1] = n_imap - 1;
                }
            }

            mocha.ldev = get_median(lrr + a, b + 1 - a, NULL);
            get_mocha_stats(pos, lrr, baf, gt_phase, n, a, b, cen_beg, cen_end, length, self->stats.baf_conc,
                            model->genome_rules->is_short_arm[model->rid], &mocha);

            double bdev_lrr_baf;
            if (model->flags & WGS_DATA) {
                minus_lrr_ad_lod_t data = {lrr + a,
                                           ad0 + a,
                                           ad1 + a,
                                           mocha.n_sites,
                                           NULL,
                                           NAN,
                                           model->lrr_bias,
                                           model->lrr_hap2dip,
                                           self->adjlrr_sd,
                                           self->stats.dispersion};
                kmin_brent(minus_lrr_ad_lod, -0.15, 0.15, (void *)&data, KMIN_EPS, &bdev_lrr_baf);
                data.err_log_prb = model->err_log_prb;
                mocha.lod_lrr_baf = -minus_lrr_ad_lod(bdev_lrr_baf, (void *)&data);
            } else {
                minus_lrr_baf_lod_t data = {lrr + a,
                                            baf + a,
                                            mocha.n_sites,
                                            NULL,
                                            NAN,
                                            model->lrr_bias,
                                            model->lrr_hap2dip,
                                            self->adjlrr_sd,
                                            self->stats.dispersion};
                kmin_brent(minus_lrr_baf_lod, -0.15, 0.15, (void *)&data, KMIN_EPS, &bdev_lrr_baf);
                data.err_log_prb = model->err_log_prb;
                mocha.lod_lrr_baf = -minus_lrr_baf_lod(bdev_lrr_baf, (void *)&data);
            }

            if (hmm_model == LRR_BAF) {
                // here you need to check whether the call would have been
                // better with the phased HMM model
                int n_hets_imap = 0;
                for (int j = beg[i]; j <= end[i]; j++)
                    if (!isnan(baf[imap_arr[j]])) {
                        n_hets_imap++;
                        hets_imap_arr[n_hets_imap - 1] = imap_arr[j];
                    }
                // TODO here it needs to pass information about the centromeres
                if (model->flags & WGS_DATA) {
                    mocha.lod_baf_phase =
                        compare_wgs_models(ad0, ad1, gt_phase, n_hets_imap, hets_imap_arr,
                                           (model->xy_major_log_prb + model->xy_minor_log_prb) / 2.0f,
                                           model->err_log_prb, model->flip_log_prb, tel_log_prb, self->stats.dispersion,
                                           model->bdev_baf_phase, model->bdev_baf_phase_n);
                } else {
                    mocha.lod_baf_phase =
                        compare_models(baf, gt_phase, n_hets_imap, hets_imap_arr,
                                       (model->xy_major_log_prb + model->xy_minor_log_prb) / 2.0f, model->err_log_prb,
                                       model->flip_log_prb, tel_log_prb, self->stats.dispersion, model->bdev_baf_phase,
                                       model->bdev_baf_phase_n);
                }
                if (mocha.lod_baf_phase > 0.8 * mocha.lod_lrr_baf) continue;

                // compute bdev, if possible
                if (n_hets_imap > 0) {
                    mocha.bdev = model->flags & WGS_DATA
                                     ? get_ad_bdev(ad0, ad1, n_hets_imap, hets_imap_arr, self->stats.dispersion)
                                     : get_baf_bdev(baf, n_hets_imap, hets_imap_arr, self->stats.dispersion);
                } else {
                    mocha.bdev = NAN;
                }
                mocha.bdev_se = NAN;
                for (int j = 0; j < n_hets_imap; j++) as[hets_imap_arr[j]] = (int8_t)SIGN(baf[hets_imap_arr[j]] - 0.5f);
            } else {
                // penalizes the LOD by the number of phase flips
                mocha.n_flips = 0;
                for (int j = beg[i]; j < end[i]; j++)
                    if (path[j] != path[j + 1]) mocha.n_flips++;

                if (model->flags & WGS_DATA) {
                    minus_ad_phase_lod_t data = {ad0,
                                                 ad1,
                                                 gt_phase,
                                                 mocha.n_hets,
                                                 imap_arr + beg[i],
                                                 path + beg[i],
                                                 NAN,
                                                 self->stats.dispersion};
                    double bdev;
                    kmin_brent(minus_ad_phase_lod, 0.1, 0.2, (void *)&data, KMIN_EPS, &bdev);
                    mocha.bdev = fabsf((float)bdev);
                    data.err_log_prb = model->err_log_prb;
                    mocha.lod_baf_phase = -minus_ad_phase_lod(mocha.bdev, (void *)&data);
                } else {
                    double bdev;
                    minus_baf_phase_lod_t data = {baf,           gt_phase, mocha.n_hets,          imap_arr + beg[i],
                                                  path + beg[i], NAN,      self->stats.dispersion};
                    kmin_brent(minus_baf_phase_lod, 0.1, 0.2, (void *)&data, KMIN_EPS, &bdev);
                    mocha.bdev = fabsf((float)bdev);
                    data.err_log_prb = model->err_log_prb;
                    mocha.lod_baf_phase = -minus_baf_phase_lod(mocha.bdev, (void *)&data);
                    // for larger bdev estimates, a model with phase might underestimate the deviation due to switch
                    // errors
                    if (mocha.bdev > self->stats.dispersion) {
                        float bdev = get_baf_bdev(baf, mocha.n_hets, imap_arr + beg[i], self->stats.dispersion);
                        if (!isnan(bdev)) mocha.bdev = bdev;
                    }
                }
                mocha.lod_baf_phase += (float)mocha.n_flips * model->flip_log_prb * (float)M_LOG10E;

                for (int j = 0; j < mocha.n_hets; j++)
                    pbaf_arr[j] = (baf[imap_arr[beg[i] + j]] - 0.5f) * (float)SIGN(path[beg[i] + j]);
                mocha.bdev_se = get_se_mean(pbaf_arr, mocha.n_hets, NULL);
                for (int j = beg[i]; j <= end[i]; j++) as[imap_arr[j]] = (int8_t)SIGN(path[j]) * gt_phase[imap_arr[j]];
            }

            mocha.type =
                mocha_type(mocha.ldev, mocha.ldev_se, mocha.bdev, mocha.bdev_se, mocha.n_hets, model->lrr_hap2dip,
                           mocha.p_arm, mocha.q_arm, model->genome_rules->is_short_arm[model->rid]);
            mocha.cf = mocha_cell_fraction(mocha.ldev, mocha.ldev_se, mocha.bdev, mocha.n_hets, mocha.type,
                                           model->lrr_hap2dip);
            mocha_table->n++;
            hts_expand(mocha_t, mocha_table->n, mocha_table->m, mocha_table->a);
            mocha_table->a[mocha_table->n - 1] = mocha;

            // update information that will be stored in the output VCF and make
            // remaining sites missing
            for (int j = a; j <= b; j++) {
                if (ldev[j] == 0) ldev[j] = float_to_int16(mocha.ldev);
                if (bdev[j] == 0) bdev[j] = float_to_int16(mocha.bdev);
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
    memcpy(self->phase_arr, as, n * sizeof(int8_t));
    free(ldev);
    free(bdev);
    free(as);
}

// computes the medoid contig for LRR regression
// TODO weight the coefficients appropriately
static int get_medoid(const float *coeffs, int n, int order) {
    int medoid_idx = -1;
    float prev = INFINITY;
    for (int i = 0; i < n; i++) {
        float next = 0.0f;
        for (int j = 0; j < n; j++) {
            if (i != j) {
                for (int k = 0; k <= order; k++) {
                    next += fabsf(coeffs[i * (order + 1) + k] - coeffs[j * (order + 1) + k]);
                }
            }
        }
        if (next < prev) {
            prev = next;
            medoid_idx = i;
        }
    }
    return medoid_idx;
}

// groups numbers in two separate distributions
static float get_lrr_cutoff(const float *v, int n) {
    if (n <= 1) return NAN;
    float *w = (float *)malloc(n * sizeof(float));
    int j = 0;
    for (int i = 0; i < n; i++) {
        if (!isnan(v[i])) w[j++] = v[i];
    }
    if (j <= 1) {
        free(w);
        return NAN;
    }
    ks_introsort_float((size_t)j, w);

    // identify a reasonable initial split allowing for some outliers
    int k = j / 2;
    int d = (int)sqrtf((float)j) - 1;
    while (k > 1 && w[k - 1] - w[d] > w[j - 1 - d] - w[k - 1]) k--;
    while (k < j && w[k] - w[d] < w[j - 1 - d] - w[k]) k++;

    // run k-means clustering EM
    while (k > 0
           && w[k - 1] - (w[(k - 1) / 2] + w[k / 2]) * 0.5f > (w[(j + k - 1) / 2] + w[(j + k) / 2]) * 0.5f - w[k - 1])
        k--;
    while (k < j && w[k] - (w[(k - 1) / 2] + w[k / 2]) * 0.5f < (w[(j + k - 1) / 2] + w[(j + k) / 2]) * 0.5f - w[k])
        k++;

    float cutoff = (k > 0 && k < j) ? (w[k - 1] + w[k]) * 0.5f : NAN;
    free(w);
    return cutoff;
}

// this function computes several contig stats and then clears the contig data from the sample
static void sample_stats(sample_t *self, const model_t *model) {
    int n = self->n;
    if (n == 0) return;
    self->n_sites += n;
    for (int i = 0; i < n; i++)
        if (self->phase_arr[i] == bcf_int8_missing) self->n_missing_gts++;

    int16_t *ad0 = self->data_arr[AD0];
    int16_t *ad1 = self->data_arr[AD1];
    float *lrr = (float *)malloc(n * sizeof(float));
    float *baf = (float *)malloc(n * sizeof(float));
    if (model->flags & WGS_DATA) {
        ad_to_lrr_baf(ad0, ad1, lrr, baf, n);
    } else {
        for (int i = 0; i < n; i++) {
            lrr[i] = int16_to_float(self->data_arr[LRR][i]);
            baf[i] = int16_to_float(self->data_arr[BAF][i]);
        }
    }
    int *imap_arr = (int *)malloc(n * sizeof(int));

    if (model->rid == model->genome_rules->x_rid) {
        int n_imap = 0;
        for (int i = 0; i < n; i++) {
            if (!isnan(baf[i])) self->n_hets++;
            int pos = model->locus_arr[self->vcf_imap_arr[i]].pos;
            if (pos > model->genome_rules->x_nonpar_beg && pos < model->genome_rules->x_nonpar_end
                && (pos < model->genome_rules->x_xtr_beg || pos > model->genome_rules->x_xtr_end)) {
                if (!isnan(baf[i])) self->x_nonpar_n_hets++;
                n_imap++;
                imap_arr[n_imap - 1] = i;
            } else if (!isnan(baf[i])) {
                if (pos <= model->genome_rules->x_nonpar_beg) {
                    self->par1_n_hets++;
                } else if (pos >= model->genome_rules->x_xtr_beg || pos <= model->genome_rules->x_xtr_end) {
                    self->xtr_n_hets++;
                } else if (pos >= model->genome_rules->x_nonpar_end) {
                    self->par2_n_hets++;
                }
            }
        }
        self->x_nonpar_lrr_median = get_median(lrr, n_imap, imap_arr);

        if (model->flags & WGS_DATA) {
            minus_lod_lkl_beta_binomial_t data = {ad0, ad1, n_imap, imap_arr};
            double x;
            kmin_brent(minus_lod_lkl_beta_binomial, 0.1, 0.2, (void *)&data, KMIN_EPS,
                       &x); // dispersions above 0.5 are not allowed
            self->x_nonpar_dispersion = (float)x;
        } else {
            self->x_nonpar_dispersion = get_sample_sd(baf, n_imap, imap_arr);
        }
    } else if (model->rid == model->genome_rules->y_rid) {
        int n_imap = 0;
        for (int i = 0; i < n; i++) {
            if (!isnan(baf[i])) self->n_hets++;
            int pos = model->locus_arr[self->vcf_imap_arr[i]].pos;
            if (pos > model->genome_rules->y_nonpar_beg && pos < model->genome_rules->y_nonpar_end
                && (pos < model->genome_rules->y_xtr_beg || pos > model->genome_rules->y_xtr_end)) {
                n_imap++;
                imap_arr[n_imap - 1] = i;
            }
        }
        self->y_nonpar_lrr_median = get_median(lrr, n_imap, imap_arr);
    } else if (model->rid == model->genome_rules->mt_rid) {
        self->mt_lrr_median = get_median(lrr, n, NULL);
    } else {
        // expand arrays if necessary
        self->n_stats++;
        hts_expand(stats_t, self->n_stats, self->m_stats, self->stats_arr);

        if (model->flags & WGS_DATA) {
            minus_lod_lkl_beta_binomial_t data = {ad0, ad1, n, NULL};
            double x;
            kmin_brent(minus_lod_lkl_beta_binomial, 0.1, 0.2, (void *)&data, KMIN_EPS,
                       &x); // dispersions above 0.5 are not allowed
            self->stats_arr[self->n_stats - 1].dispersion = (float)x;
        } else {
            self->stats_arr[self->n_stats - 1].dispersion = get_sample_sd(baf, n, NULL);
        }
        for (int i = 0; i < n; i++)
            if (!isnan(baf[i])) self->n_hets++;
        self->stats_arr[self->n_stats - 1].lrr_median = get_median(lrr, n, NULL);
        self->stats_arr[self->n_stats - 1].lrr_sd = get_sample_sd(lrr, n, NULL);

        int conc, disc;
        get_baf_conc(baf, self->phase_arr, n, &conc, &disc);
        self->stats_arr[self->n_stats - 1].baf_conc = (float)conc / (float)(conc + disc);
        self->stats_arr[self->n_stats - 1].baf_auto = get_baf_auto_corr(baf, self->phase_arr, n);
        // performs polynomial regression for LRR
        if (model->lrr_gc_order > 0 && n > model->lrr_gc_order) {
            float tss = get_tss(lrr, n);
            int ret = polyfit(lrr, model->gc_arr, n, self->vcf_imap_arr, model->lrr_gc_order,
                              self->stats_arr[self->n_stats - 1].coeffs);
            if (ret < 0) error("Polynomial regression failed\n");
            adjust_lrr(lrr, model->gc_arr, n, self->vcf_imap_arr, self->stats_arr[self->n_stats - 1].coeffs,
                       model->lrr_gc_order);
            self->stats_arr[self->n_stats - 1].coeffs[0] += get_median(lrr, n, NULL); // further adjusts by median
            float rss = get_tss(lrr, n);
            self->stats_arr[self->n_stats - 1].lrr_gc_rel_ess = 1.0f - rss / tss;
        } else if (model->lrr_gc_order != -1) {
            self->stats_arr[self->n_stats - 1].coeffs[0] = get_median(lrr, n, NULL);
            self->stats_arr[self->n_stats - 1].lrr_gc_rel_ess = NAN;
        }
        // compute autocorrelation after GC correction
        self->stats_arr[self->n_stats - 1].lrr_auto = get_lrr_auto_corr(lrr, n, NULL);
    }

    free(lrr);
    free(baf);
    free(imap_arr);
}

// this function computes the median of contig stats
static void sample_summary(sample_t *self, int n, model_t *model, int compute_gender) {
    float *tmp_arr = (float *)malloc(n * sizeof(float));
    int m_tmp = n;

    for (int i = 0; i < n; i++) {
        hts_expand(float, self[i].n_stats *(model->lrr_gc_order < 0 ? 1 : model->lrr_gc_order + 1), m_tmp, tmp_arr);

        for (int j = 0; j < self[i].n_stats; j++) tmp_arr[j] = self[i].stats_arr[j].lrr_median;
        self[i].stats.lrr_median = get_median(tmp_arr, self[i].n_stats, NULL);
        for (int j = 0; j < self[i].n_stats; j++) tmp_arr[j] = self[i].stats_arr[j].lrr_sd;
        self[i].stats.lrr_sd = get_median(tmp_arr, self[i].n_stats, NULL);
        for (int j = 0; j < self[i].n_stats; j++) tmp_arr[j] = self[i].stats_arr[j].lrr_auto;
        self[i].stats.lrr_auto = get_median(tmp_arr, self[i].n_stats, NULL);

        for (int j = 0; j < self[i].n_stats; j++) tmp_arr[j] = self[i].stats_arr[j].dispersion;
        self[i].stats.dispersion = get_median(tmp_arr, self[i].n_stats, NULL);
        for (int j = 0; j < self[i].n_stats; j++) tmp_arr[j] = self[i].stats_arr[j].baf_conc;
        self[i].stats.baf_conc = get_median(tmp_arr, self[i].n_stats, NULL);
        for (int j = 0; j < self[i].n_stats; j++) tmp_arr[j] = self[i].stats_arr[j].baf_auto;
        self[i].stats.baf_auto = get_median(tmp_arr, self[i].n_stats, NULL);

        self[i].adjlrr_sd = self[i].stats.lrr_sd;
        if (model->lrr_gc_order == 0) {
            for (int j = 0; j < self[i].n_stats; j++) tmp_arr[j] = self[i].stats_arr[j].coeffs[0];
            self[i].stats.coeffs[0] = get_median(tmp_arr, self[i].n_stats, NULL);
        } else if (model->lrr_gc_order > 0 && self[i].n_stats > 0) {
            for (int j = 0; j < self[i].n_stats; j++)
                for (int k = 0; k <= model->lrr_gc_order; k++)
                    tmp_arr[j * (model->lrr_gc_order + 1) + k] = self[i].stats_arr[j].coeffs[k];
            int medoid_idx = get_medoid(tmp_arr, self[i].n_stats, model->lrr_gc_order);
            for (int k = 0; k <= model->lrr_gc_order; k++)
                self[i].stats.coeffs[k] = tmp_arr[medoid_idx * (model->lrr_gc_order + 1) + k];
            for (int j = 0; j < self[i].n_stats; j++) tmp_arr[j] = self[i].stats_arr[j].lrr_gc_rel_ess;
            self[i].stats.lrr_gc_rel_ess = get_median(tmp_arr, self[i].n_stats, NULL);
            self[i].adjlrr_sd *= sqrtf(1.0f - self[i].stats.lrr_gc_rel_ess); // not perfect, but good enough(?)
        }
        free(self[i].stats_arr);
    }

    if (compute_gender) {
        if (model->flags & WGS_DATA) {
            if (isnan(model->lrr_cutoff)) model->lrr_cutoff = -0.3f; // arbitrary cutoff between -M_LN2 and 0
        }

        // determine LRR cutoff between haploid and diploid
        if (isnan(model->lrr_cutoff)) {
            int j = 0;
            for (int i = 0; i < n; i++)
                if (!isnan(self[i].x_nonpar_lrr_median))
                    tmp_arr[j++] = isnan(self[i].stats.lrr_median)
                                       ? self[i].x_nonpar_lrr_median
                                       : self[i].x_nonpar_lrr_median - self[i].stats.lrr_median;
            model->lrr_cutoff = get_lrr_cutoff(tmp_arr, j);
        }

        // determine gender of samples
        for (int i = 0; i < n; i++) {
            float tmp = isnan(self[i].stats.lrr_median) ? self[i].x_nonpar_lrr_median
                                                        : self[i].x_nonpar_lrr_median - self[i].stats.lrr_median;
            if (tmp < model->lrr_cutoff)
                self[i].computed_gender = GENDER_MALE;
            else if (tmp > model->lrr_cutoff)
                self[i].computed_gender = GENDER_FEMALE;
        }
    }

    free(tmp_arr);
}

static void mocha_print_stats(FILE *restrict stream, const sample_t *self, int n, int lrr_gc_order,
                              const bcf_hdr_t *hdr, int flags) {
    if (stream == NULL) return;
    char gender[4];
    gender[GENDER_UNKNOWN] = 'U';
    gender[GENDER_MALE] = 'M';
    gender[GENDER_FEMALE] = 'F';
    gender[GENDER_KLINEFELTER] = 'K';
    fputs("sample_id", stream);
    fputs("\tcomputed_gender", stream);
    fputs("\tcall_rate", stream);
    fputs(flags & WGS_DATA ? "\tcov_median" : "\tlrr_median", stream);
    fputs(flags & WGS_DATA ? "\tcov_sd" : "\tlrr_sd", stream);
    fputs(flags & WGS_DATA ? "\tcov_auto" : "\tlrr_auto", stream);
    fputs(flags & WGS_DATA ? "\tbaf_corr" : "\tbaf_sd", stream);
    fputs("\tbaf_conc", stream);
    fputs("\tbaf_auto", stream);
    fputs("\tn_sites", stream);
    fputs("\tn_hets", stream);
    fputs("\tx_nonpar_n_hets", stream);
    fputs("\tpar1_n_hets", stream);
    fputs("\txtr_n_hets", stream);
    fputs("\tpar2_n_hets", stream);
    fputs("\tx_nonpar_baf_corr", stream);
    fputs(flags & WGS_DATA ? "\tx_nonpar_cov_median" : "\tx_nonpar_lrr_median", stream);
    fputs(flags & WGS_DATA ? "\ty_nonpar_cov_median" : "\ty_nonpar_lrr_median", stream);
    fputs(flags & WGS_DATA ? "\tmt_cov_median" : "\tmt_lrr_median", stream);
    fputs("\tlrr_gc_rel_ess", stream);
    for (int j = 0; j <= lrr_gc_order; j++) fprintf(stream, "\tlrr_gc_%d", j);
    fputc('\n', stream);
    for (int i = 0; i < n; i++) {
        fputs(bcf_hdr_int2id(hdr, BCF_DT_SAMPLE, self[i].idx), stream);
        fprintf(stream, "\t%c", gender[self[i].computed_gender]);
        fprintf(stream, "\t%.4f",
                self[i].stats.call_rate ? self[i].stats.call_rate
                                        : 1.0 - (float)self[i].n_missing_gts / (float)self[i].n_sites);
        fprintf(stream, "\t%.4f", flags & WGS_DATA ? expf(self[i].stats.lrr_median) : self[i].stats.lrr_median);
        fprintf(stream, "\t%.4f",
                flags & WGS_DATA ? expf(self[i].stats.lrr_median) * self[i].stats.lrr_sd : self[i].stats.lrr_sd);
        fprintf(stream, "\t%.4f", self[i].stats.lrr_auto);
        fprintf(stream, "\t%.4f", self[i].stats.dispersion);
        fprintf(stream, "\t%.4f", self[i].stats.baf_conc);
        fprintf(stream, "\t%.4f", self[i].stats.baf_auto);
        fprintf(stream, "\t%d", self[i].n_sites);
        fprintf(stream, "\t%d", self[i].n_hets);
        fprintf(stream, "\t%d", self[i].x_nonpar_n_hets);
        fprintf(stream, "\t%d", self[i].par1_n_hets);
        fprintf(stream, "\t%d", self[i].xtr_n_hets);
        fprintf(stream, "\t%d", self[i].par2_n_hets);
        fprintf(stream, "\t%.4f", self[i].x_nonpar_dispersion);
        fprintf(stream, "\t%.4f", flags & WGS_DATA ? expf(self[i].x_nonpar_lrr_median) : self[i].x_nonpar_lrr_median);
        fprintf(stream, "\t%.4f", flags & WGS_DATA ? expf(self[i].y_nonpar_lrr_median) : self[i].y_nonpar_lrr_median);
        fprintf(stream, "\t%.4f", flags & WGS_DATA ? expf(self[i].mt_lrr_median) : self[i].mt_lrr_median);
        fprintf(stream, "\t%.4f", self[i].stats.lrr_gc_rel_ess);
        for (int j = 0; j <= lrr_gc_order; j++) fprintf(stream, "\t%.4f", self[i].stats.coeffs[j]);
        fputc('\n', stream);
    }
    if (stream != stdout && stream != stderr) fclose(stream);
}

/*********************************
 * VCF READ AND WRITE METHODS    *
 *********************************/

// moves synced bcf reader to first useful line for reader0
static inline int bcf_sr_next_line_reader0(bcf_srs_t *sr) {
    int nret = bcf_sr_next_line(sr);
    while (nret > 0 && !bcf_sr_has_line(sr, 0)) nret = bcf_sr_next_line(sr);
    return nret;
}

// write one contig
static int put_contig(bcf_srs_t *sr, const sample_t *sample, const model_t *model, htsFile *out_fh,
                      bcf_hdr_t *out_hdr) {
    int rid = model->rid;
    bcf_hdr_t *hdr = bcf_sr_get_header(sr, 0);
    bcf_sr_seek(sr, bcf_hdr_id2name(hdr, rid), 0);
    int nsmpl = bcf_hdr_nsamples(out_hdr);

    int *synced_iter = (int *)calloc(nsmpl, sizeof(int));
    float *bdev = (float *)calloc(nsmpl, sizeof(float));
    float *ldev = (float *)calloc(nsmpl, sizeof(float));
    int *as = (int *)calloc(nsmpl, sizeof(int));

    int i;
    for (i = 0; bcf_sr_next_line_reader0(sr); i++) {
        bcf1_t *line = bcf_sr_get_line(sr, 0);
        if (rid != line->rid) break;

        for (int j = 0; j < nsmpl; j++) {
            while (synced_iter[j] < sample[j].n - 1 && sample[j].vcf_imap_arr[synced_iter[j]] < i) synced_iter[j]++;
            if (sample[j].vcf_imap_arr[synced_iter[j]] == i) {
                if (sample[j].data_arr[LDEV])
                    ldev[sample[j].idx] = int16_to_float(sample[j].data_arr[LDEV][synced_iter[j]]);
                if (sample[j].data_arr[BDEV])
                    bdev[sample[j].idx] = int16_to_float(sample[j].data_arr[BDEV][synced_iter[j]]);
                if (sample[j].phase_arr) as[sample[j].idx] = sample[j].phase_arr[synced_iter[j]];
            } else {
                // if no match variant found, match the end of the contig or keep conservative
                if (i == 0 && sample[j].data_arr[BDEV])
                    bdev[sample[j].idx] = int16_to_float(sample[j].data_arr[BDEV][0]);
                if (i == 0 && sample[j].data_arr[LDEV])
                    ldev[sample[j].idx] = int16_to_float(sample[j].data_arr[LDEV][0]);
                if (sample[j].data_arr[BDEV] && int16_to_float(sample[j].data_arr[BDEV][synced_iter[j]]) == 0.0f)
                    bdev[sample[j].idx] = 0.0f;
                if (sample[j].data_arr[LDEV] && int16_to_float(sample[j].data_arr[LDEV][synced_iter[j]]) == 0.0f)
                    ldev[sample[j].idx] = 0.0f;
                if (sample[j].phase_arr) as[sample[j].idx] = 0;
            }
        }
        if (!(model->flags & WGS_DATA)) {
            bcf_update_info_float(out_hdr, line, "ADJ_COEFF", model->locus_arr[i].adjust, 9);
        }
        if (!(model->flags & NO_ANNOT)) {
            bcf_update_format_float(out_hdr, line, "Ldev", ldev, (int)nsmpl);
            bcf_update_format_float(out_hdr, line, "Bdev", bdev, (int)nsmpl);
        }
        bcf_update_format_int32(out_hdr, line, "AS", as, (int)nsmpl);

        if (bcf_write(out_fh, out_hdr, line) < 0) error("Unable to write to output VCF file\n");
    }

    free(synced_iter);
    free(ldev);
    free(bdev);
    free(as);

    return i;
}

// write header
static bcf_hdr_t *print_hdr(htsFile *out_fh, bcf_hdr_t *hdr, int argc, char *argv[], int record_cmd_line, int flags) {
    bcf_hdr_t *out_hdr = bcf_hdr_dup(hdr);
    int adj_coeff_id = bcf_hdr_id2int(hdr, BCF_DT_ID, "ADJ_COEFF");
    if (!(flags & WGS_DATA) && !bcf_hdr_idinfo_exists(hdr, BCF_HL_INFO, adj_coeff_id))
        bcf_hdr_append(out_hdr,
                       "##INFO=<ID=ADJ_COEFF,Number=9,Type=Float,Description=\"Adjust coefficients "
                       "(order=AA_BAF0,AA_BAF1,AA_LRR0,AB_BAF0,AB_BAF1,AB_LRR0,BB_BAF0,BB_BAF1,BB_LRR0)"
                       "\">");
    if (!(flags & NO_ANNOT)) {
        int ldev_id = bcf_hdr_id2int(hdr, BCF_DT_ID, "Ldev");
        if (!bcf_hdr_idinfo_exists(hdr, BCF_HL_FMT, ldev_id))
            bcf_hdr_append(out_hdr,
                           "##FORMAT=<ID=Ldev,Number=1,Type=Float,Description=\"LRR deviation "
                           "due to chromosomal alteration\">");

        int bdev_id = bcf_hdr_id2int(hdr, BCF_DT_ID, "Bdev");
        if (!bcf_hdr_idinfo_exists(hdr, BCF_HL_FMT, bdev_id))
            bcf_hdr_append(out_hdr,
                           "##FORMAT=<ID=Bdev,Number=1,Type=Float,Description=\"BAF deviation "
                           "due to chromosomal alteration\">");
    }
    int as_id = bcf_hdr_id2int(hdr, BCF_DT_ID, "AS");
    if (!bcf_hdr_idinfo_exists(hdr, BCF_HL_FMT, as_id))
        bcf_hdr_append(out_hdr,
                       "##FORMAT=<ID=AS,Number=1,Type=Integer,Description=\"Allelic "
                       "shift (1/-1 if the alternate allele is over/under represented)\">");
    if (record_cmd_line) bcf_hdr_append_version(out_hdr, argc, argv, "bcftools_plugin");
    if (bcf_hdr_write(out_fh, out_hdr) < 0) error("Unable to write to output VCF file\n");
    return out_hdr;
}

// retrieve genotypes as NC, AA, AB, BB
// assumes little endian architecture
static int bcf_get_ab_genotypes(bcf_fmt_t *fmt, int8_t *gts, int nsmpl, int allele_a, int allele_b) {
    if (!fmt || fmt->n != 2) return 0;

#define BRANCH(type_t, bcf_type_vector_end)                                                                            \
    {                                                                                                                  \
        type_t *p = (type_t *)fmt->p;                                                                                  \
        for (int i = 0; i < nsmpl; i++, p += 2) {                                                                      \
            if (p[0] == bcf_type_vector_end || bcf_gt_is_missing(p[0]) || p[1] == bcf_type_vector_end                  \
                || bcf_gt_is_missing(p[1])) {                                                                          \
                gts[i] = GT_NC;                                                                                        \
            } else {                                                                                                   \
                type_t allele_0 = bcf_gt_allele(p[0]);                                                                 \
                type_t allele_1 = bcf_gt_allele(p[1]);                                                                 \
                if (allele_0 == allele_a && allele_1 == allele_a)                                                      \
                    gts[i] = GT_AA;                                                                                    \
                else if (allele_0 == allele_b && allele_1 == allele_b)                                                 \
                    gts[i] = GT_BB;                                                                                    \
                else if ((allele_0 == allele_a && allele_1 == allele_b)                                                \
                         || (allele_0 == allele_b && allele_1 == allele_a))                                            \
                    gts[i] = GT_AB;                                                                                    \
                else                                                                                                   \
                    return -1;                                                                                         \
            }                                                                                                          \
        }                                                                                                              \
    }
    switch (fmt->type) {
    case BCF_BT_INT8:
        BRANCH(int8_t, bcf_int8_vector_end);
        break;
    case BCF_BT_INT16:
        BRANCH(int16_t, bcf_int16_vector_end);
        break;
    case BCF_BT_INT32:
        BRANCH(int32_t, bcf_int32_vector_end);
        break;
    default:
        error("Unexpected type %d\n", fmt->type);
    }
#undef BRANCH

    return 1;
}

// read one contig
static int get_contig(bcf_srs_t *sr, sample_t *sample, model_t *model) {
    int rid = model->rid;
    bcf_hdr_t *hdr = bcf_sr_get_header(sr, 0);
    bcf_sr_seek(sr, bcf_hdr_id2name(hdr, rid), 0);

    bcf_fmt_t *baf_fmt = NULL, *lrr_fmt = NULL;
    bcf_info_t *info;
    int nsmpl = bcf_hdr_nsamples(hdr);

    model->n = 0;
    model->n_flipped = 0;
    for (int j = 0; j < nsmpl; j++) sample[j].n = 0;

    if (!(model->flags & USE_NO_RULES_CHRS) && model->genome_rules->cen_beg[rid] == 0
        && model->genome_rules->cen_end[rid] == 0 && rid != model->genome_rules->mt_rid)
        return 0;

    int8_t *gts = (int8_t *)malloc(nsmpl * sizeof(int8_t));
    int8_t *phase_arr = (int8_t *)malloc(nsmpl * sizeof(int8_t));
    int *imap_arr = (int *)malloc(nsmpl * sizeof(int));
    int16_t *gt0 = (int16_t *)malloc(nsmpl * sizeof(int16_t));
    int16_t *gt1 = (int16_t *)malloc(nsmpl * sizeof(int16_t));
    int16_t *ad0 = (int16_t *)malloc(nsmpl * sizeof(int16_t));
    int16_t *ad1 = (int16_t *)malloc(nsmpl * sizeof(int16_t));
    int *last_het_pos = (int *)calloc(nsmpl, sizeof(int));
    int *last_pos = (int *)calloc(nsmpl, sizeof(int));

    int i;
    for (i = 0; bcf_sr_next_line_reader0(sr); i++) {
        bcf1_t *line = bcf_sr_get_line(sr, 0);
        if (model->filter) {
            int ret = filter_test(model->filter, line, NULL);
            if ((model->filter_logic == FLT_INCLUDE && !ret) || ret) continue;
        }
        if (rid != line->rid) break;
        int pos = line->pos + 1;

        hts_expand(locus_t, i + 1, model->m_locus, model->locus_arr);

        model->locus_arr[i].pos = pos;

        hts_expand(float, i + 1, model->m_gc, model->gc_arr);
        if (model->gc_id >= 0 && (info = bcf_get_info_id(line, model->gc_id)))
            model->gc_arr[i] = info->v1.f;
        else
            model->gc_arr[i] = NAN;

        // if failing inclusion/exclusion requirement, skip line
        if ((model->flags & FLT_EXCLUDE) && bcf_sr_get_line(sr, 1)) continue;
        if ((model->flags & FLT_INCLUDE) && !bcf_sr_get_line(sr, 1)) continue;

        // if site falls in short arm or centromere regions skip line
        if (!(model->flags & USE_SHORT_ARMS) && model->genome_rules->is_short_arm[rid]
            && pos < model->genome_rules->cen_beg[rid])
            continue;
        if (!(model->flags & USE_CENTROMERES) && pos > model->genome_rules->cen_beg[rid]
            && pos < model->genome_rules->cen_end[rid])
            continue;

        if (model->mhc_idx
            && regidx_overlap(model->mhc_idx, bcf_hdr_id2name(hdr, line->rid), line->pos, line->pos + 1, NULL))
            continue;
        if (model->kir_idx
            && regidx_overlap(model->kir_idx, bcf_hdr_id2name(hdr, line->rid), line->pos, line->pos + 1, NULL))
            continue;

        bcf_fmt_t *gt_fmt = bcf_get_fmt_id(line, model->gt_id);
        if (!bcf_get_genotype_phase(gt_fmt, phase_arr, nsmpl)) continue;

        // if neither AD nor LRR and BAF formats are present, skip line
        if (model->flags & WGS_DATA) {
            if (!bcf_get_unphased_genotype_alleles(gt_fmt, gt0, gt1, nsmpl)) continue;
            if (!bcf_get_allelic_depth(bcf_get_fmt_id(line, model->ad_id), gt0, gt1, ad0, ad1, nsmpl)) continue;
        } else {
            if ((info = bcf_get_info_id(line, model->allele_a_id)))
                model->locus_arr[i].allele_a = info->v1.i;
            else
                error("Error: ALLELE_A missing at position %s:%" PRId64 "\n", bcf_hdr_id2name(hdr, line->rid),
                      line->pos + 1);
            if ((info = bcf_get_info_id(line, model->allele_b_id)))
                model->locus_arr[i].allele_b = info->v1.i;
            else
                error("Error: ALLELE_B missing at position %s:%" PRId64 "\n", bcf_hdr_id2name(hdr, line->rid),
                      line->pos + 1);
            // missing ALLELE_A and ALLELE_B information
            if (model->locus_arr[i].allele_a < 0 || model->locus_arr[i].allele_b < 0) continue;
            // if there are no genotypes, skip line
            if (bcf_get_ab_genotypes(gt_fmt, gts, nsmpl, model->locus_arr[i].allele_a, model->locus_arr[i].allele_b)
                < 0)
                error("Error: site %s:%" PRId64 " contains non-AB alleles\n", bcf_hdr_id2name(hdr, line->rid),
                      line->pos + 1);

            if (!(lrr_fmt = bcf_get_fmt_id(line, model->lrr_id)) || !(baf_fmt = bcf_get_fmt_id(line, model->baf_id)))
                continue;

            for (int j = 0; j < nsmpl; j++) {
                if (bcf_float_is_missing(((float *)lrr_fmt->p)[sample[j].idx]))
                    ((float *)lrr_fmt->p)[sample[j].idx] = NAN;
                if (bcf_float_is_missing(((float *)baf_fmt->p)[sample[j].idx]))
                    ((float *)baf_fmt->p)[sample[j].idx] = NAN;
            }

            int is_x_nonpar = rid == model->genome_rules->x_rid && pos > model->genome_rules->x_nonpar_beg
                              && pos < model->genome_rules->x_nonpar_end
                              && (pos < model->genome_rules->x_xtr_beg || pos > model->genome_rules->x_xtr_end);
            int is_y_or_mt = rid == model->genome_rules->y_rid || rid == model->genome_rules->mt_rid;

            // adjust cluster centers and slopes, inspired by
            // (i) Staaf, J. et al. Normalization of Illumina Infinium whole-genome
            // SNP data improves copy number estimates and allelic intensity ratios.
            // BMC Bioinformatics 9, 409 (2008) & (ii) Mayrhofer, M., Viklund, B. &
            // Isaksson, A. Rawcopy: Improved copy number analysis with Affymetrix
            // arrays. Sci Rep 6, 36158 (2016)
            for (int gt = GT_AA; gt <= GT_BB; gt++) {
                if (is_y_or_mt) continue;
                int k = 0;
                for (int j = 0; j < nsmpl; j++) {
                    if (is_x_nonpar && sample[j].computed_gender == GENDER_MALE) continue;
                    if (gts[sample[j].idx] == gt) imap_arr[k++] = sample[j].idx;
                }
                float baf_b = 0.0f, baf_m = 0.0f, lrr_b = 0.0f;
                if (model->regress_baf_lrr != -1 && model->regress_baf_lrr <= k) {
                    float xss = 0.0f, yss = 0.0f, xyss = 0.0f;
                    get_cov((float *)lrr_fmt->p, (float *)baf_fmt->p, k, imap_arr, &xss, &yss, &xyss);
                    baf_m = xss == 0.0f ? 0.0f : xyss / xss;
                    for (int j = 0; j < k; j++)
                        ((float *)baf_fmt->p)[imap_arr[j]] -= baf_m * ((float *)lrr_fmt->p)[imap_arr[j]];
                }
                if (model->adj_baf_lrr != -1 && model->adj_baf_lrr <= k) {
                    baf_b = get_median((float *)baf_fmt->p, k, imap_arr) - (float)(gt - 1) * 0.5f;
                    if (isnan(baf_b)) baf_b = 0.0f;
                    for (int j = 0; j < k; j++) ((float *)baf_fmt->p)[imap_arr[j]] -= baf_b;
                    lrr_b = get_median((float *)lrr_fmt->p, k, imap_arr);
                    if (isnan(lrr_b)) lrr_b = 0.0f;
                    for (int j = 0; j < k; j++) ((float *)lrr_fmt->p)[imap_arr[j]] -= lrr_b;
                }
                // corrects males on X nonPAR
                if (is_x_nonpar) {
                    for (int j = 0; j < nsmpl; j++) {
                        if (sample[j].computed_gender != GENDER_MALE) continue;
                        if (gts[sample[j].idx] == gt) {
                            ((float *)baf_fmt->p)[sample[j].idx] -=
                                baf_m * ((float *)lrr_fmt->p)[sample[j].idx] + baf_b;
                            ((float *)lrr_fmt->p)[sample[j].idx] -= lrr_b;
                        }
                    }
                }
                model->locus_arr[i].adjust[gt - 1][0] = baf_b;
                model->locus_arr[i].adjust[gt - 1][1] = baf_m;
                model->locus_arr[i].adjust[gt - 1][2] = lrr_b;
            }

            // if allele A index is bigger than allele B index flip the BAF to make
            // sure it refers to the highest allele
            if (model->locus_arr[i].allele_a > model->locus_arr[i].allele_b) {
                for (int j = 0; j < nsmpl; j++)
                    ((float *)baf_fmt->p)[sample[j].idx] = 1.0f - ((float *)baf_fmt->p)[sample[j].idx];
            }

            for (int j = 0; j < nsmpl; j++)
                if (gts[j] != GT_AB || (is_x_nonpar && sample[j].computed_gender == GENDER_MALE) || is_y_or_mt)
                    ((float *)baf_fmt->p)[sample[j].idx] = NAN;
        }

        // read line in memory
        model->n++;
        for (int j = 0; j < nsmpl; j++) {
            if (model->flags & WGS_DATA) {
                if (ad0[j] == bcf_int16_missing && ad1[j] == bcf_int16_missing) continue;

                // site too close to last het site or hom site too close to last site
                if ((pos < last_het_pos[j] + model->min_dst)
                    || ((phase_arr[j] == bcf_int8_missing || phase_arr[j] == bcf_int8_vector_end)
                        && (pos < last_pos[j] + model->min_dst)))
                    continue;

                // substitute the last hom site with the current het site
                if (pos < last_pos[j] + model->min_dst) sample[j].n--;

                if (phase_arr[j] != bcf_int8_missing || phase_arr[j] == bcf_int8_vector_end) last_het_pos[j] = pos;
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

            } else {
                float lrr = ((float *)lrr_fmt->p)[sample[j].idx];
                float baf = ((float *)baf_fmt->p)[sample[j].idx];
                if (isnan(lrr) && isnan(baf)) continue;

                sample[j].n++;
                hts_expand(int, sample[j].n, sample[j].m_vcf_imap, sample[j].vcf_imap_arr);
                hts_expand(int8_t, sample[j].n, sample[j].m_phase, sample[j].phase_arr);
                hts_expand(int16_t, sample[j].n, sample[j].m_data[LRR], sample[j].data_arr[LRR]);
                hts_expand(int16_t, sample[j].n, sample[j].m_data[BAF], sample[j].data_arr[BAF]);
                sample[j].vcf_imap_arr[sample[j].n - 1] = i;
                sample[j].phase_arr[sample[j].n - 1] = phase_arr[j];
                sample[j].data_arr[LRR][sample[j].n - 1] = float_to_int16(lrr);
                sample[j].data_arr[BAF][sample[j].n - 1] = float_to_int16(baf);
            }
        }
    }
    free(gts);
    free(phase_arr);
    free(imap_arr);
    free(gt0);
    free(gt1);
    free(ad0);
    free(ad1);
    free(last_het_pos);
    free(last_pos);

    return i;
}

static int read_stats(sample_t *samples, const bcf_hdr_t *hdr, const char *fn, int lrr_gc_order, int flags) {
    htsFile *fp = hts_open(fn, "r");
    if (fp == NULL) error("Could not open %s: %s\n", fn, strerror(errno));

    kstring_t str = {0, 0, NULL};
    if (hts_getline(fp, KS_SEP_LINE, &str) <= 0) error("Empty file: %s\n", fn);
    tsv_t *tsv = tsv_init_delimiter(str.s, '\t');
    sample_t sample;
    memset((void *)&sample, 0, sizeof(sample));
    if (tsv_register(tsv, "sample_id", tsv_read_sample_id, (void *)&sample.idx) < 0)
        error("File %s is missing the sample_id column\n", fn);
    if (tsv_register(tsv, "computed_gender", tsv_read_computed_gender, (void *)&sample.computed_gender) < 0)
        error("File %s is missing the computed_gender column\n", fn);
    if (tsv_register(tsv, "call_rate", tsv_read_float, (void *)&sample.stats.call_rate) < 0)
        error("File %s is missing the call_rate column\n", fn);
    int lrr_median = tsv_register(tsv, flags & WGS_DATA ? "cov_median" : "lrr_median", tsv_read_float,
                                  (void *)&sample.stats.lrr_median);
    int lrr_sd =
        tsv_register(tsv, flags & WGS_DATA ? "cov_sd" : "lrr_sd", tsv_read_float, (void *)&sample.stats.lrr_sd);
    if (tsv_register(tsv, flags & WGS_DATA ? "cov_auto" : "lrr_auto", tsv_read_float, (void *)&sample.stats.lrr_auto)
        < 0)
        sample.stats.lrr_auto = NAN;
    int baf_sd =
        tsv_register(tsv, flags & WGS_DATA ? "baf_corr" : "baf_sd", tsv_read_float, (void *)&sample.stats.dispersion);
    int baf_conc = tsv_register(tsv, "baf_conc", tsv_read_float, (void *)&sample.stats.baf_conc);
    if (tsv_register(tsv, "baf_auto", tsv_read_float, (void *)&sample.stats.baf_auto) < 0) sample.stats.baf_auto = NAN;
    tsv_register(tsv, "n_sites", tsv_read_integer, (void *)&sample.n_sites);
    tsv_register(tsv, "n_hets", tsv_read_integer, (void *)&sample.n_hets);
    tsv_register(tsv, "x_nonpar_n_hets", tsv_read_integer, (void *)&sample.x_nonpar_n_hets);
    tsv_register(tsv, "par1_n_hets", tsv_read_integer, (void *)&sample.par1_n_hets);
    tsv_register(tsv, "xtr_n_hets", tsv_read_integer, (void *)&sample.xtr_n_hets);
    tsv_register(tsv, "par2_n_hets", tsv_read_integer, (void *)&sample.par2_n_hets);
    if (tsv_register(tsv, "x_nonpar_baf_corr", tsv_read_float, (void *)&sample.x_nonpar_dispersion) < 0)
        sample.x_nonpar_dispersion = NAN;
    if (tsv_register(tsv, flags & WGS_DATA ? "x_nonpar_cov_median" : "x_nonpar_lrr_median", tsv_read_float,
                     (void *)&sample.x_nonpar_lrr_median)
        < 0)
        sample.x_nonpar_lrr_median = NAN;
    if (tsv_register(tsv, flags & WGS_DATA ? "y_nonpar_cov_median" : "y_nonpar_lrr_median", tsv_read_float,
                     (void *)&sample.y_nonpar_lrr_median)
        < 0)
        sample.y_nonpar_lrr_median = NAN;
    if (tsv_register(tsv, flags & WGS_DATA ? "mt_cov_median" : "mt_lrr_median", tsv_read_float,
                     (void *)&sample.mt_lrr_median)
        < 0)
        sample.mt_lrr_median = NAN;
    int lrr_gc_rel_ess = tsv_register(tsv, "lrr_gc_rel_ess", tsv_read_float, (void *)&sample.stats.lrr_gc_rel_ess);
    if (lrr_sd < 0 || lrr_gc_rel_ess < 0) sample.adjlrr_sd = NAN;
    int lrr_gc_coeffs = 0;
    for (int i = 0; i <= lrr_gc_order; i++) {
        str.l = 0;
        ksprintf(&str, "lrr_gc_%d", i);
        if (!(tsv_register(tsv, str.s, tsv_read_float, (void *)&sample.stats.coeffs[i]) < 0)) lrr_gc_coeffs++;
    }
    int i = 0;
    while (hts_getline(fp, KS_SEP_LINE, &str) > 0) {
        if (!tsv_parse_delimiter(tsv, (bcf1_t *)hdr, str.s, '\t')) {
            int idx = sample.idx;
            if (idx < 0) continue;
            if (flags & WGS_DATA) {
                sample.stats.lrr_sd = sample.stats.lrr_sd / sample.stats.lrr_median;
                sample.stats.lrr_median = logf(sample.stats.lrr_median);
                sample.x_nonpar_lrr_median = logf(sample.x_nonpar_lrr_median);
                sample.y_nonpar_lrr_median = logf(sample.y_nonpar_lrr_median);
                sample.mt_lrr_median = logf(sample.mt_lrr_median);
            }
            if (!(lrr_sd < 0)) {
                sample.adjlrr_sd = sample.stats.lrr_sd;
                if (!(lrr_gc_rel_ess < 0)) sample.adjlrr_sd *= sqrtf(1.0f - sample.stats.lrr_gc_rel_ess);
            }
            memcpy((void *)&samples[idx], (const void *)&sample, sizeof(sample_t));
            i++;
        } else {
            error("Could not parse line: %s\n", str.s);
        }
    }
    tsv_destroy(tsv);
    free(str.s);
    hts_close(fp);
    if (i < bcf_hdr_nsamples(hdr)) error("File %s is missing samples\n", fn);
    int ret = (lrr_median < 0) || (lrr_sd < 0) || (baf_sd < 0) || (baf_conc < 0) || (lrr_gc_rel_ess < 0)
              || (lrr_gc_coeffs != lrr_gc_order + 1);
    return ret;
}

/*********************************
 * PLUGIN CODE                   *
 *********************************/

const char *about(void) { return "Runs SNP-based method for detection of MOsaic CHromosomal Alterations (MoChA).\n"; }

static const char *usage_text(void) {
    return "\n"
           "About:   MOsaic CHromosomal Alterations caller, requires phased genotypes (GT)\n"
           "         and either B-allele frequency (BAF) and Log R Ratio intensity (LRR)\n"
           "         or allelic depth coverage (AD). (version " MOCHA_VERSION
           " https://github.com/freeseek/mocha)\n"
           "Usage:   bcftools +mocha [OPTIONS] <in.vcf.gz>\n"
           "\n"
           "Required options:\n"
           "    -g, --genome <assembly>[?]      predefined genome reference rules, 'list' to print available settings, "
           "append '?' for details\n"
           "    -G, --genome-file <file>        genome reference rules, space/tab-delimited CHROM:FROM-TO,TYPE\n"
           "\n"
           "General Options:\n"
           "    -v, --variants [^]<file>        tabix-indexed [compressed] VCF/BCF file containing variants\n"
           "    -f, --apply-filters <list>      require at least one of the listed FILTER strings (e.g. \"PASS,.\")\n"
           "                                    to include (or exclude with \"^\" prefix) in the analysis\n"
           "    -e, --exclude <expr>            exclude sites for which the expression is true\n"
           "    -i, --include <expr>            select sites for which the expression is true\n"
           "    -r, --regions <region>          restrict to comma-separated list of regions\n"
           "    -R, --regions-file <file>       restrict to regions listed in a file\n"
           "        --regions-overlap 0|1|2     Include if POS in the region (0), record overlaps (1), variant "
           "overlaps (2) [1]\n"
           "    -t, --targets [^]<region>       restrict to comma-separated list of regions. Exclude regions with "
           "\"^\" prefix\n"
           "    -T, --targets-file [^]<file>    restrict to regions listed in a file. Exclude regions with \"^\" "
           "prefix\n"
           "        --targets-overlap 0|1|2     Include if POS in the region (0), record overlaps (1), variant "
           "overlaps (2) [0]\n"
           "    -s, --samples [^]<list>         comma separated list of samples to include (or exclude with \"^\" "
           "prefix)\n"
           "    -S, --samples-file [^]<file>    file of samples to include (or exclude with \"^\" prefix)\n"
           "        --force-samples             only warn about unknown subset samples\n"
           "        --input-stats <file>        input samples genome-wide statistics file\n"
           "        --only-stats                compute genome-wide statistics without detecting mosaic chromosomal "
           "alterations\n"
           "    -p  --cnp <file>                list of regions to genotype in BED format\n"
           "        --mhc <region>              MHC region to exclude from analysis (will be retained in the output)\n"
           "        --kir <region>              KIR region to exclude from analysis (will be retained in the output)\n"
           "        --threads <int>             number of extra output compression threads [0]\n"
           "\n"
           "Output Options:\n"
           "    -o, --output <file>             write output to a file [no output]\n"
           "    -O, --output-type u|b|v|z[0-9]  u/b: un/compressed BCF, v/z: un/compressed VCF, 0-9: compression level "
           "[v]\n"
           "        --no-version                do not append version and command line to the header\n"
           "    -a  --no-annotations            omit Ldev and Bdev FORMAT from output VCF (requires --output)\n"
           "        --no-log                    suppress progress report on standard error\n"
           "    -l  --log <file>                write log to file [standard error]\n"
           "    -c, --calls <file>              write chromosomal alterations calls table to a file [standard output]\n"
           "    -z  --stats <file>              write samples genome-wide statistics table to a file [no output]\n"
           "    -u, --ucsc-bed <file>           write UCSC bed track to a file [no output]\n"
           "\n"
           "HMM Options:\n"
           "        --bdev-LRR-BAF <list>       comma separated list of inverse BAF deviations for LRR+BAF model "
           "[" BDEV_LRR_BAF_DFLT
           "]\n"
           "        --bdev-BAF-phase <list>     comma separated list of inverse BAF deviations for BAF+phase model\n"
           "                                    [" BDEV_BAF_PHASE_DFLT
           "]\n"
           "        --min-dist <int>            minimum base pair distance between consecutive sites for WGS data "
           "[" MIN_DST_DFLT
           "]\n"
           "        --adjust-BAF-LRR <int>      minimum number of genotypes for a cluster to median adjust BAF and LRR "
           "(-1 for no adjustment) [" ADJ_BAF_LRR_DFLT
           "]\n"
           "        --regress-BAF-LRR <int>     minimum number of genotypes for a cluster to regress BAF against LRR "
           "(-1 for no regression) [" REGRESS_BAF_LRR_DFLT
           "]\n"
           "        --LRR-GC-order <int>        order of polynomial to regress LRR against local GC content (-1 for no "
           "regression) [" LRR_GC_ORDER_DFLT
           "]\n"
           "        --xy-major-pl               major transition phred-scaled likelihood [" XY_MAJOR_PL_DFLT
           "]\n"
           "        --xy-minor-pl               minor transition phred-scaled likelihood [" XY_MINOR_PL_DFLT
           "]\n"
           "        --auto-tel-pl               autosomal telomeres phred-scaled likelihood [" AUTO_TEL_PL_DFLT
           "]\n"
           "        --chrX-tel-pl               chromosome X telomeres phred-scaled likelihood [" CHRX_TEL_PL_DFLT
           "]\n"
           "        --chrY-tel-pl               chromosome Y telomeres phred-scaled likelihood [" CHRY_TEL_PL_DFLT
           "]\n"
           "        --error-pl                  uniform error phred-scaled likelihood [" ERR_PL_DFLT
           "]\n"
           "        --flip-pl                   phase flip phred-scaled likelihood [" FLIP_PL_DFLT
           "]\n"
           "        --short-arm-chrs <list>     list of chromosomes with short arms [" SHORT_ARM_CHRS_DFLT
           "]\n"
           "        --use-short-arms            use variants in short arms [FALSE]\n"
           "        --use-centromeres           use variants in centromeres [FALSE]\n"
           "        --use-males-xtr             use variants in XTR region for males [FALSE]\n"
           "        --use-males-par2            use variants in PAR2 region for males [FALSE]\n"
           "        --use-no-rules-chrs         use chromosomes without centromere rules  [FALSE]\n"
           "        --LRR-weight <float>        relative contribution from LRR for LRR+BAF  model [" LRR_BIAS_DFLT
           "]\n"
           "        --LRR-hap2dip <float>       difference between LRR for haploid and diploid [" LRR_HAP2DIP_DFLT
           "]\n"
           "        --LRR-cutoff <float>        cutoff between LRR for haploid and diploid used to infer gender "
           "[estimated from X nonPAR]\n"
           "\n"
           "Examples:\n"
           "    bcftools +mocha -g GRCh37 -v ^exclude.bcf -p cnps.bed -c calls.tsv -z stats.tsv input.bcf\n"
           "    bcftools +mocha -g GRCh38 -o output.bcf -Ob -c calls.tsv -z stats.tsv --LRR-weight 0.5 input.bcf\n"
           "\n";
}

static float *read_list_invf(const char *str, int *n, float min, float max) {
    char *tmp, **list = hts_readlist(str, 0, n);
    if (*n >= 128) error("Cannot handle list of 128 or more parameters: %s\n", str);
    float *ret = (float *)malloc(*n * sizeof(float));
    for (int i = 0; i < *n; i++) {
        ret[i] = 1.0f / strtof(list[i], &tmp);
        if (*tmp) error("Could not parse: %s\n", list[i]);
        if (ret[i] < min || ret[i] > max)
            error("Expected values from the interval [%f,%f], found 1/%s\n", min, max, list[i]);
        free(list[i]);
    }
    free(list);
    ks_introsort_float((size_t)*n, ret);
    return ret;
}

static int cnp_parse(const char *line, char **chr_beg, char **chr_end, uint32_t *beg, uint32_t *end, void *payload,
                     void *usr) {
    // Use the standard parser for CHROM,FROM,TO
    int ret = regidx_parse_bed(line, chr_beg, chr_end, beg, end, NULL, NULL);
    if (ret != 0) return ret;

    // Skip the fields that were parsed above
    char *ss = (char *)line;
    while (*ss && isspace(*ss)) ss++;
    for (int i = 0; i < 3; i++) {
        while (*ss && !isspace(*ss)) ss++;
        if (!*ss) return -2; // wrong number of fields
        while (*ss && isspace(*ss)) ss++;
    }
    if (!*ss) return -2;

    // Parse the payload
    int *dat = (int *)payload;
    if (strncmp(ss, "DEL", 3) == 0)
        *dat = MOCHA_CNP_LOSS;
    else if (strncmp(ss, "DUP", 3) == 0)
        *dat = MOCHA_CNP_GAIN;
    else if (strncmp(ss, "CNV", 3) == 0)
        *dat = MOCHA_CNP_CNV;
    else
        *dat = MOCHA_UNDET;
    return 0;
}

static FILE *get_file_handle(const char *str) {
    FILE *ret;
    if (strcmp(str, "-") == 0)
        ret = stdout;
    else {
        ret = fopen(str, "w");
        if (!ret) error("Failed to open %s: %s\n", str, strerror(errno));
    }
    return ret;
}

int run(int argc, char *argv[]) {
    beta_binom_null = beta_binom_init();
    beta_binom_alt = beta_binom_init();

    // program options
    char *tmp = NULL;
    int rules_is_file = 0;
    int only_stats = 0;
    int regions_is_file = 0;
    int targets_is_file = 0;
    int sample_is_file = 0;
    int force_samples = 0;
    int output_type = FT_VCF;
    int regions_overlap = 1;
    int targets_overlap = 0;
    int clevel = -1;
    int n_threads = 0;
    int record_cmd_line = 1;
    const char *filter_fname = NULL;
    const char *filter_str = NULL;
    const char *regions_list = NULL;
    const char *targets_list = NULL;
    const char *sample_names = NULL;
    const char *output_fname = NULL;
    const char *cnp_fname = NULL;
    const char *stats_fname = NULL;
    const char *mhc_reg = NULL;
    const char *kir_reg = NULL;
    char *rules = NULL;
    sample_t *sample = NULL;
    FILE *log_file = stderr;
    FILE *out_fc = stdout;
    FILE *out_fz = NULL;
    FILE *out_fu = NULL;
    bcf_hdr_t *hdr = NULL;
    bcf_hdr_t *out_hdr = NULL;
    htsFile *out_fh = NULL;
    mocha_table_t mocha_table = {0, 0, NULL};
    const char *short_arm_chrs = SHORT_ARM_CHRS_DFLT;
    const char *bdev_lrr_baf = BDEV_LRR_BAF_DFLT;
    const char *bdev_baf_phase = BDEV_BAF_PHASE_DFLT;

    // model parameters
    model_t model;
    memset(&model, 0, sizeof(model_t));
    model.xy_major_log_prb = -strtof(XY_MAJOR_PL_DFLT, NULL) / 10.0f * M_LN10;
    model.xy_minor_log_prb = -strtof(XY_MINOR_PL_DFLT, NULL) / 10.0f * M_LN10;
    model.auto_tel_log_prb = -strtof(AUTO_TEL_PL_DFLT, NULL) / 10.0f * M_LN10;
    model.chrX_tel_log_prb = -strtof(CHRX_TEL_PL_DFLT, NULL) / 10.0f * M_LN10;
    model.chrY_tel_log_prb = -strtof(CHRY_TEL_PL_DFLT, NULL) / 10.0f * M_LN10;
    model.err_log_prb = -strtof(ERR_PL_DFLT, NULL) / 10.0f * M_LN10;
    model.flip_log_prb = -strtof(FLIP_PL_DFLT, NULL) / 10.0f * M_LN10;
    model.min_dst = (int)strtol(MIN_DST_DFLT, NULL, 0);
    model.lrr_bias = strtof(LRR_BIAS_DFLT, NULL);
    model.lrr_hap2dip = strtof(LRR_HAP2DIP_DFLT, NULL);
    model.lrr_cutoff = NAN;
    model.adj_baf_lrr = (int)strtol(ADJ_BAF_LRR_DFLT, NULL, 0);
    model.regress_baf_lrr = (int)strtol(REGRESS_BAF_LRR_DFLT, NULL, 0);
    model.lrr_gc_order = (int)strtol(LRR_GC_ORDER_DFLT, NULL, 0);

    // create synced reader object
    bcf_srs_t *sr = bcf_sr_init();
    bcf_sr_set_opt(sr, BCF_SR_REQUIRE_IDX);

    static struct option loptions[] = {{"genome", required_argument, NULL, 'g'},
                                       {"genome-file", required_argument, NULL, 'G'},
                                       {"variants", required_argument, NULL, 'v'},
                                       {"apply-filters", required_argument, NULL, 'f'},
                                       {"exclude", required_argument, NULL, 'e'},
                                       {"include", required_argument, NULL, 'i'},
                                       {"regions", required_argument, NULL, 'r'},
                                       {"regions-file", required_argument, NULL, 'R'},
                                       {"regions-overlap", required_argument, NULL, 1},
                                       {"targets", required_argument, NULL, 't'},
                                       {"targets-file", required_argument, NULL, 'T'},
                                       {"targets-overlap", required_argument, NULL, 2},
                                       {"samples", required_argument, NULL, 's'},
                                       {"samples-file", required_argument, NULL, 'S'},
                                       {"force-samples", no_argument, NULL, 3},
                                       {"cnp", required_argument, NULL, 'p'},
                                       {"input-stats", required_argument, NULL, 4},
                                       {"only-stats", no_argument, NULL, 5},
                                       {"mhc", required_argument, NULL, 6},
                                       {"kir", required_argument, NULL, 7},
                                       {"threads", required_argument, NULL, 9},
                                       {"output", required_argument, NULL, 'o'},
                                       {"output-type", required_argument, NULL, 'O'},
                                       {"no-version", no_argument, NULL, 8},
                                       {"no-annotations", no_argument, NULL, 'a'},
                                       {"calls", required_argument, NULL, 'c'},
                                       {"stats", required_argument, NULL, 'z'},
                                       {"ucsc-bed", required_argument, NULL, 'u'},
                                       {"log", required_argument, NULL, 'l'},
                                       {"no-log", no_argument, NULL, 10},
                                       {"bdev-LRR-BAF", required_argument, NULL, 11},
                                       {"bdev-BAF-phase", required_argument, NULL, 12},
                                       {"min-dist", required_argument, NULL, 13},
                                       {"adjust-BAF-LRR", required_argument, NULL, 14},
                                       {"regress-BAF-LRR", required_argument, NULL, 15},
                                       {"LRR-GC-order", required_argument, NULL, 16},
                                       {"xy-major-pl", required_argument, NULL, 17},
                                       {"xy-minor-pl", required_argument, NULL, 18},
                                       {"auto-tel-pl", required_argument, NULL, 19},
                                       {"chrX-tel-pl", required_argument, NULL, 20},
                                       {"chrY-tel-pl", required_argument, NULL, 21},
                                       {"error-pl", required_argument, NULL, 22},
                                       {"flip-pl", required_argument, NULL, 23},
                                       {"short-arm-chrs", required_argument, NULL, 24},
                                       {"use-short-arms", no_argument, NULL, 25},
                                       {"use-centromeres", no_argument, NULL, 26},
                                       {"use-males-xtr", no_argument, NULL, 27},
                                       {"use-males-par2", no_argument, NULL, 28},
                                       {"use-no-rules-chrs", no_argument, NULL, 29},
                                       {"LRR-weight", required_argument, NULL, 30},
                                       {"LRR-hap2dip", required_argument, NULL, 31},
                                       {"LRR-cutoff", required_argument, NULL, 32},
                                       {NULL, 0, NULL, 0}};
    int c;
    while ((c = getopt_long(argc, argv, "h?g:G:v:f:e:i:r:R:t:T:s:S:p:o:O:al:c:z:u:", loptions, NULL)) >= 0) {
        switch (c) {
        case 'g':
            rules = optarg;
            break;
        case 'G':
            rules = optarg;
            rules_is_file = 1;
            break;
        case 'v':
            if (optarg[0] == '^') {
                filter_fname = optarg + 1;
                model.flags |= FLT_EXCLUDE;
            } else {
                filter_fname = optarg;
                model.flags |= FLT_INCLUDE;
            }
            break;
        case 'f':
            sr->apply_filters = optarg;
            break;
        case 'e':
            filter_str = optarg;
            model.filter_logic |= FLT_EXCLUDE;
            break;
        case 'i':
            filter_str = optarg;
            model.filter_logic |= FLT_INCLUDE;
            break;
        case 'r':
            regions_list = optarg;
            break;
        case 'R':
            regions_list = optarg;
            regions_is_file = 1;
            break;
        case 1:
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
        case 2:
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
        case 3:
            force_samples = 1;
            break;
        case 4:
            stats_fname = optarg;
            break;
        case 5:
            only_stats = 1;
            break;
        case 'p':
            cnp_fname = optarg;
            break;
        case 6:
            mhc_reg = optarg;
            break;
        case 7:
            kir_reg = optarg;
            break;
        case 9:
            n_threads = (int)strtol(optarg, &tmp, 0);
            if (*tmp) error("Could not parse: --threads %s\n", optarg);
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
        case 8:
            record_cmd_line = 0;
            break;
        case 'a':
            model.flags |= NO_ANNOT;
            break;
        case 'l':
            log_file = get_file_handle(optarg);
            break;
        case 10:
            model.flags |= NO_LOG;
            break;
        case 'c':
            out_fc = get_file_handle(optarg);
            break;
        case 'z':
            out_fz = get_file_handle(optarg);
            break;
        case 'u':
            out_fu = get_file_handle(optarg);
            break;
        case 11:
            bdev_lrr_baf = optarg;
            break;
        case 12:
            bdev_baf_phase = optarg;
            break;
        case 13:
            model.min_dst = (int)strtol(optarg, &tmp, 0);
            if (*tmp) error("Could not parse: --min-dist %s\n", optarg);
            break;
        case 14:
            model.adj_baf_lrr = (int)strtol(optarg, &tmp, 0);
            if (*tmp) error("Could not parse: --adjust-BAF-LRR %s\n", optarg);
            break;
        case 15:
            model.regress_baf_lrr = (int)strtol(optarg, &tmp, 0);
            if (*tmp) error("Could not parse: --regress-BAF-LRR %s\n", optarg);
            break;
        case 16:
            model.lrr_gc_order = (int)strtol(optarg, &tmp, 0);
            if (*tmp) error("Could not parse: --LRR-GC-order %s\n", optarg);
            break;
        case 17:
            model.xy_major_log_prb = -strtof(optarg, &tmp) / 10.0f * M_LN10;
            if (*tmp) error("Could not parse: --xy-major-pl %s\n", optarg);
            break;
        case 18:
            model.xy_minor_log_prb = -strtof(optarg, &tmp) / 10.0f * M_LN10;
            if (*tmp) error("Could not parse: --xy-minor-pl %s\n", optarg);
            break;
        case 19:
            model.auto_tel_log_prb = -strtof(optarg, &tmp) / 10.0f * M_LN10;
            if (*tmp) error("Could not parse: --auto-tel-pl %s\n", optarg);
            break;
        case 20:
            model.chrX_tel_log_prb = -strtof(optarg, &tmp) / 10.0f * M_LN10;
            if (*tmp) error("Could not parse: --chrX-tel-pl %s\n", optarg);
            break;
        case 21:
            model.chrY_tel_log_prb = -strtof(optarg, &tmp) / 10.0f * M_LN10;
            if (*tmp) error("Could not parse: --chrY-tel-pl %s\n", optarg);
            break;
        case 22:
            model.err_log_prb = -strtof(optarg, &tmp) / 10.0f * M_LN10;
            if (*tmp) error("Could not parse: --error-pl %s\n", optarg);
            break;
        case 23:
            model.flip_log_prb = -strtof(optarg, &tmp) / 10.0f * M_LN10;
            if (*tmp) error("Could not parse: --flip-pl %s\n", optarg);
            break;
        case 24:
            short_arm_chrs = optarg;
            break;
        case 25:
            model.flags |= USE_SHORT_ARMS;
            break;
        case 26:
            model.flags |= USE_CENTROMERES;
            break;
        case 27:
            model.flags |= USE_MALES_XTR;
            break;
        case 28:
            model.flags |= USE_MALES_PAR2;
            break;
        case 29:
            model.flags |= USE_NO_RULES_CHRS;
            break;
        case 30:
            model.lrr_bias = strtof(optarg, &tmp);
            if (*tmp) error("Could not parse: --LRR-weight %s\n", optarg);
            break;
        case 31:
            model.lrr_hap2dip = strtof(optarg, &tmp);
            if (*tmp) error("Could not parse: --LRR-hap2dip %s\n", optarg);
            break;
        case 32:
            model.lrr_cutoff = strtof(optarg, &tmp);
            if (*tmp) error("Could not parse: --LRR-cutoff %s\n", optarg);
            break;
        case 'h':
        case '?':
            error("%s", usage_text());
            break;
        default:
            error("Unknown argument: %s\n", optarg);
        }
    }

    if (!rules) {
        fprintf(log_file, "Genome reference assembly was not specified with --genome or --genome-file\n");
        error("%s", usage_text());
    }
    int len = strlen(rules);
    if (!rules_is_file && (strncmp(rules, "list", 4) == 0 || rules[len - 1] == '?'))
        genome_init_alias(log_file, rules, NULL);

    if (argc - optind > 1) error("Only one input VCF is allowed\n");

    if (model.filter_logic == (FLT_EXCLUDE | FLT_INCLUDE)) error("Only one of --include or --exclude can be given.\n");

    if (only_stats) {
        if (!out_fz) {
            fprintf(log_file, "Option --only-stats requires output option --stats\n");
            error("%s", usage_text());
        }
        if (output_fname || out_fc != stdout || out_fu) {
            fprintf(log_file, "Cannot use --only-stats with output options --output, --calls, or --ucsc-bed\n");
            error("%s", usage_text());
        }
    }

    if (output_fname == NULL && (model.flags & NO_ANNOT)) {
        fprintf(log_file, "Option --no-annotations requires option --output\n");
        error("%s", usage_text());
    }

    if (model.lrr_gc_order > MAX_ORDER) {
        fprintf(log_file, "Polynomial order must not be greater than %d: --LRR-GC-order %d\n", MAX_ORDER,
                model.lrr_gc_order);
        error("%s", usage_text());
    }

    if (model.xy_major_log_prb > model.xy_minor_log_prb) {
        fprintf(log_file,
                "Major transition phred-scaled likelihood should be greater than minor transition phred-scaled "
                "likelihood: --xy-major-pl %.2f --xy-minor-pl %.2f\n",
                -10.0f * M_LOG10E * model.xy_major_log_prb, -10.0f * M_LOG10E * model.xy_minor_log_prb);
        error("%s", usage_text());
    }

    if (model.xy_minor_log_prb > model.auto_tel_log_prb) {
        fprintf(log_file,
                "Minor transition phred-scaled likelihood should be greater than autosomal telomeres phred-scaled "
                "likelihood: --xy-minor-pl %.2f --auto-tel-pl %.2f\n",
                -10.0f * M_LOG10E * model.xy_minor_log_prb, -10.0f * M_LOG10E * model.auto_tel_log_prb);
        error("%s", usage_text());
    }

    if (model.xy_minor_log_prb > model.chrX_tel_log_prb) {
        fprintf(log_file,
                "Minor transition phred-scaled likelihood should be greater than chromosome X telomeres phred-scaled "
                "likelihood: --xy-minor-pl %.2f --chrX-tel-pl %.2f\n",
                -10.0f * M_LOG10E * model.xy_minor_log_prb, -10.0f * M_LOG10E * model.chrX_tel_log_prb);
        error("%s", usage_text());
    }

    if (model.xy_minor_log_prb > model.chrY_tel_log_prb) {
        fprintf(log_file,
                "Minor transition phred-scaled likelihood should be greater than chromosome Y telomeres phred-scaled "
                "likelihood: --xy-minor-pl %.2f --chrY-tel-pl %.2f\n",
                -10.0f * M_LOG10E * model.xy_minor_log_prb, -10.0f * M_LOG10E * model.chrY_tel_log_prb);
        error("%s", usage_text());
    }

    if (stats_fname && !isnan(model.lrr_cutoff)) {
        fprintf(log_file, "Cannot use options --input-stats and --LRR-cutoff together\n");
        error("%s", usage_text());
    }

    // parse parameters defining hidden states
    model.bdev_lrr_baf = read_list_invf(bdev_lrr_baf, &model.bdev_lrr_baf_n, -0.5f, 0.25f);
    model.bdev_baf_phase = read_list_invf(bdev_baf_phase, &model.bdev_baf_phase_n, 0.0f, 0.5f);

    // read list of regions to genotype
    if (cnp_fname) {
        model.cnp_idx = regidx_init(cnp_fname, cnp_parse, NULL, sizeof(int), NULL);
        if (!model.cnp_idx) error("Error: failed to initialize CNP regions: --cnp %s\n", cnp_fname);
        model.cnp_itr = regitr_init(model.cnp_idx);
    }
    if (mhc_reg) model.mhc_idx = regidx_init_string(mhc_reg, regidx_parse_reg, NULL, 0, NULL);
    if (kir_reg) model.kir_idx = regidx_init_string(kir_reg, regidx_parse_reg, NULL, 0, NULL);

    // input VCF
    char *input_fname = NULL;
    if (optind >= argc) {
        if (!isatty(fileno((FILE *)stdin))) input_fname = "-";
    } else
        input_fname = argv[optind];
    if (!input_fname) error("%s", usage_text());

    // read in the regions from the command line
    if (regions_list) {
        bcf_sr_set_opt(sr, BCF_SR_REGIONS_OVERLAP, regions_overlap);
        if (bcf_sr_set_regions(sr, regions_list, regions_is_file) < 0)
            error("Failed to read the regions: %s\n", regions_list);
    }
    if (targets_list) {
        bcf_sr_set_opt(sr, BCF_SR_TARGETS_OVERLAP, targets_overlap);
        if (bcf_sr_set_targets(sr, targets_list, targets_is_file, 0) < 0)
            error("Failed to read the targets: %s\n", targets_list);
    }

    if (bcf_sr_set_threads(sr, n_threads) < 0) error("Failed to create threads\n");
    if (!bcf_sr_add_reader(sr, input_fname)) error("Failed to open %s: %s\n", input_fname, bcf_sr_strerror(sr->errnum));
    if (filter_fname)
        if (!bcf_sr_add_reader(sr, filter_fname))
            error("Failed to open %s: %s\n", filter_fname, bcf_sr_strerror(sr->errnum));

    // check whether the necessary information has been included in the VCF
    hdr = bcf_sr_get_header(sr, 0);
    if (bcf_hdr_nsamples(hdr) == 0) error("Error: input VCF file has no samples\n");
    if (filter_str) model.filter = filter_init(hdr, filter_str);

    model.allele_a_id = bcf_hdr_id2int(hdr, BCF_DT_ID, "ALLELE_A");
    if (!bcf_hdr_idinfo_exists(hdr, BCF_HL_INFO, model.allele_a_id)) model.allele_a_id = -1;
    if (model.allele_a_id >= 0 && bcf_hdr_id2type(hdr, BCF_HL_INFO, model.allele_a_id) != BCF_HT_INT)
        error("Error: input VCF file ALLELE_A info field is not of integer type\n");

    model.allele_b_id = bcf_hdr_id2int(hdr, BCF_DT_ID, "ALLELE_B");
    if (!bcf_hdr_idinfo_exists(hdr, BCF_HL_INFO, model.allele_a_id)) model.allele_b_id = -1;
    if (model.allele_b_id >= 0 && bcf_hdr_id2type(hdr, BCF_HL_INFO, model.allele_b_id) != BCF_HT_INT)
        error("Error: input VCF file ALLELE_B info field is not of float type\n");

    model.gc_id = bcf_hdr_id2int(hdr, BCF_DT_ID, "GC");
    if (!bcf_hdr_idinfo_exists(hdr, BCF_HL_INFO, model.gc_id)) model.gc_id = -1;
    if (model.lrr_gc_order > 0 && model.gc_id < 0)
        error(
            "Error: input VCF has no GC info field: use \"--LRR-GC-order 0/-1\" to disable LRR adjustment through GC "
            "correction\nor use bcftools +mochatools -- -t GC to infer GC content\n");

    model.gt_id = bcf_hdr_id2int(hdr, BCF_DT_ID, "GT");
    if (!bcf_hdr_idinfo_exists(hdr, BCF_HL_FMT, model.gt_id)) error("Error: input VCF file has no GT format field\n");

    model.ad_id = bcf_hdr_id2int(hdr, BCF_DT_ID, "AD");
    if (!bcf_hdr_idinfo_exists(hdr, BCF_HL_FMT, model.ad_id)) model.ad_id = -1;
    if (model.ad_id >= 0 && bcf_hdr_id2type(hdr, BCF_HL_FMT, model.ad_id) != BCF_HT_INT)
        error("Error: input VCF file AD format field is not of integer type\n");
    if (model.ad_id >= 0 && bcf_hdr_id2length(hdr, BCF_HL_FMT, model.ad_id) != BCF_VL_R)
        error("Error: input VCF file AD format field has wrong number of values\n");

    model.baf_id = bcf_hdr_id2int(hdr, BCF_DT_ID, "BAF");
    if (!bcf_hdr_idinfo_exists(hdr, BCF_HL_FMT, model.baf_id)) model.baf_id = -1;
    if (model.baf_id >= 0 && bcf_hdr_id2type(hdr, BCF_HL_FMT, model.baf_id) != BCF_HT_REAL)
        error("Error: input VCF file BAF format field is not of float type\n");

    model.lrr_id = bcf_hdr_id2int(hdr, BCF_DT_ID, "LRR");
    if (!bcf_hdr_idinfo_exists(hdr, BCF_HL_FMT, model.lrr_id)) model.lrr_id = -1;
    if (model.lrr_id >= 0 && bcf_hdr_id2type(hdr, BCF_HL_FMT, model.lrr_id) != BCF_HT_REAL)
        error("Error: input VCF file LRR format field is not of float type\n");

    if (model.ad_id >= 0) {
        model.flags |= WGS_DATA;
        model.lrr_hap2dip = (float)M_LN2;
    } else if (model.baf_id >= 0 && model.lrr_id >= 0) {
        if (model.allele_a_id < 0 || model.allele_b_id < 0)
            error(
                "Error: input VCF file has no ALLELE_A or ALLELE_B info fields\nUse bcftools +mochatools -- -t "
                "ALLELE_A,ALLELE_B to infer ALLELE_A and ALLELE_B\n");
    } else {
        error("Error: input VCF file must contain either the AD format field or the BAF and LRR format fields\n");
    }

    // median adjustment is necessary with sequencing counts data
    if ((model.lrr_gc_order < 0) && (model.flags & WGS_DATA)) model.lrr_gc_order = 0;

    // beginning of plugin run
    fprintf(log_file, "MoChA " MOCHA_VERSION " https://github.com/freeseek/mocha\n");
    fprintf(log_file, "Genome reference: %s\n", rules);
    if (sample_names) fprintf(log_file, "Samples: %s\n", sample_names);
    if (targets_list) fprintf(log_file, "Targets: %s\n", targets_list);
    if (sr->apply_filters) fprintf(log_file, "Filters: %s\n", sr->apply_filters);
    if (model.flags & FLT_EXCLUDE) fprintf(log_file, "Variants to exclude: %s\n", filter_fname);
    if (model.flags & FLT_INCLUDE) fprintf(log_file, "Variants to include: %s\n", filter_fname);
    if (cnp_fname) fprintf(log_file, "Regions to genotype: %s\n", cnp_fname);
    if (stats_fname) fprintf(log_file, "Genome-wide statistics: %s\n", stats_fname);
    if (mhc_reg) fprintf(log_file, "MHC region to exclude from analysis: %s\n", mhc_reg);
    if (kir_reg) fprintf(log_file, "KIR region to exclude from analysis: %s\n", kir_reg);
    fprintf(log_file, "BAF deviations for LRR+BAF model: %s\n", bdev_lrr_baf);
    fprintf(log_file, "BAF deviations for BAF+phase model: %s\n", bdev_baf_phase);
    if (model.flags & WGS_DATA) {
        fprintf(log_file, "Minimum base pair distance between consecutive sites: %d\n", model.min_dst);
    } else {
        fprintf(log_file, "Minimum number of genotypes for a cluster to median adjust BAF and LRR: %d\n",
                model.adj_baf_lrr);
        fprintf(log_file, "Minimum number of genotypes for a cluster to regress BAF against LRR: %d\n",
                model.regress_baf_lrr);
    }
    fprintf(log_file, "Order of polynomial in local GC content to be used to regress LRR against GC: %d\n",
            model.lrr_gc_order);
    fprintf(log_file, "Major transition phred-scaled likelihood: %.2f\n", -10.0f * M_LOG10E * model.xy_major_log_prb);
    fprintf(log_file, "Minor transition phred-scaled likelihood: %.2f\n", -10.0f * M_LOG10E * model.xy_minor_log_prb);
    fprintf(log_file, "Autosomal telomeres phred-scaled likelihood: %.2f\n",
            -10.0f * M_LOG10E * model.auto_tel_log_prb);
    fprintf(log_file, "Chromosome X telomeres phred-scaled likelihood: %.2f\n",
            -10.0f * M_LOG10E * model.chrX_tel_log_prb);
    fprintf(log_file, "Chromosome Y telomeres phred-scaled likelihood: %.2f\n",
            -10.0f * M_LOG10E * model.chrY_tel_log_prb);
    fprintf(log_file, "Uniform error phred-scaled likelihood: %.2f\n", -10.0f * M_LOG10E * model.err_log_prb);
    fprintf(log_file, "Phase flip phred-scaled likelihood: %.2f\n", -10.0f * M_LOG10E * model.flip_log_prb);
    fprintf(log_file, "List of short arms: %s\n", short_arm_chrs);
    fprintf(log_file, "Use variants in short arms: %s\n", model.flags & USE_SHORT_ARMS ? "TRUE" : "FALSE");
    fprintf(log_file, "Use variants in centromeres: %s\n", model.flags & USE_CENTROMERES ? "TRUE" : "FALSE");
    fprintf(log_file, "Use chromosomes without centromere rules: %s\n",
            model.flags & USE_NO_RULES_CHRS ? "TRUE" : "FALSE");
    fprintf(log_file, "Relative contribution from LRR for LRR+BAF model: %.2f\n", model.lrr_bias);
    fprintf(log_file, "Difference between LRR for haploid and diploid: %.2f\n", model.lrr_hap2dip);

    // initialize genome parameters
    if (rules_is_file)
        model.genome_rules = genome_init_file(rules, hdr);
    else
        model.genome_rules = genome_init_alias(log_file, rules, hdr);
    if (!(model.flags & NO_LOG)) fprintf(log_file, "Using genome assembly from %s\n", rules);
    readlist_short_arms(model.genome_rules, short_arm_chrs, hdr);

    // subset VCF file
    if (sample_names) {
        int ret = bcf_hdr_set_samples(hdr, sample_names, sample_is_file);
        if (ret < 0)
            error("Error parsing the sample list\n");
        else if (ret > 0) {
            if (force_samples)
                fprintf(log_file, "Warn: sample #%d not found in the header... skipping\n", ret);
            else
                error(
                    "Error: sample #%d not found in the header. Use \"--force-samples\" to "
                    "ignore this error\n",
                    ret);
        }
        if (bcf_hdr_nsamples(hdr) == 0) error("Error: subsetting has removed all samples\n");
    }

    // output VCF
    if (output_fname) {
        char wmode[8];
        set_wmode(wmode, output_type, output_fname, clevel);
        out_fh = hts_open(output_fname, wmode);
        if (out_fh == NULL) error("[%s] Error: cannot write to \"%s\": %s\n", __func__, output_fname, strerror(errno));
        if (n_threads) hts_set_opt(out_fh, HTS_OPT_THREAD_POOL, sr->p);
        out_hdr = print_hdr(out_fh, hdr, argc, argv, record_cmd_line, model.flags);
    }

    int nsmpl = bcf_hdr_nsamples(hdr);
    if (nsmpl == 0) error("No samples in the VCF?\n");
    if (!(model.flags & NO_LOG)) fprintf(log_file, "Loading %d sample(s) from the VCF file\n", nsmpl);
    if (!(model.flags & WGS_DATA)) {
        if (model.regress_baf_lrr != -1 && (model.adj_baf_lrr == -1 || model.regress_baf_lrr < model.adj_baf_lrr))
            error(
                "Error: median adjustment must be performed whenever regression adjustment is "
                "performed\n--adjust-BAF-LRR %d --regress-BAF-LRR %d\n",
                model.adj_baf_lrr, model.regress_baf_lrr);
        if (nsmpl < 3 * model.adj_baf_lrr)
            fprintf(log_file,
                    "Warning: --adjust-BAF-LRR %d median adjustment will not be effective with "
                    "only %d sample(s)\n",
                    model.adj_baf_lrr, nsmpl);
        if (nsmpl < 3 * model.regress_baf_lrr)
            fprintf(log_file,
                    "Warning: --regress-BAF-LRR %d regression will not be effective with only "
                    "%d sample(s)\n",
                    model.regress_baf_lrr, nsmpl);
    }

    sample = (sample_t *)calloc(nsmpl, sizeof(sample_t));
    for (int i = 0; i < nsmpl; i++) {
        sample[i].idx = i;
        sample[i].x_nonpar_lrr_median = NAN;
        sample[i].y_nonpar_lrr_median = NAN;
        sample[i].mt_lrr_median = NAN;
    }
    if (stats_fname ? read_stats(sample, hdr, stats_fname, model.lrr_gc_order, model.flags) : 1) {
        for (int rid = 0; rid < hdr->n[BCF_DT_CTG]; rid++) {
            model.rid = rid;
            int nret = get_contig(sr, sample, &model);
            if (model.n <= 0) continue;
            if (!(model.flags & NO_LOG))
                fprintf(log_file, "Using %d out of %d read variants from contig %s\n", model.n, nret,
                        bcf_hdr_id2name(hdr, rid));
            if (model.genome_rules->length[rid] < model.locus_arr[model.n - 1].pos)
                model.genome_rules->length[rid] = model.locus_arr[model.n - 1].pos;
            for (int j = 0; j < nsmpl; j++) sample_stats(sample + j, &model);
        }

        sample_summary(sample, nsmpl, &model, stats_fname == NULL);
    }

    int cnt[4] = {0, 0, 0, 0};
    for (int i = 0; i < nsmpl; i++) cnt[sample[i].computed_gender]++;
    if (!(model.flags & NO_LOG)) {
        if (cnt[GENDER_KLINEFELTER] > 0)
            fprintf(log_file,
                    "Estimated %d sample(s) of unknown gender, %d male(s), %d female(s), and %d Klinefelter(s)\n",
                    cnt[GENDER_UNKNOWN], cnt[GENDER_MALE], cnt[GENDER_FEMALE], cnt[GENDER_KLINEFELTER]);
        else
            fprintf(log_file, "Estimated %d sample(s) of unknown gender, %d male(s), and %d female(s)\n",
                    cnt[GENDER_UNKNOWN], cnt[GENDER_MALE], cnt[GENDER_FEMALE]);
    }
    mocha_print_stats(out_fz, sample, nsmpl, model.lrr_gc_order, hdr, model.flags);

    if (!only_stats) {
        for (int rid = 0; rid < hdr->n[BCF_DT_CTG]; rid++) {
            model.rid = rid;
            int nret = get_contig(sr, sample, &model);
            if (model.n <= 0) continue;
            if (!(model.flags & NO_LOG))
                fprintf(log_file, "Using %d out of %d read variants from contig %s\n", model.n, nret,
                        bcf_hdr_id2name(hdr, rid));
            for (int j = 0; j < nsmpl; j++) {
                if (model.cnp_idx)
                    regidx_overlap(model.cnp_idx, bcf_hdr_id2name(hdr, rid), 0, model.genome_rules->length[rid],
                                   model.cnp_itr);
                sample_run(sample + j, &mocha_table, &model);
            }

            if (output_fname) {
                nret = put_contig(sr, sample, &model, out_fh, out_hdr);
                if (!(model.flags & NO_LOG))
                    fprintf(log_file, "Written %d variants for contig %s\n", nret, bcf_hdr_id2name(hdr, rid));
            }
        }

        // estimate LRR at common autosomal losses and gains
        if (!(model.flags & NO_LOG) && model.cnp_itr) {
            int n_cnp_loss = 0, n_cnp_gain = 0, n_rare_gain = 0, n_rare_loss = 0;
            float *cnp_ldev = (float *)malloc(mocha_table.n * sizeof(float));
            float *rare_ldev = (float *)malloc(mocha_table.n * sizeof(float));
            for (int i = 0; i < mocha_table.n; i++) {
                if (mocha_table.a[i].rid == model.genome_rules->x_rid) continue;
                if (mocha_table.a[i].type == MOCHA_CNP_LOSS && mocha_table.a[i].n_hets == 0)
                    cnp_ldev[n_cnp_loss++] = mocha_table.a[i].ldev;
                else if (mocha_table.a[i].type == MOCHA_CNP_GAIN && mocha_table.a[i].n_hets >= 5
                         && mocha_table.a[i].bdev >= 0.1f)
                    cnp_ldev[mocha_table.n - (++n_cnp_gain)] = mocha_table.a[i].ldev;
                else if (mocha_table.a[i].type == MOCHA_LOSS && mocha_table.a[i].n_hets == 0)
                    rare_ldev[n_rare_loss++] = mocha_table.a[i].ldev;
                else if (mocha_table.a[i].type == MOCHA_GAIN && mocha_table.a[i].n_hets >= 5
                         && mocha_table.a[i].bdev >= 0.1f)
                    rare_ldev[mocha_table.n - (++n_rare_gain)] = mocha_table.a[i].ldev;
            }
            float cnp_loss_ldev = get_median(cnp_ldev, n_cnp_loss, NULL);
            float cnp_gain_ldev = get_median(&cnp_ldev[mocha_table.n - n_cnp_gain], n_cnp_gain, NULL);
            float rare_loss_ldev = get_median(rare_ldev, n_rare_loss, NULL);
            float rare_gain_ldev = get_median(&rare_ldev[mocha_table.n - n_rare_gain], n_rare_gain, NULL);
            fprintf(log_file, "Adjusted LRR at CNP losses: %.4f\n", cnp_loss_ldev);
            fprintf(log_file, "Adjusted LRR at CNP gains: %.4f\n", cnp_gain_ldev);
            fprintf(log_file, "Adjusted LRR at rare losses: %.4f\n", rare_loss_ldev);
            fprintf(log_file, "Adjusted LRR at rare gains: %.4f\n", rare_gain_ldev);
            free(cnp_ldev);
            free(rare_ldev);
        }
    }

    // clear sample data
    for (int j = 0; j < nsmpl; j++) {
        free(sample[j].vcf_imap_arr);
        free(sample[j].data_arr[BDEV]);
        free(sample[j].data_arr[LDEV]);
        free(sample[j].phase_arr);
    }

    // clear model data
    free(model.locus_arr);
    free(model.gc_arr);
    free(model.bdev_lrr_baf);
    free(model.bdev_baf_phase);
    genome_destroy(model.genome_rules);

    // free precomputed tables
    ad_to_lrr_baf(NULL, NULL, NULL, NULL, 0);

    // write table with mosaic chromosomal alterations (and UCSC bed track)
    if (!only_stats) {
        mocha_print_calls(out_fc, mocha_table.a, mocha_table.n, hdr, model.flags, rules, model.lrr_hap2dip);
        mocha_print_ucsc(out_fu, mocha_table.a, mocha_table.n, hdr);
    }
    free(mocha_table.a);

    // close output VCF
    if (output_fname) {
        bcf_hdr_destroy(out_hdr);
        hts_close(out_fh);
    }

    if (log_file != stdout && log_file != stderr) fclose(log_file);

    // clean up
    if (model.filter) filter_destroy(model.filter);
    if (model.cnp_idx) regidx_destroy(model.cnp_idx);
    if (model.cnp_itr) regitr_destroy(model.cnp_itr);
    bcf_sr_destroy(sr);
    free(sample);

    beta_binom_destroy(beta_binom_null);
    beta_binom_destroy(beta_binom_alt);
    return 0;
}
