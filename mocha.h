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

/*
    Computes beta-binomial likelihoods using recursively pre-calculated tables
*/

#ifndef __MOCHA_H__
#define __MOCHA_H__

#include <htslib/kseq.h>
#include <htslib/ksort.h>
#include <htslib/vcf.h>
#include "bcftools.h"

// this macro from ksort.h defines the function
// float ks_ksmall_float(size_t n, float arr[], size_t kk);
KSORT_INIT_GENERIC(float)

static inline int *parse_gender(bcf_hdr_t *hdr, char *fname) {
    htsFile *fp = hts_open(fname, "r");
    if (!fp) error("Could not read: %s\n", fname);

    kstring_t str = {0, 0, NULL};
    if (hts_getline(fp, KS_SEP_LINE, &str) <= 0) error("Empty file: %s\n", fname);

    int *gender = (int *)calloc((size_t)bcf_hdr_nsamples(hdr), sizeof(int));

    int moff = 0, *off = NULL;
    char *tmp = NULL;
    do {
        int ncols = ksplit_core(str.s, 0, &moff, &off);
        if (ncols < 2) error("Could not parse the gender file: %s\n", str.s);

        int sample = bcf_hdr_id2int(hdr, BCF_DT_SAMPLE, &str.s[off[0]]);
        if (sample >= 0) {
            if (toupper(str.s[off[1]]) == 'M') {
                gender[sample] = 1;
            } else if (toupper(str.s[off[1]]) == 'F') {
                gender[sample] = 2;
            } else if (toupper(str.s[off[1]]) == 'U') {
                gender[sample] = 0;
            } else {
                gender[sample] = (int)strtol(&str.s[off[1]], &tmp, 0);
                if (*tmp) error("Could not parse gender %s from file %s\n", &str.s[off[1]], fname);
            }
        }
    } while (hts_getline(fp, KS_SEP_LINE, &str) >= 0);

    free(str.s);
    free(off);
    hts_close(fp);
    return gender;
}

// compute the median of a vector using the ksort library (with iterator)
float get_median(const float *v, int n, const int *imap) {
    if (n == 0) return NAN;
    float tmp, *w = (float *)malloc(n * sizeof(float));
    int j = 0;
    for (int i = 0; i < n; i++) {
        tmp = imap ? v[imap[i]] : v[i];
        if (!isnan(tmp)) w[j++] = tmp;
    }
    if (j == 0) {
        free(w);
        return NAN;
    }
    float ret = ks_ksmall_float((size_t)j, w, (size_t)j / 2);
    if (j % 2 == 0) ret = (ret + w[j / 2 - 1]) * 0.5f;
    free(w);
    return ret;
}

// rho =  xyss / sqrtf(xss * yss)
// m = xyss / xss
// b = ym - m * xm;
static inline int get_cov(const float *x, const float *y, int n, const int *imap, float *xss, float *yss, float *xyss) {
    if (n < 2) return -1;
    float xm = 0.0f;
    float ym = 0.0f;
    for (int i = 0; i < n; i++) {
        int idx = imap ? imap[i] : i;
        xm += x[idx];
        ym += y[idx];
    }
    xm /= n;
    ym /= n;
    *xss = 0.0f;
    *yss = 0.0f;
    *xyss = 0.0f;
    for (int i = 0; i < n; i++) {
        int idx = imap ? imap[i] : i;
        float xd = x[idx] - xm;
        float yd = y[idx] - ym;
        *xss += xd * xd;
        *yss += yd * yd;
        *xyss += xd * yd;
    }
    return 0;
}

// retrieve genotype alleles information from BCF record
// assumes little endian architecture
int bcf_get_genotype_alleles(const bcf_fmt_t *fmt, int16_t *gt0_arr, int16_t *gt1_arr, int nsmpl) {
    if (!fmt || fmt->n != 2) return 0;

#define BRANCH(type_t, bcf_type_vector_end)                                                                            \
    {                                                                                                                  \
        type_t *p = (type_t *)fmt->p;                                                                                  \
        for (int i = 0; i < nsmpl; i++, p += 2) {                                                                      \
            if (p[0] == bcf_type_vector_end || bcf_gt_is_missing(p[0]) || p[1] == bcf_type_vector_end                  \
                || bcf_gt_is_missing(p[1])) {                                                                          \
                gt0_arr[i] = bcf_int16_missing;                                                                        \
                gt1_arr[i] = bcf_int16_missing;                                                                        \
            } else {                                                                                                   \
                gt0_arr[i] = (int16_t)bcf_gt_allele(p[0]);                                                             \
                gt1_arr[i] = (int16_t)bcf_gt_allele(p[1]);                                                             \
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

// retrieve phase information from BCF record
// bcf_int8_missing if phase does not apply
// 0 if phase is not available
// 1 if higher number allele received from the mother
// -1 if higher number allele received from the father
// assumes little endian architecture
int bcf_get_genotype_phase(const bcf_fmt_t *fmt, int8_t *gt_phase_arr, int nsmpl) {
    // bcf_fmt_t *fmt = bcf_get_fmt_id(line, id);
    if (!fmt || fmt->n != 2) return 0;

#define BRANCH(type_t, bcf_type_vector_end)                                                                            \
    {                                                                                                                  \
        type_t *p = (type_t *)fmt->p;                                                                                  \
        for (int i = 0; i < nsmpl; i++, p += 2) {                                                                      \
            if (p[0] == bcf_type_vector_end || bcf_gt_is_missing(p[0]) || p[1] == bcf_type_vector_end                  \
                || bcf_gt_is_missing(p[1])) {                                                                          \
                gt_phase_arr[i] = bcf_int8_missing;                                                                    \
            } else {                                                                                                   \
                type_t gt0 = bcf_gt_allele(p[0]) > 0;                                                                  \
                type_t gt1 = bcf_gt_allele(p[1]) > 0;                                                                  \
                if (gt0 == gt1)                                                                                        \
                    gt_phase_arr[i] = bcf_int8_missing;                                                                \
                else if (!bcf_gt_is_phased(p[1]))                                                                      \
                    gt_phase_arr[i] = 0;                                                                               \
                else if (gt1 > gt0)                                                                                    \
                    gt_phase_arr[i] = 1;                                                                               \
                else                                                                                                   \
                    gt_phase_arr[i] = -1;                                                                              \
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

// retrive allelic depth information from BCF record
// assumes little endian architecture
int bcf_get_allelic_depth(const bcf_fmt_t *fmt, const int16_t *gt0_arr, const int16_t *gt1_arr, int16_t *ad0_arr,
                          int16_t *ad1_arr, int nsmpl) {
    if (!fmt) return 0;
    int nalleles = fmt->n;

#define BRANCH(type_t, bcf_type_vector_end, bcf_type_missing)                                                          \
    {                                                                                                                  \
        type_t *p = (type_t *)fmt->p;                                                                                  \
        for (int i = 0; i < nsmpl; i++, p += nalleles) {                                                               \
            if ((gt0_arr[i] != bcf_int16_missing && (int)gt0_arr[i] >= nalleles)                                       \
                || (gt1_arr[i] != bcf_int16_missing && (int)gt1_arr[i] >= nalleles))                                   \
                error(                                                                                                 \
                    "Error: found VCF record with GT alleles %d and %d and %d number of "                              \
                    "alleles\n",                                                                                       \
                    gt0_arr[i], gt1_arr[i], nalleles);                                                                 \
            if (gt0_arr[i] == bcf_int16_missing || gt1_arr[i] == bcf_int16_missing                                     \
                || p[gt0_arr[i]] == bcf_type_vector_end || p[gt0_arr[i]] == bcf_type_missing                           \
                || p[gt1_arr[i]] == bcf_type_vector_end || p[gt1_arr[i]] == bcf_type_missing) {                        \
                ad0_arr[i] = bcf_int16_missing;                                                                        \
                ad1_arr[i] = bcf_int16_missing;                                                                        \
            } else {                                                                                                   \
                type_t gt0 = gt0_arr[i] > 0;                                                                           \
                type_t gt1 = gt1_arr[i] > 0;                                                                           \
                if (gt0 == gt1) {                                                                                      \
                    ad0_arr[i] = (int16_t)(gt0_arr[i] == gt1_arr[i] ? p[gt0_arr[i]] : p[gt0_arr[i]] + p[gt1_arr[i]]);  \
                    ad1_arr[i] = bcf_int16_missing;                                                                    \
                } else {                                                                                               \
                    ad0_arr[i] = (int16_t)p[gt0_arr[i]];                                                               \
                    ad1_arr[i] = (int16_t)p[gt1_arr[i]];                                                               \
                }                                                                                                      \
            }                                                                                                          \
        }                                                                                                              \
    }
    switch (fmt->type) {
    case BCF_BT_INT8:
        BRANCH(int8_t, bcf_int8_vector_end, bcf_int8_missing);
        break;
    case BCF_BT_INT16:
        BRANCH(int16_t, bcf_int16_vector_end, bcf_int16_missing);
        break;
    case BCF_BT_INT32:
        BRANCH(int32_t, bcf_int32_vector_end, bcf_int32_missing);
        break;
    default:
        error("Unexpected type %d\n", fmt->type);
    }
#undef BRANCH

    return 1;
}

#endif
