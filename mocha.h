/* The MIT License

   Copyright (C) 2018-2022 Giulio Genovese

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
#include "tsv2vcf.h"
#include "bcftools.h"

#define GENDER_UNKNOWN 0
#define GENDER_MALE 1
#define GENDER_FEMALE 2
#define GENDER_KLINEFELTER 3

/****************************************
 * TSV FUNCTIONS                        *
 ****************************************/

// adapted from Petr Danecek's implementation of tsv_init() in bcftools/tsv2vcf.c
static inline tsv_t *tsv_init_delimiter(const char *str, char delimiter) {
    tsv_t *tsv = (tsv_t *)calloc(1, sizeof(tsv_t));
    kstring_t tmp = {0, 0, 0};
    const char *ss = str, *se = ss;
    tsv->ncols = 0;
    while (*ss) {
        if (delimiter == '\0')
            while (*se && !isspace(*se)) se++;
        else
            while (*se && *se != delimiter) se++;
        tsv->ncols++;
        tsv->cols = (tsv_col_t *)realloc(tsv->cols, sizeof(tsv_col_t) * tsv->ncols);
        tsv->cols[tsv->ncols - 1].name = NULL;
        tsv->cols[tsv->ncols - 1].setter = NULL;
        tmp.l = 0;
        kputsn(ss, se - ss, &tmp);
        if (strcasecmp("-", tmp.s)) tsv->cols[tsv->ncols - 1].name = strdup(tmp.s);
        if (!*se) break;
        se++;
        if (delimiter == '\0')
            while (*se && isspace(*se)) se++;
        ss = se;
    }
    free(tmp.s);
    return tsv;
}

static inline int tsv_parse_delimiter(tsv_t *tsv, bcf1_t *rec, char *str, char delimiter) {
    int status = 0;
    tsv->icol = 0;
    tsv->ss = tsv->se = str;
    while (*tsv->ss && tsv->icol < tsv->ncols) {
        if (delimiter == '\0')
            while (*tsv->se && !isspace(*tsv->se)) tsv->se++;
        else
            while (*tsv->se && *tsv->se != delimiter) tsv->se++;
        if (tsv->cols[tsv->icol].setter) {
            int ret = tsv->cols[tsv->icol].setter(tsv, rec, tsv->cols[tsv->icol].usr);
            if (ret < 0) return -1;
            status++;
        }
        if (*tsv->se) {
            tsv->se++;
            if (delimiter == '\0')
                while (*tsv->se && isspace(*tsv->se)) tsv->se++;
        }
        tsv->ss = tsv->se;
        tsv->icol++;
    }
    return status ? 0 : -1;
}

int tsv_read_float(tsv_t *tsv, bcf1_t *rec, void *usr) {
    float *single = (float *)usr;
    char tmp = *tsv->se;
    *tsv->se = 0;
    char *endptr;
    *single = (float)strtof(tsv->ss, &endptr);
    *tsv->se = tmp;
    return 0;
}

int tsv_read_integer(tsv_t *tsv, bcf1_t *rec, void *usr) {
    int *integer = (int *)usr;
    char tmp = *tsv->se;
    *tsv->se = 0;
    char *endptr;
    *integer = (int)strtol(tsv->ss, &endptr, 0);
    if (*endptr) error("Could not parse integer %s\n", tsv->ss);
    *tsv->se = tmp;
    return 0;
}

int tsv_read_string(tsv_t *tsv, bcf1_t *rec, void *usr) {
    char **str = (char **)usr;
    if (tsv->se == tsv->ss) {
        *str = NULL;
    } else {
        char tmp = *tsv->se;
        *tsv->se = 0;
        *str = strdup(tsv->ss);
        *tsv->se = tmp;
    }
    return 0;
}

int tsv_read_sample_id(tsv_t *tsv, bcf1_t *rec, void *usr) {
    bcf_hdr_t *hdr = (bcf_hdr_t *)rec;
    int *idx = (int *)usr;
    char tmp = *tsv->se;
    *tsv->se = 0;
    *idx = bcf_hdr_id2int(hdr, BCF_DT_SAMPLE, tsv->ss);
    *tsv->se = tmp;
    return 0;
}

int tsv_read_computed_gender(tsv_t *tsv, bcf1_t *rec, void *usr) {
    int *computed_gender = (int *)usr;
    if (toupper(*tsv->ss) == 'M') {
        *computed_gender = GENDER_MALE;
    } else if (toupper(*tsv->ss) == 'F') {
        *computed_gender = GENDER_FEMALE;
    } else if (toupper(*tsv->ss) == 'K') {
        *computed_gender = GENDER_KLINEFELTER;
    } else if (toupper(*tsv->ss) == 'U' || toupper(*tsv->ss) == 'N') {
        *computed_gender = GENDER_UNKNOWN;
    } else {
        char *endptr;
        *computed_gender = (int)strtol(tsv->ss, &endptr, 0);
        if (*endptr)
            error("Could not parse gender %s\n(Acceptable values: 1/M/m = male, 2/F/f = female, 0/U/u/N/n = missing)\n",
                  tsv->ss);
    }
    return 0;
}

/****************************************
 * BASIC STATISTICS FUNCTIONS           *
 ****************************************/

// this macro from ksort.h defines the function
// float ks_ksmall_float(size_t n, float arr[], size_t kk);
KSORT_INIT_GENERIC(float)

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

/****************************************
 * EXTRACT DATA FROM VCF FUNCTIONS      *
 ****************************************/

// retrieve unphased genotype alleles information from BCF record
// assumes little endian architecture
static inline int bcf_get_unphased_genotype_alleles(const bcf_fmt_t *fmt, int16_t *gt0_arr, int16_t *gt1_arr,
                                                    int nsmpl) {
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
                int swap = bcf_gt_allele(p[0]) > bcf_gt_allele(p[1]);                                                  \
                gt0_arr[i] = (int16_t)bcf_gt_allele(p[swap]);                                                          \
                gt1_arr[i] = (int16_t)bcf_gt_allele(p[!swap]);                                                         \
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
// bcf_int8_missing if genotype is missing
// bcf_int8_vector_end if phase does not apply
// 0 if phase is not available
// 1 if higher number allele received from the mother
// -1 if higher number allele received from the father
// assumes little endian architecture
static inline int bcf_get_genotype_phase(const bcf_fmt_t *fmt, int8_t *gt_phase_arr, int nsmpl) {
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
                type_t gt0 = bcf_gt_allele(p[0]);                                                                      \
                type_t gt1 = bcf_gt_allele(p[1]);                                                                      \
                if (gt0 == gt1)                                                                                        \
                    gt_phase_arr[i] = bcf_int8_vector_end;                                                             \
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
static inline int bcf_get_allelic_depth(const bcf_fmt_t *fmt, const int16_t *gt0_arr, const int16_t *gt1_arr,
                                        int16_t *ad0_arr, int16_t *ad1_arr, int nsmpl) {
    if (!fmt) return 0;
    int nalleles = fmt->n;

#define BRANCH(type_t, bcf_type_vector_end, bcf_type_missing)                                                          \
    {                                                                                                                  \
        type_t *p = (type_t *)fmt->p;                                                                                  \
        for (int i = 0; i < nsmpl; i++, p += nalleles) {                                                               \
            if ((gt0_arr[i] != bcf_int16_missing && (int)gt0_arr[i] >= nalleles)                                       \
                || (gt1_arr[i] != bcf_int16_missing && (int)gt1_arr[i] >= nalleles))                                   \
                error("Error: found VCF record with GT alleles %d and %d and %d number of alleles\n", gt0_arr[i],      \
                      gt1_arr[i], nalleles);                                                                           \
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
