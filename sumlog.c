/* sumlog - fast sum(log(.)), i.e. log(prod(.)) without over- or under-flow
   Copyright (C) 2008 Iain Murray

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.   */

/* Written by Iain Murray <i.murray@ed.ac.uk> */

#include <math.h>

/* This code uses ideas from:
 * http://graphics.stanford.edu/~seander/bithacks.html#IntegerLogIEEE64Float */

#ifdef __LITTLE_ENDIAN__
#define EXP_PART 1
#elifdef __BIG_ENDIAN__
#define EXP_PART 0
/* see http://gcc.gnu.org/onlinedocs/cpp/Common-Predefined-Macros.html */
#elif __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
#define EXP_PART 1
#else
#define EXP_PART 0
#endif

#ifndef min
#define min(x,y) (((x) > (y)) ? (y) : (x))
#endif
#ifndef INFINITY
#define INFINITY HUGE_VAL
#endif
#ifndef NAN
#define NAN (INFINITY - INFINITY)
#endif

int naive_sumlog(double x[], int size, double* answer)
{
    int i;
    double sum = 0;
    for (i = 0; i < size; ++i) {
        sum += log(x[i]);
    }
    *answer = sum;
    return 0;
}

int sumlog(double x[], int size, double* answer)
/* Gets each item in x as m * 2^e where m is between 1 and 2.
 * Multiplies all the m's and adds all the e's to find log(prod(x)) without
 * over- or under-flowing. */
{
    int exponent, i;
    int start = 0;
    int end;
    int type;
    int nan_flag = 0;
    int denormal_nan_inf_neg_or_zero;
    union { unsigned int u[2]; double d; } tmp;
    double prod = 1.0;
    double sum = 0.0;
    const int BATCH_SIZE = 1023; /* (2-eps)^1024 doesn't overflow */

    while (1) {
        end = min(start + BATCH_SIZE, size);
        /* Main Loop */
        for (i = start; i < end; ++i) {
            tmp.d = x[i];
            exponent = (tmp.u[EXP_PART] >> 20) - 1023;
            denormal_nan_inf_neg_or_zero = ((exponent == -1023) || (exponent > 1023));
            if (denormal_nan_inf_neg_or_zero)
                break;
            tmp.u[EXP_PART] &= 0x800FFFFF; /* 0b10000000000011111111111111111111 */
            tmp.u[EXP_PART] |= 0x3FF00000; /* 0b00111111111100000000000000000000 */
            prod *= tmp.d;
            sum += exponent;
        }

        /* Check if loop broke out because of special cases and deal with them.
         * This isn't optimized. For maximum performance, one should ensure that
         * all of the inputs to sumlog() are nice. */
        if (i < end) {
            if (tmp.d < 0) {
                *answer = NAN;
                return -1;
            }
            type = fpclassify(tmp.d);
            switch (type) {
                case FP_ZERO:
                    prod = 0.0;
                    break;
                case FP_INFINITE:
                    sum = INFINITY;
                    break;
                case FP_NAN:
                    nan_flag = 1;
                    break;
                case FP_SUBNORMAL:
                    /* denormalized numbers...YUCK! */
                    tmp.d *= (1 << 30);
                    tmp.d *= (1 << 30);
                    sum -= 60;
                    exponent = (tmp.u[EXP_PART] >> 20) - 1023;
                    tmp.u[EXP_PART] &= 0x800FFFFF; /* 0b10000000000011111111111111111111 */
                    tmp.u[EXP_PART] |= 0x3FF00000; /* 0b00111111111100000000000000000000 */
                    prod *= tmp.d;
                    sum += exponent;
                    break;
                case FP_NORMAL:
                    return -2;
                    break;
                default:
                    return -3;
            }
            /* Correct where we actually got to for next bit of code. */
            end = i + 1;
        }

        /* Stop prod from over-flowing every so often */
        if (end < size) {
            if (prod > 0) {
                tmp.d = prod;
                sum += (tmp.u[EXP_PART] >> 20) - 1023;
                tmp.u[EXP_PART] &= 0x800FFFFF; /* 0b10000000000011111111111111111111 */
                tmp.u[EXP_PART] |= 0x3FF00000; /* 0b00111111111100000000000000000000 */
                prod = tmp.d;
            }
            start = end;
        } else {
            break;
        }
    }
    if (nan_flag)
        *answer = NAN;
    else
        *answer = log(prod) + sum*M_LN2;
    return 0;
}

int sumlogf(float x[], int size, float* answer)
/* Gets each item in x as m * 2^e where m is between 1 and 2.
 * Multiplies all the m's and adds all the e's to find log(prod(x)) without
 * over- or under-flowing. */
{
    int exponent, i;
    int start = 0;
    int end;
    int type;
    int nan_flag = 0;
    int denormal_nan_inf_neg_or_zero;
    union { unsigned int u; float f; } tmp;
    float prod = 1.0;
    float sum = 0.0;
    const int BATCH_SIZE = 127; /* (2-eps)^128 doesn't overflow (for singles) */
    const int FLOAT_EXP_BIAS = 127;

    while (1) {
        end = min(start + BATCH_SIZE, size);
        /* Main Loop */
        for (i = start; i < end; ++i) {
            tmp.f = x[i];
            exponent = (tmp.u >> 23) - FLOAT_EXP_BIAS;
            denormal_nan_inf_neg_or_zero = ((exponent == -FLOAT_EXP_BIAS) || (exponent > FLOAT_EXP_BIAS));
            if (denormal_nan_inf_neg_or_zero)
                break;
            tmp.u &= 0x807FFFFF; /* 0b10000000011111111111111111111111 */
            tmp.u |= 0x3F800000; /* 0b00111111100000000000000000000000 */
            prod *= tmp.f;
            sum += exponent;
        }

        /* Check if loop broke out because of special cases and deal with them.
         * This isn't optimized. For maximum performance, one should ensure that
         * all of the inputs to sumlog() are nice. */
        if (i < end) {
            if (tmp.f < 0) {
                *answer = NAN;
                return -1;
            }
            type = fpclassify(tmp.f);
            switch (type) {
                case FP_ZERO:
                    prod = 0.0;
                    break;
                case FP_INFINITE:
                    sum = INFINITY;
                    break;
                case FP_NAN:
                    nan_flag = 1;
                    break;
                case FP_SUBNORMAL:
                    /* denormalized numbers...YUCK! */
                    tmp.f *= (1 << 30);
                    sum -= 30;
                    exponent = (tmp.u >> 23) - FLOAT_EXP_BIAS;
                    tmp.u &= 0x807FFFFF; /* 0b10000000011111111111111111111111 */
                    tmp.u |= 0x3F800000; /* 0b00111111100000000000000000000000 */
                    prod *= tmp.f;
                    sum += exponent;
                    break;
                case FP_NORMAL:
                    return -2;
                    break;
                default:
                    return -3;
            }
            /* Correct where we actually got to for next bit of code. */
            end = i + 1;
        }

        /* Stop prod from over-flowing every so often */
        if (end < size) {
            if (prod > 0) {
                tmp.f = prod;
                sum += (tmp.u >> 23) - FLOAT_EXP_BIAS;
                tmp.u &= 0x807FFFFF; /* 0b10000000011111111111111111111111 */
                tmp.u |= 0x3F800000; /* 0b00111111100000000000000000000000 */
                prod = tmp.f;
            }
            start = end;
        } else {
            break;
        }
    }
    if (nan_flag)
        *answer = NAN;
    else
        *answer = logf(prod) + sum*M_LN2;
    return 0;
}
