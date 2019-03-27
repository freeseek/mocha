/* The MIT License

   Copyright (C) 2018-2019 Giulio Genovese

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

#ifndef __BETA_BINOM_H__
#define __BETA_BINOM_H__

#include <math.h>

typedef struct _beta_binom_t beta_binom_t;

struct _beta_binom_t
{
    double p;
    double rho;
    int n1;
    int n2;
    double *log_gamma_alpha;
    double *log_gamma_beta;
    double *log_gamma_alpha_beta;
    int m_log_gamma_alpha;
    int m_log_gamma_beta;
    int m_log_gamma_alpha_beta;
};

beta_binom_t *beta_binom_init();
void beta_binom_destroy(beta_binom_t *self);

/**
 *  beta_binom_update() - updates beta binomial parameters and pre-calculated tables
 *  @p: probability of success
 *  @rho: "intra class" or "intra cluster" correlation
 *  @n1: maximum value for a, b that can be passed to beta_binom_log()
 *  @n2: maximum value for a+b that can be passed to beta_binom_log()
 */
void beta_binom_update(beta_binom_t *self, double p, double rho, int n1, int n2);

/**
 *  beta_binom_log_unsafe() - density function of the beta binomial distribution
 *  Returns the equivalent of dbeta_binom(a, a+b, p, (1 - rho) / rho, log=TRUE) from R package rmutil
 */
static inline double beta_binom_log_unsafe(const beta_binom_t *self, int a, int b)
{
    return self->log_gamma_alpha[a] + self->log_gamma_beta[b] - self->log_gamma_alpha_beta[a + b];
}

/**
 *  Same as before but it performs boundary checking before computing the log likelihood
 */
static inline double beta_binom_log(beta_binom_t *self, int a, int b)
{
    if ( a < 0 || b < 0 ) return NAN;
    if ( a > self->n1 || b > self->n1 || a+b > self->n2 )
        beta_binom_update(self, self->p, self->rho, a > b ? a : b, a + b);
    return beta_binom_log_unsafe(self, a, b);
}

#endif
