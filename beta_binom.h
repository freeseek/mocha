/* The MIT License

   Copyright (C) 2018-2021 Giulio Genovese

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

#include <stdlib.h>
#include <math.h>
#include <htslib/hts.h>

typedef struct _beta_binom_t beta_binom_t;

struct _beta_binom_t {
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

beta_binom_t *beta_binom_init() {
    beta_binom_t *self = (beta_binom_t *)calloc(1, sizeof(beta_binom_t));
    hts_expand0(double, 1, self->m_log_gamma_alpha, self->log_gamma_alpha);
    hts_expand0(double, 1, self->m_log_gamma_beta, self->log_gamma_beta);
    hts_expand0(double, 1, self->m_log_gamma_alpha_beta, self->log_gamma_alpha_beta);
    self->p = NAN;
    self->rho = NAN;
    return self;
}

void beta_binom_destroy(beta_binom_t *self) {
    free(self->log_gamma_alpha);
    free(self->log_gamma_beta);
    free(self->log_gamma_alpha_beta);
    free(self);
}

/**
 *  beta_binom_update() - updates beta binomial parameters and pre-calculated tables
 *  @p: probability of success
 *  @rho: "intra class" or "intra cluster" correlation
 *  @n1: maximum value for a, b that can be passed to beta_binom_log()
 *  @n2: maximum value for a+b that can be passed to beta_binom_log()
 */

// f(n, x) = log( \gamma(n+x) / \gamma(x) / n! )
// log_gamma_alpha[n] = f(n, \alpha)
// log_gamma_beta[n] = f(n, \beta)
// log_gamma_alpha_beta[n] = f(n, \alpha + \beta)
// see https://en.wikipedia.org/wiki/Beta-binomial_distribution#As_a_compound_distribution
// p is the probability of success
// rho is the "intra class" or "intra cluster" correlation
// in Artieri et al. 2017, overdispersion is exactly "(1 - rho) / rho"
void beta_binom_update(beta_binom_t *self, double p, double rho, int n1, int n2) {
    if (self->p != p || self->rho != rho) {
        self->p = p;
        self->rho = rho;
        self->n1 = 0;
        self->n2 = 0;
    }

    hts_expand(double, n1 + 1, self->m_log_gamma_alpha, self->log_gamma_alpha);
    hts_expand(double, n1 + 1, self->m_log_gamma_beta, self->log_gamma_beta);
    hts_expand(double, n2 + 1, self->m_log_gamma_alpha_beta, self->log_gamma_alpha_beta);

    if (rho == 0) // binomial distribution case (no overdispersion)
    {
        double log_alpha = log(p);
        double log_beta = log(1.0 - p);

        while (self->n1 < n1) {
            self->n1++;
            double log_n1 = log(self->n1);
            self->log_gamma_alpha[self->n1] = self->log_gamma_alpha[self->n1 - 1] + log_alpha - log_n1;
            self->log_gamma_beta[self->n1] = self->log_gamma_beta[self->n1 - 1] + log_beta - log_n1;
        }

        while (self->n2 < n2) {
            self->n2++;
            self->log_gamma_alpha_beta[self->n2] = self->log_gamma_alpha_beta[self->n2 - 1] - log(self->n2);
        }
    } else {
        double s = (1.0 - rho) / rho;
        double alpha = p * s;
        double beta = (1.0 - p) * s;

        while (self->n1 < n1) {
            self->n1++;
            self->log_gamma_alpha[self->n1] =
                self->log_gamma_alpha[self->n1 - 1] + log((alpha + (double)self->n1 - 1.0) / (double)self->n1);
            self->log_gamma_beta[self->n1] =
                self->log_gamma_beta[self->n1 - 1] + log((beta + (double)self->n1 - 1.0) / (double)self->n1);
        }

        while (self->n2 < n2) {
            self->n2++;
            self->log_gamma_alpha_beta[self->n2] = self->log_gamma_alpha_beta[self->n2 - 1]
                                                   + log((alpha + beta + (double)self->n2 - 1.0) / (double)self->n2);
        }
    }
}

/**
 *  beta_binom_log_unsafe() - density function of the beta binomial distribution
 *  Returns the equivalent of dbeta_binom(a, a+b, p, (1 - rho) / rho, log=TRUE) from R package
 * rmutil
 */
inline double beta_binom_log_unsafe(const beta_binom_t *self, int a, int b) {
    return self->log_gamma_alpha[a] + self->log_gamma_beta[b] - self->log_gamma_alpha_beta[a + b];
}

/**
 *  Same as before but it performs boundary checking before computing the log likelihood
 */
inline double beta_binom_log(beta_binom_t *self, int a, int b) {
    if (a < 0 || b < 0) return NAN;
    if (a > self->n1 || b > self->n1 || a + b > self->n2)
        beta_binom_update(self, self->p, self->rho, a > b ? a : b, a + b);
    return beta_binom_log_unsafe(self, a, b);
}

#endif
