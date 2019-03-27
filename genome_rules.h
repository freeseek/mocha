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
    Structure to store genome information about centromeres and sex chromosomes
*/

#ifndef __GENOME_RULES_H__
#define __GENOME_RULES_H__

#include <htslib/vcf.h>

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
} genome_rules_t;

genome_rules_t *genome_init(const bcf_hdr_t *hdr);
void genome_destroy(genome_rules_t *self);

/**
 *  genome_init_file() - initialize a genome structure from a file
 */
genome_rules_t *genome_init_file(const char *fname, const bcf_hdr_t *hdr);

/**
 *  genome_init_alias() - initialize a genome structure from an alias such as NCBI36, GRCh37, or GRCh38
 */
genome_rules_t *genome_init_alias(FILE *restrict stream, char *alias, const bcf_hdr_t *hdr);

/**
 *  readlist_short_arms() - initialize flag indicating which chromosome have short arms from a comma separated list
 */
int readlist_short_arms(genome_rules_t *self, const char *str, const bcf_hdr_t *hdr);

#endif
