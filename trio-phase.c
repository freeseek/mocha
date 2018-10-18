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

/*
    Known issues:
    - It only propagates current relative phase information forward, not backward
    An approach like duoHMM would be preferable (http://doi.org/10.1371/journal.pgen.1004234)
*/

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <htslib/vcf.h>
#include <htslib/kseq.h>
#include "bcftools.h"

#define IBD_A 0
#define IBD_B 1

typedef struct
{
    int *a;
    int n;
    int m;
}
set_t;

typedef struct
{
    set_t fathers;
    set_t mothers;
    int32_t *gt;
    int is_child;
    int flip_vote;
    int is_prev_het_flipped;
    int ibd_father;
    int ibd_mother;
}
ind_t;

typedef struct
{

    bcf_hdr_t *in_hdr, *out_hdr;
    ind_t *inds;
    set_t children;
    int prev_rid, ibd;
    int32_t *arr;
    int m_arr;
}
args_t;

args_t *args;

static int add_id(set_t *set, int id)
{
    for (int i=0; i<set->n; i++) if (id == set->a[i]) return -1;
    set->n++;
    hts_expand(int, set->n, set->m, set->a);
    set->a[set->n-1] = id;
    return 0;
}

static int inline gt_is_missing(int32_t *gt) { return gt[0]==bcf_gt_missing || gt[1]==bcf_gt_missing || gt[0]==bcf_int32_vector_end; }
static int inline gt_is_het(int32_t *gt) { return bcf_gt_allele(gt[0])!=bcf_gt_allele(gt[1]); }
static int inline gt_is_phased(int32_t *gt) { return gt_is_het(gt) ? (bcf_gt_is_phased(gt[1]) ? 1 : 0) : 1; }
static int32_t inline gt_pat(int32_t *gt) { return bcf_gt_allele(gt[0]); }
static int32_t inline gt_mat(int32_t *gt) { return bcf_gt_allele(gt[1]); }

// flips genotypes and set genotype as phased
static void inline flip_gt(int32_t *gt)
{
    *gt ^= *(gt+1);
    *(gt+1) ^= *gt;
    *gt ^= *(gt+1);
    *gt &= 0xFFFFFFFE;
    *(gt+1) |= 1;
    return;
}

// propagate information across individuals
static void ind_vote(ind_t *ind)
{
    if ( !gt_is_het(ind->gt) ) return;
    if ( ind->flip_vote > 0 )
    {
        if ( gt_is_phased(ind->gt) ) ind->is_prev_het_flipped = 1;
        flip_gt(ind->gt);
    }
    else if ( ind->flip_vote < 0 )
    {
        if ( gt_is_phased(ind->gt) ) ind->is_prev_het_flipped = 0;
        *(ind->gt+1) |= 1;
    }
    else
    {
        if( gt_is_phased(ind->gt) && ind->is_prev_het_flipped ) flip_gt(ind->gt);
    }
}

static void ind_destroy(ind_t *ind)
{
    free(ind->fathers.a);
    free(ind->mothers.a);
}

const char *about(void)
{
    return "Phase genotypes in trios.\n";
}

const char *usage(void)
{
    return
"\n"
"About: Phase genotypes in trios. (version 2018-10-18)\n"
"\n"
"Usage: bcftools +trio-phase [General Options] -- [Plugin Options]\n"
"Options:\n"
"   run \"bcftools plugin\" for a list of common options\n"
"\n"
"Plugin options:\n"
"   -p, --ped <file>            PED file\n"
"   -i, --ibd                   Whether to add IBD state for duos\n"
"\n"
"Example:\n"
"   bcftools +trio-phase file.bcf -- --ped file.ped\n"
"\n";
}

static void parse_ped(args_t *args, char *fname)
{
    htsFile *fp = hts_open(fname, "r");
    if ( !fp ) error("Could not read: %s\n", fname);

    kstring_t str = {0,0,0};
    if ( hts_getline(fp, KS_SEP_LINE, &str) <= 0 ) error("Empty file: %s\n", fname);

    int moff = 0, *off = NULL;
    do
    {
        int ncols = ksplit_core(str.s,0,&moff,&off);
        if ( ncols<4 ) error("Could not parse the ped file: %s\n", str.s);

        int father_id = bcf_hdr_id2int(args->in_hdr,BCF_DT_SAMPLE,&str.s[off[2]]);
        int mother_id = bcf_hdr_id2int(args->in_hdr,BCF_DT_SAMPLE,&str.s[off[3]]);
        int child_id = bcf_hdr_id2int(args->in_hdr,BCF_DT_SAMPLE,&str.s[off[1]]);
        if ( ( father_id < 0 && mother_id < 0 ) || child_id < 0 ) continue;
        add_id(&args->children, child_id);
        ind_t *child_ind = &args->inds[child_id];
        child_ind->is_child = 1;
        if ( father_id >= 0 ) add_id(&child_ind->fathers, father_id);
        if ( mother_id >= 0 ) add_id(&child_ind->mothers, mother_id);
    } while ( hts_getline(fp, KS_SEP_LINE, &str)>=0 );

    free(str.s);
    free(off);
    hts_close(fp);
}

int init(int argc,
         char **argv,
         bcf_hdr_t *in,
         bcf_hdr_t *out)
{
    args = (args_t *)calloc(1, sizeof(args_t));
    args->prev_rid = -1;
    args->in_hdr = in;
    args->out_hdr = out;
    args->inds = (ind_t *)calloc(bcf_hdr_nsamples(in), sizeof(ind_t));
    char *ped_fname = NULL;
    static struct option loptions[] =
    {
        {"ped", required_argument, NULL, 'p'},
        {"ibd", no_argument, NULL, 'i'},
        {0,0,0,0}
    };
    int c;
    while ((c = getopt_long(argc, argv, "?hp:i",loptions,NULL)) >= 0)
    {
        switch (c)
        {
            case 'p': ped_fname = optarg; break;
            case 'i': args->ibd = 1; break;
            case 'h':
            case '?':
            default: error("%s", usage()); break;
        }
    }
    if ( !ped_fname ) error("Expected the -p option\n");
    parse_ped(args, ped_fname);
    if ( args->ibd )
    {
        bcf_hdr_append(args->out_hdr, "##FORMAT=<ID=IBD_F,Number=1,Type=Integer,Description=\"IBD state with the father\">");
        bcf_hdr_append(args->out_hdr, "##FORMAT=<ID=IBD_M,Number=1,Type=Integer,Description=\"IBD state with the mother\">");
    }
    return 0;
}

bcf1_t *process(bcf1_t *rec)
{
    // extract genotypes from record and checks whether they follow a diploid model
    int nsmpl = bcf_hdr_nsamples(args->in_hdr);
    int ngt = bcf_get_genotypes(args->in_hdr, rec, &args->arr, &args->m_arr);
    if ( ngt<0 ) return rec;
    ngt /= nsmpl;
    if ( ngt!=2 ) return rec;
    for (int i=0; i<nsmpl; i++) args->inds[i].gt = &args->arr[ngt * i];

    // check whether the chromosome in the record has changed
    if ( rec->rid!=args->prev_rid )
    {
        args->prev_rid = rec->rid;
        for (int i=0; i<nsmpl; i++)
        {
            args->inds[i].is_prev_het_flipped = 0;
        }
    }

    // propagate information from parents to children
    for (int i=0; i<args->children.n; i++)
    {
        int child_id = args->children.a[i];
        ind_t *child_ind = &args->inds[child_id];
        if ( gt_is_missing(child_ind->gt) || !gt_is_het(child_ind->gt) ) continue;

        int vote = 0;
        for (int j=0; j<child_ind->fathers.n; j++)
        {
            int father_id = child_ind->fathers.a[j];
            int32_t *father_gt = args->arr + ngt * father_id;
            if ( gt_pat(child_ind->gt) == gt_pat(father_gt) && gt_pat(child_ind->gt) == gt_mat(father_gt) ) vote--;
            else if ( gt_mat(child_ind->gt) == gt_pat(father_gt) && gt_mat(child_ind->gt) == gt_mat(father_gt) ) vote++;
        }
        if ( vote > 0 ) child_ind->flip_vote++;
        else if ( vote < 0 ) child_ind->flip_vote--;

        vote = 0;
        for (int j=0; j<child_ind->mothers.n; j++)
        {
            int mother_id = child_ind->mothers.a[j];
            int32_t *mother_gt = args->arr + ngt * mother_id;
            if ( gt_pat(child_ind->gt) == gt_pat(mother_gt) && gt_pat(child_ind->gt) == gt_mat(mother_gt) ) vote++;
            else if ( gt_mat(child_ind->gt) == gt_pat(mother_gt) && gt_mat(child_ind->gt) == gt_mat(mother_gt) ) vote--;
        }
        if ( vote > 0 ) child_ind->flip_vote++;
        else if ( vote < 0 ) child_ind->flip_vote--;

        ind_vote(child_ind);
    }

    // propagate information from children to parents
    for (int i=0; i<args->children.n; i++)
    {
        int child_id = args->children.a[i];
        ind_t *child_ind = &args->inds[child_id];
        if ( gt_is_missing(child_ind->gt) || ( gt_is_het(child_ind->gt) && !gt_is_phased(child_ind->gt) ) ) continue;

        for (int j=0; j<child_ind->fathers.n; j++)
        {
            int father_id = child_ind->fathers.a[j];
            int32_t *father_gt = args->arr + ngt * father_id;
            ind_t *father_ind = &args->inds[father_id];
            if ( gt_pat(child_ind->gt) == gt_pat(father_gt) && child_ind->ibd_father == IBD_A ) father_ind->flip_vote--;
            if ( gt_pat(child_ind->gt) == gt_pat(father_gt) && child_ind->ibd_father == IBD_B ) father_ind->flip_vote++;
            if ( gt_pat(child_ind->gt) == gt_mat(father_gt) && child_ind->ibd_father == IBD_A ) father_ind->flip_vote++;
            if ( gt_pat(child_ind->gt) == gt_mat(father_gt) && child_ind->ibd_father == IBD_B ) father_ind->flip_vote--;
        }

        for (int j=0; j<child_ind->mothers.n; j++)
        {
            int mother_id = child_ind->mothers.a[j];
            int32_t *mother_gt = args->arr + ngt * mother_id;
            ind_t *mother_ind = &args->inds[mother_id];
            if ( gt_mat(child_ind->gt) == gt_pat(mother_gt) && child_ind->ibd_mother == IBD_A ) mother_ind->flip_vote--;
            if ( gt_mat(child_ind->gt) == gt_pat(mother_gt) && child_ind->ibd_mother == IBD_B ) mother_ind->flip_vote++;
            if ( gt_mat(child_ind->gt) == gt_mat(mother_gt) && child_ind->ibd_mother == IBD_A ) mother_ind->flip_vote++;
            if ( gt_mat(child_ind->gt) == gt_mat(mother_gt) && child_ind->ibd_mother == IBD_B ) mother_ind->flip_vote--;
        }
    }

    // propagate information from children to parents or across sites
    for (int i=0; i<nsmpl; i++)
    {
        ind_t *parent_ind = &args->inds[i];
        // make sure children don't get double voted
        if ( !parent_ind->is_child ) ind_vote(parent_ind);
        parent_ind->flip_vote = 0;
    }

    // update IBD states in children
    for (int i=0; i<args->children.n; i++)
    {
        int child_id = args->children.a[i];
        ind_t *child_ind = &args->inds[child_id];

        if ( !gt_is_missing(child_ind->gt) && gt_is_phased(child_ind->gt) )
        {
            int vote = 0;
            for (int j=0; j<child_ind->fathers.n; j++)
            {
                int father_id = child_ind->fathers.a[j];
                int32_t *father_gt = args->arr + ngt * father_id;
                if ( !gt_is_missing(father_gt) && gt_is_het(father_gt) && gt_is_phased(father_gt) )
                {
                    if ( gt_pat(child_ind->gt) == gt_pat(father_gt) ) vote++;
                    else if ( gt_pat(child_ind->gt) == gt_mat(father_gt) ) vote--;
                }
            }
            if ( vote > 0 ) child_ind->ibd_father = IBD_A;
            else if ( vote < 0 ) child_ind->ibd_father = IBD_B;

            vote = 0;
            for (int j=0; j<child_ind->mothers.n; j++)
            {
                int mother_id = child_ind->mothers.a[j];
                int32_t *mother_gt = args->arr + ngt * mother_id;
                if ( !gt_is_missing(mother_gt) && gt_is_het(mother_gt) && gt_is_phased(mother_gt) )
                {
                    if ( gt_mat(child_ind->gt) == gt_pat(mother_gt) ) vote++;
                    else if ( gt_mat(child_ind->gt) == gt_mat(mother_gt) ) vote--;
                }
            }
            if ( vote > 0 ) child_ind->ibd_mother = IBD_A;
            else if ( vote < 0 ) child_ind->ibd_mother = IBD_B;
        }
    }

    bcf_update_genotypes(args->out_hdr, rec, args->arr, args->m_arr);
    if ( args->ibd )
    {
        for (int i=0; i<nsmpl; i++) args->arr[i] = args->inds[i].ibd_father;
        bcf_update_format_int32(args->out_hdr, rec, "IBD_F", args->arr, nsmpl);
        for (int i=0; i<nsmpl; i++) args->arr[i] = args->inds[i].ibd_mother;
        bcf_update_format_int32(args->out_hdr, rec, "IBD_M", args->arr, nsmpl);
    }
    return rec;
}

void destroy(void)
{
    for (int i=0; i<bcf_hdr_nsamples(args->in_hdr); i++) ind_destroy(&args->inds[i]);
    free(args->inds);
    free(args->children.a);
    free(args->arr);
    free(args);
}
