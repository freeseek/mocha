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
    int father, mother, child;      // VCF sample index
}
trio_t;

typedef struct
{
    int *is_parent;
    int *is_child;
    int *prev_flip;
    bcf_hdr_t *in_hdr, *out_hdr;
    int ibd;
    int *ibd_father;
    int *ibd_mother;
    trio_t *trio;
    int ntrio, mtrio;
    int32_t *gt_arr;
    int mgt_arr, prev_rid;
}
args_t;

args_t *args;

const char *about(void)
{
    return "Phase genotypes in trios.\n";
}

const char *usage(void)
{
    return
"\n"
"About: Phase genotypes in trios.\n"
"Usage: bcftools +trio-phase [General Options] -- [Plugin Options]\n"
"Options:\n"
"   run \"bcftools plugin\" for a list of common options\n"
"\n"
"Plugin options:\n"
"   -p, --ped <file>            PED file\n"
"   -i, --ibd                   Whether to add IBD state for duos\n"
"\n"
"Example:\n"
"   bcftools +trio-phase file.bcf -- -p file.ped\n"
"\n";
}

// adapted from Petr Danecek's parse_ped in bcftools/plugins/trio-switch-rate.c
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

        int father = bcf_hdr_id2int(args->in_hdr,BCF_DT_SAMPLE,&str.s[off[2]]);
        int mother = bcf_hdr_id2int(args->in_hdr,BCF_DT_SAMPLE,&str.s[off[3]]);
        int child = bcf_hdr_id2int(args->in_hdr,BCF_DT_SAMPLE,&str.s[off[1]]);
        if ( ( father<0 && mother<0 ) || child<0 ) continue;
        if ( father!=-1 ) args->is_parent[father] = 1;
        if ( mother!=-1 ) args->is_parent[mother] = 1;
        args->is_child[child] = 1;

        args->ntrio++;
        hts_expand0(trio_t,args->ntrio,args->mtrio,args->trio);
        trio_t *trio = &args->trio[args->ntrio-1];
        trio->father = father;
        trio->mother = mother;
        trio->child  = child;

    } while ( hts_getline(fp, KS_SEP_LINE, &str)>=0 );

    free(str.s);
    free(off);
    hts_close(fp);
}

int init(int argc, char **argv, bcf_hdr_t *in, bcf_hdr_t *out)
{
    args = (args_t *)calloc(1, sizeof(args_t));
    args->prev_rid = -1;
    int nsmpl = bcf_hdr_nsamples(in);
    args->is_parent = (int *)calloc((size_t)nsmpl, sizeof(int));
    args->is_child = (int *)calloc((size_t)nsmpl, sizeof(int));
    args->prev_flip = (int *)calloc((size_t)nsmpl, sizeof(int));
    args->ibd_father = (int *)calloc((size_t)nsmpl, sizeof(int));
    args->ibd_mother = (int *)calloc((size_t)nsmpl, sizeof(int));
    args->in_hdr = in;
    args->out_hdr = out;
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

typedef struct
{
    int a, b, is_het, is_phased;
}
gt_t;

// extract alleles and phase from diploid genotype
static int parse_genotype(gt_t *gt, int32_t *gt_arr, int id)
{
    if ( id<0 ) return 0;
    int32_t *ptr = gt_arr + id;
    if ( ptr[0]==bcf_gt_missing ) return 0;
    if ( ptr[1]==bcf_gt_missing ) return 0;
    if ( ptr[1]==bcf_int32_vector_end ) return 0;
    gt->a = bcf_gt_allele(ptr[0]);
    gt->b = bcf_gt_allele(ptr[1]);
    gt->is_het = gt->a != gt->b;
    gt->is_phased = gt->is_het ? (bcf_gt_is_phased(ptr[1]) ? 1 : 0) : 1;
    return 1;
}

// flips genotypes and set genotype as phased
static void flip_gt(int32_t *gt)
{
    *gt ^= *(gt+1);
    *(gt+1) ^= *gt;
    *gt ^= *(gt+1);
    *gt &= 0xFFFFFFFE;
    *(gt+1) |= 1;
    return;
}

// make a decision whether to flip a genotype
static void cast_vote(int32_t *gt, int vote_flip, int *prev_flip)
{
    int is_phased = bcf_gt_is_phased(*(gt+1));
    if ( vote_flip )
    {
        if ( is_phased ) *prev_flip = vote_flip > 0;
        if( vote_flip > 0 ) flip_gt(gt);
        else *(gt+1) |= 1;
    }
    else
    {
        if( is_phased && *prev_flip ) flip_gt(gt);
    }
}

bcf1_t *process(bcf1_t *rec)
{
    // extract genotypes from record and checks whether they follow a diploid model
    int nsmpl = bcf_hdr_nsamples(args->in_hdr);
    int ngt = bcf_get_genotypes(args->in_hdr, rec, &args->gt_arr, &args->mgt_arr);
    if ( ngt<0 ) return rec;
    ngt /= nsmpl;
    if ( ngt!=2 ) return rec;

    // check whether the chromosome in the record has changed
    int i;
    if ( rec->rid!=args->prev_rid )
    {
        args->prev_rid = rec->rid;
        for (i=0; i<nsmpl; i++)
        {
            args->prev_flip[i] = 0;
        }
    }

    // this will determine whether heterozygous genotypes need to be flipped
    trio_t *trio;
    gt_t child, father, mother, parent;
    int is_father, is_mother, *vote, *vote_arr;
    vote_arr = (int *)calloc((size_t)nsmpl, sizeof(int));

    for (i=0; i<args->ntrio; i++)
    {
        trio = args->trio + i;

        // parse genotypes in the trio
        if ( !parse_genotype(&child, args->gt_arr, ngt*trio->child) ) continue;
        is_father = parse_genotype(&father, args->gt_arr, ngt*trio->father);
        is_mother = parse_genotype(&mother, args->gt_arr, ngt*trio->mother);

        // propagate information from parents to child
        if ( child.is_het )
        {
            vote = vote_arr + trio->child;
            if ( is_father )
            {
                if ( child.a == father.a || child.a == father.b ) (*vote)--;
                if ( child.b == father.a || child.b == father.b ) (*vote)++;
            }
            if ( is_mother )
            {
                if ( child.a == mother.a || child.a == mother.b ) (*vote)++;
                if ( child.b == mother.a || child.b == mother.b ) (*vote)--;
            }
            cast_vote(args->gt_arr + ngt*trio->child, *vote, args->prev_flip + trio->child);
        }

        // parse genotype of the child
        if ( !parse_genotype(&child, args->gt_arr, ngt*trio->child) ) continue;

        // propagate information from child to father
        if ( is_father && father.is_het && child.is_phased )
        {
            vote = vote_arr + trio->father;
            if ( !(*vote) && child.a == father.a && args->ibd_father[trio->child] == IBD_A ) (*vote)--;
            if ( !(*vote) && child.a == father.a && args->ibd_father[trio->child] == IBD_B ) (*vote)++;
            if ( !(*vote) && child.a == father.b && args->ibd_father[trio->child] == IBD_A ) (*vote)++;
            if ( !(*vote) && child.a == father.b && args->ibd_father[trio->child] == IBD_B ) (*vote)--;
        }

        // propagate information from child to mother
        if ( is_mother && mother.is_het && child.is_phased )
        {
            vote = vote_arr + trio->mother;
            if ( !(*vote) && child.b == mother.a && args->ibd_mother[trio->child] == IBD_A ) (*vote)--;
            if ( !(*vote) && child.b == mother.a && args->ibd_mother[trio->child] == IBD_B ) (*vote)++;
            if ( !(*vote) && child.b == mother.b && args->ibd_mother[trio->child] == IBD_A ) (*vote)++;
            if ( !(*vote) && child.b == mother.b && args->ibd_mother[trio->child] == IBD_B ) (*vote)--;
        }
    }

    // propagate information from children to parents
    for (i=0; i<nsmpl; i++)
    {
        if ( !args->is_parent[i] ) continue;
        if ( !parse_genotype(&parent, args->gt_arr, ngt*i) ) continue;
        if ( parent.is_phased && args->is_child[i] ) continue;
        if ( parent.is_het ) cast_vote(args->gt_arr + ngt*i, vote_arr[i], args->prev_flip + i);
    }

    // update IBD states in trio
    for (i=0; i<args->ntrio; i++)
    {
        trio = args->trio + i;

        // parse genotypes in the trio
        if ( !parse_genotype(&child, args->gt_arr, ngt*trio->child) ) continue;
        is_father = parse_genotype(&father, args->gt_arr, ngt*trio->father);
        is_mother = parse_genotype(&mother, args->gt_arr, ngt*trio->mother);

        if ( is_father && father.is_het && father.is_phased && child.is_phased )
        {
            if ( child.a == father.a ) args->ibd_father[trio->child] = IBD_A;
            else if ( child.a == father.b ) args->ibd_father[trio->child] = IBD_B;
        }

        if ( is_mother && mother.is_het && mother.is_phased && child.is_phased )
        {
            if ( child.b == mother.a ) args->ibd_mother[trio->child] = IBD_A;
            else if ( child.b == mother.b ) args->ibd_mother[trio->child] = IBD_B;
        }
    }

    free(vote_arr);
    bcf_update_genotypes(args->out_hdr, rec, args->gt_arr, args->mgt_arr);
    if ( args->ibd )
    {
        bcf_update_format_int32(args->out_hdr, rec, "IBD_F", args->ibd_father, nsmpl);
        bcf_update_format_int32(args->out_hdr, rec, "IBD_M", args->ibd_mother, nsmpl);
    }
    return rec;
}

void destroy(void)
{
    free(args->is_parent);
    free(args->is_child);
    free(args->prev_flip);
    free(args->trio);
    free(args->gt_arr);
    free(args);
}
