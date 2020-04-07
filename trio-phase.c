/* The MIT License

   Copyright (C) 2017-2020 Giulio Genovese

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

#define TRIO_PHASE_VERSION "2020-04-07"

#define ABSOLUTE (1 << 24)
#define TRANSMITTED (1 << 16)
#define RELATIVE (1 << 8)

typedef struct {
	int *a;
	int n;
	int m;
} set_t;

typedef struct {
	set_t fathers;
	set_t mothers;
	int32_t gt[2];
	int is_child;
	int flip_vote;
	int is_prev_het_flipped;
	int32_t ibd_father;
	int32_t ibd_mother;
} ind_t;

typedef struct {
	bcf_hdr_t *in_hdr, *out_hdr;
	ind_t *inds;
	set_t children;
	int prev_rid, ibd;
	int32_t *arr;
	int m_arr;
	int niters;
} args_t;

args_t *args;

static int add_id(set_t *set, int id)
{
	for (int i = 0; i < set->n; i++)
		if (id == set->a[i])
			return -1;
	set->n++;
	hts_expand(int, set->n, set->m, set->a);
	set->a[set->n - 1] = id;
	return 0;
}

static int inline gt_is_missing(int32_t *gt)
{
	return gt[0] == bcf_gt_missing || gt[1] == bcf_gt_missing
	       || gt[0] == bcf_int32_vector_end;
}
static int inline gt_is_het(int32_t *gt)
{
	return bcf_gt_allele(gt[0]) != bcf_gt_allele(gt[1]);
}
static int inline gt_is_phased(int32_t *gt)
{
	return gt_is_het(gt) ? (bcf_gt_is_phased(gt[1]) ? 1 : 0) : 1;
}
static int32_t inline gt_pat(int32_t *gt)
{
	return bcf_gt_allele(gt[0]);
}
static int32_t inline gt_mat(int32_t *gt)
{
	return bcf_gt_allele(gt[1]);
}

// flips genotypes and set genotype as phased
static void inline flip_gt(int32_t *gt)
{
	*gt ^= *(gt + 1);
	*(gt + 1) ^= *gt;
	*gt ^= *(gt + 1);
	*gt &= 0xFFFFFFFE;
	*(gt + 1) |= 1;
	return;
}

// propagate information across individuals
static void ind_vote(ind_t *ind)
{
	if (!gt_is_het(ind->gt))
		return;
	if (ind->flip_vote > 0) {
		flip_gt(ind->gt);
		ind->flip_vote = -ind->flip_vote;
	} else if (ind->flip_vote < 0) {
		*(ind->gt + 1) |= 1;
	}
}

// adds one to break the ties
static void inline inc_vote(int *vote, int value)
{
	if (*vote == 0)
		*vote += value > 0 ? 1 : -1;
	*vote += value;
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
	return "\n"
	       "About: Phase genotypes in trios. (version " TRIO_PHASE_VERSION
	       " https://github.com/freeseek/mocha)\n"
	       "Usage: bcftools +trio-phase [General Options] -- [Plugin Options]\n"
	       "\n"
	       "Options:\n"
	       "   run \"bcftools plugin\" for a list of common options\n"
	       "\n"
	       "Plugin options:\n"
	       "   -p, --ped <file>            PED file\n"
	       "   -i, --ibd                   Whether to add IBD state for duos\n"
	       "   -n, --niters                Number of iterations when propagating information [3]\n"
	       "\n"
	       "Example:\n"
	       "   bcftools +trio-phase file.bcf -- --ped file.ped\n"
	       "\n";
}

static void parse_ped(args_t *args, char *fname)
{
	htsFile *fp = hts_open(fname, "r");
	if (!fp)
		error("Could not read: %s\n", fname);

	kstring_t str = {0, 0, 0};
	if (hts_getline(fp, KS_SEP_LINE, &str) <= 0)
		error("Empty file: %s\n", fname);

	int moff = 0, *off = NULL;
	do {
		int ncols = ksplit_core(str.s, 0, &moff, &off);
		if (ncols < 4)
			error("Could not parse the ped file: %s\n", str.s);

		int father_id = bcf_hdr_id2int(args->in_hdr, BCF_DT_SAMPLE, &str.s[off[2]]);
		int mother_id = bcf_hdr_id2int(args->in_hdr, BCF_DT_SAMPLE, &str.s[off[3]]);
		int child_id = bcf_hdr_id2int(args->in_hdr, BCF_DT_SAMPLE, &str.s[off[1]]);
		if ((father_id < 0 && mother_id < 0) || child_id < 0)
			continue;
		add_id(&args->children, child_id);
		ind_t *child_ind = &args->inds[child_id];
		child_ind->is_child = 1;
		if (father_id >= 0)
			add_id(&child_ind->fathers, father_id);
		if (mother_id >= 0)
			add_id(&child_ind->mothers, mother_id);
	} while (hts_getline(fp, KS_SEP_LINE, &str) >= 0);

	free(str.s);
	free(off);
	hts_close(fp);
}

int init(int argc, char **argv, bcf_hdr_t *in, bcf_hdr_t *out)
{
	args = (args_t *)calloc(1, sizeof(args_t));
	args->prev_rid = -1;
	args->in_hdr = in;
	args->out_hdr = out;
	args->inds = (ind_t *)calloc(bcf_hdr_nsamples(in), sizeof(ind_t));
	for (int i = 0; i < bcf_hdr_nsamples(in); i++) {
		args->inds[i].ibd_father = bcf_int32_missing;
		args->inds[i].ibd_mother = bcf_int32_missing;
	}
	args->niters = 3;
	char *ped_fname = NULL;

	static struct option loptions[] = {{"ped", required_argument, NULL, 'p'},
					   {"ibd", no_argument, NULL, 'i'},
					   {"niters", required_argument, NULL, 'n'},
					   {0, 0, 0, 0}};
	int c;
	char *tmp = NULL;
	while ((c = getopt_long(argc, argv, "?hp:in:", loptions, NULL)) >= 0) {
		switch (c) {
		case 'p':
			ped_fname = optarg;
			break;
		case 'i':
			args->ibd = 1;
			break;
		case 'n':
			args->niters = (int)strtol(optarg, &tmp, 0);
			if (*tmp)
				error("Could not parse: --threads %s\n", optarg);
			break;
		case 'h':
		case '?':
		default:
			error("%s", usage());
			break;
		}
	}
	if (!ped_fname)
		error("Expected the -p option\n");
	parse_ped(args, ped_fname);
	if (args->ibd) {
		bcf_hdr_append(
			args->out_hdr,
			"##FORMAT=<ID=IBD_F,Number=1,Type=Integer,Description=\"IBD state with the father\">");
		bcf_hdr_append(
			args->out_hdr,
			"##FORMAT=<ID=IBD_M,Number=1,Type=Integer,Description=\"IBD state with the mother\">");
	}
	return 0;
}

bcf1_t *process(bcf1_t *rec)
{
	// extract genotypes from record and checks whether they follow a diploid model
	int nsmpl = bcf_hdr_nsamples(args->in_hdr);
	int ngt = bcf_get_genotypes(args->in_hdr, rec, &args->arr, &args->m_arr);
	if (ngt < 0)
		return rec;
	int ploidy = ngt / nsmpl;
	if (ploidy != 2)
		return rec;

	// check whether the chromosome in the record has changed
	if (rec->rid != args->prev_rid) {
		args->prev_rid = rec->rid;
		for (int i = 0; i < nsmpl; i++) {
			ind_t *sample_ind = &args->inds[i];
			sample_ind->is_prev_het_flipped = 0;
		}
		for (int i = 0; i < args->children.n; i++) {
			int child_id = args->children.a[i];
			ind_t *child_ind = &args->inds[child_id];
			child_ind->ibd_father = 0;
			child_ind->ibd_mother = 0;
		}
	}

	// copy genotypes inside individual structure while propagating information across
	// consecutive heterozygous sites
	for (int i = 0; i < nsmpl; i++) {
		ind_t *sample_ind = &args->inds[i];
		sample_ind->gt[0] = args->arr[ploidy * i];
		sample_ind->gt[1] = args->arr[ploidy * i + 1];
		if (gt_is_het(sample_ind->gt) && gt_is_phased(sample_ind->gt)
		    && sample_ind->is_prev_het_flipped)
			flip_gt(sample_ind->gt);
	}

	// run multiple times to ensure convergence
	for (int j = 0; j < args->niters; j++) {
		for (int i = 0; i < nsmpl; i++)
			args->inds[i].flip_vote = 0;

		// propagate absolute and relative phase information from parents to children
		for (int i = 0; i < args->children.n; i++) {
			int child_id = args->children.a[i];
			ind_t *child_ind = &args->inds[child_id];
			if (gt_is_missing(child_ind->gt) || !gt_is_het(child_ind->gt))
				continue;

			int absolute_vote = 0;
			int relative_vote = 0;
			for (int j = 0; j < child_ind->fathers.n; j++) {
				int father_id = child_ind->fathers.a[j];
				int32_t *father_gt = args->inds[father_id].gt;
				if (gt_pat(child_ind->gt) == gt_pat(father_gt)
				    && gt_pat(child_ind->gt) == gt_mat(father_gt))
					absolute_vote--;
				else if (gt_mat(child_ind->gt) == gt_pat(father_gt)
					 && gt_mat(child_ind->gt) == gt_mat(father_gt))
					absolute_vote++;
				else if (gt_is_phased(father_gt)) {
					if (gt_pat(child_ind->gt)
					    == bcf_gt_allele(father_gt[child_ind->ibd_father]))
						relative_vote--;
					else if (gt_mat(child_ind->gt)
						 == bcf_gt_allele(
							 father_gt[child_ind->ibd_father]))
						relative_vote++;
				}
			}
			if (absolute_vote > 0)
				inc_vote(&child_ind->flip_vote, ABSOLUTE);
			else if (absolute_vote < 0)
				inc_vote(&child_ind->flip_vote, -ABSOLUTE);
			else if (relative_vote > 0)
				inc_vote(&child_ind->flip_vote, RELATIVE);
			else if (relative_vote < 0)
				inc_vote(&child_ind->flip_vote, -RELATIVE);

			absolute_vote = 0;
			relative_vote = 0;
			for (int j = 0; j < child_ind->mothers.n; j++) {
				int mother_id = child_ind->mothers.a[j];
				int32_t *mother_gt = args->inds[mother_id].gt;
				if (gt_pat(child_ind->gt) == gt_pat(mother_gt)
				    && gt_pat(child_ind->gt) == gt_mat(mother_gt)) {
					absolute_vote++;
				} else if (gt_mat(child_ind->gt) == gt_pat(mother_gt)
					   && gt_mat(child_ind->gt) == gt_mat(mother_gt)) {
					absolute_vote--;
				} else if (gt_is_phased(mother_gt)) {
					if (gt_pat(child_ind->gt)
					    == bcf_gt_allele(mother_gt[child_ind->ibd_mother]))
						relative_vote++;
					else if (gt_mat(child_ind->gt)
						 == bcf_gt_allele(
							 mother_gt[child_ind->ibd_mother]))
						relative_vote--;
				}
			}
			if (absolute_vote > 0)
				inc_vote(&child_ind->flip_vote, ABSOLUTE);
			else if (absolute_vote < 0)
				inc_vote(&child_ind->flip_vote, -ABSOLUTE);
			else if (relative_vote > 0)
				inc_vote(&child_ind->flip_vote, RELATIVE);
			else if (relative_vote < 0)
				inc_vote(&child_ind->flip_vote, -RELATIVE);

			ind_vote(child_ind);
		}

		// propagate relative phase information from children to parents
		for (int i = 0; i < args->children.n; i++) {
			int child_id = args->children.a[i];
			ind_t *child_ind = &args->inds[child_id];

			if (gt_is_missing(child_ind->gt)
			    || (gt_is_het(child_ind->gt) && !gt_is_phased(child_ind->gt)))
				continue;

			for (int j = 0; j < child_ind->fathers.n; j++) {
				int father_id = child_ind->fathers.a[j];
				ind_t *father_ind = &args->inds[father_id];

				if (gt_pat(child_ind->gt)
				    == bcf_gt_allele(father_ind->gt[child_ind->ibd_father]))
					inc_vote(&father_ind->flip_vote, -TRANSMITTED);
				if (gt_pat(child_ind->gt)
				    == bcf_gt_allele(father_ind->gt[1 - child_ind->ibd_father]))
					inc_vote(&father_ind->flip_vote, TRANSMITTED);
			}

			for (int j = 0; j < child_ind->mothers.n; j++) {
				int mother_id = child_ind->mothers.a[j];
				ind_t *mother_ind = &args->inds[mother_id];

				if (gt_mat(child_ind->gt)
				    == bcf_gt_allele(mother_ind->gt[child_ind->ibd_mother]))
					inc_vote(&mother_ind->flip_vote, -TRANSMITTED);
				if (gt_mat(child_ind->gt)
				    == bcf_gt_allele(mother_ind->gt[1 - child_ind->ibd_mother]))
					inc_vote(&mother_ind->flip_vote, TRANSMITTED);
			}
		}

		// propagate relative phase information from children to parents
		for (int i = 0; i < nsmpl; i++)
			ind_vote(&args->inds[i]);
	}

	// update IBD states in children (if both child and parent are heterozygous and phased)
	for (int i = 0; i < args->children.n; i++) {
		int child_id = args->children.a[i];
		ind_t *child_ind = &args->inds[child_id];
		if (gt_is_missing(child_ind->gt) || !gt_is_het(child_ind->gt)
		    || !gt_is_phased(child_ind->gt))
			continue;

		int vote = 0;
		for (int j = 0; j < child_ind->fathers.n; j++) {
			int father_id = child_ind->fathers.a[j];
			int32_t *father_gt = args->inds[father_id].gt;
			if (gt_is_missing(father_gt) || !gt_is_het(father_gt)
			    || !gt_is_phased(father_gt))
				continue;

			if (gt_pat(child_ind->gt) == bcf_gt_allele(father_gt[0]))
				vote++;
			else if (gt_pat(child_ind->gt) == bcf_gt_allele(father_gt[1]))
				vote--;
		}
		if (vote > 0)
			child_ind->ibd_father = 0;
		else if (vote < 0)
			child_ind->ibd_father = 1;

		vote = 0;
		for (int j = 0; j < child_ind->mothers.n; j++) {
			int mother_id = child_ind->mothers.a[j];
			int32_t *mother_gt = args->inds[mother_id].gt;
			if (gt_is_missing(mother_gt) || !gt_is_het(mother_gt)
			    || !gt_is_phased(mother_gt))
				continue;

			if (gt_mat(child_ind->gt) == bcf_gt_allele(mother_gt[0]))
				vote++;
			else if (gt_mat(child_ind->gt) == bcf_gt_allele(mother_gt[1]))
				vote--;
		}
		if (vote > 0)
			child_ind->ibd_mother = 0;
		else if (vote < 0)
			child_ind->ibd_mother = 1;
	}

	// copy individual structure back to genotypes and update is_prev_het_flipped values
	for (int i = 0; i < nsmpl; i++) {
		ind_t *sample_ind = &args->inds[i];
		if (gt_is_het(sample_ind->gt) && gt_is_phased(sample_ind->gt)
		    && gt_is_het(&args->arr[ploidy * i])
		    && gt_is_phased(&args->arr[ploidy * i])) {
			if (bcf_gt_allele(sample_ind->gt[0])
				    == bcf_gt_allele(args->arr[ploidy * i])
			    && bcf_gt_allele(sample_ind->gt[1])
				       == bcf_gt_allele(args->arr[ploidy * i + 1]))
				sample_ind->is_prev_het_flipped = 0;
			else if (bcf_gt_allele(sample_ind->gt[0])
					 == bcf_gt_allele(args->arr[ploidy * i + 1])
				 && bcf_gt_allele(sample_ind->gt[1])
					    == bcf_gt_allele(args->arr[ploidy * i]))
				sample_ind->is_prev_het_flipped = 1;
		}
		args->arr[ploidy * i] = sample_ind->gt[0];
		args->arr[ploidy * i + 1] = sample_ind->gt[1];
	}
	bcf_update_genotypes(args->out_hdr, rec, args->arr, args->m_arr);

	if (args->ibd) {
		for (int i = 0; i < nsmpl; i++)
			args->arr[i] = args->inds[i].ibd_father;
		bcf_update_format_int32(args->out_hdr, rec, "IBD_F", args->arr, nsmpl);
		for (int i = 0; i < nsmpl; i++)
			args->arr[i] = args->inds[i].ibd_mother;
		bcf_update_format_int32(args->out_hdr, rec, "IBD_M", args->arr, nsmpl);
	}
	return rec;
}

void destroy(void)
{
	for (int i = 0; i < bcf_hdr_nsamples(args->in_hdr); i++)
		ind_destroy(&args->inds[i]);
	free(args->inds);
	free(args->children.a);
	free(args->arr);
	free(args);
}
