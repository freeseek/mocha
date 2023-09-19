/* The MIT License

   Copyright (C) 2018-2023 Giulio Genovese

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

#include <htslib/kseq.h>
#include "bcftools.h"

typedef struct {
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

typedef struct {
    const char *alias, *about, *rules;
} rules_predef_t;

// the following definitions were derived from these files:
// http://hgdownload.cse.ucsc.edu/goldenPath/hg18/database/chr{{1..22},{X,Y}}_gap.txt.gz
// http://hgdownload.cse.ucsc.edu/goldenPath/hg18/database/genomicSuperDups.txt.gz
// http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/gap.txt.gz
// http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/genomicSuperDups.txt.gz
// http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/centromeres.txt.gz
// http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/genomicSuperDups.txt.gz
static rules_predef_t rules_predefs[] = {
    {.alias = "NCBI36",
     .about = "Human Genome reference assembly NCBI36 / hg18, both chr naming conventions",
     .rules = "   1:121236957-123476957     centromere\n"
              "   2:91689898-94689898       centromere\n"
              "   3:90587544-93487544       centromere\n"
              "   4:49354874-52354874       centromere\n"
              "   5:46441398-49441398       centromere\n"
              "   6:58938125-61938125       centromere\n"
              "   7:58058273-61058273       centromere\n"
              "   8:43958052-46958052       centromere\n"
              "   9:47107499-50107499       centromere\n"
              "   10:39244941-41624941      centromere\n"
              "   11:51450781-54450781      centromere\n"
              "   12:34747961-36142961      centromere\n"
              "   13:16000000-17868000      centromere\n"
              "   14:15070000-18070000      centromere\n"
              "   15:15260000-18260000      centromere\n"
              "   16:35143302-36943302      centromere\n"
              "   17:22187133-22287133      centromere\n"
              "   18:15400898-16764896      centromere\n"
              "   19:26923622-29923622      centromere\n"
              "   20:26267569-28033230      centromere\n"
              "   21:10260000-13260000      centromere\n"
              "   22:11330000-14330000      centromere\n"
              "   X:58598737-61598737       centromere\n"
              "   Y:11253954-12308578       centromere\n"
              "   X:2709520-154583483       X_nonpar\n"
              "   X:88343457-92262165       X_xtr\n"
              "   Y:2709520-57442674        Y_nonpar\n"
              "   Y:2977958-6676600         Y_xtr\n"
              "   MT:1-16571                mitochondria\n"
              "\n"
              "   chr1:121236957-123476957  centromere\n"
              "   chr2:91689898-94689898    centromere\n"
              "   chr3:90587544-93487544    centromere\n"
              "   chr4:49354874-52354874    centromere\n"
              "   chr5:46441398-49441398    centromere\n"
              "   chr6:58938125-61938125    centromere\n"
              "   chr7:58058273-61058273    centromere\n"
              "   chr8:43958052-46958052    centromere\n"
              "   chr9:47107499-50107499    centromere\n"
              "   chr10:39244941-41624941   centromere\n"
              "   chr11:51450781-54450781   centromere\n"
              "   chr12:34747961-36142961   centromere\n"
              "   chr13:16000000-17868000   centromere\n"
              "   chr14:15070000-18070000   centromere\n"
              "   chr15:15260000-18260000   centromere\n"
              "   chr16:35143302-36943302   centromere\n"
              "   chr17:22187133-22287133   centromere\n"
              "   chr18:15400898-16764896   centromere\n"
              "   chr19:26923622-29923622   centromere\n"
              "   chr20:26267569-28033230   centromere\n"
              "   chr21:10260000-13260000   centromere\n"
              "   chr22:11330000-14330000   centromere\n"
              "   chrX:58598737-61598737    centromere\n"
              "   chrY:11253954-12308578    centromere\n"
              "   chrX:2709520-154583483    X_nonpar\n"
              "   chrX:88343457-92262165    X_xtr\n"
              "   chrY:2709520-57442674     Y_nonpar\n"
              "   chrY:2977958-6676600      Y_xtr\n"
              "   chrM:1-16571              mitochondria\n"},
    {.alias = "GRCh37",
     .about = "Human Genome reference assembly GRCh37 / hg19, both chr naming conventions",
     .rules = "   1:121535434-124535434     centromere\n"
              "   2:92326171-95326171       centromere\n"
              "   3:90504854-93504854       centromere\n"
              "   4:49660117-52660117       centromere\n"
              "   5:46405641-49405641       centromere\n"
              "   6:58830166-61830166       centromere\n"
              "   7:58054331-61054331       centromere\n"
              "   8:43838887-46838887       centromere\n"
              "   9:47367679-50367679       centromere\n"
              "   10:39254935-42254935      centromere\n"
              "   11:51644205-54644205      centromere\n"
              "   12:34856694-37856694      centromere\n"
              "   13:16000000-19000000      centromere\n"
              "   14:16000000-19000000      centromere\n"
              "   15:17000000-20000000      centromere\n"
              "   16:35335801-38335801      centromere\n"
              "   17:22263006-25263006      centromere\n"
              "   18:15460898-18460898      centromere\n"
              "   19:24681782-27681782      centromere\n"
              "   20:26369569-29369569      centromere\n"
              "   21:11288129-14288129      centromere\n"
              "   22:13000000-16000000      centromere\n"
              "   X:58632012-61632012       centromere\n"
              "   Y:10104553-13104553       centromere\n"
              "   X:2699520-154930289       X_nonpar\n"
              "   X:88456801-92375509       X_xtr\n"
              "   Y:2649520-59033286        Y_nonpar\n"
              "   Y:2917958-6616600         Y_xtr\n"
              "   MT:1-16569                mitochondria\n"
              "\n"
              "   chr1:121535434-124535434  centromere\n"
              "   chr2:92326171-95326171    centromere\n"
              "   chr3:90504854-93504854    centromere\n"
              "   chr4:49660117-52660117    centromere\n"
              "   chr5:46405641-49405641    centromere\n"
              "   chr6:58830166-61830166    centromere\n"
              "   chr7:58054331-61054331    centromere\n"
              "   chr8:43838887-46838887    centromere\n"
              "   chr9:47367679-50367679    centromere\n"
              "   chr10:39254935-42254935   centromere\n"
              "   chr11:51644205-54644205   centromere\n"
              "   chr12:34856694-37856694   centromere\n"
              "   chr13:16000000-19000000   centromere\n"
              "   chr14:16000000-19000000   centromere\n"
              "   chr15:17000000-20000000   centromere\n"
              "   chr16:35335801-38335801   centromere\n"
              "   chr17:22263006-25263006   centromere\n"
              "   chr18:15460898-18460898   centromere\n"
              "   chr19:24681782-27681782   centromere\n"
              "   chr20:26369569-29369569   centromere\n"
              "   chr21:11288129-14288129   centromere\n"
              "   chr22:13000000-16000000   centromere\n"
              "   chrX:58632012-61632012    centromere\n"
              "   chrY:10104553-13104553    centromere\n"
              "   chrX:2699520-154930289    X_nonpar\n"
              "   chrX:88456801-92375509    X_xtr\n"
              "   chrY:2649520-59033286     Y_nonpar\n"
              "   chrY:2917958-6616600      Y_xtr\n"
              "   chrM:1-16571              mitochondria\n"},
    {.alias = "GRCh38",
     .about = "Human Genome reference assembly GRCh38 / hg38, both chr naming conventions",
     .rules = "   1:122026459-124932724     centromere\n"
              "   2:92188145-94090557       centromere\n"
              "   3:90772458-93655574       centromere\n"
              "   4:49712061-51743951       centromere\n"
              "   5:46485900-50059807       centromere\n"
              "   6:58553888-59829934       centromere\n"
              "   7:58169653-61528020       centromere\n"
              "   8:44033744-45877265       centromere\n"
              "   9:43389635-45518558       centromere\n"
              "   10:39686682-41593521      centromere\n"
              "   11:51078348-54425074      centromere\n"
              "   12:34769407-37185252      centromere\n"
              "   13:16000000-18051248      centromere\n"
              "   14:16000000-18173523      centromere\n"
              "   15:17083673-19725254      centromere\n"
              "   16:36311158-38265669      centromere\n"
              "   17:22813679-26616164      centromere\n"
              "   18:15460899-20861206      centromere\n"
              "   19:24498980-27190874      centromere\n"
              "   20:26436232-30038348      centromere\n"
              "   21:10864560-12915808      centromere\n"
              "   22:12954788-15054318      centromere\n"
              "   X:58605579-62412542       centromere\n"
              "   Y:10316944-10544039       centromere\n"
              "   X:2781479-155700628       X_nonpar\n"
              "   X:89201802-93120510       X_xtr\n"
              "   Y:2781479-56887139        Y_nonpar\n"
              "   Y:3049917-6748559         Y_xtr\n"
              "   MT:1-16569                mitochondria\n"
              "\n"
              "   chr1:122026459-124932724  centromere\n"
              "   chr2:92188145-94090557    centromere\n"
              "   chr3:90772458-93655574    centromere\n"
              "   chr4:49712061-51743951    centromere\n"
              "   chr5:46485900-50059807    centromere\n"
              "   chr6:58553888-59829934    centromere\n"
              "   chr7:58169653-61528020    centromere\n"
              "   chr8:44033744-45877265    centromere\n"
              "   chr9:43389635-45518558    centromere\n"
              "   chr10:39686682-41593521   centromere\n"
              "   chr11:51078348-54425074   centromere\n"
              "   chr12:34769407-37185252   centromere\n"
              "   chr13:16000000-18051248   centromere\n"
              "   chr14:16000000-18173523   centromere\n"
              "   chr15:17083673-19725254   centromere\n"
              "   chr16:36311158-38265669   centromere\n"
              "   chr17:22813679-26616164   centromere\n"
              "   chr18:15460899-20861206   centromere\n"
              "   chr19:24498980-27190874   centromere\n"
              "   chr20:26436232-30038348   centromere\n"
              "   chr21:10864560-12915808   centromere\n"
              "   chr22:12954788-15054318   centromere\n"
              "   chrX:58605579-62412542    centromere\n"
              "   chrY:10316944-10544039    centromere\n"
              "   chrX:2781479-155700628    X_nonpar\n"
              "   chrX:89201802-93120510    X_xtr\n"
              "   chrY:2781479-56887139     Y_nonpar\n"
              "   chrY:3049917-6748559      Y_xtr\n"
              "   chrM:1-16569              mitochondria\n"},
    {
        .alias = NULL,
        .about = NULL,
        .rules = NULL,
    }};

genome_rules_t *genome_init(const bcf_hdr_t *hdr) {
    genome_rules_t *self = (genome_rules_t *)calloc(1, sizeof(genome_rules_t));
    int n = hdr->n[BCF_DT_CTG];
    self->length = (int *)calloc(n, sizeof(int));
    for (int rid = 0; rid < n; rid++) self->length[rid] = hdr->id[BCF_DT_CTG][rid].val->info[0];
    self->cen_beg = (int *)calloc(n, sizeof(int));
    self->cen_end = (int *)calloc(n, sizeof(int));
    self->is_short_arm = (int *)calloc(n, sizeof(int));
    self->x_rid = -1;
    self->y_rid = -1;
    self->mt_rid = -1;
    return self;
}

void genome_destroy(genome_rules_t *self) {
    free(self->length);
    free(self->cen_beg);
    free(self->cen_end);
    free(self->is_short_arm);
    free(self);
}

// adapted from Petr Danecek's implementation of parse_rules() in bcftools/plugins/mendelian.c
static void push_rule(genome_rules_t *self, char *line, const bcf_hdr_t *hdr) {
    // eat any leading spaces
    char *ss = (char *)line;
    while (*ss && isspace(*ss)) ss++;
    if (!*ss) return; // skip empty lines

    // chromosome name, beg, end
    char *tmp, *se = ss;
    while (*se && !isspace(*se) && *se != ':') se++;
    if (*se != ':') error("Could not parse the region: %s\n", line);
    *se = '\0'; // terminates the chromosome string
    int rid = bcf_hdr_name2id(hdr, ss);
    if (rid < 0) return;
    *se = ':'; // restores the separator
    ss = ++se;
    while (*se && !isspace(*se) && *se != '-') se++;
    if (*se != '-') error("Could not parse the region: %s\n", line);
    int beg = (int)strtol(ss, &tmp, 0);
    if (tmp == ss) error("Could not parse the region: %s\n", line);
    ss = ++se;
    int end = (int)strtol(ss, &tmp, 0);
    if (tmp == ss || beg > end) error("Could not parse the region: %s\n", line);

    // skip region
    while (*ss && !isspace(*ss)) ss++;
    while (*ss && isspace(*ss)) ss++;
    if (strncmp(ss, "centromere", 10) == 0) {
        if (self->cen_beg[rid] != 0 || self->cen_end[rid] != 0) error("Second centromere rule %s\n", line);
        self->cen_beg[rid] = beg;
        self->cen_end[rid] = end;
    } else if (strncmp(ss, "X_nonpar", 8) == 0) {
        if ((self->x_xtr_beg != 0 || self->x_xtr_end != 0) && self->x_rid != rid)
            error("Chromosome X XTR and nonPAR regions declared on different contigs: %s\n", line);
        self->x_rid = rid;
        if (self->x_nonpar_beg != 0 || self->x_nonpar_end != 0) error("Second chromosome X nonPAR rule: %s\n", line);
        self->x_nonpar_beg = beg;
        self->x_nonpar_end = end;
    } else if (strncmp(ss, "X_xtr", 5) == 0) {
        if ((self->x_nonpar_beg != 0 || self->x_nonpar_end != 0) && self->x_rid != rid)
            error("Chromosome X nonPAR and XTR regions declared on different contigs: %s\n", line);
        self->x_rid = rid;
        if (self->x_xtr_beg != 0 || self->x_xtr_end != 0) error("Second chromosome X XTR rule: %s\n", line);
        self->x_xtr_beg = beg;
        self->x_xtr_end = end;
    } else if (strncmp(ss, "Y_nonpar", 8) == 0) {
        if ((self->y_xtr_beg != 0 || self->y_xtr_end != 0) && self->y_rid != rid)
            error("Chromosome Y XTR and nonPAR regions declared on different contigs: %s\n", line);
        self->y_rid = rid;
        if (self->y_nonpar_beg != 0 || self->y_nonpar_end != 0) error("Second chromosome Y nonPAR rule: %s\n", line);
        self->y_nonpar_beg = beg;
        self->y_nonpar_end = end;
    } else if (strncmp(ss, "Y_xtr", 5) == 0) {
        if ((self->y_nonpar_beg != 0 || self->y_nonpar_end != 0) && self->y_rid != rid)
            error("Chromosome Y nonPAR and XTR regions declared on different contigs: %s\n", line);
        self->y_rid = rid;
        if (self->y_xtr_beg != 0 || self->y_xtr_end != 0) error("Second chromosome Y XTR rule: %s\n", line);
        self->y_xtr_beg = beg;
        self->y_xtr_end = end;
    } else if (strncmp(ss, "mitochondria", 12) == 0) {
        if (self->mt_rid != -1) error("Second mitochondria rule %s\n", line);
        self->mt_rid = rid;
    }
}

/**
 *  genome_init_file() - initialize a genome structure from a file
 */
// split file into lines
// adapted from Petr Danecek's implementation of regidx_init_string() in bcftools/regidx.c
genome_rules_t *genome_init_file(const char *fname, const bcf_hdr_t *hdr) {
    if (!fname) return NULL;
    kstring_t tmp = {0, 0, NULL};
    htsFile *fp = hts_open(fname, "r");
    if (!fp) error("Failed to open %s: %s\n", fname, strerror(errno));
    genome_rules_t *self = genome_init(hdr);
    while (hts_getline(fp, KS_SEP_LINE, &tmp) >= 0) push_rule(self, tmp.s, hdr);
    free(tmp.s);
    hts_close(fp);
    return self;
}

// split string into lines
// adapted from Petr Danecek's implementation of regidx_init_string() in bcftools/regidx.c
static genome_rules_t *genome_init_string(const char *str, const bcf_hdr_t *hdr) {
    if (!str) return NULL;
    genome_rules_t *self = genome_init(hdr);
    kstring_t tmp = {0, 0, NULL};
    const char *ss = str;
    while (*ss) {
        while (*ss && isspace(*ss)) ss++;
        const char *se = ss;
        while (*se && *se != '\r' && *se != '\n') se++; // equivalent to KS_SEP_LINE
        tmp.l = 0;
        kputsn(ss, (int)(se - ss), &tmp);
        push_rule(self, tmp.s, hdr);
        while (*se && isspace(*se)) se++;
        ss = se;
    }
    free(tmp.s);
    return self;
}

/**
 *  genome_init_alias() - initialize a genome structure from an alias such as NCBI36, GRCh37, or
 * GRCh38
 */
// adapted from Petr Danecek's implementation of init_rules() in bcftools/regidx.c
genome_rules_t *genome_init_alias(FILE *restrict stream, char *alias, const bcf_hdr_t *hdr) {
    const rules_predef_t *rules = rules_predefs;

    int detailed = 0, len = strlen(alias);
    if (alias[len - 1] == '?') {
        detailed = 1;
        alias[len - 1] = '\0';
    }

    while (rules->alias && strcasecmp(alias, rules->alias)) rules++;

    if (!rules->alias) {
        fprintf(stream, "\nPRE-DEFINED REFERENCE GENOME RULES\n");
        fprintf(stream, "\b * Columns are: CHROM:BEG-END centromere/nonpar\n");
        fprintf(stream, " * Coordinates are 1-based inclusive.\n");
        for (rules = rules_predefs; rules->alias; rules++) {
            fprintf(stream, "\n%s\n   .. %s\n\n", rules->alias, rules->about);
            if (detailed) fprintf(stream, "%s\n", rules->rules);
        }
        if (!detailed) {
            fprintf(stream, "\nRun as --genome <assembly> (e.g. --genome GRCh37).\n");
            fprintf(stream,
                    "To see the detailed rules definition, append a question mark (e.g. "
                    "--genome GRCh37?).\n");
            fprintf(stream, "\n");
        }
        exit(1);
    } else if (detailed) {
        fprintf(stream, "\n%s\n   .. %s\n\n", rules->alias, rules->about);
        fprintf(stream, "%s", rules->rules);
        fprintf(stream, "\n");
        exit(1);
    }
    return genome_init_string(rules->rules, hdr);
}

/**
 *  readlist_short_arms() - initialize flag indicating which chromosome have short arms from a
 * comma separated list
 */
int readlist_short_arms(genome_rules_t *self, const char *str, const bcf_hdr_t *hdr) {
    int n;
    char **list = hts_readlist(str, 0, &n);
    if (!list) return 0;
    for (int i = 0; i < n; i++) {
        int rid = bcf_hdr_name2id(hdr, list[i]);
        free(list[i]);
        if (rid < 0) continue;
        self->is_short_arm[rid] = 1;
    }
    free(list);
    return 1;
}

#endif
