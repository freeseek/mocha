ukb2txt
=======

In this tutorial we will show you how to reconstruct raw Affymetrix files for each of the 106 UK biobank batches for use with the MoChA WDL pipeline

<!--ts-->
  * [Introduction](#introduction)
  * [Download Resources](#download-resources)
  * [Build Auxiliary Tools](#build-auxiliary-tools)
  * [Split by Batches](#split-by-batches)
  * [Input Files for MoChA WDL](#input-files-for-mocha-wdl)
<!--te-->

Introduction
============

While most biobanks based on Affymetrix Axiom DNA microarrays, such as the FINNGEN and MVP biobanks, provide raw Affymetrix files split by batches, such as
```
###.AxiomGT1.report.txt
###.AxiomGT1.snp-posteriors.txt
###.AxiomGT1.calls.txt
###.AxiomGT1.confidences.txt
###.AxiomGT1.summary.txt
```
the UK biobank delivers genotypes and intensities (see <a href="https://biobank.ctsu.ox.ac.uk/crystal/label.cgi?id=263">here</a>) in non-standard file formats merged across all batches (see <a href="https://biobank.ctsu.ox.ac.uk/crystal/refer.cgi?id=531">here</a>). Furthermore, while the two Affymetrix arrays used to genotype UK biobank participants contain probesets for two recurrently somatic variants
```
AX-83208650 Affx-37797994 rs77375493 JAK2 V617F
AX-86708948 Affx-80252081 rs121913502 IDH2 R140Q
```
the JAK2 V617F marker is not present in the output files as it was removed during quality control and so the intensities from the probeset have not been made available to researchers. A <a href="http://doi.org/10.1182/blood-2015-06-652941">study</a> has shown intensities from this probeset to be a proxy for JAK2 V617F carriers. The intensities for this probeset are instead available in the FINNGEN and MVP biobanks

Finally, the provided PLINK files are encoded in such a way that the A1 and A2 alleles are, respectively, the reference and alternate alleles with respect to the GRCh37 human genome reference (unless neither of the two alleles matches the reference allele) rather than being the A and B alleles with respect to the Affymetrix probesets design. As MoChA requires knowledge of which allele is the B allele to interpret the B-allele frequency correctly, we need to recover the B-allele designation from the Affymetrix array manifest files

One advantage of reverting the UK biobank intensities and genotype files to the raw Affymetrix files split by batches (with the exception of the markers that have been removed by quality control procedures) is that the MoChA WDL pipeline can process these directly and convert the data to VCF against GRCh38 if needed

If you want to run this conversion through Cromwell, you can skip the sections below and run the <a href="ukb2txt.wdl">ukb2txt.wdl</a> with a JSON input file like the following
```json
{
  "ukb2txt.sqc_path": "gs://fc-9a7c5487-04c9-4182-b3ec-13de7f6b409b/imputed",
  "ukb2txt.cal_path": "gs://fc-9a7c5487-04c9-4182-b3ec-13de7f6b409b/genotype",
  "ukb2txt.int_path": "gs://fc-9a7c5487-04c9-4182-b3ec-13de7f6b409b/genotype"
}
```
Where `sqc_path`, `cal_path`, and `int_path` are, respectively, the paths where the private resources `ukb_sqc_v2.txt`, `ukb_cal_chr{{1..22},X,Y,XY,MT}_v2.bed`, and `ukb_int_chr{{1..22},X,Y,XY,MT}_v2.bin` are localized

To minimize network operations between storage nodes and compute nodes, chromosome and batch files for genotypes and intensities are processed on separate nodes, with the intensity file for chromosome 1 being the largest file localized to a single node (231 GiB). We advise to make sure you run the computation in the same computing location where the anonymized data files are stored to optimize latency and network bandwidth and avoid <a href="https://cloud.google.com/storage/pricing">egress costs</a>, as genotype and intensity files comprise 3.0 TiB of data. If this recommendation is followed, then running this conversion pipeline should take approximately half a day and cost $5-10 (mostly due to gzip compression). The output will consist of 2.4 TiB of gzip-compressed files split by batches

Download Resources
==================

A compute node with 7 TiB of free disk space would be able to download and process the UK biobank intensity files following the steps described below. However, if you are planning to use a cloud virtual machine, do notice that <a href="https://cloud.google.com/compute/disks-image-pricing#disk">disk storage costs</a> are higher than <a href="https://cloud.google.com/storage/pricing">bucket storage costs</a>. If you have access to the genotype intensity binary files in the cloud, you could request a virtual machine in the same cloud location as this can significantly speed the time needed to copy the intensity files into the node

You will first need the following anonymized but private resources (see <a href="https://biobank.ctsu.ox.ac.uk/crystal/refer.cgi?id=668">here</a> for instructions on how to get these files)
```
ukb_sqc_v2.txt # 252 MiB
ukb_cal_chr{{1..22},X,Y,XY,MT}_v2.bed # 92 GiB
ukb_int_chr{{1..22},X,Y,XY,MT}_v2.bin # 2.9 TiB
```

Download public resources
```
wget -nd --no-check-certificate biobank.ctsu.ox.ac.uk/crystal/crystal/auxdata/ukb_snp_posterior_chrX_haploid.bim # 1.1 MiB
wget -nd --no-check-certificate biobank.ctsu.ox.ac.uk/crystal/crystal/auxdata/ukb_snp_posterior.tar # 11 GiB
wget -nd --no-check-certificate biobank.ctsu.ox.ac.uk/crystal/crystal/auxdata/ukb_snp_bim.tar # 23 MiB
wget -nd --no-check-certificate biobank.ctsu.ox.ac.uk/crystal/crystal/auxdata/ukb_snp_posterior.batch # 1.2 KiB
wget -nd --no-check-certificate biobank.ctsu.ox.ac.uk/crystal/crystal/auxdata/ukb_snp_qc.txt # 355 MiB
wget -nd --no-check-certificate biobank.ctsu.ox.ac.uk/crystal/crystal/docs/ukb_genetic_data_description.txt # 17 KiB
wget -nd --no-check-certificate biobank.ctsu.ox.ac.uk/crystal/crystal/docs/Array_BIL_34.zip # 148 MiB
wget -nd --no-check-certificate biobank.ctsu.ox.ac.uk/crystal/crystal/docs/Array_UKB_34.zip # 149 MiB
```

For some reasons, the length of the flanking sequence field in the manifest files is capped at 250 characters. This causes the flanking sequences of three Affymetrix indels (Affx-89015252, Affx-92046163, Affx-92047500) to be truncated. The following commands can fix this issue
```
zcat Array_BIL_34.zip | \
  sed -e 's/ACTGAAGCATGACTCTTCTAAAACAAGTTAATAAT\[-\/GCATTCTACAGAAAACCTTTTAGACATTACAACCTCTAACCTAAAAAATGGAACTCTAAATACCTTTTCACCATCTCAGTAGTAATTTCTGACACATCACTTTCAGAAGGCTCATAAGCCAACTCAGATACAACCGTAATATTAATGAGGAAAGCAGTTAGTGTCAAAGGAAGAGGATTGAAAGTAGAAGTGTTATGAAACTCCTCTTTTAC/ACTGAAGCATGACTCTTCTAAAACAAGTTAATAAT\[-\/GCATTCTACAGAAAACCTTTTAGACATTACAACCTCTAACCTAAAAAATGGAACTCTAAATACCTTTTCACCATCTCAGTAGTAATTTCTGACACATCACTTTCAGAAGGCTCATAAGCCAACTCAGATACAACCGTAATATTAATGAGGAAAGCAGTTAGTGTCAAAGGAAGAGGATTGAAAGTAGAAGTGTTATGAAACTCCTCTTTTACTTTTTTGCATGCTATCAAAATCACTGCCAATAATATGAATGTCAGGTCAATTTCCATAGGTAAATCCGTTACCTTTTACCTCTTTAAAAGAAAAGTTTTGCAGAAGAGGGCTGAAAATTTCTTGAGCCATTTCAGCACAAAGAATGGAAGTTCATTTCTCACCATGATACAACTCTACCCTGCTGTCATCTTCATTGTGATGGTGGCAGAAGTTTAGCAGGGTGCAAGTGACCACTA\]AATGACATCTTTTCATGAACAATTGATAAATCTTT/' | \
  gzip > Axiom_UKBiLEVE.na34.annot.csv.gz
zcat Array_UKB_34.zip | \
  sed -e 's/ACTGAAGCATGACTCTTCTAAAACAAGTTAATAAT\[-\/GCATTCTACAGAAAACCTTTTAGACATTACAACCTCTAACCTAAAAAATGGAACTCTAAATACCTTTTCACCATCTCAGTAGTAATTTCTGACACATCACTTTCAGAAGGCTCATAAGCCAACTCAGATACAACCGTAATATTAATGAGGAAAGCAGTTAGTGTCAAAGGAAGAGGATTGAAAGTAGAAGTGTTATGAAACTCCTCTTTTAC/ACTGAAGCATGACTCTTCTAAAACAAGTTAATAAT\[-\/GCATTCTACAGAAAACCTTTTAGACATTACAACCTCTAACCTAAAAAATGGAACTCTAAATACCTTTTCACCATCTCAGTAGTAATTTCTGACACATCACTTTCAGAAGGCTCATAAGCCAACTCAGATACAACCGTAATATTAATGAGGAAAGCAGTTAGTGTCAAAGGAAGAGGATTGAAAGTAGAAGTGTTATGAAACTCCTCTTTTACTTTTTTGCATGCTATCAAAATCACTGCCAATAATATGAATGTCAGGTCAATTTCCATAGGTAAATCCGTTACCTTTTACCTCTTTAAAAGAAAAGTTTTGCAGAAGAGGGCTGAAAATTTCTTGAGCCATTTCAGCACAAAGAATGGAAGTTCATTTCTCACCATGATACAACTCTACCCTGCTGTCATCTTCATTGTGATGGTGGCAGAAGTTTAGCAGGGTGCAAGTGACCACTA\]AATGACATCTTTTCATGAACAATTGATAAATCTTT/' \
      -e 's/TATTGGCTCAATTTCTTGGGGAGGGGGTGCTGTCA\[-\/GAGATTGTTATCTGAGGATGTGACATAGATTTCTCAGGGCACAATTTCAACTACTTTTTCAGCTTTAGGGTTTTTAGATACGTTTGTACCACAATTGAGCATGGGAGGGAGAGGGGTGAGCCTAAGCAGTGATGGCTGATTTCTGTCATGTCTGTCATGTGTCCCCCAGTACCTCCAGAGGTAACTGTGCTCACAAACAGCCCTGTGGAACT/TATTGGCTCAATTTCTTGGGGAGGGGGTGCTGTCA\[-\/GAGATTGTTATCTGAGGATGTGACATAGATTTCTCAGGGCACAATTTCAACTACTTTTTCAGCTTTAGGGTTTTTAGATACGTTTGTACCACAATTGAGCATGGGAGGGAGAGGGGTGAGCCTAAGCAGTGATGGCTGATTTCTGTCATGTCTGTCATGTGTCCCCCAGTACCTCCAGAGGTAACTGTGCTCACAAACAGCCCTGTGGAACTGAGAGAGCCCAACGTCCTCATCTGTTTCATA\]GACAAGTTCACCCCACCAGTGGTCAATGTCACGTG/' \
      -e 's/GTAATAGGACTCCTTGAAGAGACTTTCGGCGGGGC\[-\/CTGGGGAAGACAGAGGGAGACACACTGAACAGTCTGGCTGTGACCTCCAGGGCCACTGTCACCCACACAGCTACCAGGTGGCCAGAGCCAGGATCTGAACCCAGGTCTGTGGGGGATCCACACCTGAATCCCATTCTTGGGGGAGTCTCATTGGCACCACAGCAGAGGAACCTCTAACCTAGGCCTTCGTTCAAGACTAGAACCTGCCCCCA/GTAATAGGACTCCTTGAAGAGACTTTCGGCGGGGC\[-\/CTGGGGAAGACAGAGGGAGACACACTGAACAGTCTGGCTGTGACCTCCAGGGCCACTGTCACCCACACAGCTACCAGGTGGCCAGAGCCAGGATCTGAACCCAGGTCTGTGGGGGATCCACACCTGAATCCCATTCTTGGGGGAGTCTCATTGGCACCACAGCAGAGGAACCTCTAACCTAGGCCTTCGTTCAAGACTAGAACCTGCCCCCACTCCAGTGCTGACCACTTCAGAGCAGAGGGGAGGCTGAAGAGGACACAGGGTCCTCAGTGTCCCAATGCCAGATCCCCACTCTCCTTGGTCACCAGCTTGTGAATCTGGGCAGTCGCCTGGCTCCTGCCTACTGTCCTGAGCCATGTTTCAGAGGGCAGGTAACAAATGAGAAGGGAAAAGTACAGCTCTAGTTCGGGGGGTGGGAGGCCGCTCTATCCTTTACTCTGAAGGCCTGGGGGAGGCTGACCTCCAGACCTGCAGCTGCCAGAAAACCCTGGGGCCCATCCACTGCTTAC\]CCATGGGGTCTGAGGAGTCAGTGATGATCACGTCG/' | \
  gzip > Axiom_UKB_WCSG.na34.annot.csv.gz
```

While Affymetrix orders marker\'s alleles according to an internal designation of alleles as A and B, the UK biobank has reordered the alleles for each marker in the genotype and SNP posterior (but not intensity) files as reference and alternate with respect to GRCh37. We will therefore recover information about which markers have A and B alleles swapped when compared to reference and alternate alleles and use this information to recover the original genotypes. The same list of markers were swapped for SNP posterior files, with the exception of ten Affymetrix indels (AX-82920347, AX-82999350, AX-83057578, AX-83106285, AX-83149193, AX-83166262, AX-83197070, AX-83253475, AX-83575267, AX-83587806) which, for unclear reasons, were inconsistently swapped between genotype files and SNP posterior files and incorrectly so in the SNP posterior files

Build Auxiliary Tools
=====================

Make sure the GNU C compiler and C Library are available to compile three auxiliary tools (Debian/Ubuntu specific if you have admin privileges)
```
sudo apt install --no-install-recommends gcc libc6-dev 1>&2
```

Generate a small ELF binary file that we will use to efficiently unpack PLINK binary files
```
echo '#include <stdio.h>
#include <stdlib.h>
#include <string.h>
int main(int argc, char **argv) {
  if (argc != 3)
    return 1;
  int i, n, nread, swap;
  if (sscanf(argv[1], "%d", &n) != 1)
    return 1;
  char *in_buffer = (char *)malloc((n + 3) / 4 * sizeof(char));
  char *out_buffer = (char *)malloc((n + 3) / 4 * 4 * sizeof(char));
  FILE *in = fopen(argv[2], "r");
  while ((nread = fread(in_buffer, 1, (n + 3) / 4, stdin))) {
    if (nread != (n + 3) / 4)
      return 1;
    if (fscanf(in, "%d\n", &swap) != 1)
      return 1;
    for (i = 0; i < (n + 3) / 4; i++) {
      char x = in_buffer[i];
      char y = (swap ? (~x & 0xAA) : (x & 0x55) << 1) ^
               ((x & 0x55) ^ ((x & 0xAA) >> 1));
      int z = (int)(y & 0x03) ^ ((int)(y & 0x0C) << 6) ^
              ((int)(y & 0x30) << 12) ^ ((int)(y & 0xC0) << 18) ^ 0x30303030;
      memcpy(&out_buffer[4 * i], &z, 4);
    }
    fwrite(out_buffer, 1, n, stdout);
  }
  free(in_buffer);
  free(out_buffer);
  return 0;
}' > unpack.c
gcc -O3 -o unpack unpack.c
```

Generate a small ELF binary file that we will use to efficiently split binary resources by batches
```
echo '#include <stdio.h>
#include <stdlib.h>
#include <string.h>
int main(int argc, char **argv) {
  if (argc != 3)
    return 1;
  int i, nread, n = 0, m = 4;
  int *offsets = (int *)malloc(m * sizeof(int));
  FILE **batches = (FILE **)malloc(m * sizeof(FILE *));
  FILE *in = fopen(argv[2], "r");
  char buffer[BUFSIZ];
  offsets[0] = 0;
  while (fscanf(in, "%d %s\n", &i, buffer) == 2) {
    sprintf(buffer + strlen(buffer), "%s", argv[1]);
    batches[n] = fopen(buffer, "w");
    n++;
    if (n == m) {
      m <<= 1;
      offsets = (int *)realloc(offsets, m * sizeof(int));
      batches = (FILE **)realloc(batches, m * sizeof(FILE *));
    }
    offsets[n] = offsets[n - 1] + i;
  }
  fclose(in);
  char *line = (char *)malloc(offsets[n] * sizeof(char *));
  while ((nread = fread(line, 1, offsets[n], stdin))) {
    if (nread != offsets[n])
      return 1;
    for (i = 0; i < n; i++)
      fwrite(line + offsets[i], 1, offsets[i + 1] - offsets[i], batches[i]);
  }
  for (i = 0; i < n; i++)
    fclose(batches[i]);
  free(offsets);
  free(batches);
  free(line);
  return 0;
}' > split.c
gcc -O3 -o split split.c
```

Generate a small ELF binary file that we will use to efficiently convert binary floats
```
echo '#include <stdio.h>
#include <stdlib.h>
#include <string.h>
int main(int argc, char **argv) {
  if (argc != 3)
    return 1;
  int i, j, n, m, nread;
  if (sscanf(argv[1], "%d", &n) != 1)
    return 1;
  if (sscanf(argv[2], "%d", &m) != 1)
    return 1;
  char *buffer = (char *)malloc(n * m * sizeof(float));
  while ((nread = fread(buffer, sizeof(float), n * m, stdin))) {
    if (nread != n * m)
      return 1;
    for (i = 0; i < n; i++)
      for (j = 0; j < m; j++)
        fprintf(stdout, "%g%c", *(float *)&buffer[4 * (n * j + i)],
                j + 1 == m ? 0x0A : 0x09);
  }
  free(buffer);
  return 0;
}' > dump.c
gcc -O3 -o dump dump.c
```

Split by Batches
================

Split sample report file into batch report files (~2 seconds)
```
sed 's/ UKBiLEVEAX_b\([1-9]\) / UKBiLEVEAX_b0\1 /' ukb_sqc_v2.txt | \
  awk 'NR==FNR {print "cel_files\tcomputed_gender\tcall_rate\tcn-probe-chrXY-ratio_gender_meanX\tcn-probe-chrXY-ratio_gender_meanY" > $1".AxiomGT1.report.txt"}
  NR>FNR {printf "%s\t%s\t%s\t%s\t%s\n",$1,$11,$7,$12,$13 > $4".AxiomGT1.report.txt"}' ukb_snp_posterior.batch -
sed 's/$/.AxiomGT1.report.txt/' ukb_snp_posterior.batch | tr '\n' '\0' | xargs -0 gzip --force --no-name
```

Split chromosome posterior files into batch posterior files (~1 minute)
```
tar xvf ukb_snp_posterior.tar 1>&2
cat ukb_snp_posterior_chr{{1..22},X,Y,XY,MT}.bin | ./split .snp-posteriors.bin <(sed 's/^/132 /' ukb_snp_posterior.batch)
```

Reconstruct Affymetrix SNP posterior files (~30 seconds per batch)
```
awk 'NR==FNR {x[$2]++} NR>FNR && FNR>1 {array=$9==0 || $9==2; print array"\t"$3;
  if ($1 in x) print array"\t"$3":1"}' ukb_snp_posterior_chrX_haploid.bim ukb_snp_qc.txt > snp-posteriors.UKBL.tsv
awk 'NR==FNR {x[$2]++} NR>FNR && FNR>1 {array=$9==1 || $9==2; print array"\t"$3;
  if ($1 in x) print array"\t"$3":1"}' ukb_snp_posterior_chrX_haploid.bim ukb_snp_qc.txt > snp-posteriors.UKBB.tsv
sed 's/ UKBiLEVEAX_b\([1-9]\) / UKBiLEVEAX_b0\1 /' ukb_sqc_v2.txt | \
  cut -d" " -f3,4 | uniq | while read array batch; do
  (echo -e "#%SnpPosteriorFormatVer=1\n#%data-order=meanX,varX,nObsMean,nObsVar,meanY,varY,covarXY\nid\tBB\tAB\tAA\tCV";
  cat $batch.snp-posteriors.bin | ./dump 1 33 | paste snp-posteriors.$array.tsv - | sed '/^0/d;s/^..//' | \
  awk '{fmt="%s\t%s,%s,%s,%s,%s,%s,%s\t%s,%s,%s,%s,%s,%s,%s\t%s,%s,%s,%s,%s,%s,%s\t%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n"
  if ($2<$16) printf fmt,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,$29,$30,$31,$32,$33,$34;
         else printf fmt,$1,$16,$17,$18,$19,$20,$21,$22,$9,$10,$11,$12,$13,$14,$15,$2,$3,$4,$5,$6,$7,$8,$23,$24,$25,$26,$27,$28,$29,$30,$31,$32,$33,$34}') | \
    gzip > $batch.AxiomGT1.snp-posteriors.txt.gz && /bin/rm $batch.snp-posteriors.bin
done
```

Split chromosome genotype calls files into batch genotype calls files (~30 minutes, could be batched by chromosomes)
```
(zcat Array_BIL_34.zip | grep -v ^# | tail -n+2;
 zcat Array_UKB_34.zip | grep -v ^# | tail -n+2) | tr -d '"' | \
  awk -F, '{if ($12==$14 && $13==$15) x=0; else if ($12==$15 && $13==$14) x=1; else x=0; print $1,x}' | \
  awk 'NR==FNR {x[$1]=$2} NR>FNR && FNR>1 {print x[$3]}' - ukb_snp_qc.txt > swap.lines
tail -qc+4 ukb_cal_chr{{1..22},X,Y,XY,MT}_v2.bed | ./unpack $(cat ukb_sqc_v2.txt | wc -l) swap.lines | \
  ./split .calls.bin <(sed 's/ UKBiLEVEAX_b\([1-9]\) / UKBiLEVEAX_b0\1 /' ukb_sqc_v2.txt | cut -d" " -f4 | uniq -c | sed 's/^ *//')
```

Reconstruct Affymetrix genotype calls files (~10 minutes per batch)
```
awk 'NR>1 {array=$9==0 || $9==2; print array"\t"$3}' ukb_snp_qc.txt > calls.UKBL.tsv
awk 'NR>1 {array=$9==1 || $9==2; print array"\t"$3}' ukb_snp_qc.txt > calls.UKBB.tsv
sed 's/ UKBiLEVEAX_b\([1-9]\) / UKBiLEVEAX_b0\1 /' ukb_sqc_v2.txt | \
  cut -d" " -f3,4 | uniq -c | while read n array batch; do
  (echo -en "#Calls: -1=NN, 0=AA, 1=AB, 2=BB\nprobeset_id\t";
  sed 's/ UKBiLEVEAX_b\([1-9]\) / UKBiLEVEAX_b0\1 /' ukb_sqc_v2.txt | \
    awk -v batch=$batch '$4==batch {print $1}' | tr '\n' '\t' | sed 's/\t$/\n/';
  fold -w$n $batch.calls.bin | tr '\n' '\r' | fold -w1 | tr '\n\r' '\t\n' | sed 's/3/-1/g' | \
  paste calls.$array.tsv - | sed '/^0/d;s/^..//') | \
    gzip > $batch.AxiomGT1.calls.txt.gz && /bin/rm $batch.calls.bin
done
```

Split chromosome intensity files into batch intensity files (~4 hours, could be batched by chromosomes)
```
cat ukb_int_chr{{1..22},X,Y,XY,MT}_v2.bin | ./split .summary.bin <(sed 's/ UKBiLEVEAX_b\([1-9]\) / UKBiLEVEAX_b0\1 /' ukb_sqc_v2.txt | cut -d" " -f4 | uniq -c | awk '{print 8*$1,$2}')
```

Reconstruct Affymetrix intensities summary files (~3.5 hours per batch)
```
awk 'NR>1 {array=$9==0 || $9==2; print array"\t"$3"-A"; print array"\t"$3"-B"}' ukb_snp_qc.txt > summary.UKBL.tsv
awk 'NR>1 {array=$9==1 || $9==2; print array"\t"$3"-A"; print array"\t"$3"-B"}' ukb_snp_qc.txt > summary.UKBB.tsv
sed 's/ UKBiLEVEAX_b\([1-9]\) / UKBiLEVEAX_b0\1 /' ukb_sqc_v2.txt | \
  cut -d" " -f3,4 | uniq -c | while read n array batch; do
  (echo -en "probeset_id\t";
  sed 's/ UKBiLEVEAX_b\([1-9]\) / UKBiLEVEAX_b0\1 /' ukb_sqc_v2.txt | \
    awk -v batch=$batch '$4==batch {print $1}' | tr '\n' '\t' | sed 's/\t$/\n/';
  cat $batch.summary.bin | ./dump 2 $n | paste summary.$array.tsv - | sed '/^0/d;s/^..//') | \
    gzip > $batch.AxiomGT1.summary.txt.gz && /bin/rm $batch.summary.bin
done
```

Input Files for MoChA WDL
=========================

Create sample table
```
sed 's/ UKBiLEVEAX_b\([1-9]\) / UKBiLEVEAX_b0\1 /' ukb_sqc_v2.txt | \
  awk 'BEGIN {print "sample_id\tbatch_id\tcel\tcomputed_gender\tcall_rate"}
  {printf "%s\t%s\t%s\t%s\t%s\n",$1,$4,$1,$11,$7/100}' > ukb.sample.tsv
```
If you wanted to use the sample IDs associated with your application ID, you should instead include in the first column of the table the IDs from the `ukbA_cal_v2_sP.fam` file which includes the IDs associated to your UK Biobank application ID. Remember, when handling retracted samples, that each ID in the first column needs to be unique for the MoChA WDL pipeline to run

Create batch table
```
sed 's/ UKBiLEVEAX_b\([1-9]\) / UKBiLEVEAX_b0\1 /' ukb_sqc_v2.txt | \
  cut -d" " -f3,4 | uniq -c | \
  awk 'BEGIN {csv["UKBB"]="Axiom_UKB_WCSG.na34.annot.csv.gz"; csv["UKBL"]="Axiom_UKBiLEVE.na34.annot.csv.gz"
  print "batch_id\tn_smpls\tcsv\tsnp\treport\tcalls\tsummary"}
  {printf "%s\t%s\t%s\t%s.AxiomGT1.snp-posteriors.txt.gz\t%s.AxiomGT1.report.txt.gz\t%s.AxiomGT1.calls.txt.gz\t%s.AxiomGT1.summary.txt.gz\n",
  $3,$1,csv[$2],$3,$3,$3,$3}' > ukb.batch.tsv
```

For the MoChA WDL to access these files in Google cloud, you can copy tables `ukb.sample.tsv` and `ukb.batch.tsv` in a Google bucket in locations such as `gs://{google-bucket}/tsvs/`, the manifest files `Axiom_UKBiLEVE.na34.annot.csv.gz` and `Axiom_UKB_WCSG.na34.annot.csv.gz` in a location such as `gs://{google-bucket}/manifests/`, and all the gzipped files created in the previous section in a location such as `gs://{google-bucket}/txts/`

Create JSON input file
```json
{
  "mocha.sample_set_id": "ukb",
  "mocha.mode": "txt",
  "mocha.target": "calls",
  "mocha.realign": true,
  "mocha.max_win_size_cm": 10.0,
  "mocha.overlap_size_cm": 2.0,
  "mocha.ref_name": "GRCh38",
  "mocha.ref_path": "gs://{google-bucket}/GRCh38",
  "mocha.manifest_path": "gs://{google-bucket}/manifests",
  "mocha.data_path": "gs://{google-bucket}/txts",
  "mocha.sample_tsv_file": "gs://{google-bucket}/tsvs/ukb.sample.tsv",
  "mocha.batch_tsv_file": "gs://{google-bucket}/tsvs/ukb.batch.tsv",
  "mocha.docker_repository": "us.gcr.io/mccarroll-mocha"
}
```

Additionally, you can also create an additional table for principal components
```
(echo sample_id pc{01..40}; cut -d" " -f1,26-65 ukb_sqc_v2.txt) | tr ' ' '\t' > ukb.pcs.tsv
```
