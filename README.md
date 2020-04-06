![](mocha_logo.png)

A bcftools extension to call mosaic chromosomal alterations starting from phased VCF files with either B Allele Frequency (BAF) and Log R Ratio (LRR) or allelic depth (AD). If you use this tool in your publication, please cite the following <a href="http://doi.org/10.1038/s41586-018-0321-x">paper</a> and <a href="http://doi.org/10.1101/653691">preprint</a>:
```
Loh P., Genovese G., McCarroll S., Price A. et al. Insights about clonal expansions from 8,342 mosaic
chromosomal alterations. Nature 559, 350–355 (2018). [PMID: 29995854] [DOI: 10.1038/s41586-018-0321-x]

Loh P., Genovese G., McCarroll S., Monogenic and polygenic inheritance become
instruments for clonal selection (2019). [DOI: 10.1101/653691]
```
and this website. For any feedback, send an email to giulio.genovese@gmail.com

<!--ts-->
   * [Usage](#usage)
   * [Installation](#installation)
   * [Download resources for GRCh37](#download-resources-for-grch37)
   * [Download resources for GRCh38](#download-resources-for-grch38)
   * [Data preparation](#data-preparation)
   * [Phasing pipeline](#phasing-pipeline)
   * [Chromosomal alterations pipeline](#chromosomal-alterations-pipeline)
   * [Allelic imbalance pipeline](#allelic-imbalance-pipeline)
   * [Plot results](#plot-results)
   * [Acknowledgements](#acknowledgements)
<!--te-->

Usage
=====

```
Usage:   bcftools mocha [OPTIONS] <in.vcf>

Required options:
    -r, --rules <assembly>[?]         predefined genome reference rules, 'list' to print available settings, append '?' for details
    -R, --rules-file <file>           genome reference rules, space/tab-delimited CHROM:FROM-TO,TYPE

General Options:
    -s, --samples [^]<list>           comma separated list of samples to include (or exclude with "^" prefix)
    -S, --samples-file [^]<file>      file of samples to include (or exclude with "^" prefix)
        --force-samples               only warn about unknown subset samples
    -t, --targets [^]<region>         restrict to comma-separated list of regions. Exclude regions with "^" prefix
    -T, --targets-file [^]<file>      restrict to regions listed in a file. Exclude regions with "^" prefix
    -f, --apply-filters <list>        require at least one of the listed FILTER strings (e.g. "PASS,.")
    -v, --variants [^]<file>          tabix-indexed [compressed] VCF/BCF file containing variants
                                      to include (or exclude with "^" prefix) in the analysis
        --threads <int>               number of extra output compression threads [0]

Output Options:
    -o, --output <file>               write output to a file [no output]
    -O, --output-type <b|u|z|v>       b: compressed BCF, u: uncompressed BCF, z: compressed VCF, v: uncompressed VCF [v]
        --no-version                  do not append version and command line to the header
    -a  --no-annotations              omit Ldev and Bdev FORMAT from output VCF (requires --output)
        --no-log                      suppress progress report on standard error
    -l  --log <file>                  write log to file [standard error]
    -m, --mosaic-calls <file>         write mosaic chromosomal alterations to a file [standard output]
    -g, --genome-stats <file>         write sample genome-wide statistics to a file [no output]
    -u, --ucsc-bed <file>             write UCSC bed track to a file [no output]

HMM Options:
    -p  --cnp <file>                  list of regions to genotype in BED format
        --bdev-LRR-BAF <list>         comma separated list of inverse BAF deviations for LRR+BAF model [-2.0,-4.0,-6.0,10.0,6.0,4.0]
        --bdev-BAF-phase <list>       comma separated list of inverse BAF deviations for BAF+phase model
                                      [6.0,8.0,10.0,15.0,20.0,30.0,50.0,80.0,100.0,150.0,200.0]
    -d, --min-dist <int>              minimum base pair distance between consecutive sites for WGS data [400]
        --no-BAF-flip                 do not correct BAF at flipped sites
        --median-BAF-adjust <int>     minimum number of heterozygous genotypes required to perform
                                      median BAF adjustment (-1 for no BAF adjustment) [5]
        --order-LRR-GC <int>          order of polynomial in local GC content to be used for polynomial
                                      regression of LRR (-1 for no LRR adjustment, 5 maximum) [2]
        --xy-prob <float>             transition probability [1e-09]
        --err-prob <float>            uniform error probability [1e-04]
        --flip-prob <float>           phase flip probability [1e-02]
        --telomere-advantage <float>  telomere advantage [1e-02]
        --centromere-penalty <float>  centromere penalty [1e-04]
        --short_arm_chrs <list>       list of chromosomes with short arms [13,14,15,21,22,chr13,chr14,chr15,chr21,chr22]
        --use_short_arms              use variants in short arms
        --use_centromeres             use variants in centromeres
        --LRR-cutoff <float>          LRR cutoff between haploid and diploid [estimated from X nonPAR]
        --LRR-hap2dip <float>         LRR difference between haploid and diploid [0.45]
        --LRR-auto2sex <float>        LRR difference between autosomes and diploid sex chromosomes [estimated from X nonPAR]
        --LRR-weight <float>          relative contribution from LRR for LRR+BAF model [0.2]

Examples:
    bcftools mocha -r GRCh37 input.bcf -v ^exclude.bcf -g stats.tsv -m mocha.tsv -p cnp.grch37.bed
    bcftools mocha -r GRCh38 input.bcf -Ob -o output.bcf -g stats.tsv -m mocha.tsv -c 1.0 --LRR-weight 0.5
```

Installation
============

Install basic tools (Debian/Ubuntu specific if you have admin privileges):
```
sudo apt install wget autoconf zlib1g-dev gzip unzip samtools bedtools
```

Optionally, you can install these libraries to activate further HTSlib features:
```
sudo apt install libbz2-dev libssl-dev liblzma-dev libgsl0-dev
```

Preparation steps
```
mkdir -p $HOME/bin $HOME/res $HOME/res/kgp && cd /tmp
```

Download latest version of <a href="https://github.com/samtools/htslib">HTSlib</a> and <a href="https://github.com/samtools/bcftools">BCFtools</a> (if not downloaded already)
```
git clone --branch=develop git://github.com/samtools/htslib.git
git clone --branch=develop git://github.com/samtools/bcftools.git
```

Add patches and code for plugin
```
/bin/rm -f bcftools/{{Makefile,main}.patch,vcfmocha.c,{beta_binom,genome_rules}.{c,h}} bcftools/plugins/{trio-phase,mochatools,importFMT,extendFMT}.c
wget -P bcftools https://raw.githubusercontent.com/freeseek/mocha/master/{{Makefile,main}.patch,vcfmocha.c,{beta_binom,genome_rules}.{c,h}}
wget -P bcftools/plugins https://raw.githubusercontent.com/freeseek/mocha/master/{trio-phase,mochatools,importFMT,extendFMT}.c
cd bcftools && patch < Makefile.patch && patch < main.patch && cd ..
```
If for any reason the patches fail with an error message, contact the <a href="mailto:giulio.genovese@gmail.com">author</a> for a fix

Compile latest version of HTSlib (optionally disable bz2, gcs, and lzma) and BCFtools (make sure you are using gcc version 5 or newer)
```
cd htslib && autoheader && (autoconf || autoconf) && ./configure --disable-bz2 --disable-gcs --disable-lzma && make && cd ..
cd bcftools && make && cd ..
/bin/cp bcftools/{bcftools,plugins/{fill-tags,fixploidy,split,trio-phase,mochatools,importFMT,extendFMT}.so} $HOME/bin/
```
Notice that you will need some functionalities missing from the base version of bcftools to run the pipeline

Make sure the directory with the plugins is available to bcftools
```
export PATH="$HOME/bin:$PATH"
export BCFTOOLS_PLUGINS="$HOME/bin"
```

Install latest version of <a href="https://data.broadinstitute.org/alkesgroup/Eagle/">Eagle</a>
```
wget -O $HOME/bin/eagle https://data.broadinstitute.org/alkesgroup/Eagle/downloads/dev/eagle_v2.4.1
chmod a+x $HOME/bin/eagle
```

Install Minimac3 (optional for array data)
```
git clone git://github.com/statgen/Minimac3.git
sed -i 's/USER_WARNINGS ?= -Werror/USER_WARNINGS ?= -Wno-format-truncation/' Minimac3/Library/libStatGenForMinimac3/general/Makefile
sed -i 's/bool legacy_count = 0/int legacy_count = 0/' Minimac3/Library/libStatGenForMinimac3/general/Parameters.cpp
sed -i 's/"\([0-9][0-9]*\)"/"\1","chr\1"/g;s/,"X"/,"X","chrX"/;s/finChromosome=="X"/(finChromosome=="X" || finChromosome=="chrX")/;s/finChromosome!="X"/(finChromosome!="X" \&\& finChromosome!="chrX")/' Minimac3/src/HaplotypeSet.cpp
sed -i 's/rHap.finChromosome!="X"/rHap.finChromosome!="X" \&\& rHap.finChromosome!="chrX"/' Minimac3/src/Imputation.cpp
cd Minimac3 && make && cd ..
/bin/cp Minimac3/bin/Minimac3{,-omp} $HOME/bin/
```

Download resources for GRCh37
=============================

Human genome reference
```
wget -O- ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz | \
  gzip -d > $HOME/res/human_g1k_v37.fasta
samtools faidx $HOME/res/human_g1k_v37.fasta
```

Genetic map
```
wget -P $HOME/res https://data.broadinstitute.org/alkesgroup/Eagle/downloads/tables/genetic_map_hg19_withX.txt.gz
```

1000 Genomes project phase 3
```
cd $HOME/res/kgp
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr{{1..22}.phase3_shapeit2_mvncall_integrated_v5a,X.phase3_shapeit2_mvncall_integrated_v1b,Y.phase3_integrated_v2a}.20130502.genotypes.vcf.gz{,.tbi}
for chr in {1..22} X Y; do
  bcftools view --no-version -Ou -c 2 ALL.chr${chr}.phase3*integrated_v[125][ab].20130502.genotypes.vcf.gz | \
  bcftools norm --no-version -Ou -m -any | \
  bcftools norm --no-version -Ob -o ALL.chr${chr}.phase3_integrated.20130502.genotypes.bcf -d none -f $HOME/res/human_g1k_v37.fasta && \
  bcftools index -f ALL.chr${chr}.phase3_integrated.20130502.genotypes.bcf
done
```

List of common germline duplications and deletions
```
wget -P $HOME/res ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/integrated_sv_map/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz{,.tbi}
bcftools query -i 'AC>1 && END-POS+1>10000 && SVTYPE!="INDEL" && (SVTYPE=="CNV" || SVTYPE=="DEL" || SVTYPE=="DUP")' \
  -f "%CHROM\t%POS0\t%END\t%SVTYPE\n" $HOME/res/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz > $HOME/res/cnp.grch37.bed
```

Minimal divergence intervals from segmental duplications (make sure your bedtools version is not affected by the groupby <a href="https://github.com/arq5x/bedtools2/issues/418">bug</a>)
```
wget -O- http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/genomicSuperDups.txt.gz | gzip -d |
  awk '!($2=="chrX" && $8=="chrY" || $2=="chrY" && $8=="chrX") {print $2"\t"$3"\t"$4"\t"$30}' > genomicSuperDups.bed

awk '{print $1,$2; print $1,$3}' genomicSuperDups.bed | \
  sort -k1,1 -k2,2n | uniq | \
  awk 'chrom==$1 {print chrom"\t"pos"\t"$2} {chrom=$1; pos=$2}' | \
  bedtools intersect -a genomicSuperDups.bed -b - | \
  bedtools sort | \
  bedtools groupby -c 4 -o min | \
  awk 'BEGIN {i=0; s[0]="+"; s[1]="-"} {if ($4!=x) i=(i+1)%2; x=$4; print $0"\t0\t"s[i]}' | \
  bedtools merge -s -c 4 -o distinct | \
  sed 's/^chr//' | grep -v gl | bgzip > $HOME/res/dup.grch37.bed.gz && \
  tabix -f -p bed $HOME/res/dup.grch37.bed.gz
```

1000 Genomes project phase 3 imputation panel for Minimac3 (optional for array data)
```
cd $HOME/res/kgp
for chr in {1..22} X; do
  if [ $chr == "X" ]; then
    bcftools +fixploidy --no-version ALL.chr$chr.phase3_integrated.20130502.genotypes.bcf -- -f 2 | \
      sed -e 's/\t0\/0/\t0|0/' -e 's/\t1\/1/\t1|1/' > tmp.chr$chr.GRCh37.vcf
  else
    bcftools view --no-version ALL.chr$chr.phase3_integrated.20130502.genotypes.bcf -o tmp.chr$chr.GRCh37.vcf
  fi
  Minimac3-omp \
    --refHaps tmp.chr$chr.GRCh37.vcf \
    --processReference \
    --myChromosome $chr \
    --prefix ALL.chr$chr.phase3_integrated.20130502.genotypes \
    --noPhoneHome
  /bin/rm tmp.chr$chr.GRCh37.vcf
done
```

Download cytoband file
```
wget -O $HOME/res/cytoBand.hg19.txt.gz http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz
```

Setup variables
```
ref="$HOME/res/human_g1k_v37.fasta"
map="$HOME/res/genetic_map_hg19_withX.txt.gz"
kgp_pfx="$HOME/res/kgp/ALL.chr"
kgp_sfx=".phase3_integrated.20130502.genotypes"
rule="GRCh37"
cnp="$HOME/res/cnp.grch37.bed"
dup="$HOME/res/dup.grch37.bed.gz"
cyto="$HOME/res/cytoBand.hg19.txt.gz"
```

Download resources for GRCh38
=============================

Human genome reference (following the suggestion from <a href="http://lh3.github.io/2017/11/13/which-human-reference-genome-to-use">Heng Li</a>)
```
wget -O- ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz | \
  gzip -d > $HOME/res/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
samtools faidx $HOME/res/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
```

Genetic map
```
wget -P $HOME/res https://data.broadinstitute.org/alkesgroup/Eagle/downloads/tables/genetic_map_hg38_withX.txt.gz
```

1000 Genomes project phase 3 (fixing contig names, removing duplicate variants, removing incomplete variants)
```
cd $HOME/res/kgp
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/ALL.chr{{1..22},X,Y}_GRCh38.genotypes.20170504.vcf.gz{,.tbi}
ref="$HOME/res/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
for chr in {1..22} X Y; do
  (bcftools view --no-version -h ALL.chr${chr}_GRCh38.genotypes.20170504.vcf.gz | \
    grep -v "^##contig=<ID=[GNh]" | sed 's/^##contig=<ID=MT/##contig=<ID=chrM/;s/^##contig=<ID=\([0-9XY]\)/##contig=<ID=chr\1/'; \
  bcftools annotate --no-version -x INFO/END ALL.chr${chr}_GRCh38.genotypes.20170504.vcf.gz | \
  bcftools view --no-version -H -c 2 | \
  grep -v "[0-9]|\.\|\.|[0-9]" | sed 's/^/chr/') | \
  bcftools norm --no-version -Ou -m -any | \
  bcftools norm --no-version -Ob -o ALL.chr${chr}_GRCh38.genotypes.20170504.bcf -d none -f $ref && \
  bcftools index -f ALL.chr${chr}_GRCh38.genotypes.20170504.bcf
done
```
Do notice though that the 1000 Genomes project team incorrectly lifted over chromosome X genotypes over PAR1 and PAR2 regions, despite this issue not being reported in the <a href="http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/README_GRCh38_liftover_20170504.txt">README</a>. This will directly affect the ability to detect chromosome Y loss events

List of common germline duplications and deletions
```
wget -P $HOME/res ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/integrated_sv_map/supporting/GRCh38_positions/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.GRCh38.vcf.gz{,.tbi}
bcftools query -i 'AC>1 && END-POS+1>10000 && SVTYPE!="INDEL" && (SVTYPE=="CNV" || SVTYPE=="DEL" || SVTYPE=="DUP")' \
  -f "chr%CHROM\t%POS0\t%END\t%SVTYPE\n" $HOME/res/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.GRCh38.vcf.gz > $HOME/res/cnp.grch38.bed
```

Minimal divergence intervals from segmental duplications (make sure your bedtools version is not affected by the groupby <a href="https://github.com/arq5x/bedtools2/issues/418">bug</a>)
```
wget -O- http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/genomicSuperDups.txt.gz | gzip -d |
  awk '!($2=="chrX" && $8=="chrY" || $2=="chrY" && $8=="chrX") {print $2"\t"$3"\t"$4"\t"$30}' > genomicSuperDups.bed

awk '{print $1,$2; print $1,$3}' genomicSuperDups.bed | \
  sort -k1,1 -k2,2n | uniq | \
  awk 'chrom==$1 {print chrom"\t"pos"\t"$2} {chrom=$1; pos=$2}' | \
  bedtools intersect -a genomicSuperDups.bed -b - | \
  bedtools sort | \
  bedtools groupby -c 4 -o min | \
  awk 'BEGIN {i=0; s[0]="+"; s[1]="-"} {if ($4!=x) i=(i+1)%2; x=$4; print $0"\t0\t"s[i]}' | \
  bedtools merge -s -c 4 -o distinct | \
  grep -v "GL\|KI" | bgzip > $HOME/res/dup.grch38.bed.gz && \
  tabix -f -p bed $HOME/res/dup.grch38.bed.gz
```

1000 Genomes project phase 3 imputation panel for Minimac3 (optional for array data)
```
cd $HOME/res/kgp
for chr in {1..22} X; do
  if [ $chr == "X" ]; then
    bcftools +fixploidy --no-version ALL.chr${chr}_GRCh38.genotypes.20170504.bcf -- -f 2 | \
      sed -e 's/\t0\/0/\t0|0/' -e 's/\t1\/1/\t1|1/' > tmp.chr$chr.GRCh38.vcf
  else
    bcftools view --no-version ALL.chr${chr}_GRCh38.genotypes.20170504.bcf -o tmp.chr$chr.GRCh38.vcf
  fi
  Minimac3-omp \
    --refHaps tmp.chr$chr.GRCh38.vcf \
    --processReference \
    --myChromosome chr$chr \
    --prefix ALL.chr${chr}_GRCh38.genotypes.20170504 \
    --noPhoneHome
  /bin/rm tmp.chr$chr.GRCh38.vcf
done
```

Download cytoband file
```
wget -O $HOME/res/cytoBand.hg38.txt.gz http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz
```

Setup variables
```
ref="$HOME/res/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
map="$HOME/res/genetic_map_hg38_withX.txt.gz"
kgp_pfx="$HOME/res/kgp/ALL.chr"
kgp_sfx="_GRCh38.genotypes.20170504"
rule="GRCh38"
cnp="$HOME/res/cnp.grch38.bed"
dup="$HOME/res/dup.grch38.bed.gz"
cyto="$HOME/res/cytoBand.hg38.txt.gz"
```

Data preparation
================

Preparation steps
```
vcf="..." # input VCF file with phased GT, LRR, and BAF
pfx="..." # output prefix
thr="..." # number of threads to use
sex="..." # file with sex information (first column sample ID, second column sex: 1=male; 2=female)
xcl="..." # VCF file with additional list of variants to exclude (optional)
ped="..." # pedigree file to use if parent child duos are present
dir="..." # directory where output files will be generated
mkdir -p $dir
```

If you want to process <b>genotype array</b> data you need a VCF file with ALLELE_A, ALLELE_B, GT, BAF, and LRR information:
```
##fileformat=VCFv4.2
##INFO=<ID=ALLELE_A,Number=1,Type=Integer,Description="A allele">
##INFO=<ID=ALLELE_B,Number=1,Type=Integer,Description="B allele">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=BAF,Number=1,Type=Float,Description="B Allele Frequency">
##FORMAT=<ID=LRR,Number=1,Type=Float,Description="Log R Ratio">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA12878
1	752566	rs3094315	G	A	.	.	ALLELE_A=1;ALLELE_B=0	GT:BAF:LRR	1|1:0.0111:-0.0798
1	776546	rs12124819	A	G	.	.	ALLELE_A=0;ALLELE_B=1	GT:BAF:LRR	0|1:0.5441:0.4959
1	798959	rs11240777	G	A	.	.	ALLELE_A=1;ALLELE_B=0	GT:BAF:LRR	0|0:0.9712:0.2276
1	932457	rs1891910	G	A	.	.	ALLELE_A=1;ALLELE_B=0	GT:BAF:LRR	1|0:0.5460:-0.1653
```
Making sure that BAF refers to the allele frequency of the reference allele if ALLELE_B=0 and of the alternate allele if ALLELE_B=1

If you do not already have a VCF file but you have Illumina or Affymetrix genotype array data, you can use the <a href="https://github.com/freeseek/gtc2vcf">gtc2vcf</a> tools to convert the data to VCF. Alternatively you can use your own scripts

Create a minimal binary VCF
```
bcftools annotate --no-version -Ob -o $dir/$pfx.unphased.bcf $vcf \
  -x ID,QUAL,^INFO/ALLELE_A,^INFO/ALLELE_B,^FMT/GT,^FMT/BAF,^FMT/LRR && \
  bcftools index -f $dir/$pfx.unphased.bcf
```

If you want to process <b>whole-genome sequence</b> data you need a VCF file with GT and AD information:
```
##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA12878
1	752566	rs3094315	G	A	.	.	.	GT:AD	1|1:0,31
1	776546	rs12124819	A	G	.	.	.	GT:AD	0|1:21,23
1	798959	rs11240777	G	A	.	.	.	GT:AD	0|0:31,0
1	932457	rs1891910	G	A	.	.	.	GT:AD	1|0:18,14
```
Make sure that AD is a "Number=R" format field (this was introduced in version <a href="https://samtools.github.io/hts-specs/VCFv4.2.pdf">4.2</a> of the VCF) or multi-allelic variants will not <a href="https://github.com/samtools/bcftools/issues/360">split properly</a>

Create a minimal binary VCF (notice that you will need a version newer than BCFtools 1.10.2 with implemented the <a href="https://github.com/samtools/bcftools/issues/360">--keep-sum</a> option)
```
bcftools view --no-version -h $vcf | sed 's/^\(##FORMAT=<ID=AD,Number=\)\./\1R/' | \
  bcftools reheader -h /dev/stdin $vcf | \
  bcftools filter --no-version -Ou -e "FMT/DP<10 | FMT/GQ<20" --set-GT . | \
  bcftools annotate --no-version -Ou -x ID,QUAL,INFO,^FMT/GT,^FMT/AD | \
  bcftools norm --no-version -Ou -m -any --keep-sum AD | \
  bcftools norm --no-version -Ob -o $dir/$pfx.unphased.bcf -f $ref && \
  bcftools index $dir/$pfx.unphased.bcf
```
This will set to missing all genotypes that have low coverage or low genotyping quality, as these can cause issues

Perform basic quality control (the generated list of variants will be excluded from modeling by both eagle and mocha)
```
n=$(bcftools query -l $dir/$pfx.unphased.bcf|wc -l); \
ns=$((n*98/100)); \
echo '##INFO=<ID=JK,Number=1,Type=Float,Description="Jukes Cantor">' | \
  bcftools annotate --no-version -Ou -a $dup -c CHROM,FROM,TO,JK -h /dev/stdin $dir/$pfx.unphased.bcf | \
  bcftools +fill-tags --no-version -Ou -- -t NS,ExcHet | \
  bcftools +mochatools --no-version -Ou -- -x $sex -G | \
  bcftools annotate --no-version -Ob -o $dir/$pfx.xcl.bcf \
    -i 'FILTER!="." && FILTER!="PASS" || JK<.02 || NS<'$ns' || ExcHet<1e-6 || AC_Sex_Test>6' \
    -x FILTER,^INFO/JK,^INFO/NS,^INFO/ExcHet,^INFO/AC_Sex_Test && \
  bcftools index -f $dir/$pfx.xcl.bcf
```
This command will create a list of variants falling within segmental duplications with low divergence (<2%), high levels of missingness (>2%), variants with excess heterozygosity (p<1e-6), and variants that correlate with sex in an unexpected way (p<1e-6). If you are using WGS data and you don't have a file with sex information, you can skip the quality control line using this information. When later running MoChA, sex will be imputed and a sex file can be computed from MoChA's output

If a file with additional variants to be excluded is available, further merge it with the generated list
```
/bin/mv $dir/$pfx.xcl.bcf $dir/$pfx.xcl.tmp.bcf && \
/bin/mv $dir/$pfx.xcl.bcf.csi $dir/$pfx.xcl.tmp.bcf.csi && \
bcftools merge --no-version -Ob -o $dir/$pfx.xcl.bcf -m none $dir/$pfx.xcl.tmp.bcf $xcl && \
bcftools index -f $dir/$pfx.xcl.bcf
```

Phasing pipeline
================

Phase VCF file by chromosome with Eagle
```
for chr in {1..22} X; do
  eagle \
    --geneticMapFile $map \
    --outPrefix $dir/$pfx.chr$chr \
    --numThreads $thr \
    --vcfRef $kgp_pfx${chr}$kgp_sfx.bcf \
    --vcfTarget $dir/$pfx.unphased.bcf \
    --vcfOutFormat b \
    --noImpMissing \
    --outputUnphased \
    --vcfExclude $dir/$pfx.xcl.bcf \
    --chrom $chr \
    --pbwtIters 3 && \
  bcftools index -f $dir/$pfx.chr$chr.bcf
done
```
Notice that you can also use alternative phasing methods that might be more effective, such as using <a href="http://www.haplotype-reference-consortium.org/">HRC</a> (use the Sanger Imputation Service, as the Michigan Imputations Server does not work with binary VCFs, does not work with VCFs with multiple chromosomes, does not work with chromosome X, and has no option for phasing without imputation). This might provide better phasing and therefore better ability to detect large events at lower cell fractions

Extract chromosomes that do not require phasing
```
bcftools view --no-version -Ob -o $dir/$pfx.other.bcf $dir/$pfx.unphased.bcf \
  -t ^$(seq -s, 1 22),X,$(seq -f chr%.0f -s, 1 22),chrX && \
bcftools index -f $dir/$pfx.other.bcf
```

Concatenate eagle output into a single VCF file and add GC/CpG content information
```
bcftools concat --no-version -Ou --threads $thr $dir/$pfx.{chr{{1..22},X},other}.bcf | \
bcftools +mochatools --no-version -Ob -o $dir/$pfx.bcf -- -f $ref && \
bcftools index -f $dir/$pfx.bcf
```

If pedigree information with duos or trios is available, you can improve the phased haplotypes from `eagle` by running the following command instead of the previous one:
```
bcftools concat --no-version -Ou --threads $thr $dir/$pfx.{chr{{1..22},X},other}.bcf | \
bcftools +mochatools --no-version -Ou -- -f $ref | \
bcftools +trio-phase --no-version -Ob -o $dir/$pfx.bcf --threads $thr -- -p $ped && \
bcftools index -f $dir/$pfx.bcf
```
(it requires a ped file)

Impute variants using Minimac3 (optional for array data)
```
for chr in {1..22} X; do
  bcftools annotate --no-version -Ou $dir/$pfx.bcf -a $dir/$pfx.xcl.bcf -m +XCL -r $chr,chr$chr | \
    bcftools view --no-version -Ov -o $dir/$pfx.chr$chr.vcf -e "XCL==1"
  Minimac3-omp \
    --cpus $thr \
    --refHaps $kgp_pfx${chr}$kgp_sfx.m3vcf.gz \
    --format GT,GP \
    --haps $dir/$pfx.chr$chr.vcf \
    --prefix $dir/$pfx.chr$chr \
    --lowMemory \
    --noPhoneHome
  bcftools query -l $dir/$pfx.bcf | \
    bcftools view --no-version -Ob -o $dir/$pfx.chr$chr.dose.bcf -S /dev/stdin $dir/$pfx.chr$chr.dose.vcf.gz
  /bin/rm $dir/$pfx.chr$chr.vcf $dir/$pfx.chr$chr.dose.vcf.gz
done
```

Concatenate imputed genotypes into a single VCF file (optional for array data)
```
bcftools concat --no-version -Ob -o $dir/$pfx.dose.bcf --threads $thr $dir/$pfx.chr{{1..22},X}.dose.vcf.gz && \
bcftools index -f $dir/$pfx.dose.bcf
```

Remove unphased VCF and single chromosome files (optional)
```
/bin/rm $dir/$pfx.{unphased,chr{{1..22},X},other}.bcf{,.csi} $dir/$pfx.chr{{1..22},X}.dose.bcf
```

Chromosomal alterations pipeline
================================

Preparation steps
```
pfx="..." # output prefix
thr="..." # number of extra threads to use
lst="..." # file with list of samples to analyze for asymmetries (e.g. samples with 1p CN-LOH)
```

Call mosaic chromosomal alterations with MoChA
```
bcftools mocha \
  --rules $rule \
  --no-version \
  --output-type b \
  --output $dir/$pfx.mocha.bcf \
  --threads $thr \
  --variants ^$dir/$pfx.xcl.bcf \
  --mosaic-calls $dir/$pfx.mocha.tsv \
  --genome-stats $dir/$pfx.stats.tsv \
  --ucsc-bed $dir/$pfx.ucsc.bed \
  --cnp $cnp \
  $dir/$pfx.bcf && \
bcftools index -f $dir/$pfx.mocha.bcf
```

The genome statistics file contains information for each sample analyzed in the VCF and it includes the following columns:
```
               SAMPLE - sample ID
           XXX_MEDIAN - median LRR or sequencing coverage across autosomes
               XXX_SD - standard deviation for LRR or sequencing coverage
             XXX_AUTO - auto correlation across consecutive sites for LRR or sequencing coverage (after GC correction)
         BAF_SD/_CORR - BAF standard deviation or beta-binomial overdispersion for read counts
             BAF_CONC - BAF phase concordance across phased heterozygous sites (see Vattathil et al. 2012)
             BAF_AUTO - phased BAF auto correlation across consecutive phased heterozygous sites
               NSITES - number of sites across the genome for model based on LRR and BAF
                NHETS - number of heterozygous sites across the genome for model based on BAF and genotype phase
       X_NONPAR_NHETS - number of heterozygous sites in the X nonPAR region
X_NONPAR_BAF_SD/_CORR - BAF standard deviation or beta-binomial overdispersion for read counts in the X nonPAR region
  X_NONPAR_XXX_MEDIAN - median LRR or sequencing coverage over the X nonPAR region
  Y_NONPAR_XXX_MEDIAN - median LRR or sequencing coverage over the Y nonPAR region
        MT_XXX_MEDIAN - median LRR or sequencing coverage over the mitochondrial genome
                  SEX - estimated sample sex from X nonPAR region (not heterozygous sites count)
              REL_ESS - LRR or sequencing coverage explained sum of squares fraction using local GC content
```

The mosaic calls file contains information about each mosaic and germline chromosomal alteration called and it includes the following columns:
```
       SAMPLE - sample ID
          SEX - inferred sample sex
        CHROM - chromosome
   BEG_XXXXXX - beginning base pair position for the call (according to XXXXXX genome reference)
   END_XXXXXX - end base pair position for the call (according to XXXXXX genome reference)
       LENGTH - base pair length of the call
        P_ARM - whether the call extends to the small arm and whether it reaches the telomere
        Q_ARM - whether the call extends to the long arm and whether it reaches the telomere
       NSITES - number of sites used for the call
        NHETS - number of heterozygous sites used for the call
     N50_HETS - N50 value for consecutive heterozygous sites distances
         BDEV - BAF deviation estimate from 0.5
      BDEV_SE - standard deviation estimate for BAF deviation
      REL_COV - relative coverage estimate from LRR or sequencing coverage
   REL_COV_SE - standard deviation estimate for relative coverage
  LOD_LRR_BAF - LOD score for model based on LRR and BAF
LOD_BAF_PHASE - LOD score for model based on BAF and genotype phase
       NFLIPS - number of phase flips for calls based on BAF and genotype phase model (-1 if LRR and BAF model used)
     BAF_CONC - BAF phase concordance across phased heterozygous sites underlying the call (see Vattathil et al. 2012)
 LOD_BAF_CONC - LOD score for model based on BAF phase concordance (genome-wide corrected)
         TYPE - Type of call based on LRR / relative coverage
           CF - estimated cell fraction based on BDEV and TYPE, or LDEV and TYPE if either BDEV or BDEV_SE are missing
```

The output VCF will contain the following extra FORMAT fields:
```
      Ldev - LRR deviation estimate
      Bdev - BAF deviation estimate
Bdev_Phase - for heterozygous calls: 1/-1 if the alternate allele is over/under represented
```

For array data, MoChA's memory requirements will depend on the number of samples (N) and the number of variants (M) in the largest contig and will amount to 9NM bytes. For example, if you are running 4,000 samples and chromosome 1 has ~80K variants, you will need approximately 2-3GB to run MoChA. For whole genome sequence data, MoChA's memory requirements will depend on the number of samples (N), the --min-dist parameter (D, 400 by default) and the length of the longest contig (L) and will amount to no more than 9NL/D, but could be significantly less, depending on how many variants you have in the VCF. If you are running 1,000 samples with default parameter --min-dist 400 and chromosome 1 is ~250Mbp long, you might need up to 5-6GB to run MoChA. Notice that for whole genome sequence data there is no need to batch too many samples together, as batching will not affect the calls made by MoChA (it will for array data unless you use option --median-BAF-adjust -1)

Notice that, depending on your application, you might want to filter the calls from MoChA. For example, the following code:
```
awk 'NR==FNR && FNR>1 && $6>.51 {x[$1]++}
  NR>FNR && (FNR==1 || !($1 in x) && $6>1e5 && $17>10 && $21!~"CNP" && $22<.5) {print}' \
  $pfx.stats.tsv $pfx.mocha.tsv > $pfx.mocha.filter.tsv
```
will generate a new table after removing samples with BAF_CONC greater than 0.51, removing calls smaller than 100kbp, removing calls with less than a LOD score of 10 for the model based on BAF and genotype phase, removing calls flagged as germline copy number polymorphisms (CNPs), and removing calls with an estimated cell fraction larger than 50%

Allelic imbalance pipeline
==========================

import results from MoChA into VCF file with imputed genotypes (optional for array data)
```
bcftools +importFMT \
  --no-version -Ob \
  --formats Ldev,Bdev,Bdev_Phase \
  $dir/$pfx.dose.bcf \
  $dir/$pfx.mocha.bcf \
  -o $dir/$pfx.dose.mocha.bcf && \
bcftools index $dir/$pfx.dose.mocha.bcf
```

Run asymmetry analyses (subset cohort, run binomial test, discard genotypes)
```
bcftools +extendFMT \
  --no-version -Ou \
  --format Bdev_Phase \
  --phase \
  --dist 500000 \
  -r $reg \
  -s $lst \
  $dir/$pfx.dose.mocha.bcf | \
bcftools +mochatools \
  --no-version -Ob \
  -o $dir/$pfx.bal.bcf \
  -- -b Bdev_Phase -G && \
bcftools index -f $dir/$pfx.bal.bcf
```

Observe results for asymmetry analyses in table format
```
fmt="%CHROM\t%POS\t%ID\t%Bal{0}\t%Bal{1}\t%Bal_Test\n"
bcftools query \
  -i "Bal_Test>6" \
  -f "$fmt" \
  $dir/$pfx.bal.bcf | \
  column -ts $'\t'
```

Split output VCF file by samples (optional)
```
bcftools annotate --no-version -Ou -x INFO $dir/$pfx.mocha.bcf | \
  bcftools +split -Ob -o $dir/
```

Plot results
============

Install basic tools (Debian/Ubuntu specific if you have admin privileges):
```
sudo apt install r-cran-ggplot2 r-cran-data.table
```

Download R scripts
```
/bin/rm -f $HOME/bin/{plot_summary,mocha_plot}.R
wget -P $HOME/bin https://raw.githubusercontent.com/freeseek/mocha/master/{plot_summary,mocha_plot}.R
chmod a+x $HOME/bin/{plot_summary,mocha_plot}.R
```

Generate summary plot
```
plot_summary.R \
  --pdf $dir/$pfx.pdf \
  --stats $dir/$pfx.stats.tsv \
  --calls $dir/$pfx.mocha.tsv
```

Plot mosaic chromosomal alterations (for array data)
```
mocha_plot.R \
  --mocha \
  --png MH0145622.png \
  --vcf $dir/$pfx.mocha.bcf \
  --samples MH0145622 \
  --regions 11:81098129-115077367 \
  --cytoband $HOME/res/cytoBand.hg19.txt.gz
```
Notice that by default MoChA will perform internal BAF (for array data) and LRR adjustments that will not be performed by the plotter so the visualized data might actually look noisier than what was actually processed

![](MH0145622.png)
Mosaic deletion from array data overlapping the ATM gene (GRCh37 coordinates). The deletion signal can be observed across LRR, BAF and phased BAF, although it is the most clear with the latter. Furthermore, evidence of three phase switch errors can be observed in the shifted phased BAF signal

Plot mosaic chromosomal alterations (for WGS data)
```
mocha_plot.R \
  --wgs \
  --mocha \
  --png CSES15_P26_140611.png \
  --vcf $dir/$pfx.mocha.bcf \
  --samples CSES15_P26_140611 \
  --regions 1:202236354-211793505 \
  --cytoband $HOME/res/cytoBand.hg19.txt.gz
```

![](CSES15_P26_140611.png)
Complex duplication overlapping the MDM4 gene (GRCh37 coordinates). Signal over heterozygous sites colored in blue shows evidence of a triplication event and signal over heterozygous sites colored in red shows evidence of a duplication event. Multiple phase switch errors can be observed in the shifted phased BAF signal

Acknowledgements
================

This work is supported by NIH grant <a href="https://projectreporter.nih.gov/project_info_description.cfm?aid=8852155">R01 HG006855</a> and the Stanley Center for Psychiatric Research and by US Department of Defense Breast Cancer Research Breakthrough Award W81XWH-16-1-0316 (project BC151244). This work would have not been possible without the efforts of Heng Li <lh3@sanger.ac.uk>, Petr Danecek <pd3@sanger.ac.uk>, John Marshall <jm18@sanger.ac.uk>, James Bonfield <jkb@sanger.ac.uk>, and Shane McCarthy <sm15@sanger.ac.uk> in building HTSlib and BCFtools and Po-Ru Loh <poruloh@broadinstitute.org> in building the Eagle phasing software
