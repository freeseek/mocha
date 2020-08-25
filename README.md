![](mocha_logo.png)

A bcftools extension to call mosaic chromosomal alterations starting from phased VCF files with either B Allele Frequency (BAF) and Log R Ratio (LRR) or allelic depth (AD). If you use this tool in your publication, please cite the following papers from <a href="http://doi.org/10.1038/s41586-018-0321-x">2018</a> and <a href="http://doi.org/10.1038/s41586-020-2430-6">2020</a>:
```
Loh P., Genovese G., McCarroll S., Price A. et al. Insights about clonal expansions from 8,342 mosaic
chromosomal alterations. Nature 559, 350–355 (2018). [PMID: 29995854] [DOI: 10.1038/s41586-018-0321-x]

Loh P., Genovese G., McCarroll S., Monogenic and polygenic inheritance become
instruments for clonal selection (2020). [PMID: 32581363] [DOI: 10.1038/s41586-020-2430-6]
```
and this website. For any feedback or questions, contact the <a href="mailto:giulio.genovese@gmail.com">author</a>

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

NOTICE: Starting from July 2020 a <a href="wdl">WDL</a> pipeline is available to run the entire MoChA pipeline from raw intensity files to final calls

```
Usage:   bcftools +mocha [OPTIONS] <in.vcf>

Required options:
    -r, --rules <assembly>[?]      predefined genome reference rules, 'list' to print available settings, append '?' for details
    -R, --rules-file <file>        genome reference rules, space/tab-delimited CHROM:FROM-TO,TYPE

General Options:
    -x, --sex <file>               file including information about the gender of the samples
        --call-rate <file>         file including information about the call_rate of the samples
    -s, --samples [^]<list>        comma separated list of samples to include (or exclude with "^" prefix)
    -S, --samples-file [^]<file>   file of samples to include (or exclude with "^" prefix)
        --force-samples            only warn about unknown subset samples
    -v, --variants [^]<file>       tabix-indexed [compressed] VCF/BCF file containing variants
    -t, --targets [^]<region>      restrict to comma-separated list of regions. Exclude regions with "^" prefix
    -T, --targets-file [^]<file>   restrict to regions listed in a file. Exclude regions with "^" prefix
    -f, --apply-filters <list>     require at least one of the listed FILTER strings (e.g. "PASS,.")
                                   to include (or exclude with "^" prefix) in the analysis
    -p  --cnp <file>               list of regions to genotype in BED format
        --mhc <region>             MHC region to exclude from analysis (will be retained in the output)
        --kir <region>             KIR region to exclude from analysis (will be retained in the output)
        --threads <int>            number of extra output compression threads [0]

Output Options:
    -o, --output <file>            write output to a file [no output]
    -O, --output-type <b|u|z|v>    b: compressed BCF, u: uncompressed BCF, z: compressed VCF, v: uncompressed VCF [v]
        --no-version               do not append version and command line to the header
    -a  --no-annotations           omit Ldev and Bdev FORMAT from output VCF (requires --output)
        --no-log                   suppress progress report on standard error
    -l  --log <file>               write log to file [standard error]
    -m, --mosaic-calls <file>      write mosaic chromosomal alterations to a file [standard output]
    -g, --genome-stats <file>      write sample genome-wide statistics to a file [no output]
    -u, --ucsc-bed <file>          write UCSC bed track to a file [no output]

HMM Options:
        --bdev-LRR-BAF <list>      comma separated list of inverse BAF deviations for LRR+BAF model [-2.0,-4.0,-6.0,10.0,6.0,4.0]
        --bdev-BAF-phase <list>    comma separated list of inverse BAF deviations for BAF+phase model
                                   [6.0,8.0,10.0,15.0,20.0,30.0,50.0,80.0,100.0,150.0,200.0]
        --min-dist <int>           minimum base pair distance between consecutive sites for WGS data [400]
        --adjust-BAF-LRR <int>     minimum number of genotypes for a cluster to median adjust BAF and LRR (-1 for no adjustment) [5]
        --regress-BAF-LRR <int>    minimum number of genotypes for a cluster to regress BAF against LRR (-1 for no regression) [15]
        --LRR-GC-order <int>       order of polynomial to regress LRR against local GC content (-1 for no regression) [2]
        --xy-prob <float>          transition probability [1e-06]
        --err-prob <float>         uniform error probability [1e-02]
        --flip-prob <float>        phase flip probability [1e-02]
        --centromere-loss <float>  penalty to avoid calls spanning centromeres [1e-04]
        --telomere-gain <float>    telomere advantage to prioritize CN-LOHs [1e-02]
        --x-telomere-gain <float>  X telomere advantage to prioritize mLOX [1e-03]
        --y-telomere-gain <float>  Y telomere advantage to prioritize mLOY [1e-04]
        --short-arm-chrs <list>    list of chromosomes with short arms [13,14,15,21,22,chr13,chr14,chr15,chr21,chr22]
        --use-short-arms           use variants in short arms [FALSE]
        --use-centromeres          use variants in centromeres [FALSE]
        --use-no-rules-chrs        use chromosomes without centromere rules  [FALSE]
        --LRR-weight <float>       relative contribution from LRR for LRR+BAF  model [0.2]
        --LRR-hap2dip <float>      difference between LRR for haploid and diploid [0.45]
        --LRR-cutoff <float>       cutoff between LRR for haploid and diploid used to infer gender [estimated from X nonPAR]

Examples:
    bcftools +mocha -r GRCh37 input.bcf -v ^exclude.bcf -g stats.tsv -m mocha.tsv -p cnp.grch37.bed
    bcftools +mocha -r GRCh38 input.bcf -Ob -o output.bcf -g stats.tsv -m mocha.tsv --LRR-weight 0.5
```

Installation
============

Install basic tools (Debian/Ubuntu specific if you have admin privileges):
```
sudo apt install wget git g++ zlib1g-dev unzip samtools bedtools bcftools libgomp1
```

Optionally, you can install these libraries to activate further HTSlib features:
```
sudo apt install libbz2-dev libssl-dev liblzma-dev libgsl0-dev
```

Preparation steps
```
mkdir -p $HOME/bin $HOME/res && cd /tmp
```

We recommend compiling the source code but, wherever this is not possible, Linux x86_64 pre-compiled binaries are available for download <a href="http://software.broadinstitute.org/software/mocha">here</a>. However, notice that you will require a copy of BCFtools 1.10 or newer (available with Ubuntu 20.04)

Download latest version of <a href="https://github.com/samtools/htslib">HTSlib</a> and <a href="https://github.com/samtools/bcftools">BCFtools</a> (if not downloaded already)
```
git clone --branch=develop git://github.com/samtools/htslib.git
git clone --branch=develop git://github.com/samtools/bcftools.git
```

Download plugins code
```
/bin/rm -f bcftools/plugins/{{mocha,beta_binom,genome_rules}.h,{mocha,trio-phase,mochatools,extendFMT}.c}
wget -P bcftools/plugins https://raw.githubusercontent.com/freeseek/mocha/master/{{mocha,beta_binom,genome_rules}.h,{mocha,trio-phase,mochatools,extendFMT}.c}
```

Compile latest version of HTSlib (optionally disable bz2, gcs, and lzma) and BCFtools (make sure you are using gcc version 5 or newer)
```
cd htslib && autoheader && (autoconf || autoconf) && ./configure --disable-bz2 --disable-gcs --disable-lzma && make && cd ..
cd bcftools && make && cd ..
/bin/cp bcftools/{bcftools,plugins/{fill-tags,fixploidy,mocha,trio-phase,mochatools,extendFMT}.so} $HOME/bin/
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
cd $HOME/res
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
cd $HOME/res
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
Notice that this conversion takes ~26GB of RAM to complete and ~4 days of CPU time

Download cytoband file
```
wget -O $HOME/res/cytoBand.hg19.txt.gz http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz
```

Setup variables
```
ref="$HOME/res/human_g1k_v37.fasta"
mhc_reg="6:28542424-33448264"
kir_reg="19:54574747-55504099"
map="$HOME/res/genetic_map_hg19_withX.txt.gz"
kgp_pfx="$HOME/res/ALL.chr"
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
cd $HOME/res
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
cd $HOME/res
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
Notice that this conversion takes ~26GB of RAM to complete and ~4 days of CPU time

Download cytoband file
```
wget -O $HOME/res/cytoBand.hg38.txt.gz http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz
```

Setup variables
```
ref="$HOME/res/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
mhc_reg="chr6:28574647-33480487"
kir_reg="chr19:54071493-54992731"
map="$HOME/res/genetic_map_hg38_withX.txt.gz"
kgp_pfx="$HOME/res/ALL.chr"
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
crt="..." # file with call rate information (first column sample ID, second column call rate)
sex="..." # file with computed gender information (first column sample ID, second column gender: 1=male; 2=female)
xcl="..." # VCF file with additional list of variants to exclude (optional)
ped="..." # pedigree file to use if parent child duos are present
dir="..." # directory where output files will be generated
mkdir -p $dir
```

If you want to process <b>genotype array</b> data you need a VCF file with ALLELE_A, ALLELE_B, GC, GT, BAF, and LRR information:
```
##fileformat=VCFv4.2
##INFO=<ID=ALLELE_A,Number=1,Type=Integer,Description="A allele">
##INFO=<ID=ALLELE_B,Number=1,Type=Integer,Description="B allele">
##INFO=<ID=GC,Number=1,Type=Float,Description="GC ratio content around the variant">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=BAF,Number=1,Type=Float,Description="B Allele Frequency">
##FORMAT=<ID=LRR,Number=1,Type=Float,Description="Log R Ratio">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA12878
1	752566	rs3094315	G	A	.	.	ALLELE_A=1;ALLELE_B=0;GC=0.3675	GT:BAF:LRR	1|1:0.0111:-0.0798
1	776546	rs12124819	A	G	.	.	ALLELE_A=0;ALLELE_B=1;GC=0.435	GT:BAF:LRR	0|1:0.5441:0.4959
1	798959	rs11240777	G	A	.	.	ALLELE_A=1;ALLELE_B=0;GC=0.4075	GT:BAF:LRR	0|0:0.9712:0.2276
1	932457	rs1891910	G	A	.	.	ALLELE_A=1;ALLELE_B=0;GC=0.6425	GT:BAF:LRR	1|0:0.5460:-0.1653
```
Making sure that BAF refers to the allele frequency of the reference allele if ALLELE_B=0 and of the alternate allele if ALLELE_B=1

If you do not already have a VCF file but you have Illumina or Affymetrix genotype array data, you can use the <a href="https://github.com/freeseek/gtc2vcf">gtc2vcf</a> tools to convert the data to VCF and you can use the mochatools plugin to fill the ALLELE_A/ALLELE_B/GC info fields. Alternatively you can use your own scripts

Create a minimal binary VCF
```
bcftools annotate --no-version -Ob -o $dir/$pfx.unphased.bcf $vcf \
  -x ID,QUAL,^INFO/ALLELE_A,^INFO/ALLELE_B,^INFO/GC,^FMT/GT,^FMT/BAF,^FMT/LRR && \
  bcftools index -f $dir/$pfx.unphased.bcf
```

If you want to process <b>whole-genome sequence</b> data you need a VCF file with GC, GT and AD information:
```
##fileformat=VCFv4.2
##INFO=<ID=GC,Number=1,Type=Float,Description="GC ratio content around the variant">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA12878
1	752566	rs3094315	G	A	.	.	GC=0.3675	GT:AD	1|1:0,31
1	776546	rs12124819	A	G	.	.	GC=0.435	GT:AD	0|1:21,23
1	798959	rs11240777	G	A	.	.	GC=0.4075	GT:AD	0|0:31,0
1	932457	rs1891910	G	A	.	.	GC=0.6425	GT:AD	1|0:18,14
```
Make sure that AD is a "Number=R" format field (this was introduced in version <a href="https://samtools.github.io/hts-specs/VCFv4.2.pdf">4.2</a> of the VCF) or multi-allelic variants will not <a href="https://github.com/samtools/bcftools/issues/360">split properly</a>

Create a minimal binary VCF (notice that you will need a version newer than BCFtools 1.10.2 with implemented the <a href="https://github.com/samtools/bcftools/issues/360">--keep-sum</a> option)
```
bcftools view --no-version -h $vcf | sed 's/^\(##FORMAT=<ID=AD,Number=\)\./\1R/' | \
  bcftools reheader -h /dev/stdin $vcf | \
  bcftools filter --no-version -Ou -e "FMT/DP<10 | FMT/GQ<20" --set-GT . | \
  bcftools annotate --no-version -Ou -x ID,QUAL,^INFO/GC,^FMT/GT,^FMT/AD | \
  bcftools norm --no-version -Ou -m -any --keep-sum AD | \
  bcftools norm --no-version -Ob -o $dir/$pfx.unphased.bcf -f $ref && \
  bcftools index $dir/$pfx.unphased.bcf
```
This will set to missing all genotypes that have low coverage or low genotyping quality, as these can cause issues

Generate a list of variants that will be excluded from modeling by both eagle and mocha (notice that you will need a version newer than BCFtools 1.10.2 with implemented the F_MISSING option, or else you should drop that filter)
```
awk -F"\t" '$2<.97 {print $1}' $crt > samples_xcl_list.txt; \
echo '##INFO=<ID=JK,Number=1,Type=Float,Description="Jukes Cantor">' | \
  bcftools annotate --no-version -Ou -a $dup -c CHROM,FROM,TO,JK -h /dev/stdin -S ^samples_xcl_list.txt $dir/$pfx.unphased.bcf | \
  bcftools +fill-tags --no-version -Ou -t ^Y,MT,chrY,chrM -- -t ExcHet,F_MISSING | \
  bcftools +mochatools --no-version -Ou -- -x $sex -G | \
  bcftools annotate --no-version -Ob -o $dir/$pfx.xcl.bcf \
    -i 'FILTER!="." && FILTER!="PASS" || INFO/JK<.02 || INFO/ExcHet<1e-6 || INFO/F_MISSING>1-.97 ||
    INFO/AC_Sex_Test<1e-6 && CHROM!="X" && CHROM!="chrX" && CHROM!="Y" && CHROM!="chrY"' \
    -x ^INFO/JK,^INFO/ExcHet,^INFO/F_MISSING,^INFO/AC_Sex_Test && \
  bcftools index -f $dir/$pfx.xcl.bcf;
/bin/rm samples_xcl_list.txt
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
Eagle's <a href="https://data.broadinstitute.org/alkesgroup/Eagle/#x1-100003.2">memory requirements</a> will depend on the number of samples in the target (Nt) and in the reference panel (Nr=2504), and the number of variants (M) in the largest contig, and will amount to 1.5(Nt+Nr)M bytes. The <a href="https://data.broadinstitute.org/alkesgroup/Eagle/#x1-110003.3">running time</a> will be ~1 minute of CPU time per genome for <a href="https://www.nature.com/articles/ng.3679#Sec18">reference-based phasing</a> with a small target and reference panel (see Supplementary Tables 2,3) and ~5 minutes of CPU time per genome for <a href="https://www.nature.com/articles/ng.3679#Sec18">non-reference-based phasing</a> with a large cohort (see Supplementary Tables 7,8). Also, by default, if the option --pbwtIters is not used, Eagle will perform one phasing iteration if Nt<Nr/2=1252, two if 1252=Nr/2<Nt<2Nr=5008, and three if 5008=2Nr<Nt and in the second and third iterations both target and reference panel haplotypes will be used as references for phasing (see <a href="https://www.nature.com/articles/ng.3679#Sec10"here</a>).

Notice that you can also use alternative phasing methods that might be more effective, such as using <a href="http://www.haplotype-reference-consortium.org/">HRC</a> (use the Sanger Imputation Service, as the Michigan Imputations Server does not work with binary VCFs, does not work with VCFs with multiple chromosomes, does not work with chromosome X, and has no option for phasing without imputation). This might provide better phasing and therefore better ability to detect large events at lower cell fractions. Notice also that phasing can also be performed across overlapping windows rather than entire chromosomes to achieve better parallelization

Extract chromosomes that do not require phasing
```
bcftools view --no-version -Ob -o $dir/$pfx.other.bcf $dir/$pfx.unphased.bcf \
  -t ^$(seq -s, 1 22),X,$(seq -f chr%.0f -s, 1 22),chrX && \
bcftools index -f $dir/$pfx.other.bcf
```

Concatenate eagle output into a single VCF file and add GC/CpG content information
```
bcftools concat --no-version -Ob -o $dir/$pfx.bcf $dir/$pfx.{chr{{1..22},X},other}.bcf && \
bcftools index -f $dir/$pfx.bcf
```
Notice that if the phasing was made in overlapping windows rather than chromosomes, the overlapping windows should be concatenated using the --ligate option in bcftools concat

If pedigree information with duos or trios is available, you can improve the phased haplotypes from `eagle` by running the following command instead of the previous one:
```
bcftools concat --no-version -Ou $dir/$pfx.{chr{{1..22},X},other}.bcf | \
bcftools +trio-phase --no-version -Ob -o $dir/$pfx.bcf -- -p $ped && \
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
bcftools concat --no-version -Ob -o $dir/$pfx.dose.bcf $dir/$pfx.chr{{1..22},X}.dose.vcf.gz && \
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
sex="..." # file with computed gender information (first column sample ID, second column gender: 1=male; 2=female)
crt="..." # file with call rate information (first column sample ID, second column call rate)
lst="..." # file with list of samples to analyze for asymmetries (e.g. samples with 1p CN-LOH)
cnp="..." # file with list of regions to genotype in BED format
mhc_reg="..." # MHC region to skip
kir_reg="..." # KIR region to skip
```

Call mosaic chromosomal alterations with MoChA
```
bcftools +mocha \
  --rules $rule \
  --sex $sex \
  --call-rate $crt \
  --no-version \
  --output-type b \
  --output $dir/$pfx.mocha.bcf \
  --variants ^$dir/$pfx.xcl.bcf \
  --mosaic-calls $dir/$pfx.calls.tsv \
  --genome-stats $dir/$pfx.stats.tsv \
  --ucsc-bed $dir/$pfx.ucsc.bed \
  --cnp $cnp \
  --mhc $mhc_reg \
  --kir $kir_reg \
  $dir/$pfx.bcf && \
bcftools index -f $dir/$pfx.mocha.bcf
```
Notice that MoChA will read input computed gender and call rate if provided, otherwise these will be estimated from the VCF. For array data these statistics are usually available from the output of the Illumina\'s GenCall or Affymetrix\'s Axiom genotyping algorithms 

The genome statistics file contains information for each sample analyzed in the VCF and it includes the following columns:
```
            sample_id - sample ID
      computed_gender - estimated sample gender from X nonPAR region (not heterozygous sites count)
           XXX_median - median LRR or sequencing coverage across autosomes
               XXX_sd - standard deviation for LRR or sequencing coverage
             XXX_auto - auto correlation across consecutive sites for LRR or sequencing coverage (after GC correction)
          baf_sd/_cor - BAF standard deviation or beta-binomial overdispersion for read counts
             baf_conc - BAF phase concordance across phased heterozygous sites (see Vattathil et al. 2012)
             baf_auto - phased BAF auto correlation across consecutive phased heterozygous sites
              n_sites - number of sites across the genome for model based on LRR and BAF
               n_hets - number of heterozygous sites across the genome for model based on BAF and genotype phase
      x_nonpar_n_hets - number of heterozygous sites in the X nonPAR region
x_nonpar_baf_sd/_corr - BAF standard deviation or beta-binomial overdispersion for read counts in the X nonPAR region
  x_nonpar_XXX_median - median LRR or sequencing coverage over the X nonPAR region
  y_nonpar_XXX_median - median LRR or sequencing coverage over the Y nonPAR region
        mt_XXX_median - median LRR or sequencing coverage over the mitochondrial genome
       lrr_gc_rel_ess - LRR or sequencing coverage explained sum of squares fraction using local GC content
             lrr_gc_X - coefficient X for polynomial in GC content fitting LRR estimates
```

The mosaic calls file contains information about each mosaic and germline chromosomal alteration called and it includes the following columns:
```
      sample_id - sample ID
computed_gender - inferred sample gender
          chrom - chromosome
     beg_XXXXXX - beginning base pair position for the call (according to XXXXXX genome reference)
     end_XXXXXX - end base pair position for the call (according to XXXXXX genome reference)
         length - base pair length of the call
          p_arm - whether the call extends to the small arm and whether it reaches the telomere
          q_arm - whether the call extends to the long arm and whether it reaches the telomere
        n_sites - number of sites used for the call
         n_hets - number of heterozygous sites used for the call
       n50_gets - N50 value for consecutive heterozygous sites distances
           bdev - BAF deviation estimate from 0.5
        bdev_se - standard deviation estimate for BAF deviation
        rel_cov - relative coverage estimate from LRR or sequencing coverage
     rel_cov_se - standard deviation estimate for relative coverage
    lod_lrr_baf - LOD score for model based on LRR and BAF
  lod_baf_phase - LOD score for model based on BAF and genotype phase
        n_flips - number of phase flips for calls based on BAF and genotype phase model (-1 if LRR and BAF model used)
       baf_conc - BAF phase concordance across phased heterozygous sites underlying the call (see Vattathil et al. 2012)
   lod_baf_conc - LOD score for model based on BAF phase concordance (genome-wide corrected)
           type - Type of call based on LRR / relative coverage
             cf - estimated cell fraction based on BDEV and TYPE, or LDEV and TYPE if either BDEV or BDEV_SE are missing
```

The output VCF will contain the following extra FORMAT fields:
```
      Ldev - LRR deviation estimate
      Bdev - BAF deviation estimate
Bdev_Phase - for heterozygous calls: 1/-1 if the alternate allele is over/under represented
```

For array data, MoChA's memory requirements will depend on the number of samples (N) and the number of variants (M) in the largest contig and will amount to 9NM bytes. For example, if you are running 4,000 samples and chromosome 1 has ~80K variants, you will need approximately 2-3GB to run MoChA. It will take ~1/3 second of CPU time per genome to process samples genotyped on the Illumina GSA DNA microarray. For whole genome sequence data, MoChA's memory requirements will depend on the number of samples (N), the --min-dist parameter (D, 400 by default) and the length of the longest contig (L) and will amount to no more than 9NL/D, but could be significantly less, depending on how many variants you have in the VCF. If you are running 1,000 samples with default parameter --min-dist 400 and chromosome 1 is ~250Mbp long, you might need up to 5-6GB to run MoChA. For whole genome sequence data there is no need to batch too many samples together, as batching will not affect the calls made by MoChA (it will for array data unless you use options --adjust-BAF-LRR -1 and --regress-BAF-LRR -1). Notice that the CPU requirements for MoChA will be negligible compared to the CPU requirements for phasing with Eagle

Depending on your application, you might want to filter the calls from MoChA. For example, the following code:
```
awk -F "\t" 'NR==FNR && FNR==1 {for (i=1; i<=NF; i++) f[$i] = i}
  NR==FNR && FNR>1 {sample_id=$(f["sample_id"]); call_rate=$(f["call_rate"]); baf_auto=$(f["baf_auto"])}
  NR==FNR && FNR>1 && (call_rate<.97 || baf_auto>.03) {xcl[sample_id]++}
  NR>FNR && FNR==1 {for (i=1; i<=NF; i++) g[$i] = i; print}
  NR>FNR && FNR>1 {sample_id=$(g["sample_id"]); len=$(g["length"]); p_arm=$(g["p_arm"]); q_arm=$(g["q_arm"]);
    bdev=$(g["bdev"]); rel_cov=$(g["rel_cov"]); lod_baf_phase=$(g["lod_baf_phase"]); type=$(g["type"]);
    if (lod_baf_phase=="nan") lod_baf_phase=0}
  NR>FNR && FNR>1 && !(sample_id in xcl) && rel_cov>0.5 && type!~"^CNP" &&
    ( len>5e6 + 5e6 * (p_arm!="N" && q_arm!="N") ||
      len>5e5 && bdev<1/10 && rel_cov<2.5 && lod_baf_phase>10 ||
      rel_cov<2.1 && lod_baf_phase>10 )' $pfx.stats.tsv $pfx.calls.tsv > $pfx.calls.filtered.tsv

awk 'NR==FNR {x[$1"_"$3"_"$4"_"$5]++} NR>FNR && ($0~"^track" || $4"_"$1"_"$2"_"$3 in x)' \
  $pfx.calls.filtered.tsv $pfx.ucsc.bed > $pfx.ucsc.filtered.bed
```
will generate a new table after removing samples with `call_rate` lower than 0.97 `baf_auto` greater than 0.03, removing calls with less than a `lod_baf_phase` score of 10 unless they are larger than 5Mbp (or 10Mbp if they span the centromere) for the model based on BAF and genotype phase, removing calls flagged as germline copy number polymorphisms (CNPs), and removing calls that are likely germline duplications similarly to how it was done in the <a href="http://doi.org/10.1038/s41586-018-0321-x">UKBB</a>

Allelic imbalance pipeline
==========================

import results from MoChA into VCF file with imputed genotypes (optional for array data)
```
bcftools annotate \
  --no-version -Ob \
  --columns FMT/Ldev,FMT/Bdev,FMT/Bdev_Phase \
  $dir/$pfx.dose.bcf \
  --annotations $dir/$pfx.mocha.bcf \
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
  --regions $reg \
  --samples $lst \
  $dir/$pfx.dose.mocha.bcf | \
bcftools +mochatools \
  --no-version \
  --output-type b \
  --output $dir/$pfx.bal.bcf \
  -- --balance Bdev_Phase \
  --drop-genotypes && \
bcftools index -f $dir/$pfx.bal.bcf
```

Observe results for asymmetry analyses in table format
```
fmt="%CHROM\t%POS\t%ID\t%Bal{0}\t%Bal{1}\t%Bal_Test\n"
bcftools query \
  --include "Bal_Test>6" \
  --format "$fmt" \
  $dir/$pfx.bal.bcf | \
  column -ts $'\t'
```

Plot results
============

Install basic tools (Debian/Ubuntu specific if you have admin privileges):
```
sudo apt install r-cran-optparse r-cran-ggplot2 r-cran-data.table
```

Download R scripts
```
/bin/rm -f $HOME/bin/{summary,pileup,mocha}_plot.R
wget -P $HOME/bin https://raw.githubusercontent.com/freeseek/mocha/master/{summary,pileup,mocha}_plot.R
chmod a+x $HOME/bin/{summary,pileup,mocha}_plot.R
```

Generate summary plot
```
summary_plot.R \
  --pdf $dir/$pfx.pdf \
  --stats $dir/$pfx.stats.tsv \
  --calls $dir/$pfx.calls.tsv
```

Generate pileup plot
```
pileup_plot.R \
  --pdf $dir/$pfx.pdf \
  --stats $dir/$pfx.stats.tsv \
  --calls $dir/$pfx.calls.tsv \
  --cytoband $HOME/res/cytoBand.hg19.txt.gz
```

Plot mosaic chromosomal alterations (for array data)
```
mocha_plot.R \
  --mocha \
  --stats $dir/$pfx.stats.tsv \
  --vcf $dir/$pfx.mocha.bcf \
  --png MH0145622.png \
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
  --stats $dir/$pfx.stats.tsv
  --vcf $dir/$pfx.mocha.bcf \
  --png CSES15_P26_140611.png \
  --samples CSES15_P26_140611 \
  --regions 1:202236354-211793505 \
  --cytoband $HOME/res/cytoBand.hg19.txt.gz
```

![](CSES15_P26_140611.png)
Complex duplication overlapping the MDM4 gene (GRCh37 coordinates). Signal over heterozygous sites colored in blue shows evidence of a triplication event and signal over heterozygous sites colored in red shows evidence of a duplication event. Multiple phase switch errors can be observed in the shifted phased BAF signal

Acknowledgements
================

This work is supported by NIH grant <a href="http://grantome.com/grant/NIH/R01-HG006855">R01 HG006855</a>, NIH grant <a href="http://grantome.com/grant/NIH/R01-MH104964">R01 MH104964</a>, US Department of Defense Breast Cancer Research Breakthrough Award W81XWH-16-1-0316 (project BC151244), and the Stanley Center for Psychiatric Research. This work would have not been possible without the efforts of Heng Li <lh3@sanger.ac.uk>, Petr Danecek <pd3@sanger.ac.uk>, John Marshall <jm18@sanger.ac.uk>, James Bonfield <jkb@sanger.ac.uk>, and Shane McCarthy <sm15@sanger.ac.uk> in building HTSlib and BCFtools and Po-Ru Loh <poruloh@broadinstitute.org> in building the Eagle phasing software
