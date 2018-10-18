![](mocha_logo.png)

A bcftools extension to call mosaic chromosomal alterations starting from phased VCF files with either B Allele Frequency (BAF) and Log R Ratio (LRR) or allelic depth (AD). If you use this tool in your publication, please cite the following <a href="http://doi.org/10.1038/s41586-018-0321-x">paper</a>:
```
Loh P., Genovese G., McCarroll S., Price A. et al. Insights about clonal expansions from 8,342 mosaic
chromosomal alterations. Nature 559, 350–355 (2018). [PMID: 29995854] [DOI: 10.1038/s41586-018-0321-x]
```
and this website. For any feedback, send an email to giulio.genovese@gmail.com

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
    -v, --variants [^]<file>          tabix-indexed [compressed] VCF/BCF file containing variants
                                      to include (or exclude with "^" prefix) in the analysis
        --threads <int>               number of extra output compression threads [0]

Output Options:
    -o, --output <file>               write output to a file [no output]
    -O, --output-type <b|u|z|v>       b: compressed BCF, u: uncompressed BCF, z: compressed VCF, v: uncompressed VCF [v]
        --no-version                  do not append version and command line to the header
    -a  --no-annotations              omit Ldev and Bdev FORMAT from output VCF (requires --output)
    -m, --mosaic-calls <file>         write mosaic chromosomal alterations to a file [standard output]
    -g, --genome-stats <file>         write sample genome-wide statistics to a file [no output]
    -u, --bed-ucsc <file>             write UCSC bed track to a file [no output]
    -l  --no-log                      suppress progress report on standard error

HMM Options:
    -x, --xy-prob <float>             transition probability [1e-09]
    -z, --err-prob <float>            uniform error probability [1e-04]
    -f, --flip-prob <float>           phase flip probability [1e-02]
    -t, --telomere-advantage <float>  telomere advantage [1e-02]
    -c, --centromere-penalty <float>  centromere penalty [1e-04]
    -p  --cnp <file>                  list of regions to genotype in BED format
    -n, --cnf <list>                  comma separated list of copy number fractions for LRR+BAF model
                                      [1.0,1.5,3.0,4.0]
    -b, --bdev <list>                 comma separated list of inverse BAF deviations for BAF+phase model
                                      [6.0,8.0,10.0,15.0,20.0,30.0,50.0,80.0,100.0,150.0,200.0]
    -d, --min-dist <int>              minimum base pair distance between consecutive sites for WGS data [400]
        --LRR-hap2dip <float>         LRR difference between haploid and diploid [estimated from X nonPAR]
        --LRR-auto2sex <float>        LRR difference between autosomes and diploid sex chromosomes [estimated from X nonPAR]
        --LRR-weight <float>          relative contribution from LRR for LRR+BAF model [0.2]
        --median-BAF-adjust <int>     minimum number of heterozygous genotypes required to perform
                                      median BAF adjustment (-1 for no BAF adjustment) [5]
        --order-LRR-GC <int>          order of polynomial in local GC content to be used for polynomial
                                      regression of LRR (-1 for no LRR adjustment, 5 maximum) [2]

Examples:
    bcftools mocha -r GRCh37 input.bcf -v ^exclude.bcf -g stats.tsv -m mocha.tsv -p cnp.grch37.bed
    bcftools mocha -r GRCh38 input.bcf -Ob -o output.bcf -g stats.tsv -m mocha.tsv -n 1.0 --LRR-weight 0.5
```

Installation
============

Install basic tools (Debian/Ubuntu specific):
```
sudo apt install wget gzip unzip samtools bedtools
```

Optionally, you can install these libraries to activate further bcftools features:
```
sudo apt install liblzma-dev libbz2-dev libgsl0-dev
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
wget -P bcftools https://raw.githubusercontent.com/freeseek/mocha/master/{Makefile.patch,main.patch,vcfnorm.patch,vcfmocha.c,sumlog.c}
wget -P bcftools/plugins https://raw.githubusercontent.com/freeseek/mocha/master/{trio-phase,mochatools,importFMT,extendFMT}.c
cd bcftools && patch < Makefile.patch && patch < main.patch && patch < vcfnorm.patch && cd ..
```

Compile latest version of HTSlib (optionally disable bz2 and lzma) and BCFtools (make sure you are using gcc version 5 or newer or else include the -std=gnu99 compilation flag)
```
cd htslib && autoheader && (autoconf || autoconf) && ./configure --disable-bz2 --disable-lzma && make && cd ..
cd bcftools && make && cd ..
/bin/cp bcftools/{bcftools,plugins/{fill-tags,,fixploidy,split,trio-phase,mochatools,importFMT,extendFMT}.so} $HOME/bin/
```
Notice that you will need some functionalities missing from the base version of bcftools to run the pipeline

Install latest version of <a href="https://data.broadinstitute.org/alkesgroup/Eagle/">Eagle</a>
```
wget -O $HOME/bin/eagle https://data.broadinstitute.org/alkesgroup/Eagle/downloads/dev/eagle_v2.4
chmod a+x $HOME/bin/eagle
```

Install Minimac3 (optional for array data)
```
git clone git://github.com/statgen/Minimac3.git
cd Minimac3
sed -i 's/USER_WARNINGS ?= -Werror/USER_WARNINGS ?= -Wno-format-truncation/' Library/libStatGenForMinimac3/general/Makefile
sed -i 's/bool legacy_count = 0/int legacy_count = 0/' Library/libStatGenForMinimac3/general/Parameters.cpp
make
cd ..
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
  $HOME/bin/bcftools view --no-version -Ou -c 2 ALL.chr${chr}.phase3*integrated_v[125][ab].20130502.genotypes.vcf.gz | \
  $HOME/bin/bcftools norm --no-version -Ou -m -any | \
  $HOME/bin/bcftools norm --no-version -Ob -o ALL.chr${chr}.phase3_integrated.20130502.genotypes.bcf -d none -f $HOME/res/human_g1k_v37.fasta && \
  $HOME/bin/bcftools index -f ALL.chr${chr}.phase3_integrated.20130502.genotypes.bcf
done
```

List of common germline duplications and deletions
```
wget -P $HOME/res ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/integrated_sv_map/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz{,.tbi}
$HOME/bin/bcftools query -i 'AC>1 && END-POS>10000 && TYPE!="INDEL" && (SVTYPE=="CNV" || SVTYPE=="DEL" || SVTYPE=="DUP")' \
  -f "%CHROM\t%POS\t%END\t%SVTYPE\n" $HOME/res/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz > $HOME/res/cnp.grch37.bed
```

List of segmental duplications (make sure your bedtools version is not affected by the groupby <a href="https://github.com/arq5x/bedtools2/issues/418">bug</a>)
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
  cut -f1-3,5 | sed 's/^chr//' | grep -v gl | bgzip > $HOME/res/dup.grch37.bed.gz && \
  tabix -f -p bed $HOME/res/dup.grch37.bed.gz
```

1000 Genomes project phase 3 imputation panel for Minimac3 (optional for array data)
```
cd $HOME/res/kgp
for chr in {1..22} X; do
  if [ $chr == "X" ]; then
    $HOME/bin/bcftools +$HOME/bin/fixploidy.so --no-version ALL.chr$chr.phase3_integrated.20130502.genotypes.bcf -- -f 2 | \
      sed -e 's/\t0\/0/\t0|0/' -e 's/\t1\/1/\t1|1/' > tmp.chr$chr.GRCh37.vcf
  else
    $HOME/bin/bcftools view --no-version ALL.chr$chr.phase3_integrated.20130502.genotypes.bcf -o tmp.chr$chr.GRCh37.vcf
  fi
  $HOME/bin/Minimac3-omp \
    --refHaps tmp.chr$chr.GRCh37.vcf \
    --processReference \
    --myChromosome $chr \
    --prefix ALL.chr$chr.phase3_integrated.20130502.genotypes \
    --noPhoneHome
  /bin/rm tmp.chr$chr.GRCh37.vcf
done
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
  ($HOME/bin/bcftools view --no-version -h ALL.chr${chr}_GRCh38.genotypes.20170504.vcf.gz | \
    grep -v "^##contig=<ID=[GNh]" | sed 's/^##contig=<ID=MT/##contig=<ID=chrM/;s/^##contig=<ID=\([0-9XY]\)/##contig=<ID=chr\1/'; \
  $HOME/bin/bcftools view --no-version -H -c 2 ALL.chr${chr}_GRCh38.genotypes.20170504.vcf.gz | \
  grep -v "[0-9]|\.\|\.|[0-9]" | sed 's/^/chr/') | \
  $HOME/bin/bcftools norm --no-version -Ou -m -any | \
  $HOME/bin/bcftools norm --no-version -Ob -o ALL.chr${chr}_GRCh38.genotypes.20170504.bcf -d none -f $ref && \
  $HOME/bin/bcftools index -f ALL.chr${chr}_GRCh38.genotypes.20170504.bcf
done
```

List of common germline duplications and deletions
```
wget -P $HOME/res ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/integrated_sv_map/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz{,.tbi}
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/liftOver
chmod a+x liftOver
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
$HOME/bin/bcftools query -i 'AC>1 && END-POS>10000 && TYPE!="INDEL" && (SVTYPE=="CNV" || SVTYPE=="DEL" || SVTYPE=="DUP")' \
  -f "chr%CHROM\t%POS\t%END\t%SVTYPE\n" $HOME/res/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz | \
    ./liftOver \
    -minMatch=0.3 \
    /dev/stdin \
    hg19ToHg38.over.chain.gz \
    $HOME/res/cnp.grch38.bed \
    /dev/stderr
```

List of segmental duplications (make sure your bedtools version is not affected by the groupby <a href="https://github.com/arq5x/bedtools2/issues/418">bug</a>)
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
  cut -f1-3,5 | grep -v "GL\|KI" | bgzip > $HOME/res/dup.grch38.bed.gz && \
  tabix -f -p bed $HOME/res/dup.grch38.bed.gz
```

1000 Genomes project phase 3 imputation	panel for Minimac3 (optional for array data)
```
cd $HOME/res/kgp
for chr in {1..22} X; do
  if [ $chr == "X" ]; then
    $HOME/bin/bcftools +$HOME/bin/fixploidy.so --no-version ALL.chr${chr}_GRCh38.genotypes.20170504.bcf -- -f 2 | \
      sed -e 's/\t0\/0/\t0|0/' -e 's/\t1\/1/\t1|1/' > tmp.chr$chr.GRCh38.vcf
  else
    $HOME/bin/bcftools view --no-version ALL.chr${chr}_GRCh38.genotypes.20170504.bcf -o tmp.chr$chr.GRCh38.vcf
  fi
  $HOME/bin/Minimac3-omp \
    --refHaps tmp.chr$chr.GRCh38.vcf \
    --processReference \
    --myChromosome chr$chr \
    --prefix ALL.chr${chr}_GRCh38.genotypes.20170504 \
    --noPhoneHome
  /bin/rm tmp.chr$chr.GRCh38.vcf
done
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

If you want to process <b>genotype array</b> data you need a VCF file with GT, BAF, and LRR information:
```
##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=BAF,Number=1,Type=Float,Description="B Allele Frequency">
##FORMAT=<ID=LRR,Number=1,Type=Float,Description="Log R Ratio">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA12878
1	752566	rs3094315	G	A	.	.	.	GT:BAF:LRR	1|1:0.9889:-0.0798
1	776546	rs12124819	A	G	.	.	.	GT:BAF:LRR	0|1:0.5441:0.4959
1	798959	rs11240777	G	A	.	.	.	GT:BAF:LRR	0|0:0.0288:0.2276
1	932457	rs1891910	G	A	.	.	.	GT:BAF:LRR	1|0:0.4540:-0.1653
```
Making sure that BAF refers to the allele frequency of what in the VCF is indicated as the alternate allele.

If you do not already have a VCF file but you have Illumina or Affymetrix genotype array data, you can use the <a href="https://github.com/freeseek/gtc2vcf">gtc2vcf</a> tools to convert the data to VCF. Alternatively you can use your own scripts.

Create a minimal binary VCF
```
$HOME/bin/bcftools annotate --no-version -Ob -o $dir/$pfx.unphased.bcf -x ID,QUAL,INFO,^FMT/GT,^FMT/BAF,^FMT/LRR $vcf && \
  $HOME/bin/bcftools index -f $dir/$pfx.unphased.bcf
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
Make sure that AD is a "Number=R" format field (this was introduced in version <a href="https://samtools.github.io/hts-specs/VCFv4.2.pdf">4.2</a> of the VCF) or multi-allelic variants will be not <a href="https://github.com/samtools/bcftools/issues/360">split properly</a>.

Create a minimal binary VCF
```
$HOME/bin/bcftools view --no-version -h $vcf | sed 's/^\(##FORMAT=<ID=AD,Number=\)\./\1R/' | \
  $HOME/bin/bcftools reheader -h /dev/stdin $vcf | \
  $HOME/bin/bcftools filter --no-version -Ou -e "FMT/DP<10 || FMT/GQ<20" --set-GT . | \
  $HOME/bin/bcftools annotate --no-version -Ou -x ID,QUAL,INFO,^FMT/GT,^FMT/AD | \
  $HOME/bin/bcftools norm --no-version -Ou -m -any -k | \
  $HOME/bin/bcftools norm --no-version -Ob -o $dir/$pfx.unphased.bcf -f $ref && \
  $HOME/bin/bcftools index $dir/$pfx.unphased.bcf
```
This will set to missing heterozygous genotypes that have low coverage or low genotyping quality, as these can cause issues.

Perform basic quality control (the generated list of variants will be excluded from modeling by both eagle and mocha)
```
n=$($HOME/bin/bcftools query -l $dir/$pfx.unphased.bcf|wc -l); \
ns=$((n*98/100)); \
echo '##INFO=<ID=JK,Number=1,Type=Float,Description="Jukes Cantor">' | \
  $HOME/bin/bcftools annotate --no-version -Ou -a $dup -c CHROM,FROM,TO,JK -h /dev/stdin $dir/$pfx.unphased.bcf | \
  $HOME/bin/bcftools +$HOME/bin/fill-tags.so --no-version -Ou -- -t NS,ExcHet | \
  $HOME/bin/bcftools +$HOME/bin/mochatools.so --no-version -Ou -- -x $sex -G | \
  $HOME/bin/bcftools annotate --no-version -Ob -o $dir/$pfx.xcl.bcf \
    -i 'FILTER!="." && FILTER!="PASS" || JK<.02 || NS<'$ns' || ExcHet<1e-6 || AC_Sex_Test>6' \
    -x FILTER,^INFO/JK,^INFO/NS,^INFO/ExcHet,^INFO/AC_Sex_Test && \
  $HOME/bin/bcftools index -f $dir/$pfx.xcl.bcf
```
This command will create a list of variants falling within segmental duplications with low divergence (<2%), high levels of missingness (>2%), variants with excess heterozygosity (p<1e-6), and variants that correlate with sex in an unexpected way (p<1e-6).

If a file with additional variants to be excluded is available, further merge it with the generated list
```
/bin/mv $dir/$pfx.xcl.bcf $dir/$pfx.xcl.tmp.bcf && \
/bin/mv $dir/$pfx.xcl.bcf.csi $dir/$pfx.xcl.tmp.bcf.csi && \
$HOME/bin/bcftools merge --no-version -Ob -o $dir/$pfx.xcl.bcf -m none $dir/$pfx.xcl.tmp.bcf $xcl && \
$HOME/bin/bcftools index -f $dir/$pfx.xcl.bcf
```

Phasing Pipeline
================

Phase VCF file by chromosome with `eagle`
```
for chr in {1..22} X; do
  $HOME/bin/eagle \
    --geneticMapFile $map \
    --chrom $chr \
    --outPrefix $dir/$pfx.chr$chr \
    --numThreads $thr \
    --vcfRef $kgp_pfx${chr}$kgp_sfx.bcf \
    --vcfTarget $dir/$pfx.unphased.bcf \
    --vcfOutFormat b \
    --noImpMissing \
    --outputUnphased \
    --vcfExclude $dir/$pfx.xcl.bcf && \
  $HOME/bin/bcftools index -f $dir/$pfx.chr$chr.bcf
done
```

Extract chromosomes that do not require phasing
```
$HOME/bin/bcftools view --no-version -Ob -o $dir/$pfx.other.bcf $dir/$pfx.unphased.bcf \
  -t ^$(seq -s, 1 22),X,$(seq -f chr%.0f -s, 1 22),chrX && \
$HOME/bin/bcftools index -f $dir/$pfx.other.bcf
```

Concatenate eagle output into a single VCF file and add GC/CpG content information
```
$HOME/bin/bcftools concat --no-version -Ou --threads $thr $dir/$pfx.{chr{{1..22},X},other}.bcf | \
$HOME/bin/bcftools +$HOME/bin/mochatools.so --no-version -Ob -o $dir/$pfx.bcf -- -f $ref && \
$HOME/bin/bcftools index -f $dir/$pfx.bcf
```

If pedigree information with trios is available, you can improve the phasing with `eagle` by running the following command instead:
```
$HOME/bin/bcftools concat --no-version -Ou --threads $thr $dir/$pfx.{chr{{1..22},X},other}.bcf | \
$HOME/bin/bcftools +$HOME/bin/mochatools.so --no-version -Ou -- -f $ref | \
$HOME/bin/bcftools +$HOME/bin/trio-phase.so --no-version -Ob -o $dir/$pfx.bcf --threads $thr -- -p $ped && \
$HOME/bin/bcftools index -f $dir/$pfx.bcf
```
(it requires a ped file)

Impute variants using Minimac3 (optional for array data)
```
for chr in {1..22} X; do
  $HOME/bin/bcftools annotate --no-version -Ou $dir/$pfx.bcf -a $dir/$pfx.xcl.bcf -m +XCL -r \$chr,chr\$chr | \
    $HOME/bin/bcftools view --no-version -Ov -o $dir/$pfx.chr$chr.vcf -e "XCL==1"
  $HOME/bin/Minimac3-omp \
    --cpus $thr \
    --refHaps $kgp_pfx${chr}$kgp_sfx.m3vcf.gz \
    --format GT,GP \
    --haps $dir/$pfx.chr$chr.vcf \
    --prefix $dir/$pfx.chr$chr \
    --lowMemory \
    --noPhoneHome
  $HOME/bin/bcftools query -l $dir/$pfx.chr$chr.bcf | \
    $HOME/bin/bcftools view --no-version -Ob -o $dir/$pfx.chr$chr.dose.bcf -S /dev/stdin $dir/$pfx.chr$chr.dose.vcf.gz
  /bin/rm $dir/$pfx.chr$chr.vcf $dir/$pfx.chr$chr.dose.vcf.gz
done
```

Concatenate imputed genotypes into a single VCF file (optional for array data)
```
$HOME/bin/bcftools concat --no-version -Ob -o $dir/$pfx.dose.bcf --threads $thr $dir/$pfx.chr{{1..22},X}.dose.vcf.gz && \
$HOME/bin/bcftools index -f $dir/$pfx.dose.bcf
```

Remove unphased VCF and single chromosome files (optional)
```
/bin/rm $dir/$pfx.{unphased,chr{{1..22},X},other}.bcf{,.csi} $dir/$pfx.chr{{1..22},X}.dose.bcf
```

Chromosomal Alterations Pipeline
================================

Preparation steps
```
pfx="..." # output prefix
thr="..." # number of extra threads to use
lst="..." # file with list of samples to analyze for asymmetries (e.g. samples with 1p CNN-LOH)
```

Call mosaic chromosomal alterations with MoChA
```
$HOME/bin/bcftools mocha \
  -r $rule \
  --no-version -Ob \
  -o $dir/$pfx.mocha.bcf \
  --threads $thr \
  --variants ^$dir/$pfx.xcl.bcf \
  -m $dir/$pfx.mocha.tsv \
  -g $dir/$pfx.stats.tsv \
  -u $dir/$pfx.ucsc.bed \
  -p $cnp \
  $dir/$pfx.bcf && \
$HOME/bin/bcftools index -f $dir/$pfx.mocha.bcf
```

import results from MoChA into VCF file with imputed genotypes (optional for array data)
```
$HOME/bin/bcftools +$HOME/bin/importFMT.so \
  --no-version -Ob \
  --formats Ldev,Bdev,Bdev_Phase \
  $dir/$pfx.dose.bcf \
  $dir/$pfx.mocha.bcf \
  -o $dir/$pfx.dose.mocha.bcf && \
$HOME/bin/bcftools index $dir/$pfx.dose.mocha.bcf
```

Run asymmetry analyses (subset cohort, run binomial test, discard genotypes)
```
$HOME/bin/bcftools +$HOME/bin/extendFMT.so \
  --no-version -Ou \
  --format Bdev_Phase \
  --phase \
  --dist 500000 \
  -r $reg \
  -s $lst \
  $dir/$pfx.dose.mocha.bcf | \
$HOME/bin/bcftools +$HOME/bin/mochatools.so \
  --no-version -Ob \
  -o $dir/$pfx.bal.bcf \
  -- -b Bdev_Phase -G && \
$HOME/bin/bcftools index -f $dir/$pfx.bal.bcf
```

Observe results for asymmetry analyses in table format
```
fmt="%CHROM\t%POS\t%ID\t%Bal{0}\t%Bal{1}\t%Bal_Test\n"
$HOME/bin/bcftools query \
  -i "Bal_Test>6" \
  -f "$fmt" \
  $dir/$pfx.bal.bcf | \
  column -ts $'\t'
```

Split output VCF file by samples (optional)
```
$HOME/bin/bcftools +$HOME/bin/split.so -Ob -o $dir/ $dir/$pfx.mocha.bcf -k FMT
```

Plot Results
============

Install basic tools (Debian/Ubuntu specific):
```
sudo apt install r-cran-ggplot2 r-cran-data.table r-cran-gridextra
```

Download R scripts
```
wget -P $HOME/bin https://raw.githubusercontent.com/freeseek/mocha/master/plot_{summary,mocha}.R
chmod a+x $HOME/bin/plot_{summary,mocha}.R
```

Generate summary plot
```
$HOME/bin/plot_summary.R $dir/$pfx.pdf $dir/$pfx.stats.tsv $dir/$pfx.mocha.tsv
```

Plot mosaic chromosomal alterations
```
for sm in $($HOME/bin/bcftools query -l $dir/$pfx.mocha.bcf | tr ' ' '_'); do
  $HOME/bin/plot_mocha.R $dir/$sm.pdf $dir/$sm.bcf
done
```

Acknowledgements
================

This work is supported by NIH grant <a href="https://projectreporter.nih.gov/project_info_description.cfm?aid=8852155">R01 HG006855</a> and the Stanley Center for Psychiatric Research and by US Department of Defense Breast Cancer Research Breakthrough Award W81XWH-16-1-0316 (project BC151244). This work would have not been possible without the efforts of Heng Li <lh3@sanger.ac.uk>, Petr Danecek <pd3@sanger.ac.uk>, John Marshall <jm18@sanger.ac.uk>, James Bonfield <jkb@sanger.ac.uk>, and Shane McCarthy <sm15@sanger.ac.uk> in building HTSlib and BCFtools and Po-Ru Loh <poruloh@broadinstitute.org> in building the Eagle phasing software.
