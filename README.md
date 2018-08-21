MoChA (MOsaic CHromosomal Alterations caller)
=============================================

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
```

Installation
============

Install basic tools (Debian/Ubuntu specific):
```
sudo apt-get install wget liblzma-dev libbz2-dev libgsl0-dev gzip samtools unzip
```

Preparation steps
```
mkdir -p $HOME/bin $HOME/res $HOME/res/chrs && cd /tmp
```

Download latest version of `htslib` and `bcftools` (if not downloaded already)
```
git clone --branch=develop git://github.com/samtools/htslib.git
git clone --branch=develop git://github.com/samtools/bcftools.git
```

Add patches and code for plugin
```
wget -P bcftools https://raw.githubusercontent.com/freeseek/mocha/master/{Makefile.patch,main.patch,vcfnorm.patch,vcfmocha.c,sumlog.c}
wget -P bcftools/plugins https://raw.githubusercontent.com/freeseek/mocha/master/{trio-phase,mochatools}.c
cd bcftools && patch < Makefile.patch && patch < main.patch && patch < vcfnorm.patch && cd ..
```

Compile latest version of `htslib` (optionally disable `bz2` and `lzma`) and `bcftools`
```
cd htslib && autoheader && (autoconf || autoconf) && ./configure --disable-bz2 --disable-lzma && make && cd ..
cd bcftools && make && cd ..
/bin/cp bcftools/{bcftools,plugins/{fill-tags,split,trio-phase,mochatools}.so} $HOME/bin/
```
Notice that you will need some functionalities missing from the base version of bcftools to run the pipeline

Install latest version of `eagle`
```
wget -O $HOME/bin/eagle https://data.broadinstitute.org/alkesgroup/Eagle/downloads/dev/eagle_v2.4
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
cd $HOME/res/chrs
wget http://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/chr{1-4,5-9,10-15,16-22X}.1kg.phase3.v5a.{bref,vcf}.zip
for chrs in 1-4 5-9 10-15 16-22X; do unzip chr$chrs.1kg.phase3.v5a.vcf.zip; done
for chr in {1..22} X; do
  $HOME/bin/bcftools view --no-version -Ob chr$chr.1kg.phase3.v5a.vcf.gz -o chr$chr.1kg.phase3.v5a.bcf && \
  $HOME/bin/bcftools index -f chr$chr.1kg.phase3.v5a.bcf
done
```

List of common germline duplications and deletions
```
wget -P $HOME/res ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/integrated_sv_map/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz{,.tbi}
bcftools query -i 'AC>1 && END-POS>10000 && TYPE!="INDEL" && (SVTYPE=="CNV" || SVTYPE=="DEL" || SVTYPE=="DUP")' \
  -f "%CHROM\t%POS\t%END\t%SVTYPE\n" $HOME/res/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz > $HOME/res/cnp.grch37.bed
```

Setup variables
```
ref="$HOME/res/human_g1k_v37.fasta"
map="$HOME/res/genetic_map_hg19_withX.txt.gz"
kgp_pfx="$HOME/res/chrs/chr"
kgp_sfx=".1kg.phase3.v5a.bcf"
rule="GRCh37"
cnp="$HOME/res/cnp.grch37.bed"
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

1000 Genomes project phase 3 (fixing contig names and removing duplicate variants)
```
cd $HOME/res/chrs
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/ALL.chr{{1..22},X,Y}_GRCh38.genotypes.20170504.vcf.gz{,.tbi}
ref="$HOME/res/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
for chr in {1..22} X Y; do
  ($HOME/bin/bcftools view --no-version -h ALL.chr${chr}_GRCh38.genotypes.20170504.vcf.gz | \
    grep -v "^##contig=<ID=[GNh]" | sed 's/^##contig=<ID=MT/##contig=<ID=chrM/;s/^##contig=<ID=\([0-9XY]\)/##contig=<ID=chr\1/'; \
  $HOME/bin/bcftools view --no-version -H -c 2 ALL.chr${chr}_GRCh38.genotypes.20170504.vcf.gz | sed 's/^/chr/') | \
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
bcftools query -i 'AC>1 && END-POS>10000 && TYPE!="INDEL" && (SVTYPE=="CNV" || SVTYPE=="DEL" || SVTYPE=="DUP")' \
  -f "chr%CHROM\t%POS\t%END\t%SVTYPE\n" $HOME/res/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz | \
    ./liftOver \
    -minMatch=0.3 \
    /dev/stdin \
    hg19ToHg38.over.chain.gz \
    $HOME/res/cnp.grch38.bed \
    /dev/stderr
```

Setup variables
```
ref="$HOME/res/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
map="$HOME/res/genetic_map_hg38_withX.txt.gz"
kgp_pfx="$HOME/res/chrs/ALL.chr"
kgp_sfx="_GRCh38.genotypes.20170504.bcf"
rule="GRCh38"
cnp="$HOME/res/cnp.grch38.bed"
```

Phasing Pipeline
================

Preparation steps
```
tbl="..." # input Illumina GenomeStudio table
vcf="..." # input VCF file with phased GT, LRR, and BAF
pfx="..." # output prefix
thr="..." # number of threads to use
sex="..." # file with sex information (first column sample ID, second column sex: 1=male; 2=female)
ped="..." # pedigree file to use if parent child duos are present
dir="mocha"
mkdir -p $dir
```

If you already have a VCF file with BAF and LRR information, create a minimal binary VCF
```
$HOME/bin/bcftools annotate --no-version -Ob -o $dir/$pfx.unphased.bcf -x FILTER,INFO,^FMT/GT,^FMT/BAF,^FMT/LRR $vcf && \
  $HOME/bin/bcftools index -f $dir/$pfx.unphased.bcf
```
Otherwise, if you have Illumina genotype array data, use the <a href="https://github.com/freeseek/gtc2vcf">gtc2vcf</a> plugin to convert the data to VCF

If you prefer to bring your BAF and LRR data to VCF format using your own scripts, it should have the following format:
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

Perform basic quality control
```
n=$(bcftools query -l $dir/$pfx.unphased.bcf|wc -l)
ns=$((n*98/100))
$HOME/bin/bcftools +$HOME/bin/fill-tags.so --no-version -Ou $dir/$pfx.unphased.bcf -- -t NS,ExcHet | \
  $HOME/bin/bcftools +$HOME/bin/mochatools.so --no-version -Ou -- -x $sex -G | \
  $HOME/bin/bcftools annotate --no-version -Ob -o $dir/$pfx.xcl.bcf -i "NS<$ns || ExcHet<1e-6 || AC_Sex_Test>6" \
    -x FILTER,^INFO/NS,^INFO/ExcHet,^INFO/AC_Sex_Test && \
  $HOME/bin/bcftools index -f $dir/$pfx.xcl.bcf
```
(the generated list of variants will be excluded from modeling by both eagle and mocha)

Phase VCF file by chromosome with `eagle`
```
for chr in {1..22} X; do
  $HOME/bin/eagle \
    --geneticMapFile $map \
    --chrom $chr \
    --outPrefix $dir/$pfx.chr$chr \
    --numThreads $thr \
    --vcfRef $kgp_pfx$chr$kgp_sfx \
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

Remove unphased VCF and single chromosome files (optional)
```
/bin/rm $dir/$pfx.{unphased,chr{{1..22},X},other}.bcf{,.csi}
```

Chromosomal Alterations Pipeline
================================

Preparation steps
```
pfx="..." # output prefix
thr="..." # number of extra threads to use
lst="..." # file with list of samples to analyze for asymmetries (e.g. samples with 1p CNN-LOH)
```

Call mosaic chromosomal alterations with `mocha`
```
$HOME/bin/bcftools mocha \
  --no-version -Ob \
  -o $dir/$pfx.mocha.bcf \
  --threads $thr \
  --variants ^$dir/$pfx.xcl.bcf \
  -m $dir/$pfx.mocha.tsv \
  -g $dir/$pfx.stats.tsv \
  -u $dir/$pfx.ucsc.bed \
  -r $rule \
  -p $cnp \
  $dir/$pfx.bcf && \
$HOME/bin/bcftools index -f $dir/$pfx.mocha.bcf
```

Run asymmetry analyses (subset cohort, run binomial test, discard genotypes)
```
$HOME/bin/bcftools +$HOME/bin/mochatools.so \
  --no-version -Ou -r X \
  -o $dir/$pfx.ase.bcf \
  $dir/$pfx.mocha.bcf \
  -- -b -S $lst -G && \
$HOME/bin/bcftools index -f $dir/$pfx.ase.bcf
```

Observe results for asymmetry analyses in table format
```
fmt="%CHROM\t%POS\t%ID\t%ASE{0}\t%ASE{1}\t%ASE_Test\n"
bcftools query -i "ASE_Test<1e-6" -f "$fmt" $dir/$pfx.ase.bcf | tlmn
```

Split output VCF file by samples (optional)
```
$HOME/bin/bcftools +$HOME/bin/split.so -Ob -o $dir/ $dir/$pfx.mocha.bcf -k FMT
```

Plot Results
============

Install basic tools (Debian/Ubuntu specific):
```
sudo apt-get install r-cran-ggplot2 r-cran-data.table r-cran-gridextra
```

Download R scripts
```
wget -P $HOME/bin bcftools https://raw.githubusercontent.com/freeseek/mocha/master/plot_{summary,mocha}.R
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
