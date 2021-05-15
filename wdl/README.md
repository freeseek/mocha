MoChA WDL pipeline
==================

This page contains instructions for how to run the <a href="mocha.wdl">MoChA WDL</a> pipeline to detect mosaic chromsosomal alterations. The MoChA WDL pipeline first runs Illumina\'s GenCall or Affymetrix\'s Axiom genotyping algorithms, then gtc2vcf to format Illumina or Affymetrix genotype data in more compliant VCF containers, then runs <a href="https://odelaneau.github.io/shapeit4/">SHAPEIT4</a> (or <a href="https://alkesgroup.broadinstitute.org/Eagle/">Eagle</a>) to phase genotypes across overlapping genome windows, and finally it runs MoChA to detect mosaic chromosomal alterations. It can also be used to analyze whole genome sequencing data. Using a highly performant and parallelizable design, this workflow can be scaled to satisfy the needs for large biobanks. The workflow also allows for automatic realigning of manifest files to a genome reference of choice making it effortless to use GRCh38 even when GRCh38 manifest files are not available from the array manufacturer. For any feedback or questions, contact the <a href="mailto:giulio.genovese@gmail.com">author</a>

<!--ts-->
   * [Input Modes](#input-modes)
   * [Input Files](#input-files)
   * [Primary Options](#primary-options)
   * [Secondary Options](#secondary-options)
   * [Outputs](#outputs)
   * [Illumina Example](#illumina-example)
   * [Affymetrix Example](#affymetrix-example)
   * [Imputation Pipeline](#imputation-pipeline)
   * [Allelic Shift Pipeline](#allelic-shift-pipeline)
   * [Polygenic Score Pipeline](#polygenic-score-pipeline)
   * [Running with Terra](#running-with-terra)
   * [Running with Cromwell on Google Cloud Platform](#running-with-cromwell-on-google-cloud-platform)
   * [Filtering Output Calls](#filtering-output-calls)
   * [Scattered Tasks](#scattered-tasks)
   * [Troubleshooting](#troubleshooting)
   * [Implementation Notes](#implementation-notes)
   * [Dockerfiles](#dockerfiles)
   * [Acknowledgements](#acknowledgements)
<!--te-->

Input modes
===========

![](mocha_wdl.png)

Due to different scenarios for how DNA microarray data is made available, the MoChA WDL pipeline can accept seven possible different data types (visualized in the previous figure as labels in red clouds):

| mode     | description                                        |
|----------|----------------------------------------------------|
| **idat** | intensity files (two per sample)                   |
| **gtc**  | genotype/intensities files (one per sample)        |
| **cel**  | intensity files (one per sample)                   |
| **chp**  | genotype/intensities files (one per sample)        |
| **txt**  | genotype/intensities tables (two per batch)        |
| **vcf**  | unphased variant call format files (one per batch) |
| **pvcf** | phased variant call format files (one per batch)   |

If the **mode** variable is not set to be one of these seven values, the behavior will be undefined. Although we advise against it, if you have an Illumina GenomeStudio file and no other options to access the raw data, you can convert the file to a compliant VCF using <a href="https://github.com/freeseek/gtc2vcf">gtc2vcf</a> and then use the WDL pipeline through the **vcf** mode. However, in this scenario there is no way to select a genome reference different from the one chosen when the GenomeStudio file was created

Notice that if you want to use the **idat** mode, this requires to run either the <a href="https://support.illumina.com/array/array_software/illumina-array-analysis-platform/downloads.html">Illumina Array Analysis Platform</a> Genotyping Command Line Interface or the <a href="https://support.illumina.com/array/array_software/beeline/downloads.html">Illumina AutoConvert</a> Software. Neither software is free (see <a href="https://support.illumina.com/downloads/iaap-genotyping-orchestrated-workflow.html">EULA</a>) and both are covered by a <a href="https://patents.google.com/patent/US7035740">patent</a> expiring in 2024 so check first to make sure you are allowed to run this software. Similarly, while MoChA is completely open source and released under the MIT license, there is a <a href="https://patents.google.com/patent/WO2019079493">patent</a> that covers the method. However, the use of MoChA for research purposes is unrestricted

You can run the pipeline all the way to the initial conversion to unphased VCF without phasing by setting the **target** variable to **vcf** or all the way to phased VCF without running MoChA by setting the **target** variable to **pvcf**. You can then input these VCFs back into the pipeline using, respectively, modes **vcf** or **pvcf**, although you will need to provide information about **computed_gender** and **call_rate** that are generated during the conversion to VCF. If you want to run the pipeline on data from whole genome sequencing, then you can use the **vcf** mode while making sure to selecte the **wgs** boolean variable appropriately. VCFs from whole genome sequencing data require specific handling when being merged and phased and we do not advise to use the pipeline on single sample whole genomes sequencing VCFs, as they don't record any information at sites where the sample is homozygous reference

The pipeline will perform a minimal amount of quality control. It is up to the user to correctly handle failed samples and duplicates. In particular, if a large amount of duplicate samples is present, it is important to include as input a list of samples to remove during the QC, or this step could run with incorrect assumptions. Furthermore, when including the conversion to VCF step, the pipeline will output three files, **sample_id_file**, **computed_gender_file**, and **call_rate_file**. We highly recommend to check for call rate failure rate and for discrepancies between reported gender and computed gender that might be indicative of sample swaps

Currently no arrays other than the Affymetrix Genome-Wide Human SNP Array 6.0 are supported in **cel** mode. Newer Axiom arrays used in biobanks have less uniform ways to run the genotype clustering algorithms that are difficult to package into a pipeline. In these scenarios we recommend to perform the conversion from CEL files to table output separately and then use the pipeline in **txt** mode

Input Files
===========

The following diagram depicts all possible files that you might need to input into the workflow:

![](mocha_input.png)

Not all files will be needed and which ones are needed will be defined by the analysis mode selected

For Illumina data you will need to provide a CSV manifest file, a BPM manifest file, and an EGT cluster file. Illumina provides these files for their arrays <a href="https://support.illumina.com/array/downloads.html">here</a> or through their customer support. You will need to make sure you are using the correct manifest and cluster files or else the behavior will be undefined. All Illumina arrays should be supported. Notice that some custom Illumina arrays, such as the Global Screening Array + Multi Disease (GSAMD), are extensions of more mainstream Illumina arrays, such as the Global Screening Array (GSA). A manifest file for the latter will likely worker for data from the former, but you will end up missing some custom probes in doing so. Illumina only provides manifest files for custom arrays through their customer support

For Affymetrix data you will need to provide a CSV manifest file. If you are providing CEL files you will need to provide an XML option file and then a ZIP file including all the files listed within the XML option file. You will need to manually create the ZIP file. Alternatively, you can also provide a mix of SNP posterior files containing cluster information and either CHP files or calls TXT and summary TXT files containing genotype and intensities information, all files that can be generated with apt-probeset-genotype or the apt-genotype-axiom <a href="https://www.thermofisher.com/us/en/home/life-science/microarray-analysis/microarray-analysis-partners-programs/affymetrix-developers-network/affymetrix-power-tools.html">Affymetrix Power Tools</a> softwares. However, only analysis types AxiomGT1 and birdseed are supported. In particular, Affymetrix arrays before the Genome-Wide Human SNP Array 6.0 are not supported

The **vcf** and **pvcf** modes require the user to provide a set of VCF files together with a master tracker table including sample IDs, computed gender, and call rate information. If the data is from DNA microarray cohorts, then the VCFs provided must include the INFO fields ALLELE_A, ALLELE_B, and GC and FORMAT fields GT, BAF, and LRR and we recommend using <a href="https://github.com/freeseek/gtc2vcf">gtc2vcf</a> to generate such a VCF. If the data is from WGS cohorts, then the VCFs provided must include the INFO field GC and the FORMAT fields GT and AD. In the **pvcf** mode the VCF file must be provided already phased and VCF files with variants to exclude from further analysis should also be provided. Notice that you can define a unique set of variants to exclude in the analysis across all batches through the optional file **extra_xcl_vcf_file**. This can be useful for WGS cohorts to exclude variants that are known from other analyses as having unreliable allelic balance due to misalignment issues. Notice also that for WGS cohorts multi-allelic variants should be split as explained <a href="https://github.com/freeseek/mocha#data-preparation">here</a>

Allowed columns in the sample table:

| mode            | gzipped | idat | gtc  | cel  | chp  | txt  | vcf  | pvcf |
|-----------------|---------|------|------|------|------|------|------|------|
| sample_id       |         | req. | req. | req. | req  | req. | req. | req. |
| batch_id        |         | req. | req. | req. | req  | req. |      |      |
| computed_gender |         |      |      |      |      |      | req. | req. |
| call_rate       |         |      |      |      |      |      | req. | reg. |
| green_idat      | allowed | req. |      |      |      |      |      |      |
| red_idat        | allowed | req. |      |      |      |      |      |      |
| gtc             | allowed |      | req. |      |      |      |      |      |
| cel             | allowed |      |      | req  |      | req. |      |      |
| chp             | allowed |      |      |      | req. |      |      |      |

All files in this table can be provided in gzip format. The pipeline will process these seamlessly

It is extremely important that the **sample_id** column contains unique IDs. Repeated IDs will cause undefined behavior. You can use `awk -F"\t" 'x[$1]++' sample.tsv` to verify that the first column contains unique IDs. If this is not the case, you can use something like
```
awk -F"\t" -v OFS="\t" '{x[$1]++; if (x[$1]>1) $1=$1"-"(x[$1]-1); print}' sample.tsv
```
to generate a new table with additional suffixes in the sample IDs to make the IDs unique. It is instead allowed to have the same file names in other columns as long as the **batch_id** for these duplicates is different

When using the pipeline in **idat**, **gtc**, **cel**, or **chp** mode, it is extremely important that IDAT, GTC, CEL, or CHP files within the same batch come from the same DNA microarray model. Discrepancies across the input files within the same batch will cause undefined behavior. It is okay to use files from different DNA microarrays across batches but we highly recommend that these being different versions of the same array (e.g\. Global Screening Array v1.0 vs\. the Global Screening Array v2.0)

We further recommend batches of 4,000-5,000 samples as this is small enough to make most tasks able to run on preemptible computing and large enough to make MoChA able to perform effective within-batch adjustments for BAF and LRR. For DNA microarrays such as the Global Screening Array and the Multi-Ethnic Global Array, this batch sizes will generate VCFs smaller than, respectively, 20GB and 60GB

If using the **idat** mode, the green and red IDAT file names have to match and end, respectively, with suffixes **\_Grn.idat** and **\_Red.idat**. If the file names have been changed and these suffixes have been modified, unfortunately you will have to manually change the names to first satisfy this requirement or else the conversion to GTC will fail as the Illumina software requires this naming standard

When using the **txt** mode, you don't have to provide CEL files, but you have to provide the **cel** column with a list of CEL files names grouped by batch that match the CEL file names present in the calls TXT and summary TXT files. These will be used exclusively to rename the IDs, especially in situations where different batches contain the same CEL file IDs. You will also have to make sure that **calls** and **summary** tables contain the same samples in the same order and that these are indeed present in the **report** table or else the behavior will be undefined. To verify the list of samples in each table, you can use the following commands:
```
grep -v ^# AxiomGT1.calls.txt | head -n1 | tr '\t' '\n' | tail -n+2
grep -v ^# AxiomGT1.summary.txt | head -n1 | tr '\t' '\n' | tail -n+2
grep -v ^# AxiomGT1.report.txt | cut -f1 | tail -n+2
```

If you are using the pipeline in **vcf** or **pvcf** mode you have to supply **computed_gender** and **call_rate** information. This can be extracted from a previous run of the pipeline that included the conversion to VCF using a command like the following:
```
(echo -e "sample_id\tcomputed_gender\tcall_rate"; \
paste -d $'\t' sample_id.lines computed_gender.lines | \
paste -d $'\t' - call_rate.lines) > sample.tsv
```

Allowed columns in the batch table:

| mode          | gzipped | idat | gtc  | cel  | chp  | txt  | vcf  | pvcf |
|---------------|---------|------|------|------|------|------|------|------|
| batch_id      |         | req. | req. | req. | req. | req. | req. | req. |
| csv           | allowed | req. | req. | req. | req. | req. |      |      |
| sam           |         | opt. | opt. | opt. | opt. | opt. |      |      |
| bpm           | allowed | req. | req. |      |      |      |      |      |
| egt           | allowed | req. | req. |      |      |      |      |      |
| xml           |         |      |      | req. |      |      |      |      |
| zip           |         |      |      | req. |      |      |      |      |
| path          |         | opt. | opt. | opt. | opt. | opt. | opt. | opt. |
| snp           | allowed |      |      |      | req. | req. |      |      |
| probeset_ids  | allowed |      |      | opt. | opt. | opt. |      |      |
| calls         | allowed |      |      |      |      | req. |      |      |
| confidences   | allowed |      |      |      |      | opt. |      |      |
| summary       | allowed |      |      |      |      | req. |      |      |
| report        | allowed |      |      |      |      | req. |      |      |
| vcf           |         |      |      |      |      |      | req. | req. |
| vcf_index     |         |      |      |      |      |      | req. | req. |
| xcl_vcf       |         |      |      |      |      |      |      | req. |
| xcl_vcf_index |         |      |      |      |      |      |      | req. |

Some files in this table can be provided in gzip format and the **sam** alignment file can be provided either as a SAM or as a BAM file. The pipeline will process these seamlessly

It is extremely important that the **batch_id** column contains unique IDs. Repeated IDs will cause undefined behavior. You can use `awk -F"\t" 'x[$1]++' batch.tsv` to verify that the first column contains unique IDs. It is also equally important that, if the boolean **realign** variable is left to its default **false** value, then the **csv** manifest files must be provided with respect to the same reference genome as the one selected in the **ref_name** variable. A mismatch in the two references will cause undefined behavior. For newer DNA microarrays, Illumina follows a convention of providing manifest files ending with the **1** suffix for GRCh37 (e.g\. `Multi-EthnicGlobal_D1.csv`) amd ending with the **2** suffix for GRCh38 (e.g\. `Multi-EthnicGlobal_D2.csv`). We recommend to always use GRCh38 with the latter type of manifest files and, whenever GRCh38 manifests are not available from Illumina, to still use GRCh38 by setting the boolean **realign** variable to **true**

If you are running the pipeline in **idat** mode, it is important that the **bpm** file name matches the internal descriptor file name or a safety check during the conversion from GTC to VCF will fail. You can verify the internal descriptor file name with the command `bcftools +gtc2vcf -b $bpm_file`. We have observed discrepancies only in old Illumina BPM manifest files. If necessary, you can use set the boolean flag **do_not_check_bpm** to turn off this safety check

The **snp** files are usually the ones generated by the Affymetrix Power Tools as **AxiomGT1.snp-posteriors.txt** files. If the **sam** files are provided then the genome coordinates in the **csv** manifests will be updated during conversion to **vcf** even if the boolean flag **realign** is set to false

If the CSV manifest file is provided according to a genome reference different from the one you wish to run the MoChA analysis on, you can either provide a SAM/BAM file with realignments for the surrounding flanking sequences for each SNP as described <a href="https://github.com/freeseek/gtc2vcf#using-an-alternative-genome-reference">here</a>, or you can set the **realign** boolean flag to **true** so that the SNPs in the manifest files will be automatically realigned and, when the VCF files are created, the entries will be assigned coordinates for the new reference on the fly

To reduce redundancy in your tables, you can omit the paths of the files included in each group (as depicted in the previous figure) and include them once for each group in the optional variable **manifest_path**, **ref_path**, and in the optional **path** column in the **batch_tsv_file** table

The **batch_id** column defines batch membership through a categorical variable. Any strings can be used but the user should avoid to use a string containing a substring used internally by the pipeline (modifiable by the user though the variable **delim**) to refine batches in sub batches

Primary Options
===============

The following are the primary options that you can set in the main input json file. Parameters that are indicated with ? are optional (unless required by the specific mode selected)

| key                    | type     | description                                                                                       |
|------------------------|----------|---------------------------------------------------------------------------------------------------|
| sample_set_id          | String   | cohort name that will be used as prefix for temporary and output files                            |
| mode                   | String   | pipeline mode, one of: **idat**, **gtc**, **cel**, **chp**, **txt**, **vcf**, or **pvcf**         |
| target                 | String?  | pipeline final state, one of: **vcf**, **pvcf**, **calls**, or **pngs** [**pngs**]                |
| realign                | Boolean? | whether manifest file should be realigned (not in **vcf** or **pvcf** mode) [false]               |
| wgs                    | Boolean? | whether the input VCF contains a high coverage whole genome sequencing data [false]               |
| gtc_output             | Boolean? | whether to output the conversion from IDAT to GTC (only in **idat** mode) [false]                 |
| chp_output             | Boolean? | whether to output the conversion from CEL to CHP (only in **cel** mode) [false]                   |
| idat_batch_size        | Int?     | largest batch size for conversion of IDAT files to GTC [48] (only in **idat** mode)               |
| gtc_batch_size         | Int?     | largest batch size for conversion of GTC files to VCF [1024] (only in **idat** and **gtc** mode)  |
| chp_batch_size         | Int?     | largest batch size for conversion of CHP files to VCF [1024] (only in **cel** and **chp** mode)   |
| max_win_size_cm        | Float?   | maximum windows size in cM for phasing [300.0]                                                    |
| overlap_size_cm        | Float?   | required overlap size in cM for consecutive windows [5.0]                                         |
| sample_call_rate_thr   | Float?   | Minimum sample call rate allowed for quality control and plotting purposes [0.97]                 |
| variant_call_rate_thr  | Float?   | Minimum sample call rate allowed for plotting purposes [0.97]                                     |
| baf_auto_thr           | Float?   | Minimum sample call rate allowed for plotting purposes [0.03]                                     |
| ext_string             | String?  | Extension suffix for VCF files with AS information [as]                                           |
| ref_name               | String?  | name of reference genome, with resource default files for GRCh37 and GRCh38 [GRCh38]              |
| ref_path               | String?  | path for reference genome resources (needed unless all resources are provided with full path)     |
| ref_fasta              | String?  | reference sequence [GCA_000001405.15_GRCh38_no_alt_analysis_set.fna/human_g1k_v37.fasta]          |
| n_chrs                 | Int?     | number of chromosomes [23]                                                                        |
| mhc_reg                | String?  | chromosome for MHC [chr6:27518932-33480487/6:27486711-33448264]                                   |
| kir_reg                | String?  | chromosome for KIR [chr19:54071493-54992731/19:54574747-55504099]                                 |
| dup_file               | String?  | file with location of segmental duplications [segdups.bed.gz]                                     |
| genetic_map_file       | String?  | genetic map [genetic_map_hg38_withX.txt.gz/genetic_map_hg19_withX.txt.gz]                         |
| cnp_file               | String?  | file with location of copy number polymorphism [cnps.bed]                                         |
| cyto_file              | String?  | file with location of cytoband regions [cytoBand.txt.gz]                                          |
| panel_pfx              | String?  | prefix for phasing reference [CCDG_14151_B01_GRM_WGS_2020-08-05_/ALL.chr]                         |
| panel_sfx              | String?  | suffix for phasing reference [.filtered.phased.bcf/.phase3_integrated.20130502.genotypes.bcf]     |
| panel_idx              | String?  | index extension for phasing reference [.bcf]                                                      |
| n_panel_smpls          | Int?     | number of samples in phasing reference [3202/2504]                                                |
| manifest_path          | String?  | path for manifest file resources if these are provided without path []                            |
| data_path              | String?  | path for data files (overrides **path** column in **batch_tsv_file**)                             |
| sample_tsv_file        | File     | TSV file with sample information                                                                  |
| batch_tsv_file         | File     | TSV file with batch information                                                                   |
| ped_file               | File?    | optional PED file for improved phasing with trios                                                 |
| duplicate_samples_file | File?    | optional file with list of duplicate samples that should not be use in task vcf_qc                |
| extra_xcl_vcf_file     | File?    | optional VCF file with list of additional variants to exclude from analysis, mostly for WGS data  |
| phase_extra_args       | String?  | extra arguments for SHAPEIT4/Eagle                                                                |
| mocha_extra_args       | String?  | extra arguments for the MoChA plugin                                                              |
| basic_bash_docker      | String?  | docker to run basic bash scripts [ubuntu:latest]                                                  |
| pandas_docker          | String?  | docker to run task ref_scatter [amancevice/pandas:slim]                                           |
| docker_registry        | String?  | location of docker images [us.gcr.io/mccarroll-mocha]                                             |
| bcftools_docker        | String?  | docker to run tasks requiring BCFtools [bcftools:1.11-yyyymmdd]                                   |
| iaap_cli_docker        | String?  | docker to run task idat2gtc [iaap_cli:1.11-yyyymmdd]                                              |
| autoconvert_docker     | String?  | docker to run task idat2gtc [autoconvert:1.11-yyyymmdd]                                           |
| apt_docker             | String?  | docker to run task cel2chp [apt:1.11-yyyymmdd]                                                    |
| shapeit4_docker        | String?  | docker to run task vcf_phase [shapeit4:1.11-yyyymmdd]                                             |
| eagle_docker           | String?  | docker to run task vcf_phase [eagle:1.11-yyyymmdd]                                                |
| r_mocha_docker         | String?  | docker to run tasks mocha_{plot,summary} [r_mocha:1.11-yyyymmdd]                                  |

The **ref_path** variable should contain the path to the genome reference resources. These are available for download <a href="http://software.broadinstitute.org/software/mocha">here</a> for either the GRCh37 or GRCh38 human genome reference

However, if you do not provide the **ref_path** variable, variables **ref_fasta**, **dup_file**, **genetic_map_file**, **cnp_file**, **cyto_file**, and **panel_pfx** will need to be provided with their full path. Also notice that the reference genome file should come with the fasta index file and, if you request the manifest files to be realigned, you will need to make sure it also comes with the five corresponding bwa index files. The fasta index file needs to be such that the first 23 entries correspond to the 22 human autosomes and chromosome X, in no specific order

The **manifest_path** variable should contain the path to all the CSV/BPM/EGT/XML/ZIP/SAM files necessary to run the analyses. If these manifest files are located in different directories, then provide them with their full path and leave the **manifest_path** variable empty

Docker images will be automatically pulled by the specified **docker_registry** Google bucket. As these GCR buckets are Requester pays buckets, if you run computations in locations separate from where these buckets are localized, you will incur some additional costs. Docker images are located in **us.gcr.io/mccarroll-mocha**, **eu.gcr.io/mccarroll-mocha**, and **asia.gcr.io/mccarroll-mocha** covering, respectively, US, EU, and ASIA <a href="https://cloud.google.com/storage/docs/locations">locations</a>. If you are planning to run the pipeline with resources other than Google Cloud, then you can also use your own container registry. If you are running Cromwell on a computational framework that will simply download the images once and then reuse them, then you can also use **gcr.io/mccarroll-mocha**

Secondary Options
=================

There are options specific to single tasks in the pipeline that can be used and they are mostly to handle specific corner cases. See the softwares involved to learn about what each option changes

| key                  | type          | task                   | description                                                             |
|----------------------|---------------|------------------------|-------------------------------------------------------------------------|
| autoconvert          | Boolean       | idat2gtc               | if true uses AutoConvert rather than IAAP CLI [false]                   |
| do_not_check_bpm     | Boolean       | gtc2vcf                | do not check whether BPM and GTC files match manifest file name [false] |
| table_output         | Boolean       | cel2chp                | output matrices of tab delimited genotype calls and confidences [false] |
| do_not_use_reference | Boolean       | vcf_phase              | whether to phase without using a reference panel [false]                |
| eagle                | Boolean       | vcf_phase              | if true uses Eagle rather than SHAPEIT4 to phase [false]                |
| delim                | String        | {idat,gtc,chp}_scatter | string delimiter used to define IDAT/GTC/CHP sub batches [\"~\"]        |
| chip_type            | Array[String] | cel2chp                | list of chip types to check library and CEL files against               |
| tags                 | Array[String] | {gtc,chp,txt}2vcf      | list of FORMAT tags to output in the VCF [\"GT,BAF,LRR\"]               |
| gc_window_size       | Int           | {gtc,chp,txt}2vcf      | window size in bp used to compute the GC content [200]                  |

When using older Illumina arrays IAAP CLI might not work and furthermore the BPM sanity check for consistency between BPM files and GTC files might also fail. In these cases, it is advised to turn the two flags **autoconvert** and **do_not_use_reference** to true. Similarly, when using older Affymetrix arrays, the chip type information in the CEL files might not reflect the chip type information in the XML file. You can include in the **chip_type** array the chip type labels included in the CEL files

Outputs
=======

After a pipeline run, assuming the **target** variable is set to the default **pngs** value, the following outputs will be generated:

| key                  | type         | idat | gtc  | cel  | chp  | txt  | vcf  | pvcf | description                                         |
|----------------------|--------------|------|------|------|------|------|------|------|-----------------------------------------------------|
| green_idat_tsv_file  | File?        | yes  |      |      |      |      |      |      | summary table for green IDAT files                  |
| red_idat_tsv_file    | File?        | yes  |      |      |      |      |      |      | summary table for red IDAT files                    |
| gtc_tsv_file         | File?        | yes  | yes  |      |      |      |      |      | summary table for GTC files                         |
| cel_tsv_file         | File?        |      |      | yes  |      |      |      |      | summary table for CEL files                         |
| affy_tsv_file        | File?        |      |      | yes  | yes  | yes  |      |      | summary table for CHP or TXT files                  |
| sample_id_file       | File         | yes  | yes  | yes  | yes  | yes  | yes  | yes  | sample ID list                                      |
| computed_gender_file | File         | yes  | yes  | yes  | yes  | yes  | yes  | yes  | sample computed gender list                         |
| call_rate_file       | File?        | yes  | yes  | yes  | yes  | yes  | yes  |      | sample call rate list                               |
| mocha_stats_file     | File?        | yes  | yes  | yes  | yes  | yes  | yes  | yes  | summary table for samples processed with MoChA      |
| mocha_calls_file     | File?        | yes  | yes  | yes  | yes  | yes  | yes  | yes  | table of mosaic chromosomal alterations calls       |
| mocha_ucsc_bed       | File?        | yes  | yes  | yes  | yes  | yes  | yes  | yes  | bed file with mosaic chromosomal alterations calls  |
| mocha_summary_pdf    | File?        | yes  | yes  | yes  | yes  | yes  | yes  | yes  | summary of MoChA run across all batches             |
| mocha_pileup_pdf     | File?        | yes  | yes  | yes  | yes  | yes  | yes  | yes  | pileup of calls passing quality control filters     |
| png_files            | Array[File]? | yes  | yes  | yes  | yes  | yes  | yes  | yes  | single call plots for all calls passing filters     |
| mendel_files         | Array[File]? | opt. | opt. | opt. | opt. | opt. | opt. |      | Mendelian summary tables if ped_file was provided   |
| gtc_files            | Array[File]? | yes  |      |      |      |      |      |      | Illumina GTC files (if **gtc_output** is true)      |
| chp_files            | Array[File]? |      |      | yes  |      |      |      |      | Affymetrix CHP files (if **chp_output** is true)    |
| snp_files            | Array[File]? |      |      | yes  |      |      |      |      | Affymetrix SNP posteriors files (if **chp_output**) |
| vcf_files            | Array[File]  | yes  | yes  | yes  | yes  | yes  | yes  | yes  | output VCF files                                    |
| vcf_idxs             | Array[File]  | yes  | yes  | yes  | yes  | yes  | yes  | yes  | output VCF indexes                                  |
| pgt_vcf_files        | Array[File]? | yes  | yes  | yes  | yes  | yes  | yes  |      | output VCF files with phased GTs only               |
| pgt_vcf_idxs         | Array[File]? | yes  | yes  | yes  | yes  | yes  | yes  |      | indexes for output VCF files with phased GTs only   |
| xcl_vcf_file         | File?        | yes  | yes  | yes  | yes  | yes  | yes  |      | VCF of variants excluded from analysis              |
| xcl_vcf_idx          | File?        | yes  | yes  | yes  | yes  | yes  | yes  |      | index for VCF of variants excluded from analysis    |
| output_tsv_file      | File         | yes  | yes  | yes  | yes  | yes  | yes  | yes  | table with output VCF files                         |

Illumina Example
================

The following example consists of <a href="http://doi.org/doi:10.18129/B9.bioc.hapmap370k">40 HapMap samples</a> genotyped in 2007 with the Illumina HumanCNV370-Duov1 array

Download manifest data:
```
wget ftp://webdata:webdata@ussd-ftp.illumina.com/downloads/ProductFiles/HumanCNV370/HumanCNV370-Duo/humancnv370v1_c.bpm
wget ftp://webdata2:webdata2@ussd-ftp.illumina.com/downloads/ProductFiles/HumanCNV370/HumanCNV370-Duo/HumanCNV370v1_C.egt
wget https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL6nnn/GPL6986/suppl/GPL6986_HumanCNV370v1_C.csv.gz
gunzip GPL6986_HumanCNV370v1_C.csv.gz
/bin/mv GPL6986_HumanCNV370v1_C.csv HumanCNV370v1_C.csv
```

Download IDAT files and sample trackers:
```
wget http://bioconductor.org/packages/release/data/annotation/src/contrib/hapmap370k_1.0.1.tar.gz
tar xzvf hapmap370k_1.0.1.tar.gz --strip-components=3 hapmap370k/inst/idatFiles
```
Then make sure you copy all the IDAT files in the `gs://{google-bucket}/idats` Google bucket

Define options to run the WDL:
```json
{
  "mocha.sample_set_id": "hapmap370k",
  "mocha.mode": "idat",
  "mocha.realign": true,
  "mocha.max_win_size_cm": 300.0,
  "mocha.overlap_size_cm": 5.0,
  "mocha.ref_name": "GRCh38",
  "mocha.ref_path": "gs://{google-bucket}/GRCh38",
  "mocha.manifest_path": "gs://{google-bucket}/manifests",
  "mocha.data_path": "gs://{google-bucket}/idats",
  "mocha.batch_tsv_file": "gs://{google-bucket}/tsvs/hapmap370k.batch.tsv",
  "mocha.sample_tsv_file": "gs://{google-bucket}/tsvs/hapmap370k.sample.tsv",
  "mocha.ped_file": "gs://{google-bucket}/hapmap370k.ped",
  "mocha.docker_registry": "us.gcr.io/mccarroll-mocha",
  "mocha.do_not_check_bpm": true
}
```
Notice that for this example we need the secondary boolean flag **do_not_check_bpm** checked as true due to the Illumina manifest being in an old format

The **hapmap370k.batch.tsv** table could look like this:

| batch_id | csv                 | bpm                 | egt                 |
|----------|---------------------|---------------------|---------------------|
| A        | HumanCNV370v1_C.csv | humancnv370v1_c.bpm | HumanCNV370v1_C.egt |
| B        | HumanCNV370v1_C.csv | humancnv370v1_c.bpm | HumanCNV370v1_C.egt |

The **hapmap370k.sample.tsv** table could look like this:

| sample_id | batch_id | green_idat            | red_idat              |
|-----------|----------|-----------------------|-----------------------|
| NA07029   | A        | 4019585367_A_Grn.idat | 4019585367_A_Red.idat |
| NA18502   | B        | 4019585376_B_Grn.idat | 4019585376_B_Red.idat |
| NA10851   | A        | 4019585401_B_Grn.idat | 4019585401_B_Red.idat |
| NA18912   | B        | 4019585402_A_Grn.idat | 4019585402_A_Red.idat |
| NA18500   | A        | 4019585413_A_Grn.idat | 4019585413_A_Red.idat |
| NA10859   | B        | 4019585415_B_Grn.idat | 4019585415_B_Red.idat |
| NA12057   | A        | 4019585422_A_Grn.idat | 4019585422_A_Red.idat |
| NA07000   | B        | 4019585433_A_Grn.idat | 4019585433_A_Red.idat |
| NA07034   | A        | 4019585455_A_Grn.idat | 4019585455_A_Red.idat |
| NA18501   | B        | 4019585483_B_Grn.idat | 4019585483_B_Red.idat |
| NA18912-1 | A        | 4019585498_B_Grn.idat | 4019585498_B_Red.idat |
| NA12057-1 | B        | 4019585506_B_Grn.idat | 4019585506_B_Red.idat |
| NA12056   | A        | 4019585508_B_Grn.idat | 4019585508_B_Red.idat |
| NA06991   | B        | 4019585512_B_Grn.idat | 4019585512_B_Red.idat |
| NA06994   | A        | 4019585520_B_Grn.idat | 4019585520_B_Red.idat |
| NA18501-1 | B        | 4019585575_A_Grn.idat | 4019585575_A_Red.idat |
| NA06985   | A        | 4019585596_A_Grn.idat | 4019585596_A_Red.idat |
| NA18502-1 | B        | 4019585597_A_Grn.idat | 4019585597_A_Red.idat |
| NA18858   | A        | 4030186100_B_Grn.idat | 4030186100_B_Red.idat |
| NA18505   | B        | 4030186109_A_Grn.idat | 4030186109_A_Red.idat |
| NA11882   | A        | 4030186125_A_Grn.idat | 4030186125_A_Red.idat |
| NA11994   | B        | 4030186132_A_Grn.idat | 4030186132_A_Red.idat |
| NA18504   | A        | 4030186167_B_Grn.idat | 4030186167_B_Red.idat |
| NA18503   | B        | 4030186197_A_Grn.idat | 4030186197_A_Red.idat |
| NA12155   | A        | 4030186219_B_Grn.idat | 4030186219_B_Red.idat |
| NA11993   | B        | 4030186254_B_Grn.idat | 4030186254_B_Red.idat |
| NA07000-1 | A        | 4030186263_B_Grn.idat | 4030186263_B_Red.idat |
| NA11993-1 | B        | 4030186332_A_Grn.idat | 4030186332_A_Red.idat |
| NA12156   | A        | 4030186339_A_Grn.idat | 4030186339_A_Red.idat |
| NA06991-1 | B        | 4030186347_A_Grn.idat | 4030186347_A_Red.idat |
| NA10860   | A        | 4030186396_B_Grn.idat | 4030186396_B_Red.idat |
| NA10861   | B        | 4030186415_A_Grn.idat | 4030186415_A_Red.idat |
| NA11995   | A        | 4030186434_B_Grn.idat | 4030186434_B_Red.idat |
| NA10859-1 | B        | 4030186513_A_Grn.idat | 4030186513_A_Red.idat |
| NA11881   | A        | 4031058010_B_Grn.idat | 4031058010_B_Red.idat |
| NA10851-1 | B        | 4031058082_A_Grn.idat | 4031058082_A_Red.idat |
| NA18506   | A        | 4031058102_B_Grn.idat | 4031058102_B_Red.idat |
| NA11882-1 | B        | 4031058127_B_Grn.idat | 4031058127_B_Red.idat |
| NA06993   | A        | 4031058132_A_Grn.idat | 4031058132_A_Red.idat |
| NA06993-1 | B        | 4031058211_B_Grn.idat | 4031058211_B_Red.idat |

And the **hapmap370k.ped** file could look like this (only the first four columns are required):

|      |         |           |           |   |   |     |
|------|---------|-----------|-----------|---|---|-----|
| 1340 | NA07029 | NA06994   | NA07000   | 1 | 0 | CEU |
| 1340 | NA07029 | NA06994   | NA07000-1 | 1 | 0 | CEU |
| 1341 | NA07048 | NA07034   | NA07055   | 1 | 0 | CEU |
| 1341 | NA06991 | NA06993   | NA06985   | 2 | 0 | CEU |
| 1341 | NA06991 | NA06993-1 | NA06985   | 2 | 0 | CEU |
| 1344 | NA10851 | NA12056   | NA12057   | 1 | 0 | CEU |
| 1344 | NA10851 | NA12056   | NA12057-1 | 1 | 0 | CEU |
| 1347 | NA10859 | NA11881   | NA11882   | 2 | 0 | CEU |
| 1347 | NA10859 | NA11881   | NA11882-1 | 2 | 0 | CEU |
| 1362 | NA10860 | NA11992   | NA11993   | 1 | 0 | CEU |
| 1362 | NA10860 | NA11992   | NA11993-1 | 1 | 0 | CEU |
| 1362 | NA10861 | NA11994   | NA11995   | 2 | 0 | CEU |
| 1408 | NA10831 | NA12155   | NA12156   | 2 | 0 | CEU |
| Y004 | NA18500 | NA18501   | NA18502   | 1 | 0 | YRI |
| Y004 | NA18500 | NA18501-1 | NA18502   | 1 | 0 | YRI |
| Y004 | NA18500 | NA18501   | NA18502-1 | 1 | 0 | YRI |
| Y004 | NA18500 | NA18501-1 | NA18502-1 | 1 | 0 | YRI |
| Y005 | NA18503 | NA18504   | NA18505   | 1 | 0 | YRI |
| Y012 | NA18860 | NA18859   | NA18858   | 1 | 0 | YRI |
| Y028 | NA18914 | NA18913   | NA18912   | 1 | 0 | YRI |
| Y028 | NA18914 | NA18913   | NA18912-1 | 1 | 0 | YRI |

Affymetrix Example
==================

The following example consists of <a href="http://doi.org/10.1038/ng.238">270 HapMap samples</a> including 60 trios genotyped in 2007 with the Affymetrix Genome-Wide Human SNP Array 6.0

Download manifest data:
```
wget http://www.affymetrix.com/Auth/support/downloads/library_files/genomewidesnp6_libraryfile.zip
wget http://www.affymetrix.com/Auth/analysis/downloads/lf/genotyping/GenomeWideSNP_6/SNP6_supplemental_axiom_analysis_files.zip
wget http://www.affymetrix.com/Auth/analysis/downloads/na35/genotyping/GenomeWideSNP_6.na35.annot.csv.zip
unzip -oj genomewidesnp6_libraryfile.zip CD_GenomeWideSNP_6_rev3/Full/GenomeWideSNP_6/LibFiles/GenomeWideSNP_6.{cdf,chrXprobes,chrYprobes,specialSNPs}
unzip -o SNP6_supplemental_axiom_analysis_files.zip GenomeWideSNP_6.{generic_prior.txt,apt-probeset-genotype.AxiomGT1.xml,AxiomGT1.sketch}
mv GenomeWideSNP_6.apt-probeset-genotype.AxiomGT1.xml GenomeWideSNP_6.AxiomGT1.xml
unzip -o GenomeWideSNP_6.na35.annot.csv.zip GenomeWideSNP_6.na35.annot.csv
sed -i 's/chrX-probes/special-snps"       currentValue="GenomeWideSNP_6.specialSNPs" \/>\n        <Parameter name="chrX-probes/' GenomeWideSNP_6.AxiomGT1.xml
sed -i 's/target-sketch/read-models-brlmmp" currentValue="GenomeWideSNP_6.generic_prior.txt" \/>\n        <Parameter name="target-sketch/' GenomeWideSNP_6.AxiomGT1.xml
zip GenomeWideSNP_6.AxiomGT1.zip GenomeWideSNP_6.{cdf,chrXprobes,chrYprobes,specialSNPs,generic_prior.txt,AxiomGT1.sketch}
rm GenomeWideSNP_6.{cdf,chrXprobes,chrYprobes,specialSNPs,generic_prior.txt,AxiomGT1.sketch}
```

Download CEL files:
```
wget ftp://ftp.ncbi.nlm.nih.gov/hapmap/raw_data/hapmap3_affy6.0/{Broad_hapmap3_r2_Affy6_cels_excluded,SCALE,GIGAS,SHELF}.tgz
tar xzvf Broad_hapmap3_r2_Affy6_cels_excluded.tgz --strip-components=1 \
  Broad_hapmap3_r2_Affy6_cels_excluded/SCALE_g_GAINmixHapMapAffy1_GenomeWideEx_6_{H11_30996,D01_30828,F06_30912,C08_30938,H03_30868,D12_31004,F02_30848,B08_30936,B02_30840,C10_30970}.CEL \
  Broad_hapmap3_r2_Affy6_cels_excluded/GIGAS_g_GAINmixHapMapAffy2_GenomeWideEx_6_D09_31352.CEL \
  Broad_hapmap3_r2_Affy6_cels_excluded/SHELF_g_GAINmixHapMapAffy3_GenomeWideEx_6_{B05_31476,A01_31410,E05_31482,E04_31466,G11_31582}.CEL
echo {SCALE,GIGAS,SHELF} | tr ' ' '\n' | xargs -i tar xzvf {}.tgz
```
Then make sure you copy all the CEL files in the `gs://{google-bucket}/cels` Google bucket

Define options to run the WDL:
```json
{
  "mocha.sample_set_id": "hapmapSNP6",
  "mocha.mode": "cel",
  "mocha.realign": true,
  "mocha.max_win_size_cm": 100.0,
  "mocha.overlap_size_cm": 5.0,
  "mocha.ref_name": "GRCh38",
  "mocha.ref_path": "gs://{google-bucket}/GRCh38",
  "mocha.manifest_path": "gs://{google-bucket}/manifests",
  "mocha.data_path": "gs://{google-bucket}/cels",
  "mocha.batch_tsv_file": "gs://{google-bucket}/tsvs/hapmapSNP6.batch.tsv",
  "mocha.sample_tsv_file": "gs://{google-bucket}/tsvs/hapmapSNP6.sample.tsv",
  "mocha.ped_file": "gs://{google-bucket}/hapmapSNP6.ped",
  "mocha.docker_registry": "us.gcr.io/mccarroll-mocha",
  "mocha.chip_type": ["GenomeWideEx_6"]
}
```
Notice that for this example we needed the **chip_type** secondary option set due to the CEL files containing the chip name **GenomeWideEx_6** rather than **GenomeWideSNP_6**

The **hapmapSNP6.batch.tsv** table could look like this:

| batch_id | csv                            | xml                          | zip                          |
|----------|--------------------------------|------------------------------|------------------------------|
| SCALE    | GenomeWideSNP_6.na35.annot.csv | GenomeWideSNP_6.AxiomGT1.xml | GenomeWideSNP_6.AxiomGT1.zip |
| GIGAS    | GenomeWideSNP_6.na35.annot.csv | GenomeWideSNP_6.AxiomGT1.xml | GenomeWideSNP_6.AxiomGT1.zip |
| SHELF    | GenomeWideSNP_6.na35.annot.csv | GenomeWideSNP_6.AxiomGT1.xml | GenomeWideSNP_6.AxiomGT1.zip |

The **hapmapSNP6.sample.tsv** and **hapmapSNP6.ped** files could be generated with the following commands:
```
wget ftp://ftp.ncbi.nlm.nih.gov/hapmap/raw_data/hapmap3_affy6.0/{passing,excluded}_cels_sample_map.txt
cat {passing,excluded}_cels_sample_map.txt | grep "SCALE\|GIGAS\|SHELF" | \
  awk -F"\t" -v OFS="\t" 'BEGIN {print "sample_id","batch_id","cel"}
  {sm=$1; if (sm in x) sm=sm"-"x[sm]; print sm,substr($2,1,5),$2; x[$1]++}' > hapmapSNP6.sample.tsv
wget ftp://ftp.ncbi.nlm.nih.gov/hapmap/phase_3/relationships_w_pops_051208.txt
awk -F"\t" -v OFS="\t" 'NR==FNR {x[$1]++} NR>FNR && $3 in x || $4 in x {print;
  if ($3"-1" in x) {$3=$3"-1"; print} if ($4"-1" in x) {$4=$4"-1"; print}}' \
  hapmapSNP6.batch.tsv relationships_w_pops_051208.txt > hapmapSNP6.ped
```

Filtering Output Calls
======================

Once you have completed a full run of the MoChA pipeline, you will receive two tables, **mocha_stats_tsv** and **mocha_calls_tsv**, the first one with overall statistics about the samples, and the second one with the list of chromosomal alterations calls made by the MoChA software. The columns and suggestions for filtering in the table are described in the software <a href="../mocha#chromosomal-alterations-pipeline">documentation</a>

Not all calls made by MoChA will be mosaic chromosomal alterations. These are some suggestions for how to think about filtering the calls:

* Calls with **type** CNP are inherited copy number polymorphism. To clean the data, MoChA attempts to call common deletions and duplications during the first pass through the data. These can be safely removed from your callset
* Calls with **lod_baf_phase** less than 10 are often constitutional rather than mosaic events. Constitutional events are more likely to be germline rather than somatic
* Calls with **length** less than 500 kbp and **rel_cov** greater than 2.5 are often germline duplications that can have a strong allelic imbalance signal. It might be difficult to distinguish these alterations from genuine somatic amplifications but to edge on the conservative side we recommend to filter these out unless in loci with strong prior evidence for mosaicism
* Samples that have low **call_rate** (by default Illumina and Affymetrix recommend less than 0.97) or high **baf_auto** values (by default we recommend a 0.03 threshold, although the correct threshold could depend on DNA microarray type) could be low quality or contaminated samples with widespread allelic imbalance across the genome. Plotting the values for these two statistics can help identify outliers and select the correct thresholds

After filtering out your callset, plotting the prevalence of the number of samples with mosaic chromosomal alterations as a function of age might provide the best evidence for an acceptable filtering strategy if young individuals have little to no evidence for mosaicism. Notice that there is not one unique strategy that works across all datasets. Chromosome X and 15q11 might require special care. If your samples are from cell lines for example, the inactive X and the paternal region in 15q11 tend to replicate later than their active counterparts due to epigenetic mechanisms and this can cause both regions to exhibit significant amounts of allelic imbalance that is of replication timing origin

Scattered Tasks
===============

The workflow consists of the following scattered tasks:

| task                | description                                                              |
|---------------------|--------------------------------------------------------------------------|
| csv2bam             | realigns the sequences in the manifest files to the reference genome     |
| idat2gtc/cel2chp    | converts raw intensity IDAT/CEL files to GTC/CHP files with genotypes    |
| {gtc,chp,txt}2vcf   | imports genotypes and intensities to a compliant set of VCF files        |
| {gtc,chp}2vcf_merge | merges previous outputs across sub bactches to make single batch VCFs    |
| vcf_scatter         | extracts genotypes and scatter the VCF files across the genome windows   |
| vcf_merge           | joins the VCF files into VCF slices including all samples across batches |
| vcf_qc              | identifies a subset of variants that performed poorly                    |
| vcf_phase           | phases each genome window across the whole cohort                        |
| vcf_split           | splits the phased VCF slices back into VCF shards                        |
| vcf_concat          | concatenates the VCF shards back into VCF files for each batch           |
| vcf_import          | imports the phased genotypes back into the original VCFs                 |
| vcf_mocha           | runs the MoChA framework across batches                                  |
| mocha_plot          | plots the most confident events called by MoChA across batches           |

![](mocha_scatter.png)

To make best use of the elasticity of cloud computing, the pipeline attempts to achieve parallelization in four different ways:

* by splitting the conversion from IDAT to GTC in small sub batches (controlled by **idat_batch_size**)
* by splitting the conversion from GTC to VCF in small sub batches (controlled by **gtc_batch_size**)
* by splitting the conversion from CHP to VCF in small sub batches (controlled by **chp_batch_size**)
* by splitting the genome in overlapping windows each phased independently (controlled by **max_win_size_cm** and **overlap_size_cm**)

All tasks but the phasing tasks are set by default to run on a single CPU. All tasks but cel2chp and txt2vcf are set by default to run on preemptible computing on the first run, and to then proceed to run as non-preemptible is they are preempted

Terra has a 2,300 task limit per project (called \"Hog Limits\") but a single job should work well within this limit. Tasks idat2gtc/cel2chp and tasks {gtc,chp,txt}2vcf will output tables with metadata about their input files which will be aggregated across batches. From personal experience, even running a small cohort takes a minimum of ~1 hour on Terra, as there are enough housekeeping tasks that need to be run

Imputation pipeline
===================

Once you have run the MoChA pipeline you will obtain as an output an array of `pgt_vcf_files` VCF files (one per batch). These files can be used as a basis for imputation and then to extend the mosaic chromosomal calls to imputed variants to build a map of which alleles are over-represented and under-represented. Imputation can be run with <a href="https://jmarchini.org/impute5/">IMPUTE5</a> or with <a href="http://faculty.washington.edu/browning/beagle/beagle.html">Beagle5</a>, and both require as input VCF files without missing or unphased genotypes. Notice that, while Beagle5 is completely open source and released under the GPLv3 license, the source code for IMPUTE5 is not provided and the software is intended only for academics carrying out research and not for use by consumers or commercial business. If you are interested in using this imputation pipeline commercially, you can start using Beagle5 by setting the corresponding boolean flag in the pipeline and then negotiate a license for IMPUTE5 if the cost savings apply to you. We have observed IMPUTE5 to be 2-3x cheaper to run than Beagle5

The following are the primary options that you can set in the main input json file for use with <a href="impute.wdl">impute.wdl</a>

| key               | type           | description                                                                                            |
|-------------------|----------------|--------------------------------------------------------------------------------------------------------|
| sample_set_id     | String         | cohort name that will be used as prefix for temporary and output files                                 |
| mode              | String?        | pipeline mode, one of: **pgt** or **imp** [**pgt**]                                                    |
| target            | String?        | pipeline final state, one of: **imp** or **ext** [**pngs**]                                            |
| max_win_size_cm   | Float?         | maximum windows size in cM for phasing [10.0]                                                          |
| overlap_size_cm   | Float?         | required overlap size in cM for consecutive windows [2.0]                                              |
| format_id         | String?        | format ID to extend to imputed haplotypes [AS]                                                         |
| ext_string        | String?        | extension string for the imputed VCF with the extended format [as]                                     |
| target_chrs       | Array[String]? | list of chromosomes to impute                                                                          |
| ref_name          | String?        | name of reference genome, with resource default files for GRCh37 and GRCh38 [GRCh38]                   |
| ref_path          | String?        | path for reference genome resources (needed unless all resources are provided with full path)          |
| ref_fasta_fai     | String?        | reference sequence index [GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai/human_g1k_v37.fasta.fai] |
| ref_n_chrs        | Int?           | number of chromosomes [23]                                                                             |
| genetic_map_file  | String?        | genetic map [genetic_map_hg38_withX.txt.gz/genetic_map_hg19_withX.txt.gz]                              |
| panel_pfx         | String?        | prefix for phasing reference [CCDG_14151_B01_GRM_WGS_2020-08-05_/ALL.chr]                              |
| panel_sfx         | String?        | Suffix for phasing reference [.filtered.phased.bcf/.phase3_integrated.20130502.genotypes.bcf]          |
| panel_idx         | String?        | index extension for phasing reference [.bcf]                                                           |
| n_panel_smpls     | Int?           | number of samples in phasing reference [3202/2504]                                                     |
| data_path         | String?        | path for data files (overrides **path** column in **batch_tsv_file**)                                  |
| batch_tsv_file    | File           | TSV file with batch information                                                                        |
| convert_panel     | Boolean?       | whether to convert the reference panel to IMP5 format [true]                                           |
| beagle            | Boolean?       | whether to run Beagle5 rather than IMPUTE5 [false]                                                     |
| out_ds            | Boolean?       | whether imputation VCFs should contain the FORMAT/DS field (Genotype dosages) [true]                   |
| out_gp            | Boolean?       | whether imputation VCFs should contain the FORMAT/GP field (Genotype probabilities) [false]            |
| out_ap            | Boolean?       | whether imputation VCFs should contain the FORMAT/AP field (ALT haplotype probabilities) [false]       |
| impute_extra_args | String?        | extra arguments for IMPUTE5/Beagle5                                                                    |
| basic_bash_docker | String?        | docker to run basic bash scripts [ubuntu:latest]                                                       |
| pandas_docker     | String?        | docker to run task ref_scatter [amancevice/pandas:slim]                                                |
| docker_registry   | String?        | location of docker images [us.gcr.io/mccarroll-mocha]                                                  |
| bcftools_docker   | String?        | docker to run tasks requiring BCFtools [bcftools:1.11-yyyymmdd]                                        |
| impute5_docker    | String?        | docker to run tasks requiring IMPUTE5 [impute5:1.11-yyyymmdd]                                          |
| beagle5_docker    | String?        | docker to run tasks requiring Beagle5 [beagle5:1.11-yyyymmdd]                                          |

Make sure all the `pgt_vcf_files` from the MoChA pipeline are first available in the `gs://{google-bucket}/vcfs` directory, together with the `vcf_files` including the MoChA calls. We strongly advise to use **max_win_size_cm** smaller than **30.0** as imputation of large chromosome windows can be prohibitevely memory intensive due to the high density of variants in the reference panels

Define options to run the WDL:
```json
{
  "impute.sample_set_id": "hapmap370k",
  "impute.mode": "pgt",
  "impute.target": "ext",
  "impute.batch_tsv_file": "gs://{google-bucket}/tsvs/hapmap370k.batch.tsv",
  "impute.max_win_size_cm": 50.0,
  "impute.overlap_size_cm": 5.0,
  "impute.target_chrs": ["chr12", "chrX"],
  "impute.ref_name": "GRCh38",
  "impute.ref_path": "gs://{google-bucket}/GRCh38",
  "impute.data_path": "gs://{google-bucket}/vcfs",
  "impute.beagle": false,
  "mocha.docker_registry": "us.gcr.io/mccarroll-mocha"
}
```

The `batch_tsv_file` input table should look like this:

| batch_id | n_smpls | vcf                   |  vcf_index                | pgt_vcf              | pgt_vcf_index            |
|----------|---------|-----------------------|---------------------------|----------------------|--------------------------|
| A        | 20      | hapmap370k.A.as.bcf   | hapmap370k.A.as.bcf.csi   | hapmap370k.A.pgt.bcf | hapmap370k.A.pgt.bcf.csi |
| B        | 20      | hapmap370k.B.as.bcf   | hapmap370k.B.as.bcf.csi   | hapmap370k.B.pgt.bcf | hapmap370k.B.pgt.bcf.csi |

The `output_tsv_file` from the MoChA pipeline will include these columns

Notice that the imputation pipeline acts independently on each batch and on each chromosome so that you do not have to run all batches and chromosomes at once in one job. Running imputation as multiple jobs will not affect the final result. You can control which batches are imputed through the `batch_tsv_file` input table and you can control which chromosomes are imputed through the `target_chrs` string array

The workflow consists of the following scattered tasks:

| task                | description                                                                  | run if and only if    |
|---------------------|------------------------------------------------------------------------------|-----------------------|
| get_panel_markers   | retrieve size of the reference panel from the VCF index files                | **mode** is **pgt**   |
| panel_scatter       | splits reference panel chromosome VCF files into required windows            | **mode** is **pgt**   |
| init panel          | generates site only reference panel VCF and convert to IMPUTE5/Beagle format | **mode** is **pgt**   |
| vcf_scatter         | splits phased genotype batches into required windows                         | **mode** is **pgt**   |
| vcf_impute          | imputes haplotypes across windows and batches                                | **mode** is **pgt**   |
| vcf_ligate          | concatenate imputed haplotypes into single chromosome VCFs across batches    | **mode** is **pgt**   |
| get_markers         | retrieve size of each imputed VCF from the VCF index files                   | **target** is **ext** |
| vcf_extend          | extend format field across imputed haplotypes                                | **target** is **ext** |

If you already have been provided with imputation VCFs, you might prefer to avoid running imputation again. In this case, it is possible to run only the extension part of the pipeline which will extend allelic shift information for MoChA calls from genotyped to imputed heterozygous sites. However, the pipeline expects the data to be split by chromosomes and the same batches used by the MoChA pipeline. If you have been provided with chromosome imputed VCFs for the whole cohort, you will have manually split these across batches. You might also want to only work with genotypes rather than dosages and genotype probabilities to avoid generating redundant files. If you have a `sample_tsv_file`, you can use it as follows to split a VCF in batches:

```
for chr in {1..22} X; do
  awk -F "\t" -v OFS="\t" -v chr=$chr 'NR==1 {for (i=1; i<=NF; i++) f[$i] = i}
    NR>1 {print $(f["sample_id"]),"-","hapmap370k."$(f["batch_id"])".chr"chr}' hapmap370k.sample.tsv > samples_file.txt
  bcftools annotate --no-version -Ou -x QUAL,FILTER,INFO,^FMT/GT hapmap370k.chr$chr.vcf.gz | bcftools +split -o . -Ob -G samples_file.txt
  cut -f3 samples_file.txt | sort | uniq | xargs -i bcftools index -f {}.bcf
done
```
Make sure you are using BCFtools version 1.11 or newer, as previous versions are missing the required splitting functionality. Notice that the above command will also strip the imputed VCF of all attributes except ID and FMT/GT, as this will allow the allelic shift pipeline to run faster

Allelic Shift Pipeline
======================

The following are the primary options that you can set in the main input json file for use with <a href="shift.wdl">shift.wdl</a>

| key               | type     | description                                                                                   |
|-------------------|----------|-----------------------------------------------------------------------------------------------|
| sample_set_id     | String   | cohort name that will be used as prefix for temporary and output files                        |
| region            | String   | region to extract for allelic shift analysis                                                  |
| as_id             | String?  | allelic shift ID [AS]                                                                         |
| ext_string        | String   | extension string for the imputed VCF with the extended format [as]                            |
| ref_path          | String?  | path for reference genome resources (needed unless all resources are provided with full path) |
| cyto_file         | String?  | file with location of cytoband regions [cytoBand.txt.gz]                                      |
| data_path         | String?  | path for data files (overrides **path** column in **batch_tsv_file**)                         |
| batch_tsv_file    | File     | TSV file with batch information                                                               |
| samples_file      | File?    | list of samples to include in the allelic shift analysis                                      |
| drop_genotypes    | Boolean? | whether to drop genotypes after computing allelic shift counts [true]                         |
| phred_score       | Boolean? | whether the binomial p-values should be recorded as phred scores [true]                       |
| plot              | Boolean? | whether to generate a summary plot [true]                                                     |
| basic_bash_docker | String?  | docker to run basic bash scripts [ubuntu:latest]                                              |
| docker_registry   | String?  | location of docker images [us.gcr.io/mccarroll-mocha]                                         |
| bcftools_docker   | String?  | docker to run tasks requiring BCFtools [bcftools:1.11-yyyymmdd]                               |
| r_mocha_docker    | String?  | docker to run tasks mocha_{plot,summary} [r_mocha:1.11-yyyymmdd]                              |

If you wanted to study allelic shift for mosaic loss of chromosomes X (mLOX), you could define individuals with putative mLOX as follows:
```
awk -F "\t" 'NR==FNR && FNR==1 {for (i=1; i<=NF; i++) f[$i] = i}
  NR==FNR && FNR>1 {sample_id=$(f["sample_id"]); call_rate=$(f["call_rate"]); baf_auto=$(f["baf_auto"])}
  NR==FNR && FNR>1 && (call_rate<.97 || baf_auto>.03) {xcl[sample_id]++}
  NR>FNR && FNR==1 {for (i=1; i<=NF; i++) g[$i] = i}
  NR>FNR && FNR>1 {sample_id=$(g["sample_id"]); computed_gender=$(g["computed_gender"]);
  chrom=$(g["chrom"]); sub("^chr", "", chrom); len=$(g["length"]); rel_cov=$(g["rel_cov"])}
  NR>FNR && FNR>1 && !(sample_id in xcl) && computed_gender=="F" &&
  chrom=="X" && len>1e8 && rel_cov<2.5 {print sample_id}' hapmap370k.stats.tsv hapmap370k.calls.tsv > hapmap370k.mLOX.lines
```

Define options to run the WDL:
```json
{
  "shift.sample_set_id": "hapmap370k",
  "shift.region": "chrX",
  "shift.samples_file": "gs://{google-bucket}/tsvs/hapmap370k.mLOX.lines",
  "shift.batch_tsv_file": "gs://{google-bucket}/tsvs/hapmap370k.batch.tsv",
  "shift.ref_path": "gs://mccarroll-mocha/GRCh38",
  "shift.data_path": "gs://{google-bucket}/imp_vcfs",
  "mocha.docker_registry": "us.gcr.io/mccarroll-mocha"
}
```

The `batch_tsv_file` input table should look like this:

| batch_id | chr1_imp_vcf             | chr1_imp_vcf_index           | ... | chrX_imp_vcf             | chrX_imp_vcf_index           |
|----------|--------------------------|------------------------------|-----|--------------------------|------------------------------|
| A        | hapmap370k.A.chr1.as.bcf | hapmap370k.A.chr1.as.bcf.csi | ... | hapmap370k.A.chrX.as.bcf | hapmap370k.A.chrX.as.bcf.csi |
| B        | hapmap370k.B.chr1.as.bcf | hapmap370k.B.chr1.as.bcf.csi | ... | hapmap370k.B.chrX.as.bcf | hapmap370k.B.chrX.as.bcf.csi |

The `output_tsv_file` from the imputation pipeline will include these columns

Once successfully run, the allelic shift pipeline will output a VCF named like `hapmap370k.chrX.as.bcf` with an `AS` INFO filed containing allelic shift counts. If you have two (or more) such VCFs (e.g. `hapmapSNP6.chrX.as.bcf`), these can be easily combined with the following command:
```
bcftools merge --no-version -Ou -i AS:sum -m none hapmap370k.chrX.as.bcf hapmapSNP6.chrX.as.bcf | \
  bcftools +mochatools --no-version -Ou -- --test AS | \
  bcftools query -i 'binom(INFO/AS)<1e-6' -f "%CHROM\t%POS\t%ID\t%REF\t%ALT\t%AS{0}\t%AS{1}\t%binom_AS\n"
```

Polygenic Score Pipeline
========================

The following are the primary options that you can set in the main input json file for use with <a href="score.wdl">score.wdl</a>

| key               | type           | description                                                                                            |
|-------------------|----------------|--------------------------------------------------------------------------------------------------------|
| sample_set_id     | String         | cohort name that will be used as prefix for temporary and output files                                 |
| vcf               | Boolean?       | whether the input summary statistics files are in VCF format [false]                                   |
| sample_header     | String?        | what column header to use for the column with samples IDs [sample_id]                                  |
| region            | String?        | region to extract for polygenic score analysis                                                         |
| tag               | String?        | VCF format tag to use for polygenic score computation among GP, AP, HDS, DS, GT, and AS                |
| ext_string        | String?        | extension string for the output polygenic score table [scores]                                         |
| summary_path      | String?        | path for summary statistics files (unless all files are provided with full path)                       |
| summary_files     | Array[String]  | summary statistics files                                                                               |
| summary_idxs      | Array[String]? | summary statistics files indexes (if files provided as VCF)                                            |
| q_score_thr       | Array[Float]?  | list of p-value thresholds                                                                             |
| covars_file       | File?          | covariate file to be used to computed adjusted polygenic scores                                        |
| ref_name          | String?        | name of reference genome, with resource default files for GRCh37 and GRCh38 [GRCh38]                   |
| ref_path          | String?        | path for reference genome resources (needed unless all resources are provided with full path)          |
| ref_fasta_fai     | String?        | reference sequence index [GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai/human_g1k_v37.fasta.fai] |
| ref_n_chrs        | Int?           | number of chromosomes [23]                                                                             |
| data_path         | String?        | path for data files (overrides **path** column in **batch_tsv_file**)                                  |
| batch_tsv_file    | File           | TSV file with batch information                                                                        |
| samples_file      | File?          | list of samples to include in the polygenic score analysis                                             |
| exclude_str       | String?        | exclusion criterias for variants (e.g. [INFO<0.8]) not to be used with include_str                     |
| include_str       | String?        | inclusion criterias for variants (e.g. [AF>0.01 && AF<0.99]) not to be used with exclude_str           |
| basic_bash_docker | String?        | docker to run basic bash scripts [ubuntu:latest]                                                       |
| docker_registry   | String?        | location of docker images [us.gcr.io/mccarroll-mocha]                                                  |
| bcftools_docker   | String?        | docker to run tasks requiring BCFtools [bcftools:1.11-yyyymmdd]                                        |
| r_mocha_docker    | String?        | docker to run tasks mocha_{plot,summary} [r_mocha:1.11-yyyymmdd]                                       |

Download polygenic scores summary statistics for blood cell counts and liftOver to GRCh38:
```
for i in {163..191}; do wget http://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/PGS000$i/ScoringFiles/PGS000$i.txt.gz; done
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/liftOver
chmod a+x liftOver
for i in {163..191}; do
  (zcat PGS000$i.txt.gz | head -n10;
  zcat PGS000$i.txt.gz | grep -v ^# | \
    awk -F"\t" 'NR==1 {for (i=1; i<=NF; i++) f[$i] = i}
    NR>1 {chr="chr"$(f["chr_name"]); pos=$(f["chr_position"]); ref=$(f["reference_allele"])
    gsub("\t", "@", $0); printf "%s\t%d\t%d\t%s\t0\t+\n", chr, pos-1, pos+length(ref), $0}' | \
    ./liftOver /dev/stdin hg19ToHg38.over.chain.gz /dev/stdout /dev/stderr) | \
    awk '{rev["A"]="T"; rev["C"]="G"; rev["G"]="C"; rev["T"]="A"}
    NR==10 {for (i=1; i<=NF; i++) f[$i] = i; snp=f["rsID"]; chr=f["chr_name"]; pos=f["chr_position"];
    ref=f["reference_allele"]; alt=f["effect_allele"]} NR<=10
    NR>10 {split($4,a,"@"); a[chr]=$1; a[pos]=$2+1;
    if ($6=="-") {gsub("[ATCG]","[TAGC]",a[ref]); gsub("[ATCG]","[TAGC]",a[alt]);
    x=""; for (i=length(a[ref]); i>0; i--) x=x substr(a[ref],i,1); a[ref]=x;
    x=""; for (i=length(a[alt]); i>0; i--) x=x substr(a[alt],i,1); a[alt]=x}
    for (i=1; i<length(a); i++) printf "%s\t",a[i]; printf a[length(a)]"\n"}' | gzip > PGS000$i.hg38.txt.gz
done
gsutil -m cp PGS000{163..191}.hg38.txt.gz gs://{google-bucket}/scores/
```

Define options to run the WDL to compute the polygenic scores related to mLOX:
```json
{
  "score.sample_set_id": "hapmap370k",
  "score.region": "chrX",
  "score.tag": "AS",
  "score.ext_string": "mLOX",
  "score.summary_path": "gs://{google-bucket}/GRCh38/scores",
  "score.summary_files": ["PGS000163.hg38.txt.gz", "PGS000164.hg38.txt.gz", "PGS000165.hg38.txt.gz", "PGS000166.hg38.txt.gz",
    "PGS000167.hg38.txt.gz", "PGS000168.hg38.txt.gz", "PGS000169.hg38.txt.gz", "PGS000170.hg38.txt.gz", "PGS000171.hg38.txt.gz",
    "PGS000172.hg38.txt.gz", "PGS000173.hg38.txt.gz", "PGS000174.hg38.txt.gz", "PGS000175.hg38.txt.gz", "PGS000176.hg38.txt.gz",
    "PGS000177.hg38.txt.gz", "PGS000178.hg38.txt.gz", "PGS000179.hg38.txt.gz", "PGS000180.hg38.txt.gz", "PGS000181.hg38.txt.gz",
    "PGS000182.hg38.txt.gz", "PGS000183.hg38.txt.gz", "PGS000184.hg38.txt.gz", "PGS000185.hg38.txt.gz", "PGS000186.hg38.txt.gz",
    "PGS000187.hg38.txt.gz", "PGS000188.hg38.txt.gz", "PGS000189.hg38.txt.gz", "PGS000190.hg38.txt.gz", "PGS000191.hg38.txt.gz"],
  "score.samples_file": "gs://{google-bucket}/tsvs/hapmap370k.mLOX.lines",
  "score.batch_tsv_file": "gs://{google-bucket}/tsvs/hapmap370k.batch.tsv",
  "score.ref_path": "gs://mccarroll-mocha/GRCh38",
  "score.data_path": "gs://{google-bucket}/imp_vcfs",
  "score.docker_registry": "us.gcr.io/mccarroll-mocha"
}
```

Define options to run the WDL to compute genome-wide polygenic scores for all samples including adjusted polygenic scores if a table `hapmap370k.pcs.tsv` with principal components is available:
```json
{
  "score.sample_set_id": "hapmap370k",
  "score.tag": "DS",
  "score.summary_path": "gs://{google-bucket}/scores",
  "score.summary_files": ["PGS000163.hg38.txt.gz", "PGS000164.hg38.txt.gz", "PGS000165.hg38.txt.gz", "PGS000166.hg38.txt.gz",
    "PGS000167.hg38.txt.gz", "PGS000168.hg38.txt.gz", "PGS000169.hg38.txt.gz", "PGS000170.hg38.txt.gz", "PGS000171.hg38.txt.gz",
    "PGS000172.hg38.txt.gz", "PGS000173.hg38.txt.gz", "PGS000174.hg38.txt.gz", "PGS000175.hg38.txt.gz", "PGS000176.hg38.txt.gz",
    "PGS000177.hg38.txt.gz", "PGS000178.hg38.txt.gz", "PGS000179.hg38.txt.gz", "PGS000180.hg38.txt.gz", "PGS000181.hg38.txt.gz",
    "PGS000182.hg38.txt.gz", "PGS000183.hg38.txt.gz", "PGS000184.hg38.txt.gz", "PGS000185.hg38.txt.gz", "PGS000186.hg38.txt.gz",
    "PGS000187.hg38.txt.gz", "PGS000188.hg38.txt.gz", "PGS000189.hg38.txt.gz", "PGS000190.hg38.txt.gz", "PGS000191.hg38.txt.gz"],
  "score.covars_file": "gs://{google-bucket}/tsvs/hapmap370k.pcs.tsv",
  "score.batch_tsv_file": "gs://{google-bucket}/tsvs/hapmap370k.batch.tsv",
  "score.ref_path": "gs://mccarroll-mocha/GRCh38",
  "score.data_path": "gs://{google-bucket}/imp_vcfs",
  "score.docker_registry": "us.gcr.io/mccarroll-mocha"
}
```

The `batch_tsv_file` input table should look like this:

| batch_id | chr1_imp_vcf             | chr1_imp_vcf_index           | ... | chrX_imp_vcf             | chrX_imp_vcf_index           |
|----------|--------------------------|------------------------------|-----|--------------------------|------------------------------|
| A        | hapmap370k.A.chr1.as.bcf | hapmap370k.A.chr1.as.bcf.csi | ... | hapmap370k.A.chrX.as.bcf | hapmap370k.A.chrX.as.bcf.csi |
| B        | hapmap370k.B.chr1.as.bcf | hapmap370k.B.chr1.as.bcf.csi | ... | hapmap370k.B.chrX.as.bcf | hapmap370k.B.chrX.as.bcf.csi |

The `output_tsv_file` from the imputation pipeline will include these columns

Running with Terra
==================

While the pipeline can be run with a Cromwell server alone, even on your laptop, it has been fully tested to work in the <a href="https://terra.bio/">Terra</a> environment developed at the <a href="https://www.broadinstitute.org/data-sciences-platform">Broad Institute</a>. To setup you workspace to work with the MoChA WDL on Terra you will have to:

* <a href="https://support.terra.bio/hc/en-us/articles/360034677651-Account-setup-and-exploring-Terra">Setup an account</a> with Terra associated with a billing account
* Create a new Terra workspace
* Find the <a href="https://portal.firecloud.org/?return=terra#methods/mocha/mocha/">MoChA method</a> in the Broad Methods Repository and have it exported to your workspace (you can choose **sample_set** as Root Entity Type)
* Setup resources for GRCh37 or GRCh38 (available for download <a href="http://software.broadinstitute.org/software/mocha/">here</a>) by unpacking them in a directory in your own Google bucket, such as **gs://{google-bucket}/GRCh38/** and making sure that the **ref_path** variable points to that path including the trailing **/**
* Make sure the bucket with the reference genome resources can be accessed by your Terra service account. Your Terra service account will be labeled like **pet-012345678909876543210@TERRA-PROJECT.iam.gserviceaccount.com** and you will have to add it to the list of members with permission to access the bucket with the role "Storage Object Viewer". You can use a command like the following: **gsutil iam ch serviceAccount:pet-012345678909876543210@TERRA-PROJECT.iam.gserviceaccount.com:objectViewer gs://{google-bucket}**
* Upload your CSV/BPM/EGT/XML/ZIP manifest files in the same location in your Google bucket, such as **gs://{google-bucket}/manifests/** and making sure that the **manifest_path** variable points to that path including the trailing **/**. Alternatively you can leave the **manifest_path** variable empty and include the full paths in the **batch_tsv_file** table
* Upload your IDAT/GTC/CEL/CHP/TXT/VCF files to your Google bucket, if you have uploaded CHP/TXT files, make sure you upload also the corresponding SNP **AxiomGT1.snp-posteriors.txt** and report **AxiomGT1.report.txt** files
* Format two TSV tables describing samples and batches as explained in the inputs section, upload them to your Google bucket, and make sure variables **sample_tsv_file** and **batch_tsv_file** are set to their location. If you include the file names without their absolute paths, you can include the path including the trailing **/** in the **path** column but you have to make sure all data files from the same batch are in the same directory. See below for examples for Illumina and Affymetrix arrays
* Create a JSON file with all the variables including the cohort ID (**sample_set_id**), the **mode** of the analysis, the location of your tables and resource files, parameters to select the amount of desired parallelization. See below for examples for Illumina and Affymetrix arrays
* From your workspace, go to your WORKSPACE tab, select the MoChA workflow that you had previously imported from the Broad Methods Repository
* Select option \"Run workflow with inputs defined by file paths\" and we further recommend to select options \"Use call caching\" and \"Delete intermediate outputs\" (as these can increase 4-5x the storage footprint)
* Click the on \"upload json\" and select the JSON files with the variable defining your run. Then click on \"SAVE\" and then click on \"RUN ANALYSIS\"
* A job will have spawned that you will be able to monitor through the Job Manager to check the pipeline progress. While monitoring the progress, you will be able to open some summary outputs that are generated before the pipeline fully completes

Running with Cromwell on Google Cloud Platform
==============================================

If you want to run the pipeline without Terra, you can set up your own Cromwell server instead. This is notoriously difficult. These are the steps that I personally advice to take:

* Initialize your Google Cloud configuration:
```
$ gcloud init
```
* If your institution has not provided you with one, create a billing account from the <a href="https://console.cloud.google.com/billing">Billing</a> page
* Create a project `{google-project}` in your Google Cloud Platform account from the <a href="https://console.cloud.google.com/cloud-resource-manager">Manage resources</a> page
* Find default service account for `{google-project}` from the <a href="https://console.cloud.google.com/iam-admin/iam">IAM</a> page which should be labeled as `{google-number}-compute@developer.gserviceaccount.com`
* Download a private json key for your service account through the following command:
```
$ gcloud iam service-accounts keys create {google-project}.key.json --iam-account={google-number}-compute@developer.gserviceaccount.com
```
or from the <a href="https://console.cloud.google.com/iam-admin/serviceaccounts">Service accounts</a> page by selecting the `{google-project}` first and then the `{google-number}-compute@developer.gserviceaccount.com` service account
* Enable the Cloud Life Sciences API for your Google project from the Google <a href="https://console.developers.google.com/apis/library/lifesciences.googleapis.com">console</a>
* Enable the Compute Engine API for your Google project from the Google <a href="https://console.cloud.google.com/apis/library/compute.googleapis.com">console</a>
* Enable the Google Cloud Storage JSON API for your Google project from the Google <a href="https://console.cloud.google.com/apis/library/storage-api.googleapis.com">console</a>
```
$ for api in lifesciences compute storage-api; do
  gcloud services enable $api.googleapis.com
done
```
* You need the following roles available to your service account: `Cloud Life Sciences Workflows Runner`, `Service Account User`, `Storage Object Admin`. To add these roles you need to have been assigned the rights to change roles you can use the following commands (or manually from the <a href="https://console.cloud.google.com/iam-admin/iam">IAM</a> page):
```
$ for role in lifesciences.workflowsRunner iam.serviceAccountUser storage.objectAdmin; do
  gcloud projects add-iam-policy-binding {google-project} --member serviceAccount:{google-number}-compute@developer.gserviceaccount.com --role roles/$role
done
```
* Create a Google bucket `gs://{google-bucket}` making sure request pays is not activated from the <a href="https://console.cloud.google.com/storage/browser">Storage</a> page or through the following command:
```
gsutil mb -p {google-project} -l us-central1-a -b on gs://{google-bucket}
```
* Create a Google virtual machine (VM) with name `INSTANCE-ID` from the <a href="https://console.cloud.google.com/compute/instances">VM instances</a> page. The **n1-standard-1** (1 vCPU, 3.75 GB memory) VM with Ubuntu 21.04 will be sufficient, but make sure to provide at least 100GB of space for the boot disk, rather than the default 10GB, as the mySQL database will easily fail with the default space settings

* Start the virtual machine with the following command:
```
$ gcloud compute instances start INSTANCE-ID --project {google-project} --zone us-central1-a
```
* Copy the private json key for your service account to the VM:
```
$ gcloud compute scp {google-project}.key.json INSTANCE-ID: --project mccarroll-mocha --zone us-central1-a
```
* Login to the VM with the following command (notice this command will enable ssh tunneling):
```
$ gcloud compute ssh INSTANCE-ID --project {google-project} --zone us-central1-a -- -L 8000:localhost:8000
```
* As the suggested VM costs approximately $35/month to run, when you are not running any computations, to avoid incurring unnecessary costs, you can stop the virtual machine with the following command:
```
$ gcloud compute instances stop INSTANCE-ID --project {google-project} --zone us-central1-a
```
* Once in the virtual machine, install some basic packages as well as the WOMtool and the Cromwell server (replace `XY` with the current version):
```
$ sudo apt update && sudo apt install default-jre-headless mysql-server jq
$ wget https://github.com/broadinstitute/cromwell/releases/download/XY/womtool-XY.jar
$ wget https://github.com/broadinstitute/cromwell/releases/download/XY/cromwell-XY.jar
```
* Download the <a href="https://github.com/broadinstitute/cromwell/blob/develop/cromwell.example.backends/PAPIv2.conf">PAPIv2.conf</a> configuration file with the following command:
```
$ wget https://raw.githubusercontent.com/broadinstitute/cromwell/develop/cromwell.example.backends/PAPIv2.conf
```
* This basic configuration file contains only the **backend** configuration stanza. To configure your Cromwell server you have to edit additional configuration stanzas before the **backend** one
* Add the **webservice** configuration stanza to the `PAPIv2.conf` configuration file to make sure the <a href="https://cromwell.readthedocs.io/en/stable/Configuring/#server">server</a> will only be accessible from the local machine (by default it is open to any interface):
```
webservice {
  port = 8000
  interface = 127.0.0.1
}
```
* Add the following **akka** configuration stanza to the `PAPIv2.conf` configuration file to make sure Cromwell has enough time to respond to requests for metadata even for large workflows, as it would be the case for MoChA runs on biobank-size cohorts:
```
akka {
  http {
    server {
      request-timeout = 300s
    }
  }
}
```
* Add the **google** configuration stanza to the `PAPIv2.conf` configuration file (leave `service-account` and `service_account` verbatim):
```
google {
  application-name = "cromwell"
  auths = [
    {
      name = "service-account"
      scheme = "service_account"
      json-file = "{google-project}.key.json"
    }
  ]
}
```
* Change `auth = "application-default"` to `auth = "service-account"` in the `PAPIv2.conf` configuration file in both instances where it occurs within the **backend** configuration stanza
* Change `project = "my-cromwell-workflows"` and `project = "google-billing-project"` to `project = "{google-project}"` in the `PAPIv2.conf` configuration file within the **backend** configuration stanza
* Change `root = "gs://my-cromwell-workflows-bucket"` to `root = "gs://{google-bucket}/cromwell/executions"` in the `PAPIv2.conf` configuration file within the **backend** configuration stanza
* As Cromwell will need to load some input files to properly organize the batching, it will need the <a href="https://cromwell.readthedocs.io/en/stable/filesystems/Filesystems/#engine-filesystems">engine filesystem</a> activated for reading files, Add the **engine** configuration stanza to the `PAPIv2.conf` configuration file (leave `service-account` verbatim):
```
engine {
  filesystems {
    gcs {
      auth = "service-account"
      project = "{google-project}"
    }
  }
}
```
* Add the **database** configuration stanza to the `PAPIv2.conf` configuration file (as explained <a href="https://cromwell.readthedocs.io/en/stable/Configuring/#database">here</a> and <a href="https://cromwell.readthedocs.io/en/develop/tutorials/PersistentServer/">here</a>):
```
database {
  profile = "slick.jdbc.MySQLProfile$"
  db {
    driver = "com.mysql.cj.jdbc.Driver"
    url = "jdbc:mysql://localhost/cromwell?rewriteBatchedStatements=true"
    user = "user"
    password = "pass"
    connectionTimeout = 60000
  }
}
```
* Activate <a href="https://cromwell.readthedocs.io/en/stable/cromwell_features/CallCaching">CallCaching</a> to allow to re-run jobs without re-running tasks that have already completed by adding the **call-caching** configuration stanza:
```
call-caching {
  enabled = true
  invalidate-bad-cache-results = true
}
```
* If you want to configure the directory of the mysql database to something other than `/var/lib/mysql` make sure to set the `datadir` variable in the `/etc/mysql/mysql.conf.d/mysqld.cnf` file
* Start the mySQL server and initialize the root user with the following command (use cromwell as the default root password):
```
$ sudo systemctl start mysql
$ sudo mysql_secure_installation
```
* Login into the mySQL database and run the following commands to create a database to be used by Cromwell:
```
$ sudo mysql --user=root --password=cromwell --execute "SET GLOBAL validate_password.policy=LOW"
$ sudo mysql --user=root --password=cromwell --execute "SET GLOBAL validate_password.length=4"
$ sudo mysql --user=root --password=cromwell --execute "CREATE USER 'user'@'localhost' IDENTIFIED BY 'pass'"
$ sudo mysql --user=root --password=cromwell --execute "GRANT ALL PRIVILEGES ON * . * TO 'user'@'localhost'"
$ sudo mysql --user=root --password=cromwell --execute "CREATE DATABASE cromwell"
```
* Start the Cromwell server on the Google VM with the following command:
```
$ (java -XX:MaxRAMPercentage=90 -Dconfig.file=PAPIv2.conf -jar cromwell-XY.jar server &)
```
* If you get the error `Bind failed for TCP channel on endpoint [/127.0.0.1:8000]`, it means some other service is already using port `8000`. Run `lsof -i:8000` to see which service. If it is another Cromwell server, you can stop that with:
```
$ killall java
```
* Now it would be a good time to test whether the configuration worked well. Edit and generate the following `hello.wdl` workflow file:
```
version development

workflow myWorkflow {
  call myTask
}

task myTask {
  command {
    echo "hello world"
  }
  output {
    String out = read_string(stdout())
  }
  runtime {
    docker: "ubuntu:latest"
  }
}
```
* Submit the workflow to the Cromwell server:
```
$ java -jar cromwell-XY.jar submit hello.wdl
```
* You will receive an output as follows with the `{id}` of the submitted job:
```
[yyyy-mm-dd hh:mm:ss,ss] [info] Slf4jLogger started
[yyyy-mm-dd hh:mm:ss,ss] [info] Workflow 01234567-89ab-dcef-0123-456789abcdef submitted to http://localhost:8000
```
* Create an `options.json` file for Cromwell that should look like this (additional options can be used from <a href="https://cromwell.readthedocs.io/en/stable/wf_options/Google/">here</a>):
```
{
  "delete_intermediate_output_files": true,
  "final_workflow_outputs_dir": "gs://{google-bucket}/cromwell/outputs",
  "use_relative_output_paths": true,
  "final_workflow_log_dir": "gs://{google-bucket}/cromwell/wf_logs",
  "final_call_logs_dir": "gs://{google-bucket}/cromwell/call_logs"
}
```
* Download the MoChA WDL pipeline:
```
$ curl https://raw.githubusercontent.com/freeseek/mocha/master/wdl/mocha.wdl -o mocha.wdl
```
* To verify that your input `{sample-set-id}.json` file is correcly formatted, you can use the WOMtool as follows:
```
$ java -jar womtool-XY.jar validate mocha.wdl -i {sample-set-id}.json
```
* To submit a job, use the command (you can run this straight from your computer if you logged to the VM with ssh tunneling):
```
$ java -jar cromwell-XY.jar submit mocha.wdl -i {sample-set-id}.json -o options.json
```
* To monitor the status of jobs submitted to the server, you can use:
```
$ curl -X GET http://localhost:8000/api/workflows/v1/query | jq
```
* To monitor the status of a specific job and extract the metadata:
```
$ curl -X GET http://localhost:8000/api/workflows/v1/{id}/metadata | jq
```
* To monitor a submitted job with workflow ID `{id}` just open your browser and go to URL:
```
http://localhost:8000/api/workflows/v1/{id}/timing
```
* To abort a running job, you can run:
```
$ curl -X POST http://localhost:8000/api/workflows/v1/{id}/abort | jq
```
* It is possible for the Cromwell server to crash. In this case it is okay to restart the Cromwell server. This will automatically connect with the mySQL server and resume the jobs that were running
* It is also possible for the mySQL database to grow to such an extent that it will then completely fill the VM disk. To check the size of the mySQL database, run:
```
$ sudo ls -l /var/lib/mysql/cromwell
```
* If you want to flush the metadata database to recover disk space from the VM, login into the database and run (though notice that you will have first to stop the Cromwell server and then restart it):
```
$ killall java # to stop the Cromwell server
$ sudo mysql --user=root --password=cromwell --execute "DROP DATABASE cromwell"
$ sudo mysql --user=root --password=cromwell --execute "CREATE DATABASE cromwell"
$ (java -XX:MaxRAMPercentage=90 -Dconfig.file=PAPIv2.conf -jar cromwell-XY.jar server &)
```
* If you want to clean temporary files and log files, after you have properly moved all the output files you need, you can delete the workflow executions directory defined as `root` in the `PAPIv2.conf` configuration file. Remember that failure to do proper cleanup could cause you to incur unexpected storage costs. After making sure there are no active jobs running with the Cromwell server, you can remove the executions and logs directory with the following command:
```
$ gsutil -m rm -r gs://{google-bucket}/cromwell/executions gs://{google-bucket}/cromwell/call_logs
```
Notice that either flushing the metadata from the database or removing the workflow executions files will invalidate the cache from previously run tasks

When running computations over large cohorts, it is <a href="https://cloud.google.com/compute/docs/instances/create-start-preemptible-instance#best_practices">best practice</a> to run the jobs on nights and weekends, as this will increase the chances that preemptible jobs will not be stopped and replaced by a non-preemptible version. This will make jobs run both faster and cheaper

Troubleshooting
===============

These are some of the messages that you might receive when something goes wrong:

* If you run the pipeline on Terra with the `Delete intermediate options` flag selected and your workflow keeps showing as Running even after the final outputs have been generated, it is possible that the Cromwell server behind Terra might have failed while deleting the intermediate outputs. This is an <a href="https://support.terra.bio/hc/en-us/community/posts/360071861791-Job-seems-stuck-indefinitely-at-the-delete-intermediate-files-step-and-does-not-complete">issue</a> that is being patched
* `Failed to evaluate input 'disk_size' (reason 1 of 1): ValueEvaluator[IdentifierLookup]: No suitable input for ...`: this indicates that Cromwell was unable to find the size of one of the input files for a task, most likely because the file does not exist where indicated by the user
* `The job was stopped before the command finished. PAPI error code 2. Execution failed: generic::unknown: pulling image: docker pull: running ["docker" "pull" "###"]: exit status 1 (standard error: "Error response from daemon: Get https://###: unknown: Project 'project:###' not found or deleted.\n")`: this means that one of the docker images provided does not exist
* `Job exit code 255. Check gs://###/stderr for more information. PAPI error code 9. Please check the log file for more details:`: if this is an error provided by the task cel2chp, it means that the `apt-probeset-genotype` command has encountered an error. Reading the `stderr` file should easily provide an explanation
* `The job was stopped before the command finished. PAPI error code 10. The assigned worker has failed to complete the operation`: this could mean that the job was preempted despite the fact that it was not running in preemptible computing (see <a href="https://support.terra.bio/hc/en-us/community/posts/360046714292-Are-you-experiencing-PAPI-error-code-10-Read-this-for-a-workaround-">here</a>)
* `The compute backend terminated the job. If this termination is unexpected, examine likely causes such as preemption, running out of disk or memory on the compute instance, or exceeding the backend's maximum job duration`: this could be an indication that a task was killed as it requested too much memory
* `The job was aborted from outside Cromwell`: this could be an indication that a task was killed as it requested too much memory
* idat2gtc task fails with internal `stderr` message `Normalization failed! Unable to normalize!`: this means that either the IDAT sample is of very poor quality and it cannot be processed by the GenCall algorithm or that you have matched the IDAT with the wrong Illumina BPM manifest file
* idat2gtc task fails with internall `stderr` message `System.Exception: Unrecoverable Error...Exiting! Unable to find manifest entry ######## in the cluster file!`: this means that you are using the incorrect Illumina EGT cluster file
* If when monitor the status of the job you get the error: `Job Manager is running but encountered a problem getting data from its workflow server. Click here to start over. 500: Internal Server Error` then it means that there is too much metadata input and output into the tasks for the Job Manager to handle te request. This metadata limit is a known issue currently being worked on
* When you mismatch a BPM manifest file with an IDAT in task idat2gtc, iaap-cli, while outputting an error message such as `Normalization failed for sample: ########! This is likely a BPM and IDAT mismatch. ERROR: Index was outside the bounds of the array.` will not return an error code, causing the pipeline to fail at the next task
* If you are running with your own Cromwell server using the PAPIv2 API and you get error `Error attempting to Execute cromwell.backend.google.pipelines.common.api.PipelinesApiRequestManager$UserPAPIApiException: Unable to complete PAPI request due to a problem with the request (The caller does not have permission).` whenever Cromwell tries to submit a task, the the cause is the service account that you are using to run the computations with Google Cloud does not have the <a href="https://cloud.google.com/life-sciences/docs/concepts/access-control#roles">Cloud Life Sciences</a> Workflows Runner (`lifesciences.workflowsRunner`) role set
* If you are running with your own Cromwell server using the PAPIv2 API and some of your tasks start running but they then fail with and you have the error `Error attempting to Execute cromwell.backend.google.pipelines.common.api.PipelinesApiRequestManager$UserPAPIApiException: Unable to complete PAPI request due to a problem with the request (Error: checking service account permission: caller does not have access to act as the specified service account: "MY_NUMBER-compute@developer.gserviceaccount.com").` in the workflow log, then the cause is the service account not having the <a href="https://cloud.google.com/iam/docs/service-accounts#user-role">Service Account User</a> (`roles/iam.serviceAccountUser`) role set
* If you are running with your own Cromwell server using the PAPIv2 API and all of your tasks start running but they all fail with a log file including just the line `yyyy/mm/dd hh:mm:ss Starting container setup.` the cause is the service account not having the <a href="https://cloud.google.com/storage/docs/access-control/iam-roles">Storage Object</a> Admin (`storage.objectAdmin`) role set
* When running the phasing step on a very large cohort, we have noticed that tasks including the MHC can run very slowly. This is due to the very special nature of the MHC and there is currently no solution. Currently these tasks might dominate the speed of the whole pipeline

If you receive an error that you do not understand, feel free contact the <a href="mailto:giulio.genovese@gmail.com">author</a> for troubleshooting

Implementation Notes
====================

The following softwares are used by the various steps of the pipelines:

| software                                                          | description                                    | license     |
|-------------------------------------------------------------------|------------------------------------------------|-------------|
| <a href="https://support.illumina.com/array/array_software/illumina-array-analysis-platform/downloads.html">Illumina Array Analysis Platform</a> | runs the GenCall algorithm to call genotypes   | proprietary |
| <a href="https://support.illumina.com/array/array_software/beeline/downloads.html">Illumina AutoConvert</a> | runs the GenCall algorithm to call genotypes   | proprietary |
| <a href="https://www.thermofisher.com/us/en/home/life-science/microarray-analysis/microarray-analysis-partners-programs/affymetrix-developers-network/affymetrix-power-tools.html">Affymetrix Power Tools</a> | runs the Axiom algorithm to call genotypes     | GPLv2       |
| <a href="https://github.com/lh3/bwa">BWA</a>                      | maps DNA sequences against a reference genome  | GPLv3 / MIT |
| <a href="https://github.com/samtools/htslib">HTSlib</a>           | a C library for accessing SAMs, CRAMs and VCFs | MIT / BSD   |
| <a href="https://github.com/samtools/samtools">Samtools</a>       | manipulates SAMs and BAMs                      | MIT         |
| <a href="https://github.com/samtools/bcftools">BCFtools</a>       | manipulates VCFs and BCFs                      | MIT         |
| <a href="https://github.com/freeseek/gtc2vcf">gtc2vcf</a>         | converts DNA microarray data to VCF            | MIT         |
| <a href="https://github.com/freeseek/mocha">MoChA</a>             | calls mosaic chromosomal alterations           | MIT         |
| <a href="https://github.com/odelaneau/shapeit4">SHAPEIT4</a>      | estimates haplotype phase from genotypes       | MIT         |
| <a href="https://github.com/poruloh/Eagle">Eagle</a>              | estimates haplotype phase from genotypes       | GPLv3       |
| <a href="https://jmarchini.org/impute5/">IMPUTE5</a>              | imputes genotypes from phased haplotypes       | proprietary |
| <a href="http://faculty.washington.edu/browning/beagle/beagle.html">Beagle5</a> | imputes genotypes from phased haplotypes       |  GPLv3       |
| <a href="https://github.com/broadinstitute/cromwell">Cromwell</a> | workflow management system                     | BSD         |

For users and developers that want to understand the logic and ideas behind the pipeline, here is a list of important aspects:

* The WDL is written according to the version development specification as it requires grouping by batch and sub batch through the <a href="https://github.com/openwdl/wdl/blob/main/versions/development/SPEC.md#mapxarrayy-collect_by_keyarraypairxy">collect_bt_key()</a> function together with the <a href="https://github.com/openwdl/wdl/blob/main/versions/development/SPEC.md#mapxy-as_maparraypairxy">as_map()</a> and <a href="https://github.com/openwdl/wdl/blob/main/versions/development/SPEC.md#arrayx-keysmapxy">keys()</a> functions
* To achieve parallelization both across batches and across genome windows, the <a href="https://github.com/openwdl/wdl/blob/main/versions/development/SPEC.md#arrayarrayx-transposearrayarrayx">transpose()</a> function is used to move from one model to the other
* To minimize memory load on Cromwell, computed gender and call rate maps are passed to tasks **vcf_qc** and **vcf_mocha** as files rather than WDL maps, and similarly sample IDs are passed to tasks **vcf_merge** as files rather than WDL arrays
* As the current version of Cromwell does not accept optional outputs in tasks dispatched to Google Cloud, we have to **touch** optional output files in tasks **cel2affy** and **vcf_split**
* As the current version of Cromwell does not delocalize on Google Cloud lists of output files whose names are determined during runtime, we use the trick of delocalizing a **Directory** for tasks **idat2gtc**, **cel2affy**, **vcf_scatter**, and **vcf_split**
* Terra does not allow <a href="https://support.terra.bio/hc/en-us/articles/360037128572-Scatter-Gather-Parallelism">scatters</a> with width larger than 35,000 (although Cromwell by default allows width up to 1,000,000), so to allow cohorts with more than 35,000 samples to be able to run on Terra we avoid scattering over each sample in the cohort
* As estimating the sizes of a large array of files is extremely time consuming in Google Cloud, causing tasks to spend a long time in PreparingJob state before starting to localize the files, for array of IDAT, GTC, CEL, and CHP files we estimate the size of the first file in an array and assume all the other files have similar sizes
* As the Terra job manager crashes if too much metadata is transferred between the workflow and the tasks, we try to transfer metadata to tasks through files wherever possible
* To avoid SIGPIPE error 141 when piping a stream to the UNIX command `head`, we use the `|| if [[ $? -eq 141 ]]; then true; else exit $?; fi` trick
* To avoid return error 1 when grep returns no matches, we use the `|| if [[ $? -eq 1 ]]; then true; else exit $?; fi` trick
* To avoid return error 2 when gunzip is used on a file localized as a hard link, we use `gunzip --force`
* Due to a numpy <a href="https://github.com/numpy/numpy/issues/10448">BUG in interp</a>, the option `period = np.inf` needs to be added when interpolating chromosome coordinates
* Conversion from IDAT to GTC for some arrays might fail to infer gender when using IAAP CLI. For this reason, we use by default the IAAP CLI flag **\--gender-estimate-call-rate-threshold -0.1** (similarly to how it is done by the Broad Genomics platform in <a href="https://github.com/broadinstitute/warp/blob/IlluminaGenotypingArray_v1.11.0/tasks/broad/IlluminaGenotypingArrayTasks.wdl">IlluminaGenotypingArrayTasks.wdl</a>)
* Conversion from IDAT to GTC using IAAP CLI is 2-3x faster when run on multiple IDATs at once than when run individually on each IDAT (the latter being how the Broad Genomics platform runs it)
* Conversion from IDAT to GTC for Omni5 arrays requires more than 8GB of RAM
* Conversion of many GTC or CHP files to VCF can take a long time and be difficult to run on preemptible cloud computing, so making it hierarchical, while requiring more CPU cycles, can actually be cheaper. By default, up to 1,024 GTC or CHP files are converted at once
* When converting from GTC to VCF and from CHP to VCF, heavy random access to the reference genome is needed, so it is important that enough memory to cache the reference is provided or else the task can run excruciatingly slowly
* Genotyping using the Affymetrix Power Tools is best run on many samples at once. Affymetrix recommends, in their <a href="https://assets.thermofisher.com/TFS-Assets/LSG/manuals/axiom_genotyping_solution_analysis_guide.pdf">data analysis guide</a>, to batch as large a batch size as computationally feasible, or up to 4800 samples. However, this can take up to 16 hours to execute, it is not parallelizable, and the genotyping software from Affymetrix is not multi-threaded. For this reason task **cel2chp** is set to run on non-preemptible computing by default
* Conversion to VCF for Affymetrix data can be done either from AGCC CHP files or from matrices of tab delimited genotype calls and confidences. We recommend the former as it can be more easily parallelizable by splitting the conversion into sub batches. When converting the latter, the whole batch will need to be converted into VCF in one single task. For this reason task **txt2vcf** is set to run on non-preemptible computing by default
* All tasks that output a VCF can, if requested, output either in compressed or uncompressed format. For VCFs with BAF and LRR intensities, we observed a modest ~25% space saving, most likely due to the high entropy in the intensity measurements. For VCF files with exclusively the GT format field, we observed around ~75% space saving
* Due to high sequence divergence between HLA class I and class II genes in the MHC (<a href="http://doi.org/10.1101/gr.213538.116">Norman et al. 2017</a>), heterozygous variants in the MHC region show unusual degrees of allelic imbalance. To avoid false positive mosaic chromosomal alterations we mask the ~6Mbp segment between rs9468147 and rs9366822 on chromosome 6. We further mask a ~1Mbp KIR region between rs4806703 and rs34679212 on chromosome 19
* Due to residual correlation between the BAF and LRR, observed both for Illumina arrays (<a href="http://doi.org/10.1101/gr.5686107">Oosting et al\. 2007</a> and <a href="http://doi.org/10.1186/1471-2105-9-409">Staaf et al\. 2008</a>) and Affymetrix arrays (<a href="http://doi.org/10.1038/srep36158">Mayrhofer et al\. 2016</a>), MoChA performs a BAF correction by regressing BAF values against LRR values. This improves detection of mosaic chromosomal alterations at low cell fractions. However, for this correction to be effective, batches need to include 100s of samples
* While Illumina internally computes BAF and LRR and stores these values in the GTC files, we instead recompute these values from normalized intensities and genotype cluster centers following the method first described by Illumina in <a href="http://doi.org/10.1101/gr.5402306">Peiffer et al\. 2006</a> but we do non truncate the BAF values between 0 and 1 like Illumina does. This allows to better estimate the residual correlation between the BAF and the LRR for homozygous calls
* When phasing a dataset with a reference panel only the variants present both in the reference panel and the dataset will be phased. As the 1000 Genomes project reference panel for <a href="https://www.internationalgenome.org/announcements/updated-GRCh38-liftover/">GRCh38</a> is missing the chromosome X PAR1/PAR2 regions and these are exactly the regions that need to be phased to detect mosaic loss of chromosome Y, the most common chromosomal alteration in clonal hematopoiesis, we recommend using SHAPEIT4 without a reference panel to detect mosaic loss of chromosome Y at low cell fractions. This is the default behavior for cohorts whose sample size is more than twice the size of the 1000 Genomes reference panel (2,504 samples)
* SHAPEIT4 does not have an option to disable imputation of missing target genotypes, like Eagle does. This is problematic as homozygous variants improperly imputed as heterozygous can be interpreted as huge BAF imbalances. To overcome this limitation of SHAPEIT4, we use BCFtools annotate to import the genotypes back in the target VCF using the **\--columns -FMT/GT** option to replace only existing unphased genotype values without modifying missing genotypes
* SHAPEIT4 is used to phase genotypes across genome windows with overlapping ends. To avoid phase switch errors across windows which would negatively affect the ability to detect large mosaic chromosomal alterations at low cell fractions, we use BCFtools concat with the option **\--ligate** to ligate phased VCFs by matching phase at phased haplotypes over overlapping windows
* Genotype data is transmitted across tasks in binary VCF format (BCF) rather than plain text VCF format as this is much more efficient to process especially for many samples. The relationship between BCF and VCF is similar to the one between BAM and SAM. The whole pipeline, both for Illumina and Affymetrix data, can completely avoid any intermediate conversion to non-binary format, providing significant time and cost savings (though notice that in **txt** mode the data provided by the user is not in binary format)
* Many large biobanks perform their genotyping across a span of several years. This causes sometimes to have to switch from one DNA microarray version to a newer version of the same DNA microarray (e.g\. the Global Screening Array v1.0 vs\. the Global Screening Array v2.0). While it is not possible to process samples on different arrays in the same batch, the batched design of the pipeline allows to include samples on separate arrays in the same workflow as long as they are assigned to different batches. One small drawback is that, when the whole cohort is phased across batches, only variants that are shared across all batches will be phased. Therefore we recommend to only process samples genotyped on the same array or on fairly similar arrays. As an example, we absolutely recommend to not process samples from the Global Screening Array and the Multi-Ethnic Global Array in the same workflow, although the user will not be prevented from doing so
* The pipeline is designed to be minimalistic and require the user to provide the minimal amount of information necessary to run

Dockerfiles
===========

Dockerfile for BCFtools:
```
FROM ubuntu:21.04
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get -qqy update --fix-missing && \
    apt-get -qqy install --no-install-recommends \
                 wget \
                 bwa \
                 samtools \
                 bcftools && \
    wget http://software.broadinstitute.org/software/gtc2vcf/gtc2vcf_1.11-20210315_amd64.deb && \
    dpkg -i gtc2vcf_1.11-20210315_amd64.deb && \
    wget http://software.broadinstitute.org/software/mocha/bio-mocha_1.11-20210315_amd64.deb && \
    dpkg -i bio-mocha_1.11-20210315_amd64.deb && \
    apt-get -qqy purge --auto-remove --option APT::AutoRemove::RecommendsImportant=false \
                 wget && \
    apt-get -qqy clean && \
    rm -rf gtc2vcf_1.11-20210315_amd64.deb \
           bio-mocha_1.11-20210315_amd64.deb \
           /var/lib/apt/lists/*
```

Dockerfile for Illumina AutoConvert Software:
```
FROM ubuntu:21.04
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get -qqy update --fix-missing && \
    apt-get -qqy install --no-install-recommends \
                 wget \
                 bcftools \
                 unzip \
                 msitools \
                 cabextract \
                 gcc \
                 libc6-dev \
                 libmono-system-windows-forms4.0-cil && \
    wget http://software.broadinstitute.org/software/gtc2vcf/gtc2vcf_1.11-20210315_amd64.deb && \
    dpkg -i gtc2vcf_1.11-20210315_amd64.deb && \
#   wget https://support.illumina.com/content/dam/illumina-support/documents/downloads/software/beeline/autoconvert-software-v2-0-1-installer.zip && \
    wget --no-check-certificate https://www.dropbox.com/s/tm02cu6t0ib1us7/autoconvert-software-v2-0-1-installer.zip && \
    unzip autoconvert-software-v2-0-1-installer.zip && \
    msiextract AutoConvertInstaller.msi && \
    wget ftp://webdata2:webdata2@ussd-ftp.illumina.com/downloads/software/genomestudio/genomestudio-software-v2-0-4-5-installer.zip && \
    unzip -j genomestudio-software-v2-0-4-5-installer.zip && \
    cabextract GenomeStudioInstaller.exe && \
    msiextract a0 && \
    mv Illumina/GenomeStudio\ 2.0/Heatmap.dll Illumina/AutoConvert\ 2.0 && \
    wget --no-check-certificate https://raw.githubusercontent.com/freeseek/gtc2vcf/master/nearest_neighbor.c && \
    gcc -fPIC -shared -O2 -o Illumina/AutoConvert\ 2.0/libMathRoutines.dll.so nearest_neighbor.c && \
    mv Illumina/AutoConvert\ 2.0 /opt && \
    apt-get -qqy purge --auto-remove --option APT::AutoRemove::RecommendsImportant=false \
                 wget \
                 unzip \
                 msitools \
                 cabextract \
                 gcc \
                 libc6-dev && \
    apt-get -qqy clean && \
    rm -rf gtc2vcf_1.11-20210315_amd64.deb \
           autoconvert-software-v2-0-1-installer.zip \
           AutoConvertInstaller.msi \
           genomestudio-software-v2-0-4-5-installer.zip \
           GenomeStudioInstaller.exe \
           0 u0 u1 u2 u3 u4 u5 a0 \
           nearest_neighbor.c \
           Illumina \
           /var/lib/apt/lists/*
```

Dockerfile for Illumina Array Analysis Platform Genotyping Command Line Interface:
```
FROM ubuntu:21.04
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get -qqy update --fix-missing && \
    apt-get -qqy install --no-install-recommends \
                 wget \
                 bcftools \
                 icu-devtools && \
    wget http://software.broadinstitute.org/software/gtc2vcf/gtc2vcf_1.11-20210315_amd64.deb && \
    dpkg -i gtc2vcf_1.11-20210315_amd64.deb && \
    wget ftp://webdata2:webdata2@ussd-ftp.illumina.com/downloads/software/iaap/iaap-cli-linux-x64-1.1.0.tar.gz && \
    mkdir /opt/iaap-cli && \
    tar xzvf iaap-cli-linux-x64-1.1.0.tar.gz -C /opt iaap-cli-linux-x64-1.1.0/iaap-cli --strip-components=1 && \
    ln -s /opt/iaap-cli/iaap-cli /usr/local/bin/iaap-cli && \
    apt-get -qqy purge --auto-remove --option APT::AutoRemove::RecommendsImportant=false \
                 wget && \
    apt-get -qqy clean && \
    rm -rf gtc2vcf_1.11-20210315_amd64.deb \
           iaap-cli-linux-x64-1.1.0.tar.gz \
           /var/lib/apt/lists/*
```

Dockerfile for Affymetrix Power Tools:
```
FROM ubuntu:21.04
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get -qqy update --fix-missing && \
    apt-get -qqy install --no-install-recommends \
                 wget \
                 bcftools \
                 unzip && \
    wget http://software.broadinstitute.org/software/gtc2vcf/gtc2vcf_1.11-20210315_amd64.deb && \
    dpkg -i gtc2vcf_1.11-20210315_amd64.deb && \
    wget --no-check-certificate https://downloads.thermofisher.com/APT/2.11.3/apt_2.11.3_linux_64_bit_x86_binaries.zip && \
    unzip -ojd /usr/local/bin apt_2.11.3_linux_64_bit_x86_binaries.zip apt_2.11.3_linux_64_bit_x86_binaries/bin/apt-probeset-genotype && \
    chmod a+x /usr/local/bin/apt-probeset-genotype && \
    apt-get -qqy purge --auto-remove --option APT::AutoRemove::RecommendsImportant=false \
                 wget && \
    apt-get -qqy clean && \
    rm -rf gtc2vcf_1.11-20210315_amd64.deb \
           apt_2.11.3_linux_64_bit_x86_binaries.zip \
           /var/lib/apt/lists/*
```

Dockerfile for ggplot2 and scripts to plot MoChA calls:
```
FROM ubuntu:21.04
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get -qqy update --fix-missing && \
    apt-get -qqy install --no-install-recommends \
                 wget \
                 bcftools \
                 r-cran-optparse \
                 r-cran-ggplot2 \
                 r-cran-data.table \
                 r-cran-reshape2 && \
    wget http://software.broadinstitute.org/software/mocha/bio-mocha_1.11-20210315_amd64.deb && \
    dpkg -i bio-mocha_1.11-20210315_amd64.deb && \
    apt-get -qqy purge --auto-remove --option APT::AutoRemove::RecommendsImportant=false \
                 wget && \
    apt-get -qqy clean && \
    rm -rf bio-mocha_1.11-20210315_amd64.deb \
           /var/lib/apt/lists/*
```

Dockerfile for SHAPEIT4:
```
FROM ubuntu:21.04
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get -qqy update --fix-missing && \
    apt-get -qqy install --no-install-recommends \
                 wget \
                 g++ \
                 make \
                 libboost-iostreams-dev \
                 libboost-program-options-dev \
                 libhts-dev \
                 libbz2-dev \
		 libssl-dev \
                 libboost-iostreams1.74.0 \
                 libboost-program-options1.74.0 \
                 bcftools && \
    wget --no-check-certificate https://github.com/odelaneau/shapeit4/archive/v4.2.1.tar.gz && \
    tar xzf v4.2.1.tar.gz && \
    cd shapeit4-4.2.1 && \
    sed -i 's/^HTSLIB_INC=\$(HOME)\/Tools\/htslib-1.11$/HTSLIB_INC=-Ihtslib/' makefile && \
    sed -i 's/^HTSLIB_LIB=\$(HOME)\/Tools\/htslib-1.11\/libhts.a$/HTSLIB_LIB=-lhts/' makefile && \
    sed -i 's/^BOOST_LIB_IO=\/usr\/lib\/x86_64-linux-gnu\/libboost_iostreams.a$/BOOST_LIB_IO=-lboost_iostreams/' makefile && \
    sed -i 's/^BOOST_LIB_PO=\/usr\/lib\/x86_64-linux-gnu\/libboost_program_options.a$/BOOST_LIB_PO=-lboost_program_options/' makefile && \
    make && \
    mv bin/shapeit4.2 /usr/bin/ && \
    cd .. && \
    apt-get -qqy purge --auto-remove --option APT::AutoRemove::RecommendsImportant=false \
                 wget \
                 g++ \
                 make \
                 libboost-iostreams-dev \
                 libboost-program-options-dev \
                 libhts-dev \
                 libbz2-dev \
		 libssl-dev && \
    apt-get -qqy clean && \
    rm -rf v4.2.1.tar.gz \
           shapeit4-4.2.1 \
           /var/lib/apt/lists/*
```

Dockerfile for Eagle:
```
FROM ubuntu:21.04
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get -qqy update --fix-missing && \
    apt-get -qqy install --no-install-recommends \
                 git \
                 g++ \
                 make \
                 libboost-iostreams-dev \
                 libboost-program-options-dev \
                 libopenblas-dev \
                 libhts-dev \
                 libboost-iostreams1.74.0 \
                 libboost-program-options1.74.0 \
                 libgomp1 \
                 libopenblas0 \
                 bcftools && \
    git config --global http.sslVerify false && \
    git clone https://github.com/poruloh/Eagle.git && \
    cd Eagle/src && \
    sed -i 's/^BLAS_DIR =.*/BLAS_DIR =/' Makefile && \
    sed -i 's/^BOOST_INSTALL_DIR =.*/BOOST_INSTALL_DIR =/' Makefile && \
    sed -i 's/^HTSLIB_DIR =.*/HTSLIB_DIR =/' Makefile && \
    sed -i 's/^ZLIB_STATIC_DIR =.*/ZLIB_STATIC_DIR =/' Makefile && \
    sed -i 's/^LIBSTDCXX_STATIC_DIR =.*/LIBSTDCXX_STATIC_DIR =/' Makefile && \
    sed -i 's/^GLIBC_STATIC_DIR =.*/GLIBC_STATIC_DIR =/' Makefile && \
    sed -i 's/ memcpy.o$//' Makefile && \
    make && \
    mv eagle /usr/bin && \
    cd ../.. && \
    apt-get -qqy purge --auto-remove --option APT::AutoRemove::RecommendsImportant=false \
                 git \
                 g++ \
                 make \
                 libboost-iostreams-dev \
                 libboost-program-options-dev \
                 libopenblas-dev \
                 libhts-dev && \
    apt-get -qqy clean && \
    rm -rf Eagle \
           /var/lib/apt/lists/*
```

Dockerfile for IMPUTE5:
```
FROM ubuntu:21.04
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get -qqy update --fix-missing && \
    apt-get -qqy install --no-install-recommends \
                 wget \
                 bcftools \
                 unzip && \
    wget -O impute5_v1.1.4.zip --no-check-certificate https://www.dropbox.com/sh/mwnceyhir8yze2j/AAAPJVJv3kI2glXGDdc3sSHga/impute5_v1.1.4.zip?dl=0 && \
    unzip -jd /usr/bin/ impute5_v1.1.4.zip impute5_v1.1.4/impute5_1.1.4_static impute5_v1.1.4/imp5Converter_1.1.4_static && \
    ln -s impute5_1.1.4_static /usr/bin/impute5 && \
    ln -s imp5Converter_1.1.4_static /usr/bin/imp5Converter && \
    apt-get -qqy purge --auto-remove --option APT::AutoRemove::RecommendsImportant=false \
                 wget \
                 unzip && \
    apt-get -qqy clean && \
    rm -rf impute5_v1.1.4.zip \
           /var/lib/apt/lists/*
```

Dockerfile for Beagle:
```
FROM ubuntu:21.04
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get -qqy update --fix-missing && \
    apt-get -qqy install --no-install-recommends \
                 wget \
                 bcftools \
                 default-jre-headless && \
    wget -P /usr/bin/ --no-check-certificate http://faculty.washington.edu/browning/beagle/beagle.18May20.d20.jar && \
    wget -P /usr/bin http://faculty.washington.edu/browning/beagle/bref3.18May20.d20.jar && \
    ln -s beagle.18May20.d20.jar /usr/bin/beagle.jar && \
    ln -s bref3.18May20.d20.jar /usr/bin/bref3.jar && \
    apt-get -qqy purge --auto-remove --option APT::AutoRemove::RecommendsImportant=false \
                 wget && \
    apt-get -qqy autoremove && \
    apt-get -qqy clean && \
    rm -rf /var/lib/apt/lists/*
```

Acknowledgements
================

This work is supported by NIH grant <a href="http://grantome.com/grant/NIH/R01-HG006855">R01 HG006855</a>, NIH grant <a href="http://grantome.com/grant/NIH/R01-MH104964">R01 MH104964</a>, and the Stanley Center for Psychiatric Research. This work was spurred from discussions within the COVID-19 host genetics and clonal hematopoiesis working groups: special thanks to Mark Daly, Andrea Ganna, Philip Awadalla, Pradeep Natarajan, Kelly Bolton, Alexander Bick and many others
