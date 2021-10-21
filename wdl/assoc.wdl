version development

## Copyright (c) 2021 Giulio Genovese
##
## Version 2021-10-20
##
## Contact Giulio Genovese <giulio.genovese@gmail.com>
##
## This WDL workflow runs association analyses with REGENIE and PLINK2
##
## Cromwell version support
## - Successfully tested on v70
##
## Distributed under terms of the MIT License

struct Reference {
  File fasta_fai
  Int n_chrs
  File genetic_map_file
  String pca_exclusion_regions
}

workflow assoc {
  input {
    String sample_set_id
    Float max_win_size_cm_step2 = 20.0
    File sample_tsv_file
    File? remove_samples_file
    Int min_mac = 10
    Float min_maf = 0.01
    File? covar_tsv_file
    File pheno_tsv_file
    String dosage_field = "DS"
    String? sex_specific # male female
    String space_character = "_"
    Boolean binary = true
    Int min_case_count = 20
    Int min_sex_count = 20
    Int bsize = 500
    Boolean loocv = true
    String? regenie_step0_extra_args
    String? regenie_step1_extra_args
    String? regenie_step2_extra_args
    String? plink_extra_args
    Boolean step1 = true
    Boolean pgt_output = false
    Boolean pca = false
    Boolean step2 = true
    Boolean cis = false
    Int pca_ndim = 20
    File? input_loco_lst
    String? input_loco_path
    File? input_firth_lst
    String? input_firth_path

    String ref_name = "GRCh38"
    String? ref_path
    String? ref_fasta_fai
    Int? ref_n_chrs
    String? genetic_map_file
    String? pca_exclusion_regions

    File? batch_step1_tsv_file # batch_id path pgt_vcf pgt_vcf_index
    String? data_step1_path
    File? batch_step2_tsv_file # batch_id path chr1_imp_vcf chr1_imp_vcf_index chr2_imp_vcf chr2_imp_vcf_index ...
    String? data_step2_path
    String basic_bash_docker = "debian:stable-slim"
    String pandas_docker = "amancevice/pandas:slim"
    String docker_repository = "us.gcr.io/mccarroll-mocha"
    String bcftools_docker = "bcftools:1.13-20211015"
    String regenie_docker = "regenie:1.13-20211015"
  }

  String docker_repository_with_sep = docker_repository + if docker_repository != "" && docker_repository == sub(docker_repository, "/$", "") then "/" else ""

  String ref_path_with_sep = select_first([ref_path, ""]) + if defined(ref_path) && select_first([ref_path]) == sub(select_first([ref_path]), "/$", "") then "/" else ""
  Reference ref = object {
    fasta_fai: ref_path_with_sep + select_first([ref_fasta_fai, if ref_name == "GRCh38" then "GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai" else if ref_name == "GRCh37" then "human_g1k_v37.fasta.fai" else None]),
    n_chrs: select_first([ref_n_chrs, 23]),
    genetic_map_file: ref_path_with_sep + select_first([genetic_map_file, if ref_name == "GRCh38" then "genetic_map_hg38_withX.txt.gz" else if ref_name == "GRCh37" then "genetic_map_hg19_withX.txt.gz" else None]),
    pca_exclusion_regions: select_first([pca_exclusion_regions, if ref_name == "GRCh38" then "5:43999898-52204166,6:24999772-33532223,8:8142478-12142491,11:44978449-57232526" else if ref_name == "GRCh37" then "5:44000000-51500000,6:25000000-33500000,8:8000000-12000000,11:45000000-57000000" else None]),
  }

  Array[Array[String]] ref_fasta_fai_tbl = transpose(read_tsv(ref.fasta_fai))
  scatter (idx in range(ref.n_chrs)) {
    String chrs = ref_fasta_fai_tbl[0][idx]
    Int lens = ref_fasta_fai_tbl[1][idx]
  }

  call prune_file {
    input:
      sample_tsv_file = sample_tsv_file,
      pheno_tsv_file = pheno_tsv_file,
      remove_samples_file = remove_samples_file,
      covar_tsv_file = covar_tsv_file,
      sex_specific = sex_specific,
      space_character = space_character,
      binary = binary,
      min_case_count = min_case_count,
      min_sex_count = min_sex_count,
      filebase = sample_set_id,
      docker = basic_bash_docker
  }

  # REGENIE step 1

  if (step1 || pca) {
    call ref_scatter as ref_scatter_step1 {
      input:
        chrs = chrs,
        lens = lens,
        genetic_map_file = ref.genetic_map_file,
        max_win_size_cm = 300.0, # until regenie updates, step 1 cannot be parallelized beyond the 23 chromosomes
        overlap_size_cm = 0.0,
        docker = pandas_docker
    }

    # read table with batches information (scatter could be avoided if there was a tail() function)
    Array[Array[String]] batch_step1_tsv = read_tsv(select_first([batch_step1_tsv_file]))
    Int n_batches_step1 = length(batch_step1_tsv)-1
    scatter (idx in range(n_batches_step1)) { Array[String] batch_step1_tsv_rows = batch_step1_tsv[(idx+1)] }
    Map[String, Array[String]] batch_step1_tbl = as_map(zip(batch_step1_tsv[0], transpose(batch_step1_tsv_rows)))

    # compute data paths for each batch, if available (scatter could be avoided if there was a contains_key() function)
    scatter (key in keys(batch_step1_tbl)) { Boolean? is_key_step1_equal_path = if key == "path" then true else None }
    scatter (idx in range(n_batches_step1)) {
      String data_step1_paths = select_first([data_step1_path, if length(select_all(is_key_step1_equal_path))>0 then batch_step1_tbl["path"][idx] else ""])
      String data_step1_paths_with_sep = data_step1_paths + (if data_step1_paths == "" || sub(data_step1_paths, "/$", "") != data_step1_paths then "" else "/")
    }

    scatter (idx in range(n_batches_step1)) {
      call vcf_scatter as pgt_scatter {
        input:
          vcf_file = data_step1_paths_with_sep[idx] + batch_step1_tbl["pgt_vcf"][idx],
          intervals_bed = ref_scatter_step1.intervals_bed,
          docker = docker_repository_with_sep + bcftools_docker
      }
    }

    Array[Array[File]] interval_slices = transpose(pgt_scatter.vcf_files)
    scatter (idx in range(length(interval_slices))) {
      if (length(interval_slices[idx])>1) {
        call vcf_merge as pgt_merge {
          input:
            vcf_files = interval_slices[idx],
            filebase = sample_set_id + "." + idx,
            docker = docker_repository_with_sep + bcftools_docker
        }
      }

      call pgt_prune {
        input:
          vcf_file = select_first([pgt_merge.vcf_file, interval_slices[idx][0]]),
          sample_tsv_file = sample_tsv_file,
          remove_samples_file = remove_samples_file,
          sex_specific = sex_specific,
          space_character = space_character,
          min_mac = min_mac,
          min_maf = min_maf,
          docker = docker_repository_with_sep + regenie_docker
      }

      Int n_smpls = select_first([pgt_merge.n_smpls, pgt_scatter.n_smpls])

      if (step1) {
        call regenie_step0 {
          input:
            idx = idx,
            n_phenos = length(prune_file.pheno_names),
            n_smpls = n_smpls,
            n_markers = pgt_prune.n_markers,
            bed_file = pgt_prune.bed_file,
            bim_file = pgt_prune.bim_file,
            fam_file = pgt_prune.fam_file,
            covar_file = prune_file.covar,
            pheno_file = prune_file.pheno,
            binary = binary,
            bsize = bsize,
            loocv = loocv,
            regenie_step0_extra_args = regenie_step0_extra_args,
            filebase = sample_set_id,
            docker = docker_repository_with_sep + regenie_docker
        }
      }
    }

    call pgt_concat {
      input:
        bed_files = pgt_prune.bed_file,
        bim_files = pgt_prune.bim_file,
        fam_files = pgt_prune.fam_file,
        filebase = sample_set_id + ".prune",
        docker = basic_bash_docker
    }

    if (pca) {
      call plink_pca {
        input:
          n_smpls = n_smpls[0],
          n_markers = pgt_concat.n_markers,
          pca_ndim = pca_ndim,
          bed_file = pgt_concat.bed_file,
          bim_file = pgt_concat.bim_file,
          fam_file = pgt_concat.fam_file,
          exclusion_regions = ref.pca_exclusion_regions,
          filebase = sample_set_id,
          docker = docker_repository_with_sep + regenie_docker
      }
    }

    if (step1) {
      Array[Array[File]] l0_files = transpose(select_all(regenie_step0.l0_files))
      scatter (idx in range(length(l0_files))) {
        call regenie_step1 {
          input:
            pheno_name = prune_file.pheno_names[idx],
            n_smpls = n_smpls[0],
            n_markers = pgt_concat.n_markers,
            bed_file = pgt_concat.bed_file,
            bim_file = pgt_concat.bim_file,
            fam_file = pgt_concat.fam_file,
            covar_file = prune_file.covar,
            pheno_file = prune_file.pheno,
            n_markers_array = pgt_prune.n_markers,
            l0_files = l0_files[idx],
            binary = binary,
            bsize = bsize,
            loocv = loocv,
            regenie_step1_extra_args = regenie_step1_extra_args,
            filebase = sample_set_id,
            docker = docker_repository_with_sep + regenie_docker
        }
        String? loco_lines = if regenie_step1.loco_line == "" then None else regenie_step1.loco_line
        File? loco_file = if defined(loco_lines) then regenie_step1.loco_file else None
        if (binary) {
          String? firth_lines = if select_first([regenie_step1.firth_line]) == "" then None else select_first([regenie_step1.firth_line])
          File? firth_file = if defined(firth_lines) then select_first([regenie_step1.firth_file]) else None
        }
      }
      # unnecessary task for compatibility with Terra https://support.terra.bio/hc/en-us/community/posts/360071465631-write-lines-write-map-write-tsv-write-json-fail-when-run-in-a-workflow-rather-than-in-a-task
      call serialize_lines as loco_lst { input: lines = select_all(loco_lines), filename = sample_set_id + "_pred.list", docker = basic_bash_docker }
      if (binary) {
        call serialize_lines as firth_lst { input: lines = select_all(firth_lines), filename = sample_set_id + "_firth.list", docker = basic_bash_docker }
      }
    }
  }

  if (!step1) {
    if (defined(input_loco_lst) && defined(input_loco_path)) {
      scatter (line in read_lines(select_first([input_loco_lst]))) {
        File input_loco_files = sub(select_first([input_loco_path]), "/$", "") + '/' + sub(line, "^.*[ \t]", "")
      }
    }

    if (binary && defined(input_firth_lst) && defined(input_firth_path)) {
      scatter (line in read_lines(select_first([input_firth_lst]))) {
        File input_firth_files = sub(select_first([input_firth_path]), "/$", "") + '/' + sub(line, "^.*[ \t]", "")
      }
    }
  }

  # REGENIE step 2

  if (step2 || cis) {
    call ref_scatter as ref_scatter_step2 {
      input:
        chrs = chrs,
        lens = lens,
        genetic_map_file = ref.genetic_map_file,
        max_win_size_cm = max_win_size_cm_step2,
        overlap_size_cm = 0.0,
        docker = pandas_docker
    }

    Array[Array[String]] intervals_tbl = transpose(read_tsv(ref_scatter_step2.intervals_bed))
    # this is a trick to table how many intervals you will use for each chromosome
    Map[String, Array[Int]] chr_map = collect_by_key(zip(intervals_tbl[0], range(length(intervals_tbl[0]))))

    # read table with batches information (scatter could be avoided if there was a tail() function)
    Array[Array[String]] batch_step2_tsv = read_tsv(select_first([batch_step2_tsv_file]))
    Int n_batches_step2 = length(batch_step2_tsv)-1
    scatter (idx in range(n_batches_step2)) { Array[String] batch_step2_tsv_rows = batch_step2_tsv[(idx+1)] }
    Map[String, Array[String]] batch_step2_tbl = as_map(zip(batch_step2_tsv[0], transpose(batch_step2_tsv_rows)))

    # compute data paths for each batch, if available (scatter could be avoided if there was a contains_key() function)
    scatter (key in keys(batch_step2_tbl)) { Boolean? is_key_step2_equal_path = if key == "path" then true else None }
    scatter (idx in range(n_batches_step2)) {
      String data_step2_paths = select_first([data_step2_path, if length(select_all(is_key_step2_equal_path))>0 then batch_step2_tbl["path"][idx] else ""])
      String data_step2_paths_with_sep = data_step2_paths + (if data_step2_paths == "" || sub(data_step2_paths, "/$", "") != data_step2_paths then "" else "/")
    }

    scatter (p in cross(range(n_batches_step2), range(ref.n_chrs))) {
      File imp_vcf_file = data_step2_paths_with_sep[p.left] + batch_step2_tbl[("chr" + sub(chrs[p.right], "^chr", "") + "_imp_vcf")][p.left]
      if (length(chr_map[(chrs[p.right])]) > 1) {
        call vcf_scatter {
          input:
            vcf_file = imp_vcf_file,
            intervals_bed = ref_scatter_step2.intervals_bed,
            chr = chrs[p.right],
            dosage_field = dosage_field,
            docker = docker_repository_with_sep + bcftools_docker
        }
      }
      Int cross_idx = p.right
      Array[File] scatter_vcf_files = select_first([vcf_scatter.vcf_files, [imp_vcf_file]])
    }

    Map[Int, Array[Array[File]]] idx2vcf_files = collect_by_key(zip(cross_idx, scatter_vcf_files))
    scatter (idx in range(ref.n_chrs)) { Array[Array[File]] slices_vcf_files = transpose(idx2vcf_files[idx]) }
    Array[Array[File]] matrix_vcf_files = flatten(slices_vcf_files)

    if (step2) {
      # generate list of expected output association files
      scatter (line in read_lines(select_first([firth_lst.file, loco_lst.file, input_firth_lst, input_loco_lst]))) {
        String regenie_suffix = sub(line, " .*$", "") + ".regenie.gz" # https://github.com/broadinstitute/cromwell/issues/5549
      }
    }

    # merging has to happen at the VCF level as plink2 does not currently merge pgen files
    # https://www.cog-genomics.org/plink/2.0/data#pmerge
    scatter (idx in range(length(matrix_vcf_files))) {
      if (length(matrix_vcf_files[idx])>1) {
        call vcf_merge {
          input:
            vcf_files = matrix_vcf_files[idx],
            filebase = sample_set_id + "." + idx,
            docker = docker_repository_with_sep + bcftools_docker
        }
      }

      call vcf2pgen {
        input:
          vcf_file = select_first([vcf_merge.vcf_file, matrix_vcf_files[idx][0]]),
          sample_tsv_file = sample_tsv_file,
          remove_samples_file = remove_samples_file,
          dosage_field = dosage_field,
          ref_name = ref_name,
          sex_specific = sex_specific,
          space_character = space_character,
          docker = docker_repository_with_sep + regenie_docker
      }

      if (step2) {
        call regenie_step2 {
          input:
            n_smpls = vcf2pgen.n_smpls,
            n_markers = vcf2pgen.n_markers,
            pgen_file = vcf2pgen.pgen_file,
            pvar_file = vcf2pgen.pvar_file,
            psam_file = vcf2pgen.psam_file,
            covar_file = prune_file.covar,
            pheno_file = prune_file.pheno,
            suffix = select_first([regenie_suffix]), # https://github.com/broadinstitute/cromwell/issues/5549
            binary = binary,
            bsize = bsize,
            regenie_step2_extra_args = regenie_step2_extra_args,
            loco_lst = select_first([loco_lst.file, input_loco_lst]),
            loco_files = select_first([loco_files, input_loco_files]),
            firth_lst = if binary then select_first([firth_lst.file, input_firth_lst]) else None,
            firth_files = if binary then select_first([firth_files, input_firth_files]) else None,
            docker = docker_repository_with_sep + regenie_docker
        }
      }
    }

    if (step2) {
      Array[Array[File]] regenie_matrix_files = transpose(select_all(regenie_step2.regenie_files))
      scatter (idx in range(length(regenie_matrix_files))) {
        call assoc_concat as regenie_concat {
          input:
            assoc_files = regenie_matrix_files[idx],
            filebase = sample_set_id + "_" + select_first([regenie_suffix])[idx],
            docker = docker_repository_with_sep + bcftools_docker
        }
      }
    }

    if (cis) {
      Array[String] lines = if defined(loco_lst.file) || defined(input_loco_lst) then read_lines(select_first([loco_lst.file, input_loco_lst])) else prune_file.pheno_names
      scatter (idx in range(length(lines))) {
        String plink_pheno_names = if defined(loco_lst.file) || defined(input_loco_lst) then sub(lines[idx], " .*$", "") else lines[idx]
        String plink_pheno_chrs = sub(sub(sub(plink_pheno_names, "_.*$", ""), "[pq]*$", ""), "Y", "X")
        Int? cis_idx = if sub(plink_pheno_chrs, "[0-9X]+", "") == "" then idx else None
      }
      # this map, given a chromosome (1, 2, ..., X), returns the indexes of the intervals for that chromosomes
      scatter (chr in intervals_tbl[0]) { String chr_string = sub(chr, "^chr", "") }
      Map[String, Array[Int]] chr2idx = collect_by_key(zip(chr_string, range(length(intervals_tbl[0]))))
      # the following code checks that the chromosome name is in the list of available chromosomes but is not perfect
      scatter (idx in select_all(cis_idx)) {
        Array[Pair[Int, Int]] pheno_interval_pairs = cross([idx], chr2idx[(plink_pheno_chrs[idx])])
      }
      # maybe I should test whether the interval falls under the event or not
      scatter (p in flatten(pheno_interval_pairs)) {
        Int pheno_idx = p.left
        call plink_glm {
          input:
            chr_num = sub(plink_pheno_chrs[p.left], "X", "23"),
            pheno_name = plink_pheno_names[p.left],
            n_smpls = vcf2pgen.n_smpls[p.right],
            n_markers = vcf2pgen.n_markers[p.right],
            pgen_file = vcf2pgen.pgen_file[p.right],
            pvar_file = vcf2pgen.pvar_file[p.right],
            psam_file = vcf2pgen.psam_file[p.right],
            covar_file = prune_file.covar,
            pheno_file = prune_file.pheno,
            binary = binary,
            loco_file =  if defined(loco_files) || defined(input_loco_files) then select_first([loco_files, input_loco_files])[p.left] else None,
            plink_extra_args = plink_extra_args,
            docker = docker_repository_with_sep + regenie_docker
        }
      }

      Map[Int, Array[File]] idx2assoc_files = collect_by_key(zip(pheno_idx, plink_glm.assoc_file))
      scatter (idx in select_all(cis_idx)) {
        call assoc_concat as plink_concat {
          input:
            assoc_files = idx2assoc_files[idx],
            zst = true,
            filebase = sample_set_id + "_" + plink_pheno_names[idx] + ".glm." + (if binary then "logistic.hybrid" else "linear") + ".gz",
            docker = docker_repository_with_sep + regenie_docker
        }
      }
    }
  }

  output {
    File? bed_file = if pgt_output then pgt_concat.bed_file else None
    File? bim_file = if pgt_output then pgt_concat.bim_file else None
    File? fam_file = if pgt_output then pgt_concat.fam_file else None
    File? eigenvec_file = plink_pca.eigenvec_file
    File? eigenval_file = plink_pca.eigenval_file
    File? loco_lst_file = loco_lst.file
    Array[File]? loco_files = if defined(loco_file) then select_all(select_first([loco_file])) else None
    File? firth_lst_file = firth_lst.file
    Array[File]? firth_files = if defined(firth_file) then select_all(select_first([firth_file])) else None
    Array[File]? regenie_files = regenie_concat.file
    Array[File]? regenie_indexes = regenie_concat.index
    Array[File]? plink_files = plink_concat.file
    Array[File]? plink_indexes = plink_concat.index
  }
}

task get_n {
  input {
    File file

    String docker
    Int cpu = 1
    Int disk_size = 10
    Float memory = 3.5
    Int preemptible = 1
    Int maxRetries = 0
  }

  command <<<
    set -euo pipefail
    mv "~{file}" .
    grep -v ^# "~{basename(file)}" | wc -l
    rm "~{basename(file)}"
  >>>

  output {
    Int n = read_int(stdout())
  }

  runtime {
    docker: docker
    cpu: cpu
    disks: "local-disk " + disk_size + " HDD"
    memory: memory + " GiB"
    preemptible: preemptible
    maxRetries: maxRetries
  }
}

# use of !(a!=b) due to bug Cromwell team will not fix: https://github.com/broadinstitute/cromwell/issues/5602
task prune_file {
  input {
    File sample_tsv_file
    File pheno_tsv_file
    File? remove_samples_file
    File? covar_tsv_file
    String? sex_specific
    String space_character
    Boolean binary
    Int min_case_count
    Int min_sex_count
    String filebase

    String docker
    Int cpu = 1
    Int disk_size = 10
    Float memory = 3.5
    Int preemptible = 1
    Int maxRetries = 0
  }

  command <<<
    set -euo pipefail
    echo "~{sep="\n" select_all([sample_tsv_file, pheno_tsv_file, remove_samples_file, covar_tsv_file])}" | \
      tr '\n' '\0' | xargs -0 mv -t .
    awk -F"\t" 'NR==1 {for (i=1; i<=NF; i++) f[$i] = i}
      NR>1 {sex=substr($(f["computed_gender"]),1,1); if (toupper(sex)=="M" || sex==1) printf "%s\t1\n",$(f["sample_id"])}' \
      "~{basename(sample_tsv_file)}" > "~{filebase}.male"
    awk -F"\t" 'NR==1 {for (i=1; i<=NF; i++) f[$i] = i}
      NR>1 {sex=substr($(f["computed_gender"]),1,1); if (toupper(sex)=="F" || sex==2) printf "%s\t2\n",$(f["sample_id"])}' \
      "~{basename(sample_tsv_file)}" > "~{filebase}.female"
    cat "~{filebase + "." + if defined(sex_specific) then select_first([sex_specific]) else "male\" \"" + filebase + ".female"}" | \
      ~{if defined(remove_samples_file) then
        "awk -F\"\\t\" 'NR==FNR {x[$1]++} NR>FNR && !($1 in x)' \"" + basename(select_first([remove_samples_file])) + "\" - | \\\n"
        else "" + if defined(covar_tsv_file) then
        "awk -F\"\\t\" 'NR==FNR {x[$1]++} NR>FNR && $1 in x' \"" + basename(select_first([covar_tsv_file])) + "\" - | \\\n"
        else ""}awk -F"\t" 'NR==FNR {x[$1]++} NR>FNR && (FNR==1 || $1 in x)' - "~{basename(pheno_tsv_file)}" > "~{filebase}.tmp"
    cat "~{filebase}.male" "~{filebase}.female" | \
      awk -F"\t" 'NR==FNR {sex[$1]=$2} NR>FRN && FNR==1 {for (i=2; i<=NF; i++) pheno[i] = $i}
      NR>FNR && FNR>1 {for (i=2; i<=NF; i++) {if ($i==0) ctrls[i]++; if ($i==1) cases[i]++
      if (sex[$1]==1 && $i!="NA") males[i]++; if (sex[$1]==2 && $i!="NA") females[i]++}}
      END {for (i in pheno); printf "%s\t%d\t%d\t%d\t%d\n",pheno[i],ctrls[i],cases[i],males[i],females[i]}' \
      - "~{filebase}.tmp" > "~{filebase}.cnt"
    awk -F"\t" 'NR==FNR ~{if binary then "&& $2>=" + min_case_count + " && $3>=" + min_case_count else ""} && $~{
      if defined(sex_specific) && !(select_first([sex_specific]) != "male") then "4"
      else if defined(sex_specific) && !(select_first([sex_specific]) != "female") then "5"
      else "4>=" + min_sex_count + " && $5"}>=~{min_sex_count} {keep[$1]++}
      NR>FNR {if (FNR==1) {for (i=2; i<=NF; i++) if ($i in keep) col[j++]=i; printf "FID\tIID"}
      else {gsub(" ","~{space_character}",$1); printf "0\t%s",$1} for (i=0; i<j; i++) printf "\t%s",$col[i]; printf "\n"}' \
      "~{filebase}.cnt" "~{filebase}.tmp" > "~{filebase}.phe"
    ~{if defined(covar_tsv_file) then "cat \"" + filebase + "." +
      (if defined(sex_specific) then select_first([sex_specific]) else "male\" + \"" + filebase + ".female") + "\" | \\\n" +
      "  awk -F\"\\t\" 'NR==FNR {sex[$1]=$2} NR>FNR {if (FNR==1) {col=\"sex\"; $1=\"FID\\tIID\"}\n" +
      "  if (FNR>1) {col=sex[$1]; gsub(\" \",\"" + space_character + "\",$1); $1=\"0\\t\"$1} printf \"%s\\t%s\\n\",$0,col}' \\\n" +
      "  - \"" + basename(select_first([covar_tsv_file])) + "\" > \"" + filebase + ".cov\"\n"
      else ""}head -n1 "~{filebase}.phe" | cut -f3- | tr '\t' '\n'
    rm "~{filebase}.male" "~{filebase}.female" "~{filebase}.tmp" "~{filebase}.cnt"
    echo "~{sep="\n" select_all([sample_tsv_file, pheno_tsv_file, remove_samples_file, covar_tsv_file])}" | \
      sed 's/^.*\///' | tr '\n' '\0' | xargs -0 rm
  >>>

  output {
    Array[String] pheno_names = read_lines(stdout())
    File pheno = filebase + ".phe"
    File? covar = if defined(covar_tsv_file) then filebase + ".cov" else None
  }

  runtime {
    docker: docker
    cpu: cpu
    disks: "local-disk " + disk_size + " HDD"
    memory: memory + " GiB"
    preemptible: preemptible
    maxRetries: maxRetries
  }
}

task ref_scatter {
  input {
    Array[String] chrs
    Array[String] lens
    File genetic_map_file
    Float max_win_size_cm
    Float overlap_size_cm

    String docker
    Int cpu = 1
    Int disk_size = 10
    Float memory = 3.5
    Int preemptible = 1
    Int maxRetries = 0
  }

  command <<<
    set -euo pipefail
    mv "~{genetic_map_file}" .
    chrs=~{write_lines(chrs)}
    lens=~{write_lines(lens)}
    paste -d $'\t' $chrs $lens > chr2len.tsv
    python3 <<CODE
    import sys, pandas as pd, numpy as np
    chr2len = {}
    with open('chr2len.tsv') as f:
      for line in f:
        (key, val) = line.split('\t')
        chr2len[key] = int(val)
    df_map = pd.read_csv('~{basename(genetic_map_file)}', delim_whitespace = True, header = 0, names = ['CHR', 'POS' ,'RATE', 'CM'])
    df_out = {}
    for chr, df_group in df_map.groupby('CHR'):
      fai_chr = str(chr) if str(chr) in chr2len else 'chr' + str(chr) if 'chr' + str(chr) in chr2len else 'X' if 'X' in chr2len else 'chrX' if 'chrX' in chr2len else None
      if fai_chr:
        chr_cm_len = max(df_group['CM'])
        n_win = np.ceil((chr_cm_len - ~{overlap_size_cm})/(~{max_win_size_cm} - ~{overlap_size_cm}))
        win_size = (chr_cm_len - ~{overlap_size_cm}) / n_win + ~{overlap_size_cm}
        cm_begs = (win_size - ~{overlap_size_cm}) * np.arange(1, n_win)
        cm_ends = (win_size - ~{overlap_size_cm}) * np.arange(1, n_win) + ~{overlap_size_cm}
        pos_begs = np.concatenate(([1], np.interp(cm_begs, df_group['CM'], df_group['POS'], period = np.inf).astype(int)))
        pos_ends = np.concatenate((np.interp(cm_ends, df_group['CM'], df_group['POS'], period = np.inf).astype(int), [chr2len[fai_chr]]))
        df_out[fai_chr] = pd.DataFrame.from_dict({'CHR': fai_chr, 'BEG': pos_begs, 'END': pos_ends})
    df = pd.concat([df_out[fai_chr] for fai_chr in chr2len.keys()])
    df[['CHR', 'BEG', 'END']].to_csv('ref_scatter.bed', sep='\t', header = False, index = False)
    CODE
    rm chr2len.tsv
    rm "~{basename(genetic_map_file)}"
  >>>

  output {
    File intervals_bed = "ref_scatter.bed"
  }

  runtime {
    docker: docker
    cpu: cpu
    disks: "local-disk " + disk_size + " HDD"
    memory: memory + " GiB"
    preemptible: preemptible
    maxRetries: maxRetries
  }
}

task vcf_scatter {
  input {
    File vcf_file
    File intervals_bed
    String? chr
    String? dosage_field

    String docker
    Int? cpu_override
    Int? disk_size_override
    Float? memory_override
    Int preemptible = 1
    Int maxRetries = 0
  }

  Float vcf_size = size(vcf_file, "GiB")
  Int disk_size = select_first([disk_size_override, ceil(10.0 + 3.0 * vcf_size)])
  Float memory = select_first([memory_override, 3.5])
  Int cpu = select_first([cpu_override, if memory > 6.5 then 2 * ceil(memory / 13) else 1])
  String filebase = basename(basename(vcf_file, ".bcf"), ".vcf.gz")

  command <<<
    set -euo pipefail
    mv "~{vcf_file}" .
    mv "~{intervals_bed}" .
    ~{if defined(chr) then
      "mv \"" + basename(intervals_bed) + "\" \"" + basename(intervals_bed, ".bed") + ".all.bed\"\n" +
      "awk -v chr=\"" + chr + "\" '$1==chr' \"" + basename(intervals_bed, ".bed") + ".all.bed\" > \"" + basename(intervals_bed) + "\"\n" +
      "rm \"" + basename(intervals_bed, ".bed") + ".all.bed\""
      else ""}
    bcftools query --list-samples "~{basename(vcf_file)}" | wc -l > n_smpls.int
    awk -F"\t" '{print $1":"$2"-"$3"\t"NR-1}' "~{basename(intervals_bed)}" > regions.lines
    bcftools annotate \
      --no-version \
      --output-type u \
      --remove ID,QUAL,FILTER,INFO,^FMT/GT~{if defined(dosage_field) then ",^FMT/DS" else ""} \
      ~{if cpu > 1 then "--threads " + (cpu - 1) else ""} \
      "~{basename(vcf_file)}" | \
    bcftools norm \
      --no-version \
      --rm-dup exact \
      ~{if cpu > 1 then "--threads " + (cpu - 1) else ""} | \
    bcftools +scatter \
      --no-version \
      --output-type b \
      --output vcfs \
      ~{if cpu > 1 then "--threads " + (cpu - 1) else ""} \
      --scatter-file regions.lines \
      --prefix "~{filebase}."
    while read reg i; do
      bcftools query --format "\n" "vcfs/~{filebase}.$i.bcf" | wc -l
      bcftools index --force "vcfs/~{filebase}.$i.bcf"
    done < regions.lines > n_markers.lines
    cut -f2 regions.lines | sed 's/^/vcfs\/~{filebase}./;s/$/.bcf/'
    cut -f2 regions.lines | sed 's/^/vcfs\/~{filebase}./;s/$/.bcf.csi/' > vcf_idxs.lines
    rm "~{basename(vcf_file)}"
    rm "~{basename(intervals_bed)}"
    rm regions.lines
  >>>

  output {
    Int n_smpls = read_int("n_smpls.int")
    Array[Int] n_markers = read_lines("n_markers.lines")
    Directory vcfs = "vcfs"
    Array[File] vcf_files = read_lines(stdout())
    Array[File] vcf_idxs = read_lines("vcf_idxs.lines")
  }

  runtime {
    docker: docker
    cpu: cpu
    disks: "local-disk " + disk_size + " HDD"
    memory: memory + " GiB"
    preemptible: preemptible
    maxRetries: maxRetries
  }
}

task vcf_merge {
  input {
    Array[File]+ vcf_files
    String filebase
    Boolean wgs = false

    String docker
    Int? cpu_override
    Int? disk_size_override
    Float? memory_override
    Int preemptible = 1
    Int maxRetries = 0
  }

  Float vcf_size = size(vcf_files, "GiB")
  Int disk_size = select_first([disk_size_override, ceil(10 + 2.0 * vcf_size)])
  Float memory = select_first([memory_override, 3.5])
  Int cpu = select_first([cpu_override, if memory > 6.5 then 2 * ceil(memory / 13) else 1])

  command <<<
    set -euo pipefail
    ~{if length(vcf_files) > 1 then "vcf_files=" else ""}~{if length(vcf_files) > 1 then write_lines(vcf_files) else ""}
    ~{if length(vcf_files) > 1 then "cat $vcf_files | tr '\\n' '\\0' | xargs -0 mv -t .\n" +
      "sed -i 's/^.*\\///' $vcf_files\n" +
      "bcftools merge \\\n" +
      "  --no-version \\\n" +
      "  --output-type b \\\n" +
      "  --output \"" + filebase + ".bcf\" \\\n" +
      (if wgs then "--missing-to-ref \\\n" else "") +
      "  --file-list $vcf_files \\\n" +
      "  --merge none \\\n" +
      "  --no-index \\\n" +
      (if cpu > 1 then "  --threads " + (cpu - 1) else "")
      else "mv \"" + vcf_files[0] + "\" \"" + filebase + ".bcf\""}
    bcftools index --force "~{filebase}.bcf"
    bcftools query --list-samples "~{filebase}.bcf" | wc -l
    ~{if length(vcf_files) > 1 then "cat $vcf_files | tr '\\n' '\\0' | xargs -0 rm" else ""}
  >>>

  output {
    File vcf_file = filebase + ".bcf"
    File vcf_idx = filebase + ".bcf.csi"
    Int n_smpls = read_int(stdout())
  }

  runtime {
    memory: memory + " GiB"
    disks: "local-disk " + disk_size + " HDD"
    cpu: cpu
    docker: docker
    preemptible: preemptible
    maxRetries: maxRetries
  }
}

# this command needs PLINK 1.9 as conversion from VCF using PLINK 2.0 is inefficient
# https://groups.google.com/g/plink2-users/c/hsByNOklyA0
# the U sex needs to be encoded as 0 as this is the only accepted value for PLINK 1.9
# https://groups.google.com/g/plink2-users/c/z7YJYa677NQ
# --compression-level is not working in BCFtools 1.13 https://github.com/samtools/bcftools/issues/1528
# use of !(a!=b) due to bug Cromwell team will not fix: https://github.com/broadinstitute/cromwell/issues/5602
task pgt_prune {
  input {
    File vcf_file
    File? sample_tsv_file
    File? remove_samples_file
    String? sex_specific
    String space_character
    Int min_mac
    Float min_maf

    String docker
    Int cpu = 1
    Int? disk_size_override
    Float memory = 3.5
    Int preemptible = 1
    Int maxRetries = 0
  }

  String filebase = basename(basename(vcf_file, ".bcf"), ".vcf.gz")
  Float vcf_size = size(vcf_file, "GiB")
  Int disk_size = select_first([disk_size_override, ceil(10.0 + 3.0 * vcf_size)])

  command <<<
    set -euo pipefail
    echo "~{sep="\n" select_all([vcf_file, sample_tsv_file, remove_samples_file])}" | \
      tr '\n' '\0' | xargs -0 mv -t .
    ~{if defined(sample_tsv_file) then
        "awk 'NR==1 {for (i=1; i<=NF; i++) f[$i] = i}\n" +
        "  NR>1 {id=$(f[\"sample_id\"]); gsub(\" \",\"" + space_character + "\",id);\n" +
        "  print 0,id,toupper(substr($(f[\"computed_gender\"]),1,1))}' \"" + basename(select_first([sample_tsv_file])) + "\" | \\\n" +
        "  sed 's/U$/0/' > \"" + filebase + ".sex\"\n"
      else ""}~{if defined (sex_specific) && !(select_first([sex_specific]) != "male") then
        "  awk '$3==\"M\" || $3==1' \"" + filebase + ".sex\" > " + select_first([sex_specific]) + ".fam\n"
      else if defined (sex_specific) && !(select_first([sex_specific]) != "female") then
        "  awk '$3==\"F\" || $3==2' \"" + filebase + ".sex\" > " + select_first([sex_specific]) + ".fam\n"
      else ""}~{if defined(remove_samples_file) then "bcftools view --no-version -Ou --samples-file \"^" + basename(select_first([remove_samples_file])) + "\" --force-samples \"" + basename(vcf_file) + "\" |\n  "
      else ""}bcftools +fill-tags --no-version -Ou --include 'sum(AC)>=~{min_mac} && AN-sum(AC)>=~{min_mac} && MAF>=~{min_maf}' ~{if defined(remove_samples_file) then""
      else "\"" + basename(vcf_file) + "\" "}-- --tags AC,AN,MAF | \
      bcftools +add-variantkey --no-version -Ou | \
      bcftools annotate --no-version -Ob --set-id +'%VKX' --remove FILTER,INFO,^FMT/GT | \
#      bcftools view --no-version -Ob --compression-level 0 | \
      plink1.9 \
        --bcf /dev/stdin \
        ~{if defined(sample_tsv_file) then "--update-sex \"" + filebase + ".sex\"" else ""} \
        ~{if defined(sex_specific) then "--keep " + select_first([sex_specific]) + ".fam" else ""} \
        --keep-allele-order \
        --vcf-idspace-to ~{space_character} \
        --const-fid \
        --allow-extra-chr 0 \
        --make-bed \
        --out "~{filebase}" \
        1>&2
    ~{if defined (sex_specific) then "rm " + select_first([sex_specific]) + ".fam" else ""}
    echo "~{sep="\n" select_all([vcf_file, sample_tsv_file, remove_samples_file])}" | \
      sed 's/^.*\///' | tr '\n' '\0' | xargs -0 rm
    plink1.9 \
      --bfile "~{filebase}" \
      --keep-allele-order \
      --indep 50 5 2 \
      --out "~{filebase}" \
      1>&2
    cat "~{filebase}.prune.in" | wc -l
    plink1.9 \
      --bfile "~{filebase}" \
      --keep-allele-order \
      --extract "~{filebase}.prune.in" \
      --make-bed \
      --out "~{filebase}.prune" \
      1>&2
    ~{if defined(sample_tsv_file) then "rm \"" + filebase + ".sex\"" else ""}
    rm "~{filebase}.bed" "~{filebase}.bim" "~{filebase}.fam" "~{filebase}.prune.in"
  >>>

  output {
    Int n_markers = read_int(stdout())
    File bed_file = filebase + ".prune.bed"
    File bim_file = filebase + ".prune.bim"
    File fam_file = filebase + ".prune.fam"
  }

  runtime {
    memory: memory + " GiB"
    disks: "local-disk " + disk_size + " HDD"
    cpu: cpu
    docker: docker
    preemptible: preemptible
    maxRetries: maxRetries
  }
}

task pgt_concat {
  input {
    Array[File] bed_files
    Array[File] bim_files
    Array[File] fam_files
    String filebase

    String docker
    Int cpu = 1
    Int? disk_size_override
    Float memory = 3.5
    Int preemptible = 0
    Int maxRetries = 0
  }

  Float bed_size = size(bed_files, "GiB")
  Float bim_size = size(bim_files, "GiB")
  Float fam_size = size(fam_files, "GiB")
  Int disk_size = select_first([disk_size_override, ceil(10.0 + 2.0 * (bed_size + bim_size + fam_size))])

  command <<<
    set -euo pipefail
    bed_files=~{write_lines(bed_files)}
    bim_files=~{write_lines(bim_files)}
    fam_files=~{write_lines(fam_files)}
    cat $bed_files $bim_files $fam_files | tr '\n' '\0' | xargs -0 mv -t .
    sed -i 's/^.*\///' $bed_files $bim_files $fam_files
    (echo -en "\x6C\x1B\x01"; cat $bed_files | tr '\n' '\0' | xargs -0 tail -qc+4) > "~{filebase}.bed"
    cat $bim_files | tr '\n' '\0' | xargs -0 cat > "~{filebase}.bim"
    head -n1 $fam_files | tr '\n' '\0' | xargs -0 cat > "~{filebase}.fam"
    cat "~{filebase}.bim" | wc -l
    cat $bed_files $bim_files $fam_files | tr '\n' '\0' | xargs -0 rm
  >>>

  output {
    Int n_markers = read_int(stdout())
    File bed_file = filebase + ".bed"
    File bim_file = filebase + ".bim"
    File fam_file = filebase + ".fam"
  }

  runtime {
    memory: memory + " GiB"
    disks: "local-disk " + disk_size + " HDD"
    cpu: cpu
    docker: docker
    preemptible: preemptible
    maxRetries: maxRetries
  }
}

task plink_pca {
  input {
    Int n_smpls
    Int n_markers
    Int? pca_ndim
    File bed_file
    File bim_file
    File fam_file
    String exclusion_regions
    String filebase

    String docker
    Int? cpu_override
    Int? disk_size_override
    Float? mult_override
    Float? memory_override
    Int preemptible = 0
    Int maxRetries = 0

    Float mult = 10 # TODO not sure here what I should request
  }

  Float bed_size = size(bed_file, "GiB")
  Float bim_size = size(bim_file, "GiB")
  Float fam_size = size(fam_file, "GiB")
  Int disk_size = select_first([disk_size_override, ceil(10.0 + 2.0 * (bed_size + bim_size + fam_size))])
  Float memory = select_first([memory_override, 3.5 + mult * n_markers * n_smpls / 1024 / 1024 / 1024])
  Int cpu = select_first([cpu_override, if memory > 6.5 then 2 * ceil(memory / 13) else 1])

  command <<<
    echo "~{sep="\n" select_all([bed_file, bim_file, fam_file])}" | \
      tr '\n' '\0' | xargs -0 mv -t .
    echo "~{exclusion_regions}" | \
      tr ',[:\-]' '\n ' | awk '{print $0,"r"NR}' > exclusion_regions.txt
    plink2 \
      --bed "~{basename(bed_file)}" \
      --bim "~{basename(bim_file)}" \
      --fam "~{basename(fam_file)}" \
      --exclude range exclusion_regions.txt \
      --pca ~{pca_ndim} approx \
      --out "~{filebase}"
    rm exclusion_regions.txt
    echo "~{sep="\n" select_all([bed_file, bim_file, fam_file])}" | \
      sed 's/^.*\///' | tr '\n' '\0' | xargs -0 rm
    set -euo pipefail
  >>>

  output {
    File eigenvec_file = filebase + ".eigenvec"
    File eigenval_file = filebase + ".eigenval"
  }

  runtime {
    memory: memory + " GiB"
    disks: "local-disk " + disk_size + " HDD"
    cpu: cpu
    docker: docker
    preemptible: preemptible
    maxRetries: maxRetries
  }
}

# see https://github.com/rgcgithub/regenie/wiki/Further-parallelization-for-level-0-models-in-Step-1
task regenie_step0 {
  input {
    Int idx
    Int n_phenos
    Int n_smpls
    Int n_markers
    File bed_file
    File bim_file
    File fam_file
    File? covar_file
    File pheno_file
    Boolean binary
    Int bsize
    Boolean loocv
    String? regenie_step0_extra_args
    String filebase

    String docker
    Int? cpu_override
    Int? disk_size_override
    Float? memory_override
    Int preemptible = 1
    Int maxRetries = 0

    Float mult = 6.0 # TODO this is likely incorrect ... I probably need the number of phenotypes tested
  }

  Float bed_size = size(bed_file, "GiB")
  Float bim_size = size(bim_file, "GiB")
  Float fam_size = size(fam_file, "GiB")
  Int disk_size = select_first([disk_size_override, ceil(10.0 + bed_size + bim_size + fam_size + mult * n_markers * n_smpls / 1024 / 1024 / 1024)])
  Float memory = select_first([memory_override, 3.5 + 16.0 * n_smpls * bsize / 1024 / 1024 / 1024])
  Int cpu = select_first([cpu_override, if memory > 6.5 then 2 * ceil(memory / 13) else 1])
  Int n_bins = ceil(1.0 * n_markers / bsize)

  command <<<
    set -euo pipefail
    echo "~{sep="\n" select_all([bed_file, bim_file, fam_file, covar_file, pheno_file])}" | \
      tr '\n' '\0' | xargs -0 mv -t .
    echo -e "~{n_markers} ~{bsize}\n~{filebase}_job~{idx+1} ~{n_bins} ~{n_markers}" > ~{filebase}_~{idx+1}.master
    cut -f2 "~{basename(bim_file)}" > "~{filebase}_job~{idx+1}.snplist"
    regenie \
      --step 1 \
      --bed "~{basename(bed_file, ".bed")}" \
      --phenoFile "~{basename(pheno_file)}" \
      ~{if defined(covar_file) then "--covarFile \"" + basename(select_first([covar_file])) + "\"" else ""} \
      --gz \
      ~{if binary then "--bt" else ""} \
      --bsize ~{bsize} \
      ~{if loocv then "--loocv" else ""} \
      --run-l0 "~{filebase}_~{idx+1}.master",1 \
      ~{if defined(regenie_step0_extra_args) then regenie_step0_extra_args else ""} \
      --threads ~{cpu} \
      --out "~{filebase}" \
      1>&2
    mkdir l0
    seq 1 ~{n_phenos} | sed 's/^/~{filebase}_job~{idx+1}_l0_Y/' | tr '\n' '\0' | xargs -0 mv -t l0
    seq 1 ~{n_phenos} | sed 's/^/l0\/~{filebase}_job~{idx+1}_l0_Y/'
    echo "~{sep="\n" select_all([bed_file, bim_file, fam_file, covar_file, pheno_file])}" | \
      sed 's/^.*\///' | tr '\n' '\0' | xargs -0 rm
  >>>

  output {
    File master_file = filebase + "_" + (idx + 1) + ".master"
    Directory l0_dir = "l0"
    Array[File] l0_files = read_lines(stdout())
  }

  runtime {
    memory: memory + " GiB"
    disks: "local-disk " + disk_size + " HDD"
    cpu: cpu
    docker: docker
    preemptible: preemptible
    maxRetries: maxRetries
  }
}

# see https://github.com/rgcgithub/regenie/wiki/Further-parallelization-for-level-0-models-in-Step-1
task regenie_step1 {
  input {
    String pheno_name
    Int n_smpls
    Int n_markers
    File bed_file
    File bim_file
    File fam_file
    File? covar_file
    File pheno_file
    Array[Int] n_markers_array
    Array[File] l0_files
    Boolean binary
    Int bsize
    Boolean loocv
    String? regenie_step1_extra_args
    String filebase

    String docker
    Int? cpu_override
    Int? disk_size_override
    Float? memory_override
    Int preemptible = 0
    Int maxRetries = 0

    Float mult = 6.0 # TODO this is likely incorrect ... I probably need the number of phenotypes tested
  }

  Float bed_size = size(bed_file, "GiB")
  Float bim_size = size(bim_file, "GiB")
  Float fam_size = size(fam_file, "GiB")
  Float l0_size = size(l0_files, "GiB")
  Int disk_size = select_first([disk_size_override, ceil(10.0 + bed_size + bim_size + fam_size + l0_size + mult * n_markers * n_smpls / 1024 / 1024 / 1024)])
  Float memory = select_first([memory_override, 3.5 + bed_size + bim_size + fam_size])
  Int cpu = select_first([cpu_override, 2 * ceil(memory / 13)]) # always require at least two CPUs

  command <<<
    set -euo pipefail
    markers_lines=~{write_lines(n_markers_array)}
    l0_files=~{write_lines(l0_files)}
    echo "~{sep="\n" select_all([bed_file, bim_file, fam_file, covar_file, pheno_file])}" | \
      cat - $l0_files | tr '\n' '\0' | xargs -0 mv -t .
    cat $l0_files | sed 's/^.*\///;s/.*/& &/;s/[0-9]*$/1/' | while read src dst; do mv --no-clobber $src $dst; done
    sed -i 's/^.*\///;s/[0-9]*$/1/' $l0_files
    awk -F"\t" -v OFS="\t" 'NR==1 {for (i=1; i<=NF; i++) f[$i] = i}
      {print $(f["FID"]),$(f["IID"]),$(f["~{pheno_name}"])}' "~{basename(pheno_file)}" > "~{basename(pheno_file)}.~{pheno_name}"
    awk 'BEGIN {print "~{n_markers} ~{bsize}"} {print "~{filebase}_job"NR,1+int(($0-1)/~{bsize}),$0}' \
      $markers_lines > "~{filebase}.master"
    regenie \
      --step 1 \
      --bed "~{basename(bed_file, ".bed")}" \
      --phenoFile "~{basename(pheno_file)}.~{pheno_name}" \
      ~{if defined(covar_file) then "--covarFile \"" + basename(select_first([covar_file])) + "\"" else ""} \
      ~{if binary then "--bt" else ""} \
      --bsize ~{bsize} \
      ~{if loocv then "--loocv" else ""} \
      --run-l1 "~{filebase}.master" \
      --keep-l0 \
      --gz \
      ~{if binary then "--write-null-firth" else ""} \
      ~{if defined(regenie_step1_extra_args) then regenie_step1_extra_args else ""} \
      --threads ~{cpu} \
      --out "~{filebase}" \
      1>&2
    rm "~{basename(pheno_file)}.~{pheno_name}" "~{filebase}.master"
    mkdir loco
    if [ -f "~{filebase}_1.loco.gz" ]; then
      mv "~{filebase}_1.loco.gz" "loco/~{filebase}_~{pheno_name}.loco.gz"
      sed -i 's/\([^ \t]*\) .*\//\1 /;s/_1.loco.gz/_~{pheno_name}.loco.gz/' "~{filebase}_pred.list"
    else
      touch "loco/~{filebase}_~{pheno_name}.loco.gz"
    fi
    ~{if binary then
      "mkdir firth\n" +
      "if [ -f \"" + filebase + "_1.firth.gz\" ]; then\n" +
      "  mv \"" +  filebase + "_1.firth.gz\" \"firth/" + filebase + "_" + pheno_name + ".firth.gz\"\n" +
      "  sed -i 's/\\([^ \\t]*\\) .*\\//\\1 /;s/_1.firth.gz/_" + pheno_name + ".firth.gz/' \"" + filebase + "_firth.list\"\n" +
      "else\n" +
      "  touch \"firth/" + filebase + "_" + pheno_name + ".firth.gz\"\n" +
      "fi" else ""}
    echo "~{sep="\n" select_all([bed_file, bim_file, fam_file, covar_file, pheno_file])}" | \
      sed 's/^.*\///' | cat - $l0_files | tr '\n' '\0' | xargs -0 rm
  >>>

  output {
    String loco_line = read_string(filebase + "_pred.list")
    File loco_file = "loco/" + filebase + "_" + pheno_name + ".loco.gz"
    String? firth_line = if binary then read_string(filebase + "_firth.list") else None
    File? firth_file = if binary then "firth/" + filebase + "_" + pheno_name + ".firth.gz" else None
  }

  runtime {
    memory: memory + " GiB"
    disks: "local-disk " + disk_size + " HDD"
    cpu: cpu
    docker: docker
    preemptible: preemptible
    maxRetries: maxRetries
  }
}

task serialize_lines {
  input {
    Array[String] lines
    String filename

    String docker
    Int cpu = 1
    Int disk_size = 10
    Float memory = 3.5
    Int preemptible = 1
    Int maxRetries = 0
  }

  command <<<
    set -euo pipefail
    lines=~{write_lines(lines)}
    mv $lines "~{filename}"
  >>>

  output {
    File file = filename
  }

  runtime {
    docker: docker
    cpu: cpu
    disks: "local-disk " + disk_size + " HDD"
    memory: memory + " GiB"
    preemptible: preemptible
    maxRetries: maxRetries
  }
}

# number of variants and sample should be input here
# option psam-cols=fid,sex required due to https://github.com/rgcgithub/regenie/issues/105
task vcf2pgen {
  input {
    File vcf_file
    File? sample_tsv_file
    File? remove_samples_file
    String dosage_field
    String? ref_name
    String? sex_specific
    String space_character

    String docker
    Int cpu = 1
    Int? disk_size_override
    Float memory = 3.5
    Int preemptible = 1
    Int maxRetries = 0
  }

  String filebase = basename(basename(vcf_file, ".bcf"), ".vcf.gz")
  Float vcf_size = size(vcf_file, "GiB")
  Int disk_size = select_first([disk_size_override, ceil(10.0 + 3.0 * vcf_size)])

  command <<<
    set -euo pipefail
    echo "~{sep="\n" select_all([vcf_file, sample_tsv_file, remove_samples_file])}" | \
      tr '\n' '\0' | xargs -0 mv -t .
    ~{if defined(sample_tsv_file) then
        "awk 'NR==1 {for (i=1; i<=NF; i++) f[$i] = i}\n" +
        "  NR>1 {id=$(f[\"sample_id\"]); gsub(\" \",\"" + space_character + "\",id);\n" +
        "  print 0,id,$(f[\"computed_gender\"])}' \"" + basename(select_first([sample_tsv_file])) + "\" > \"" + filebase + ".sex\""
      else ""}
    ~{if defined(remove_samples_file) then
      "sed -i 's/ /" + space_character + "/g;s/^/0\\t/' \"" + basename(select_first([remove_samples_file])) + "\""
      else ""}
    bcftools query -l "~{basename(vcf_file)}" | wc -l
    plink2 \
      --bcf "~{basename(vcf_file)}" dosage=~{dosage_field} \
      ~{if defined(sample_tsv_file) then "--update-sex \"" + filebase + ".sex\" " else ""}\
      ~{if defined(sex_specific) then "--keep-" + select_first([sex_specific]) + "s " else ""}\
      --vcf-idspace-to ~{space_character} \
      --const-fid \
      --allow-extra-chr 0 \
      --make-pgen erase-phase psam-cols=fid,sex \
      ~{if !(ref_name != "GRCh38") then "--split-par hg38 " else if !(ref_name != "GRCh37") then "--split-par hg37 " else ""}\
      --out "~{filebase}" \
      1>&2
    grep -v ^# "~{filebase}.pvar" | wc -l > "~{filebase}.nvar"
    ~{if defined(sample_tsv_file) then "rm \"" + filebase + ".sex\"" else ""}
    echo "~{sep="\n" select_all([vcf_file, sample_tsv_file, remove_samples_file])}" | \
      sed 's/^.*\///' | tr '\n' '\0' | xargs -0 rm
  >>>

  output {
    Int n_smpls = read_int(stdout())
    Int n_markers = read_int(filebase + ".nvar")
    File pgen_file = filebase + ".pgen"
    File pvar_file = filebase + ".pvar"
    File psam_file = filebase + ".psam"
  }

  runtime {
    memory: memory + " GiB"
    disks: "local-disk " + disk_size + " HDD"
    cpu: cpu
    docker: docker
    preemptible: preemptible
    maxRetries: maxRetries
  }
}

# loco file needs to be pruned due to bug https://github.com/rgcgithub/regenie/issues/199
# write_lines() hack needed due to bug https://github.com/broadinstitute/cromwell/issues/5540
task regenie_step2 {
  input {
    Int n_smpls
    Int n_markers
    File pgen_file
    File pvar_file
    File psam_file
    File? covar_file
    File pheno_file
    Array[String] suffix # suffix array passed due to bug https://github.com/broadinstitute/cromwell/issues/5549
    Boolean binary
    File loco_lst
    Array[File] loco_files
    File? firth_lst
    Array[File]? firth_files
    Int bsize
    String? regenie_step2_extra_args

    String docker
    Int? cpu_override
    Int? disk_size_override
    Float memory = 3.5
    Int preemptible = 1
    Int maxRetries = 0
  }

  String filebase = basename(pgen_file, ".pgen")
  Float pgen_size = size(pgen_file, "GiB")
  Float pvar_size = size(pvar_file, "GiB")
  Float psam_size = size(psam_file, "GiB")
  Int disk_size = select_first([disk_size_override, ceil(10.0 + pgen_size + pvar_size + psam_size)])
  Int cpu = select_first([cpu_override, 2 * ceil(memory / 13)]) # always require at least two CPUs

  command <<<
    set -euo pipefail
    loco_files=~{write_lines(loco_files)}
    firth_files=~{if defined(firth_files) then write_lines(select_first([firth_files])) else ""}
    echo "~{sep="\n" select_all([pgen_file, pvar_file, psam_file, covar_file, pheno_file, loco_lst, firth_lst])}" | \
      cat - $loco_files~{if defined(firth_files) then " $firth_files" else ""} | tr '\n' '\0' | xargs -0 mv -t .
    ~{if defined(firth_lst) then
      "awk 'NR==FNR {x[$1]++} NR>FNR && $1 in x' \"" + basename(select_first([firth_lst])) + "\" \"" + basename(loco_lst) + "\" > \"" + basename(loco_lst) + ".alt\""
      else ""}
    regenie \
      --step 2 \
      --pgen "~{basename(pgen_file, ".pgen")}" \
      --phenoFile "~{basename(pheno_file)}" \ \
      ~{if defined(covar_file) then "--covarFile \"" + basename(select_first([covar_file])) + "\"" else ""} \
      --pred "~{if defined(firth_lst) then basename(loco_lst) + ".alt" else basename(loco_lst)}" \
      ~{if binary then "--bt" else ""} \
      --bsize ~{bsize} \
      --firth --approx --firth-se \
      ~{if binary && defined(firth_lst) then "--use-null-firth \"" +  basename(select_first([firth_lst])) + "\"" else ""} \
      ~{if defined(regenie_step2_extra_args) then regenie_step2_extra_args else ""} \
      --gz \
      --threads ~{cpu} \
      --out "~{filebase}" \
      1>&2
    ~{if defined(firth_lst) then "rm \"" + basename(loco_lst) + ".alt\"" else ""}
    echo "~{sep="\n" select_all([pgen_file, pvar_file, psam_file, covar_file, pheno_file, loco_lst, firth_lst])}" | \
      cat - $loco_files~{if defined(firth_files) then " $firth_files" else ""} | sed 's/^.*\///' | tr '\n' '\0' | xargs -0 rm
  >>>

  output {
    Array[File] regenie_files = prefix(filebase + "_", suffix) # suffix array passed due to bug https://github.com/broadinstitute/cromwell/issues/5549
  }

  runtime {
    memory: memory + " GiB"
    disks: "local-disk " + disk_size + " HDD"
    cpu: cpu
    docker: docker
    preemptible: preemptible
    maxRetries: maxRetries
  }
}

task assoc_concat {
  input {
    Array[File]+ assoc_files
    Boolean zst = false
    String filebase

    String docker
    Int cpu = 1
    Int? disk_size_override
    Float memory = 3.5
    Int preemptible = 1
    Int maxRetries = 0
  }

  Float assoc_size = size(assoc_files, "GiB")
  Int disk_size = select_first([disk_size_override, ceil(10.0 + 2.0 * assoc_size)])

  command <<<
    set -euo pipefail
    assoc_files=~{write_lines(assoc_files)}
    cat $assoc_files | tr '\n' '\0' | xargs -0 mv -t .
    sed -i 's/^.*\///' $assoc_files
    cat $assoc_files | tr '\n' '\0' | xargs -0 ~{if zst then "-n1 plink2 --zst-decompress" else "zcat"} | \
      awk 'NR==1 || $0!~"^CHROM" && $0!~"^#CHROM"' | \
      sed 's/^23/X/;s/^PAR[12]/X/' | tr ' ' '\t' | \
      bgzip > "~{filebase}"
    tabix --begin 2 --end 2 --force --skip-lines 1 "~{filebase}"
    cat $assoc_files | tr '\n' '\0' | xargs -0 rm
  >>>

  output {
    File file = filebase
    File index = filebase + ".tbi"
  }

  runtime {
    docker: docker
    cpu: cpu
    disks: "local-disk " + disk_size + " HDD"
    memory: memory + " GiB"
    preemptible: preemptible
    maxRetries: maxRetries
  }
}

# use of !(a!=b) due to bug Cromwell team will not fix: https://github.com/broadinstitute/cromwell/issues/5602
# plink2 will return error code 7 when covariate-only Firth regression fails to converge
# plink2 will return error code 13 when no diploid variants remain for --glm hetonly
task plink_glm {
  input {
    String chr_num
    String pheno_name
    Int n_smpls
    Int n_markers
    File pgen_file
    File pvar_file
    File psam_file
    File? loco_file
    File? covar_file
    File pheno_file
    Boolean binary
    String? plink_extra_args

    String docker
    Int cpu = 1
    Int? disk_size_override
    Float memory = 3.5
    Int preemptible = 1
    Int maxRetries = 0
  }

  String filebase = basename(pgen_file, ".pgen")
  Float pgen_size = size(pgen_file, "GiB")
  Float pvar_size = size(pvar_file, "GiB")
  Float psam_size = size(psam_file, "GiB")
  Int disk_size = select_first([disk_size_override, ceil(10.0 + pgen_size + pvar_size + psam_size)])

  command <<<
    set -euo pipefail
    echo "~{sep="\n" select_all([pgen_file, pvar_file, psam_file, covar_file, pheno_file, loco_file])}" | \
      tr '\n' '\0' | xargs -0 mv -t .
    ~{if defined(loco_file) then
      "zcat \"" + basename(select_first([loco_file])) + "\" | \\\n" +
      "awk 'BEGIN {print \"IID\\tLOCO\"} NR==1 {for (i=2; i<=NF; i++) {sub(\"^0_\", \"\", $i); f[i] = $i}}\n" +
      "  $1==\"" + chr_num + "\" {for (i=2; i<=NF; i++) print f[i]\"\\t\"$i}' " +
      (if defined(covar_file) then "| \\\n" +
        "awk -F\"\\t\" -v OFS=\"\\t\" 'NR==FNR {x[$1]=$2} NR>FNR && $2 in x {print $0\"\\t\"x[$2]}' - \"" + basename(select_first([covar_file])) + "\" "
        else "") +
      "> \"" +  filebase + ".cov\""
    else ""}
    plink2 \
      --pgen "~{basename(pgen_file)}" \
      --pvar "~{basename(pvar_file)}" \
      --psam "~{basename(psam_file)}" \
      ~{if defined(loco_file) then "--covar \"" + filebase + ".cov\"" else
        if defined(covar_file) then "--covar \"" + basename(select_first([covar_file])) + "\"" else ""} \
      ~{if defined(loco_file) || defined(covar_file) then "--covar-variance-standardize" else ""} \
      ~{if defined(loco_file) || defined(covar_file) then "--require-covar" else ""} \
      --1 --pheno "~{basename(pheno_file)}" \
      --pheno-name ~{pheno_name} \
      --require-pheno \
      --glm zs log10 hetonly hide-covar~{if !defined(loco_file) && !defined(covar_file) then " allow-no-covars" else ""} cc-residualize firth-residualize cols=+a1freq,+machr2 \
      ~{if defined(plink_extra_args) then plink_extra_args else ""} \
      --out "~{filebase}" \
      1>&2 || if [[ $? -eq 7 || $? -eq 13 ]]; then
        echo -en "\x28\xb5\x2f\xfd\x24\x00\x01\x00\x00\x99\xe9\xd8\x51" > "~{filebase}.~{pheno_name}.glm.~{if binary then "logistic.hybrid" else "linear"}.zst"
      else exit $?; fi
    ~{if defined(loco_file) then "rm \"" + filebase + ".cov\"" else ""}
    echo "~{sep="\n" select_all([pgen_file, pvar_file, psam_file, covar_file, pheno_file, loco_file])}" | \
      sed 's/^.*\///' | tr '\n' '\0' | xargs -0 rm
  >>>

  output {
    File assoc_file = filebase + "." + pheno_name + ".glm." + (if binary then "logistic.hybrid" else "linear") + ".zst"
  }

  runtime {
    memory: memory + " GiB"
    disks: "local-disk " + disk_size + " HDD"
    cpu: cpu
    docker: docker
    preemptible: preemptible
    maxRetries: maxRetries
  }
}
