version development

## Copyright (c) 2020 Giulio Genovese
##
## Version 2020-07-22
##
## Contact Giulio Genovese <giulio.genovese@gmail.com>
##
## This WDL workflow runs MoChA on a cohort of samples genotyped with either Illumina or Affymetrix DNA microarrays
##
## Cromwell version support
## - Successfully tested on v52
##
## Distributed under terms of the MIT License

struct Reference {
  String name
  String fasta_ref
  Int n_chrs
  String mhc_reg
  String dup_file
  String genetic_map_file
  String cnp_file
  String cyto_file
  String kgp_pfx
  String kgp_sfx
  Int n_smpls
}

workflow mocha {
  input {
    String sample_set_id
    String mode # idat gtc cel chp txt vcf pvcf
    Boolean realign = false
    Boolean gtc_output = false # only for idat
    Boolean chp_output = false # only for cel
    Boolean do_not_phase = false # only for idat gtc cel chp txt mode
    Boolean do_not_mocha = false # only for idat gtc cel chp txt vcf mode
    Boolean sequencing = false
    Int idat_batch_size = 48
    Int gtc_batch_size = 1024
    Int chp_batch_size = 1024
    Int phase_threads = 4
    Float max_win_size_cm = 300.0
    Float overlap_size_cm = 5.0
    String ref_path = ""
    String manifest_path = ""
    File docker_json_file # ubuntu pandas iaap_cli autoconvert apt gtc2vcf bcftools shapeit4 eagle mocha mocha_plot
    File ref_json_file # name fasta_ref n_chrs mhc_reg dup_file genetic_map_file cnp_file cyto_file kgp_pfx kgp_sfx n_smpls
    File sample_tsv_file # sample_id batch_id green_idat red_idat gtc cel chp computed_gender call_rate
    File batch_tsv_file # batch_id path bpm csv egt sam xml zip snp report calls confidences summary vcf vcf_index xcl_vcf xcl_vcf_index
    File? ped_file
    File? extra_xcl_vcf_file
    Boolean autoconvert = false
    Boolean? table_output
    Boolean? do_not_check_bpm
    Boolean? do_not_use_reference
    Boolean eagle = false
    String delim = "~"
    Array[String]? chip_type
    Array[String] tags = ["GT", "BAF", "LRR"]
    Int? gc_window_size
  }

  Boolean mode_is_vcf = mode == "vcf" || mode == "pvcf"

  # read table with genome reference resources and compute scatter intervals
  Map[String, String] docker = read_json(docker_json_file)

  # read table with genome reference resources and compute scatter intervals
  Reference ref = read_json(ref_json_file)

  # read table with batchs information (scatter could be avoided if there was a tail() function)
  call tsv_sorted as batch_sorted_tsv { input: tsv_file = batch_tsv_file, column = "batch_id", docker = docker["ubuntu"] }
  Array[Array[String]] batch_tsv = read_tsv(batch_sorted_tsv.file)
  Int n_batches = length(batch_tsv)-1
  scatter (idx in range(n_batches)) { Array[String] batch_tsv_rows = batch_tsv[(idx+1)] }
  Map[String, Array[String]] batch_tbl = as_map(zip(batch_tsv[0], transpose(batch_tsv_rows)))
  Array[String] batches = batch_tbl["batch_id"]

  # compute data paths for each batch, if available (scatter could be avoided if there was a contains_key() function)
  scatter (key in keys(batch_tbl)) { Boolean? is_key_equal_path = if key == "path" then true else None }
  scatter (idx in range(n_batches)) { String data_paths = if length(select_all(is_key_equal_path))>0 then batch_tbl["path"][idx] else "" }

  # aligns manifest file to human genome reference if requested
  if (realign && !mode_is_vcf) {
    # hack due to lack of unique function
    Array[String] csv_files = keys(collect_by_key(zip(batch_tbl["csv"], batches)))
    scatter (csv_file in csv_files) {
      call csv2bam {
        input:
          plugin = if mode == "idat" || mode == "gtc" then "gtc2vcf" else "affy2vcf",
          csv_file = manifest_path + csv_file,
          fasta_ref = ref_path + ref["fasta_ref"],
          fasta_ref_idxs = prefix(ref_path + ref["fasta_ref"] + ".", ["amb", "ann", "bwt", "pac", "sa"]),
          docker = docker["gtc2vcf"],
      }
    }
    Map[String, File] csv2sam = as_map(zip(csv_files, csv2bam.bam_file))
  }
  # compute sam file for each batch, if available (scatter could be avoided if there was a contains_key() function)
  scatter (key in keys(batch_tbl)) { Boolean? is_key_equal_sam = if key == "sam" then true else None }
  scatter (idx in range(n_batches)) {
    File? sams = if realign && !mode_is_vcf then select_first([csv2sam])[(batch_tbl["csv"][idx])]
      else if length(select_all(is_key_equal_sam))>0 then manifest_path + batch_tbl["sam"][idx] else None
  }

  if (!mode_is_vcf) {
    # resort table with sample information and extract sample_id column
    call tsv_sorted as sample_sorted_tsv { input: tsv_file = sample_tsv_file, column = "batch_id", docker = docker["ubuntu"] }
    call tsv_column as sample_id_lines { input: tsv_file = sample_sorted_tsv.file, column = "sample_id", docker = docker["ubuntu"] }
    call tsv_column as batch_id_lines { input: tsv_file = sample_sorted_tsv.file, column = "batch_id", docker = docker["ubuntu"] }

    # process Illumina data
    if (mode == "idat") {
      call tsv_column as green_idat_lines { input: tsv_file = sample_sorted_tsv.file, column = "green_idat", docker = docker["ubuntu"] }
      call tsv_column as red_idat_lines { input: tsv_file = sample_sorted_tsv.file, column = "red_idat", docker = docker["ubuntu"] }
      # group samples by IDAT batches
      call batch_scatter as idat_scatter { input: batch_id_file = batch_id_lines.file, sub_batch_size = idat_batch_size, delim = delim, docker = docker["ubuntu"] }
      Map[String, Array[String]] idat_batch2green_idat = collect_by_key(zip(read_lines(idat_scatter.sub_batch_id), read_lines(green_idat_lines.file)))
      Map[String, Array[String]] idat_batch2red_idat = collect_by_key(zip(read_lines(idat_scatter.sub_batch_id), read_lines(red_idat_lines.file)))
      Map[String, Int] idat_sub_batch2idx = as_map(zip(idat_scatter.sub_batches, idat_scatter.idxs))
      scatter (idat_batch in idat_scatter.sub_batches) {
        Int idat_idx = idat_sub_batch2idx[idat_batch]
        call idat2gtc {
          input:
            bpm_file = manifest_path + batch_tbl["bpm"][idat_idx],
            egt_file = manifest_path + batch_tbl["egt"][idat_idx],
            green_idat_files = prefix(data_paths[idat_idx], idat_batch2green_idat[idat_batch]),
            red_idat_files = prefix(data_paths[idat_idx], idat_batch2red_idat[idat_batch]),
            autoconvert = autoconvert,
            filebase = sample_set_id + "." + idat_batch,
            docker = docker[if autoconvert then "autoconvert" else "iaap_cli"],
        }
      }
      call tsv_concat as idat_tsv { input: tsv_files = idat2gtc.idat_tsv, filebase = sample_set_id + ".idat", docker = docker["ubuntu"] }
    }

    if (mode == "gtc") {
      call tsv_column as gtc_lines { input: tsv_file = sample_sorted_tsv.file, column = "gtc", docker = docker["ubuntu"] }
    }

    if (mode == "idat" || mode == "gtc") {
      # group samples by GTC batches
      call batch_scatter as gtc_scatter { input: batch_id_file = batch_id_lines.file, sub_batch_size = gtc_batch_size, delim = delim, docker = docker["ubuntu"] }
      Array[String]+ input_gtc_files = if mode == "idat" then flatten(select_first([idat2gtc.gtc_files])) else read_lines(select_first([gtc_lines.file]))
      Map[String, Array[String]] gtc_batch2gtc_files = collect_by_key(zip(read_lines(gtc_scatter.sub_batch_id), input_gtc_files))
      scatter (file in input_gtc_files) { String gtc_barcode = basename(file, ".gtc") }
      Map[String, Array[String]] gtc_batch2barcode = collect_by_key(zip(read_lines(gtc_scatter.sub_batch_id), gtc_barcode))
      Map[String, Array[String]] gtc_batch2sample_id = collect_by_key(zip(read_lines(gtc_scatter.sub_batch_id), read_lines(sample_id_lines.file)))
      Map[String, Int] gtc_sub_batch2idx = as_map(zip(gtc_scatter.sub_batches, gtc_scatter.idxs))
      scatter (gtc_batch in gtc_scatter.sub_batches) {
        Int gtc_idx = gtc_sub_batch2idx[gtc_batch]
        call gtc2vcf {
          input:
            tags = tags,
            bpm_file = manifest_path + batch_tbl["bpm"][gtc_idx],
            csv_file = manifest_path + batch_tbl["csv"][gtc_idx],
            egt_file = manifest_path + batch_tbl["egt"][gtc_idx],
            fasta_ref = ref_path + ref["fasta_ref"],
            fasta_ref_fai = ref_path + ref["fasta_ref"] + ".fai",
            gc_window_size = gc_window_size,
            gtc_files = if mode == "idat" then gtc_batch2gtc_files[gtc_batch] else prefix(data_paths[gtc_idx], gtc_batch2gtc_files[gtc_batch]),
            do_not_check_bpm = do_not_check_bpm,
            sam_file = sams[gtc_idx],
            reheader_map = as_map(zip(gtc_batch2barcode[gtc_batch], gtc_batch2sample_id[gtc_batch])),
            filebase = sample_set_id + "." + gtc_batch,
            docker = docker["gtc2vcf"],
        }
      }
      call tsv_concat as gtc_tsv { input: tsv_files = gtc2vcf.gtc_tsv, filebase = sample_set_id + ".gtc", docker = docker["ubuntu"] }

      # this job can be long, so it is better to run as non-preemptible
      Map[Int, Array[String]] idx2gtc2vcf_files = collect_by_key(zip(gtc_scatter.idxs, gtc2vcf.vcf_file))
      Map[Int, Array[String]] idx2gtc2vcf_idxs = collect_by_key(zip(gtc_scatter.idxs, gtc2vcf.vcf_idx))
      scatter (idx in range(n_batches)) {
        if (length(idx2gtc2vcf_files[idx]) > 1) {
          call vcf_merge as gtc2vcf_merge {
            input:
              vcf_files = idx2gtc2vcf_files[idx],
              filebase = sample_set_id + "." + batches[idx],
              docker = docker["gtc2vcf"]
          }
        }
        File gtc2vcf_files = select_first([gtc2vcf_merge.vcf_file, idx2gtc2vcf_files[idx][0]])
        File gtc2vcf_idxs = select_first([gtc2vcf_merge.vcf_idx, idx2gtc2vcf_idxs[idx][0]])
      }
    }

    if (mode == "cel" || mode == "txt") {
      call tsv_column as cel_lines { input: tsv_file = sample_sorted_tsv.file, column = "cel", docker = docker["ubuntu"] }
    }

    # process Affymetrix data
    if (mode == "cel") {
      Map[String, Array[String]] batch2cel_files = collect_by_key(zip(read_lines(batch_id_lines.file), read_lines(select_first([cel_lines.file]))))
      scatter (idx in range(n_batches)) {
        call cel2affy as cel2chp {
          input:
            xml_file = manifest_path + batch_tbl["xml"][idx],
            zip_file = manifest_path + batch_tbl["zip"][idx],
            cel_files = prefix(data_paths[idx], batch2cel_files[(batches[idx])]),
            chip_type = chip_type,
            table_output = table_output,
            filebase = sample_set_id + "." + batches[idx],
            docker = docker["apt"]
        }
      }
      call tsv_concat as cel_tsv { input: tsv_files = cel2chp.cel_tsv, filebase = sample_set_id + ".cel", docker = docker["ubuntu"] }
    }

    if (mode == "chp") {
      call tsv_column as chp_lines { input: tsv_file = sample_sorted_tsv.file, column = "chp", docker = docker["ubuntu"] }
    }

    if (mode == "cel" || mode == "chp") {
      # group samples by CHP batches
      call batch_scatter as chp_scatter { input: batch_id_file = batch_id_lines.file, sub_batch_size = chp_batch_size, delim = delim, docker = docker["ubuntu"] }
      Array[String]+ input_chp_files = if mode == "cel" then flatten(select_first([cel2chp.chp_files])) else read_lines(select_first([chp_lines.file]))
      Map[String, Array[String]] chp_batch2chp_files = collect_by_key(zip(read_lines(chp_scatter.sub_batch_id), input_chp_files))
      scatter (file in input_chp_files) { String chp_barcode = basename(basename(file, ".AxiomGT1.chp"), ".birdseed-v2.chp") }
      Map[String, Array[String]] chp_batch2barcode = collect_by_key(zip(read_lines(chp_scatter.sub_batch_id), chp_barcode))
      Map[String, Array[String]] chp_batch2sample_id = collect_by_key(zip(read_lines(chp_scatter.sub_batch_id), read_lines(sample_id_lines.file)))
      Map[String, Int] chp_sub_batch2idx = as_map(zip(chp_scatter.sub_batches, chp_scatter.idxs))
      scatter (chp_batch in chp_scatter.sub_batches) {
        Int chp_idx = chp_sub_batch2idx[chp_batch]
        call affy2vcf as chp2vcf {
          input:
            tags = tags,
            csv_file = manifest_path + batch_tbl["csv"][chp_idx],
            fasta_ref = ref_path + ref["fasta_ref"],
            fasta_ref_fai = ref_path + ref["fasta_ref"] + ".fai",
            gc_window_size = gc_window_size,
            snp_file = if mode == "cel" then select_first([cel2chp.snp_file])[chp_idx] else data_paths[chp_idx] + batch_tbl["snp"][chp_idx],
            chp_files = if mode == "cel" then chp_batch2chp_files[chp_batch] else prefix(data_paths[chp_idx], chp_batch2chp_files[chp_batch]),
            sam_file = sams[chp_idx],
            reheader_map = as_map(zip(chp_batch2barcode[chp_batch], chp_batch2sample_id[chp_batch])),
            filebase = sample_set_id + "." + chp_batch,
            docker = docker["gtc2vcf"]
        }
      }
      call tsv_concat as chp_tsv { input: tsv_files = select_all(chp2vcf.chp_tsv), filebase = sample_set_id + ".chp", docker = docker["ubuntu"] }

      # this job can be long, so it is better to run as non-preemptible
      Map[Int, Array[String]] idx2chp2vcf_files = collect_by_key(zip(chp_scatter.idxs, chp2vcf.vcf_file))
      Map[Int, Array[String]] idx2chp2vcf_idxs = collect_by_key(zip(chp_scatter.idxs, chp2vcf.vcf_idx))
      scatter (idx in range(n_batches)) {
        if (length(idx2chp2vcf_files[idx]) > 1) {
          call vcf_merge as chp2vcf_merge {
            input:
              vcf_files = idx2chp2vcf_files[idx],
              filebase = sample_set_id + "." + batches[idx],
              docker = docker["gtc2vcf"]
          }
        }
        File chp2vcf_files = select_first([chp2vcf_merge.vcf_file, idx2chp2vcf_files[idx][0]])
        File chp2vcf_idxs = select_first([chp2vcf_merge.vcf_idx, idx2chp2vcf_idxs[idx][0]])
      }
    }

    if (mode == "txt") {
      scatter (key in keys(batch_tbl)) { Boolean? is_key_equal_confidences = if key == "confidences" then true else None }
      scatter (file in read_lines(select_first([cel_lines.file]))) { String cel_barcode = basename(file, ".CEL") }
      Map[String, Array[String]] batch2barcode = collect_by_key(zip(read_lines(batch_id_lines.file), cel_barcode))
      Map[String, Array[String]] batch2sample_id = collect_by_key(zip(read_lines(batch_id_lines.file), read_lines(sample_id_lines.file)))
      scatter (idx in range(n_batches)) {
        # this job can be long, so it is better to run as non-preemptible
        call affy2vcf as txt2vcf {
          input:
            tags = tags,
            csv_file = manifest_path + batch_tbl["csv"][idx],
            fasta_ref = ref_path + ref["fasta_ref"],
            fasta_ref_fai = ref_path + ref["fasta_ref"] + ".fai",
            gc_window_size = gc_window_size,
            calls_file = data_paths[idx] + batch_tbl["calls"][idx],
            confidences_file = if length(select_all(is_key_equal_confidences))>0 then data_paths[idx] + batch_tbl["confidences"][idx] else None,
            summary_file = data_paths[idx] + batch_tbl["summary"][idx],
            snp_file = data_paths[idx] + batch_tbl["snp"][idx],
            sam_file = sams[idx],
            reheader_map = as_map(zip(batch2barcode[(batches[idx])], batch2sample_id[(batches[idx])])),
            filebase = sample_set_id + "." + batches[idx],
            docker = docker["gtc2vcf"]
        }
        File reports = data_paths[idx] + batch_tbl["report"][idx]
      }
      call tsv_concat as txt_tsv { input: tsv_files = reports, filebase = sample_set_id + ".txt", docker = docker["ubuntu"] }
    }

    Array[File] unphased_vcf_files = select_first([gtc2vcf_files, chp2vcf_files, txt2vcf.vcf_file])
    Array[File] unphased_vcf_idxs = select_first([gtc2vcf_idxs, chp2vcf_idxs, txt2vcf.vcf_idx])
    call lst_flatten as flatten_sample_id_lines { input: lst_files = select_first([gtc2vcf.sample_id_lines, chp2vcf.sample_id_lines, txt2vcf.sample_id_lines]), filebase = "sample_id", docker = docker["ubuntu"] }
  }

  if (mode_is_vcf) {
    call tsv_column as vcf_sample_id_lines { input: tsv_file = sample_tsv_file, column = "sample_id", docker = docker["ubuntu"] }
  }

  call tsv_column as computed_gender_lines { input: tsv_file = select_first([gtc_tsv.file, chp_tsv.file, txt_tsv.file, sample_tsv_file]), column = "computed_gender", docker = docker["ubuntu"] }
  if (mode != "pvcf" && !do_not_phase) {
    call ref_scatter {
      input:
        n_chrs = ref["n_chrs"],
        fasta_ref_fai = ref_path + ref["fasta_ref"] + ".fai",
        genetic_map_file = ref_path + ref["genetic_map_file"],
        max_win_size_cm = max_win_size_cm,
        overlap_size_cm = overlap_size_cm,
        docker = docker["pandas"]
    }
    scatter (idx in range(n_batches)) {
      call vcf_scatter {
        input:
          vcf_file = if mode_is_vcf then data_paths[idx] + batch_tbl["vcf"][idx] else select_first([unphased_vcf_files])[idx],
          intervals_bed = ref_scatter.intervals_bed,
          filebase = sample_set_id + "." + batches[idx] + ".",
          docker = docker["bcftools"]
      }
    }
    if (defined(extra_xcl_vcf_file)) {
      call vcf_scatter as xcl_vcf_scatter {
        input:
          vcf_file = select_first([extra_xcl_vcf_file]),
          intervals_bed = ref_scatter.intervals_bed,
          filebase = basename(select_first([extra_xcl_vcf_file]), ".bcf") + ".",
          docker = docker["bcftools"]
      }
    }

    call lst_concat as sample_id_split_tsv { input: lst_files = vcf_scatter.sample_id_lines, filebase = "split_sample_id", docker = docker["ubuntu"] }
    call tsv_column as call_rate_lines { input: tsv_file = select_first([gtc_tsv.file, chp_tsv.file, txt_tsv.file, sample_tsv_file]), column = "call_rate", docker = docker["ubuntu"] }
    Int n_smpls = length(flatten(read_tsv(sample_id_split_tsv.file)))
    Boolean use_reference = !select_first([do_not_use_reference, !eagle && n_smpls > 2 * ref["n_smpls"]])
    Array[Array[File]] interval_slices = transpose(vcf_scatter.vcf_files)
    Array[String] chrs = transpose(read_tsv(ref_scatter.intervals_bed))[0]
    scatter (idx in range(length(chrs))) {
      call vcf_merge {
        input:
          vcf_files = interval_slices[idx],
          filebase = sample_set_id + "." + idx,
          docker = docker["bcftools"]
      }
      call vcf_qc {
        input:
          vcf_file = vcf_merge.vcf_file,
          vcf_idx = vcf_merge.vcf_idx,
          mhc_reg = ref["mhc_reg"],
          dup_file = ref_path + ref["dup_file"],
          sample_id_file = select_first([flatten_sample_id_lines.file, vcf_sample_id_lines.file]),
          computed_gender_file = computed_gender_lines.file,
          call_rate_file = call_rate_lines.file,
          extra_xcl_vcf_file = if defined(extra_xcl_vcf_file) then select_first([xcl_vcf_scatter.vcf_files])[idx] else None,
          docker = docker["bcftools"]
      }
      call vcf_phase {
        input:
          n_smpls = n_smpls,
          n_vars = vcf_qc.n_vars,
          unphased_vcf_file = vcf_merge.vcf_file,
          unphased_vcf_idx = vcf_merge.vcf_idx,
          genetic_map_file = ref_path + ref["genetic_map_file"],
          n_ref_smpls = if use_reference then ref["n_smpls"] else None,
          ref_vcf_file = if use_reference then ref_path + ref["kgp_pfx"] + chrs[idx] + ref["kgp_sfx"] else None,
          ref_vcf_idx = if use_reference then ref_path + ref["kgp_pfx"] + chrs[idx] + ref["kgp_sfx"] + ".csi" else None,
          xcl_vcf_file = vcf_qc.xcl_vcf_file,
          xcl_vcf_idx = vcf_qc.xcl_vcf_idx,
          chr = if eagle then None else chrs[idx],
          sequencing = sequencing,
          eagle = eagle,
          docker = docker[if eagle then "eagle" else "shapeit4"],
          cpu = phase_threads
      }
      call vcf_split {
        input:
          vcf_file = vcf_phase.phased_vcf_file,
          batches = batches,
          sample_id_file = sample_id_split_tsv.file,
          ped_file = ped_file,
          rule = if defined(ped_file) then ref["name"] else None,
          docker = docker["bcftools"]
      }
    }

    call vcf_concat as xcl_vcf_concat {
      input:
        vcf_files = vcf_qc.xcl_vcf_file,
        filebase = sample_set_id + ".xcl",
        docker = docker["bcftools"]
    }

    Array[Array[File]] batch_slices = transpose(vcf_split.vcf_files)
    scatter (idx in range(n_batches)) {
      call vcf_concat {
        input:
          vcf_files = batch_slices[idx],
          filebase = sample_set_id + "." + batches[idx] + ".GT",
          docker = docker["bcftools"]
      }
      call vcf_import {
        input:
          phased_vcf_file = vcf_concat.vcf_file,
          phased_vcf_idx = vcf_concat.vcf_idx,
          unphased_vcf_file = if mode_is_vcf then data_paths[idx] + batch_tbl["vcf"][idx] else select_first([unphased_vcf_files])[idx],
          unphased_vcf_idx = if mode_is_vcf then data_paths[idx] + batch_tbl["vcf_index"][idx] else select_first([unphased_vcf_idxs])[idx],
          docker = docker["bcftools"]
      }
    }
  }

  if (!do_not_phase && !do_not_mocha) {
    scatter (idx in range(n_batches)) {
      call vcf_mocha {
        input:
          filebase = sample_set_id + "." + batches[idx],
          rule = ref["name"],
          vcf_file = if mode == "pvcf" then data_paths[idx] + batch_tbl["vcf"][idx] else select_first([vcf_import.vcf_file])[idx],
          vcf_idx = if mode == "pvcf" then data_paths[idx] + batch_tbl["vcf_index"][idx] else select_first([vcf_import.vcf_idx])[idx],
          sample_id_file = select_first([flatten_sample_id_lines.file, vcf_sample_id_lines.file]),
          computed_gender_file = computed_gender_lines.file,
          xcl_vcf_file = if mode == "pvcf" then data_paths[idx] + batch_tbl["xcl_vcf"][idx] else xcl_vcf_concat.vcf_file,
          xcl_vcf_idx = if mode == "pvcf" then data_paths[idx] + batch_tbl["xcl_vcf_index"][idx] else xcl_vcf_concat.vcf_idx,
          cnp_file = ref_path + ref["cnp_file"],
          docker = docker["mocha"]
      }
    }
    call tsv_concat as mocha_stats_tsv { input: tsv_files = vcf_mocha.stats_tsv, filebase = sample_set_id + ".stats", docker = docker["ubuntu"] }
    call tsv_concat as mocha_calls_tsv { input: tsv_files = vcf_mocha.calls_tsv, filebase = sample_set_id + ".calls", docker = docker["ubuntu"] }
    call mocha_summary { input: calls_tsv = mocha_calls_tsv.file, stats_tsv = mocha_stats_tsv.file, ucsc_beds = vcf_mocha.ucsc_bed, filebase = sample_set_id, docker = docker["mocha_plot"] }
  }

  output {
    File? ref_intervals_bed = ref_scatter.intervals_bed
    File? idat_tsv_file = idat_tsv.file
    File? gtc_tsv_file = gtc_tsv.file
    File? cel_tsv_file = cel_tsv.file
    File? chp_tsv_file = chp_tsv.file
    File? txt_tsv_file = txt_tsv.file
    File sample_id_file = select_first([flatten_sample_id_lines.file, vcf_sample_id_lines.file])
    File computed_gender_file = computed_gender_lines.file
    File? call_rate_file = call_rate_lines.file
    File? mocha_stats_file = mocha_stats_tsv.file
    File? mocha_calls_file = mocha_calls_tsv.file
    File? mocha_summary_pdf = mocha_summary.pdf
    File? mocha_ucsc_bed = mocha_summary.ucsc_bed
    Array[File]? bam_files = csv2bam.bam_file
    Array[File]? mendel_files = if mode != "pvcf" && !do_not_phase && defined(ped_file) then select_all(select_first([vcf_split.mendel_tsv])) else None
    Array[File]? gtc_files = if mode == "idat" && gtc_output then flatten(select_first([idat2gtc.gtc_files])) else None
    Array[File]? chp_files = if mode == "cel" && chp_output then flatten(select_first([cel2chp.chp_files])) else None
    Array[File]? snp_files = if mode == "cel" && chp_output then select_first([cel2chp.snp_file]) else None
    Array[File] vcf_files = select_first([vcf_mocha.mocha_vcf_file, vcf_import.vcf_file, unphased_vcf_files])
    Array[File] vcf_idxs = select_first([vcf_mocha.mocha_vcf_idx, vcf_import.vcf_idx, unphased_vcf_idxs])
    File? xcl_vcf_file = xcl_vcf_concat.vcf_file
    File? xcl_vcf_idx = xcl_vcf_concat.vcf_idx
  }

  meta {
    author: "Giulio Genovese (with help from Chris Whelan)"
    email: "giulio.genovese@gmail.com"
    description: "See the [MoChA](https://github.com/freeseek/mocha) website for more information."
  }
}

# uses https://stackoverflow.com/questions/22464786/ignoring-bash-pipefail-for-error-code-141#comment60412687_33026977
task tsv_sorted {
  input {
    File tsv_file
    String column

    String docker
    Int cpu = 1
    Int disk_size = 10
    Float memory = 3.5
    Int preemptible = 1
    Int maxRetries = 0
  }

  String filebase = basename(tsv_file, ".tsv")

  command <<<
    set -euo pipefail
    col=$((grep -v ^# "~{tsv_file}" || if [[ $? -eq 141 ]]; then true; else exit $?; fi) | \
      head -n1 | tr '\t' '\n' | awk -F"\t" '$0=="~{column}" {print NR}')
    grep -v ^# "~{tsv_file}" | \
      (read -r; printf "%s\n" "$REPLY"; sort -k $col,$col -s -t $'\t') > "~{filebase}.sorted.tsv"
    rm "~{tsv_file}"
  >>>

  output {
    File file = filebase + ".sorted.tsv"
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

task tsv_column {
  input {
    File tsv_file
    String column

    String docker
    Int cpu = 1
    Int disk_size = 10
    Float memory = 3.5
    Int preemptible = 1
    Int maxRetries = 0
  }

  command <<<
    set -euo pipefail
    grep -v ^# "~{tsv_file}" | \
      awk -F"\t" 'NR==1 {for (i=1; i<=NF; i++) if ($i=="~{column}") col=i} NR>1 {print $col}' > "~{column}.lines"
    rm "~{tsv_file}"
  >>>

  output {
    File file = column + ".lines"
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


# uses https://stackoverflow.com/questions/22464786/ignoring-bash-pipefail-for-error-code-141#comment60412687_33026977
# code to avoid https://github.com/broadinstitute/cromwell/issues/5540
task tsv_concat {
  input {
    Array[File]+ tsv_files
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
    tsv_files=~{write_lines(tsv_files)}
    ((grep -v ^# "~{tsv_files[0]}" || if [[ $? -eq 141 ]]; then true; else exit $?; fi) | head -n1
    while IFS= read -r file
    do
      grep -v ^# "$file" | tail -n+2
    done < $tsv_files) > "~{filebase}.tsv"
    cat $tsv_files | tr '\n' '\0' | xargs -0 rm
  >>>

  output {
    File file = filebase + ".tsv"
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

task lst_flatten {
  input {
    Array[File] lst_files
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
    lst_files=~{write_lines(lst_files)}
    cat $lst_files | tr '\n' '\0' | xargs -0 awk 1 > "~{filebase}.lines"
    cat $lst_files | tr '\n' '\0' | xargs -0 rm
  >>>

  output {
    File file = filebase + ".lines"
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

task lst_concat {
  input {
    Array[File] lst_files
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
    lst_files=~{write_lines(lst_files)}
    cat $lst_files | tr '\n' '\0' | xargs -0 -n 1 awk '{printf $0"\t"} END {printf "\n"}' | sed 's/\t$//' > "~{filebase}.tsv"
    cat $lst_files | tr '\n' '\0' | xargs -0 rm
  >>>

  output {
    File file = filebase + ".tsv"
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

# this task generates batches from batchs and then returns them in the same order as batches
task batch_scatter {
  input {
    File batch_id_file
    Int sub_batch_size
    String delim

    String docker
    Int cpu = 1
    Int disk_size = 10
    Float memory = 3.5
    Int preemptible = 1
    Int maxRetries = 0
  }

  command <<<
    set -euo pipefail
    awk -F"\t" 'NR==FNR {x[$0]++} NR>FNR {n=x[$0]/int((x[$0]-1)/~{sub_batch_size}+1); print $0"~{delim}"int(y[$0]/n); y[$0]++}' \
      "~{batch_id_file}" \
      "~{batch_id_file}" \
      > sub_batch_ids.lines
    uniq sub_batch_ids.lines
    uniq sub_batch_ids.lines | awk -F"~{delim}" -v OFS="~{delim}" '{NF--} !x[$0]++ {idx++} {print idx-1}' 1>&2
    rm "~{batch_id_file}"
  >>>

  output {
    File sub_batch_id = "sub_batch_ids.lines"
    Array[String] sub_batches = read_lines(stdout())
    Array[Int] idxs = read_lines(stderr())
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

# once String sep(String, Array[String]) is implemented code can be modified, as the current syntax is non-compliant
task csv2bam {
  input {
    String plugin = "gtc2vcf"
    File csv_file
    File fasta_ref
    Array[File]+ fasta_ref_idxs

    String docker
    Int cpu = 1
    Int? disk_size_override
    Float? memory_override
    Int preemptible = 1
    Int maxRetries = 0
  }

  Float csv_size = size(csv_file, "GiB")
  Float ref_size = size(fasta_ref, "GiB")
  Float index_size = size(fasta_ref_idxs, "GiB")
  Int disk_size = select_first([disk_size_override, ceil(10.0 + 2.0 * csv_size + ref_size + index_size)])
  # if gtc2vcf was memory efficient this requirement could be eased
  Float memory = select_first([memory_override, ceil(7.25 + 2.5 * csv_size)])

  command <<<
    set -euo pipefail
    echo "~{sep="\n" flatten([[csv_file, fasta_ref], fasta_ref_idxs])}" | \
      tr '\n' '\0' | xargs -0 -n 1 ln -s
    bcftools +~{plugin} \
      --csv "~{basename(csv_file)}" \
      --fasta-flank | \
      bwa mem~{if cpu > 1 then " -t " + cpu else ""} -M "~{basename(fasta_ref)}" - | \
      samtools view -bS -o "~{basename(csv_file, ".csv")}.bam"
    echo "~{sep="\n" flatten([[csv_file, fasta_ref], fasta_ref_idxs])}" | \
      awk -F/ '{print $0"\n"$NF}' | tr '\n' '\0' | xargs -0 rm
  >>>

  output {
    File bam_file = basename(csv_file, ".csv") + ".bam"
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

# https://support.terra.bio/hc/en-us/community/posts/360071476431-Terra-fails-to-delocalize-files-listed-through-read-lines-
task idat2gtc {
  input {
    File bpm_file
    File egt_file
    Array[File]+ green_idat_files
    Array[File]+ red_idat_files
    Boolean is_gender_autocall = true # from dsde-pipelines/tasks/IlluminaGenotypingArrayTasks.wdl
    Boolean autoconvert = false
    String filebase

    String docker
    Int cpu = 1
    Int? disk_size_override
    Float? memory_override
    Int preemptible = 1
    Int maxRetries = 0
  }

  Float bpm_size = size(bpm_file, "GiB")
  Float egt_size = size(egt_file, "GiB")
  Float green_idat_size = length(green_idat_files) * size(green_idat_files[0], "GiB")
  Float red_idat_size = length(red_idat_files) * size(red_idat_files[0], "GiB")
  Int disk_size = select_first([disk_size_override, ceil(10.0 + bpm_size + egt_size + 2.0 * green_idat_size + 2.0 * red_idat_size)])
  Float memory = select_first([memory_override, 3.5 + 5.0 * bpm_size + 5.0 * egt_size]) # will request more than 8GB for Omni5

  command <<<
    set -euo pipefail
    green_idat_files=~{write_lines(green_idat_files)}
    red_idat_files=~{write_lines(red_idat_files)}
    echo "~{sep="\n" flatten([[bpm_file, egt_file]])}" | \
      cat - $green_idat_files $red_idat_files | tr '\n' '\0' | xargs -0 -n 1 ln -s
    awk -F/ '{print $NF}' $green_idat_files $red_idat_files > idat_file.lines
    bcftools +gtc2vcf --idat --gtcs idat_file.lines --output "~{filebase}.idat.tsv"
    ~{if autoconvert then
      "sed -i 's/^   <NumberOfThreads>4<\\/NumberOfThreads>\\r$/   <NumberOfThreads>" + cpu + "<\\/NumberOfThreads>\\r/' \\\n" +
      "  \"/opt/AutoConvert 2.0/AutoCallConfig.xml\"\n" +
      "mkdir gtcs\n" +
      "cat $green_idat_files | tr '\\n' '\\0' | \\\n" +
      "  xargs -0 -i mono \\\n" +
      "  \"/opt/AutoConvert 2.0/AutoConvert.exe\" \\\n" +
      "  \"{}\" \\\n" +
      "  gtcs \\\n" +
      "  \"" + basename(bpm_file) + "\" \\\n" +
      "  \"" + basename(egt_file) + "\" \\\n"
    else
      "iaap-cli \\\n" +
      "  gencall \\\n" +
      "  \"" + basename(bpm_file) + "\" \\\n" +
      "  \"" + basename(egt_file) + "\" \\\n" +
      "  gtcs \\\n" +
      "  --idat-folder . \\\n" +
      "  --output-gtc \\\n" +
      (if is_gender_autocall then "  --gender-estimate-call-rate-threshold -0.1 \\\n" else "") +
      (if cpu > 1 then "  --num-threads " + cpu + " \\\n" else "")}  1>&2
    awk -F/ '{print "gtcs/"$NF}' $green_idat_files | sed 's/_Grn\.idat$/.gtc/'
    echo "~{sep="\n" flatten([[bpm_file, egt_file]])}" | \
      cat - $green_idat_files $red_idat_files | awk -F/ '{print $0"\n"$NF}' | tr '\n' '\0' | xargs -0 rm
    rm idat_file.lines
  >>>

  output {
    File idat_tsv = filebase + ".idat.tsv"
    Directory gtcs = "gtcs"
    Array[File] gtc_files = read_lines(stdout())
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

# once String sep(String, Array[String]) is implemented --chip-type can be simplified, as the current syntax is non-compliant
# https://support.terra.bio/hc/en-us/community/posts/360071476431-Terra-fails-to-delocalize-files-listed-through-read-lines-
task cel2affy {
  input {
    File? xml_file
    File? zip_file
    Array[File]+ cel_files
    File? cdf
    File? chrX_snps
    File? special_snps
    File? chrX_probes
    File? chrY_probes
    Array[String]? chip_type
    Boolean table_output = false
    String? qmethod_spec
    File? read_models_brlmmp
    String? set_analysis_name
    File? target_sketch
    String? set_gender_method
    String? em_gender
    Float? female_thresh
    Float? male_thresh
    String filebase

    String docker
    Int cpu = 1 # unfortunately this task is not parallelizable
    Int? disk_size_override
    Float? memory_override
    Int preemptible = 0
    Int maxRetries = 0
  }

  Float xml_size = size(xml_file, "GiB")
  Float zip_size = size(zip_file, "GiB")
  Float cel_size = length(cel_files) * size(cel_files[0], "GiB")
  Int disk_size = select_first([disk_size_override, ceil(10.0 + xml_size + zip_size + (if table_output then 4.0 else 3.0) * cel_size)])
  Float memory = select_first([memory_override, 2 * 7.25])

  String? zip_dir = if defined(zip_file) then basename(select_first([zip_file]), ".zip") else None
  String analysis_name = if defined(set_analysis_name) then set_analysis_name else "AxiomGT1"

  command <<<
    set -euo pipefail
    cel_files=~{write_lines(cel_files)}
    echo "~{sep="\n" select_all([xml_file, zip_file, cdf, chrX_snps, special_snps, chrX_probes, chrY_probes, read_models_brlmmp, target_sketch])}" | \
      cat - $cel_files | tr '\n' '\0' | xargs -0 -n 1 ln -s
    ~{if defined(zip_file) then "unzip -jd \"" + zip_dir + "\" \"" + basename(select_first([zip_file])) + "\" 1>&2" else ""}
    awk -F/ 'BEGIN {print "cel_files"} {print $NF}' $cel_files > cel_file.lines
    bcftools +affy2vcf --cel --chps $cel_files --output "~{filebase}.cel.tsv"
    apt-probeset-genotype \
      ~{if defined(zip_dir) then "--analysis-files-path \"" + zip_dir + "\"" else ""} \
      ~{if defined(xml_file) then "--xml-file \"" + basename(select_first([xml_file])) + "\"" else ""} \
      --cel-files cel_file.lines \
      ~{if defined(cdf) then "--cdf-file \"" + basename(select_first([cdf])) + "\"" else ""} \
      ~{if defined(chrX_snps) then "--chrX-snps \"" + basename(select_first([chrX_snps])) + "\"" else ""} \
      ~{if defined(special_snps) then "--special-snps \"" + basename(select_first([special_snps])) + "\"" else ""} \
      ~{if defined(chrX_probes) then "--chrX-probes \"" + basename(select_first([chrX_probes])) + "\"" else ""} \
      ~{if defined(chrY_probes) then "--chrY-probes \"" + basename(select_first([chrY_probes])) + "\"" else ""} \
      ~{if defined(chip_type) then "--chip-type " else ""}~{sep=" --chip-type " chip_type} \
      --table-output ~{table_output} \
      --cc-chp-output \
      ~{if table_output then "--summaries" else ""} \
      ~{if defined(qmethod_spec) then "--qmethod-spec " + qmethod_spec else ""} \
      ~{if defined(read_models_brlmmp) then "--read-models-brlmmp \"" + basename(select_first([read_models_brlmmp])) + "\"" else ""} \
      --write-models \
      ~{if defined(target_sketch) then "--target-sketch \"" + basename(select_first([target_sketch])) + "\"" else ""} \
      ~{if defined(set_analysis_name) then "--set-analysis-name " + set_analysis_name else ""} \
      ~{if defined(set_gender_method) then "--set-gender-method " + set_gender_method else ""} \
      ~{if defined(em_gender) then "--em-gender " + em_gender else ""} \
      ~{if defined(female_thresh) then "--female-thresh " + female_thresh else ""} \
      ~{if defined(male_thresh) then "--male-thresh " + male_thresh else ""}
    for sfx in snp-posteriors report~{if table_output then " calls confidences summary normalized-summary" else ""}; do
      mv "~{analysis_name}.$sfx.txt" "~{filebase}.$sfx.txt"
    done
    awk -F/ '{print "cc-chp/"$NF}' $cel_files | sed 's/\.CEL$/.~{analysis_name}.chp/'
    echo "~{sep="\n" select_all([xml_file, zip_file, cdf, chrX_snps, special_snps, chrX_probes, chrY_probes, read_models_brlmmp, target_sketch])}" | \
      cat - $cel_files | awk -F/ '{print $0"\n"$NF}' | tr '\n' '\0' | xargs -0 rm
    ~{if defined(zip_file) then "rm -r \"" + zip_dir + "\"" else ""}
    rm cel_file.lines
    touch "~{filebase}.calls.txt" \
      "~{filebase}.confidences.txt" \
      "~{filebase}.summary.txt" \
      "~{filebase}.normalized-summary.txt"
  >>>

  output {
    File cel_tsv = filebase + ".cel.tsv"
    Directory chps = "cc-chp"
    Array[File] chp_files = read_lines(stdout())
    File log_file = "apt-probeset-genotype.log"
    File snp_file = filebase + ".snp-posteriors.txt"
    File report_file = filebase + ".report.txt"
    File? calls_file = filebase + ".calls.txt"
    File? confidences_file = filebase + ".confidences.txt"
    File? summary_file = filebase + ".summary.txt"
    File? normalized_summary_file = filebase + ".normalized-summary.txt"
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

# current sep() syntax is non-compliant
task gtc2vcf {
  input {
    Array[String]? tags
    File? bpm_file
    File? csv_file
    File? egt_file
    File fasta_ref
    File fasta_ref_fai
    Int? gc_window_size
    Array[File]+ gtc_files
    Int capacity = 32768
    Boolean use_gtc_sample_names = false
    Boolean do_not_check_bpm = false
    File? sam_file
    Map[String, String]? reheader_map
    String filebase
    Boolean uncompressed = false

    String docker
    Int cpu = 1
    Int? disk_size_override
    Float? memory_override
    Int preemptible = 1
    Int maxRetries = 0
  }

  Float bpm_size = size(bpm_file, "GiB")
  Float csv_size = size(csv_file, "GiB")
  Float egt_size = size(egt_file, "GiB")
  Float ref_size = size(fasta_ref, "GiB")
  Float gtc_size = length(gtc_files) * size(gtc_files[0], "GiB")
  Float sam_size = size(sam_file, "GiB")
  Int disk_size = select_first([disk_size_override, ceil(10.0 + bpm_size + csv_size + egt_size + ref_size + 4.0 * gtc_size + sam_size)])
  # due to heavy random access to the reference genome, it is important here that enough memory to cache the reference is provided
  Float memory = select_first([memory_override, 3.5 + 2.0 * csv_size + 2.0 * egt_size + ref_size +
    length(gtc_files) * capacity * 19 / 1024 / 1024 / 1024])

  command <<<
    set -euo pipefail
    gtc_files=~{write_lines(gtc_files)}
    ~{if defined(reheader_map) then "reheader_file=" else ""}~{if defined(reheader_map) then write_map(select_first([reheader_map])) else ""}
    echo "~{sep="\n" select_all([bpm_file, csv_file, egt_file, fasta_ref, fasta_ref_fai, sam_file])}" | \
      cat - $gtc_files | tr '\n' '\0' | xargs -0 -n 1 ln -s
    awk -F/ '{print $NF}' $gtc_files > gtc_file.lines
    bcftools +gtc2vcf \
      --no-version \
      --output-type u \
      ~{if defined(tags) then "--tags " else ""}~{sep="," tags} \
      ~{if defined(bpm_file) then "--bpm \"" + basename(select_first([bpm_file])) + "\"" else ""} \
      ~{if defined(csv_file) then "--csv \"" + basename(select_first([csv_file])) + "\"" else ""} \
      ~{if defined(egt_file) then "--egt \"" + basename(select_first([egt_file])) + "\"" else ""} \
      --fasta-ref "~{basename(fasta_ref)}" \
      ~{if defined(gc_window_size) then "--gc-window-size " + select_first([gc_window_size]) else ""} \
      --gtcs gtc_file.lines \
      ~{if capacity != 32768 then "--capacity " + capacity else ""} \
      ~{if use_gtc_sample_names then "--use-gtc-sample-names" else ""} \
      ~{if do_not_check_bpm then "--do-not-check-bpm" else ""} \
      --extra "~{filebase}.gtc.tsv" \
      ~{if cpu > 1 then "--threads " + (cpu - 1) else ""} \
      ~{if defined(sam_file) then "--sam-flank \"" + basename(select_first([sam_file])) + "\"" else ""} | \
      ~{if defined(reheader_map) then "bcftools reheader --samples $reheader_file |" else ""} \
      bcftools sort --output-type u --temp-dir ./bcftools-sort.XXXXXX | \
      bcftools norm --no-version --output-type ~{if uncompressed then "u" else "b"} --output "~{filebase}.bcf" --check-ref x --fasta-ref "~{basename(fasta_ref)}"~{if cpu > 1 then " --threads " + (cpu - 1) else ""}
    bcftools index --force "~{filebase}.bcf"
    bcftools query --list-samples "~{filebase}.bcf" > "~{filebase}.sample_id.lines"
    echo "~{sep="\n" select_all([bpm_file, csv_file, egt_file, fasta_ref, fasta_ref_fai, sam_file])}" | \
      cat - $gtc_files | awk -F/ '{print $0"\n"$NF}' | tr '\n' '\0' | xargs -0 rm
    rm gtc_file.lines
  >>>

  output {
    File gtc_tsv = filebase + ".gtc.tsv"
    File vcf_file = filebase + ".bcf"
    File vcf_idx = filebase + ".bcf.csi"
    File sample_id_lines = filebase + ".sample_id.lines"
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

# the current sep() syntax is non-compliant and https://github.com/broadinstitute/cromwell/issues/5540
task affy2vcf {
  input {
    Array[String]? tags
    File csv_file
    File fasta_ref
    File fasta_ref_fai
    Int? gc_window_size
    File? calls_file
    File? confidences_file
    File? summary_file
    File? snp_file
    Array[File]? chp_files
    File? sam_file
    Map[String, String]? reheader_map
    String filebase
    Boolean uncompressed = false

    String docker
    Int cpu = 1
    Int? disk_size_override
    Float? memory_override
    Int? preemptible_override
    Int maxRetries = 0
  }

  Float csv_size = size(csv_file, "GiB")
  Float ref_size = size(fasta_ref, "GiB")
  Float calls_size = size(calls_file, "GiB")
  Float confidences_size = size(confidences_file, "GiB")
  Float summary_size = size(summary_file, "GiB")
  Float snp_size = size(snp_file, "GiB")
  Float chp_size = if defined(chp_files) then length(select_first([chp_files])) * size(select_first([chp_files])[0], "GiB") else 0.0

  Float sam_size = size(sam_file, "GiB")
  Int disk_size = select_first([disk_size_override, ceil(10.0 + csv_size + ref_size + 4.0 * calls_size + 4.0 * confidences_size + 4.0 * summary_size + snp_size + 4.0 * chp_size + sam_size)])
  # due to heavy random access to the reference genome, it is important here that enough memory to cache the reference is provided
  Float memory = select_first([memory_override, 3.5 + 2.0 * csv_size + 2.0 * snp_size + ref_size +
    (if defined(chp_files) then length(select_first([chp_files])) * 32768 / 1024 / 1024 / 1024 else 0.0)])
  Int preemptible = select_first([preemptible_override, if defined(chp_files) then 1 else 0])

  command <<<
    set -euo pipefail
    ~{if defined(chp_files) then "chp_files=" else ""}~{if defined(chp_files) then write_lines(select_first([chp_files])) else ""}
    ~{if defined(reheader_map) then "reheader_file=" else ""}~{if defined(reheader_map) then write_map(select_first([reheader_map])) else ""}
    echo "~{sep="\n" select_all([csv_file, fasta_ref, fasta_ref_fai, calls_file, confidences_file, summary_file, snp_file, sam_file])}" | \
      ~{if defined(chp_files) then "cat - $chp_files | " else ""}tr '\n' '\0' | xargs -0 -n 1 ln -s
    ~{if defined(chp_files) then "awk -F/ '{print $NF}' $chp_files > chp_file.lines" else ""}
    bcftools +affy2vcf \
      --no-version \
      --output-type u \
      ~{if defined(tags) then "--tags " else ""}~{sep="," tags} \
      --csv "~{basename(csv_file)}" \
      --fasta-ref "~{basename(fasta_ref)}" \
      ~{if defined(gc_window_size) then "--gc-window-size " + select_first([gc_window_size]) else ""} \
      ~{if defined(calls_file) then "--calls \"" + basename(select_first([calls_file])) + "\"" else ""} \
      ~{if defined(confidences_file) then "--confidences \"" + basename(select_first([confidences_file])) + "\"" else ""} \
      ~{if defined(summary_file) then "--summary \"" + basename(select_first([summary_file])) + "\"" else ""} \
      ~{if defined(snp_file) then "--snp \"" + basename(select_first([snp_file])) + "\"" else ""} \
      ~{if defined(chp_files) then "--chps chp_file.lines" else ""} \
      ~{if defined(chp_files) then "--extra " + filebase + ".chp.tsv" else ""} \
      ~{if cpu > 1 then "--threads " + (cpu - 1) else ""} \
      ~{if defined(sam_file) then "--sam-flank \"" + basename(select_first([sam_file])) + "\"" else ""} | \
      ~{if defined(reheader_map) then "bcftools reheader --samples $reheader_file |" else ""} \
      bcftools sort --output-type u --temp-dir ./bcftools-sort.XXXXXX | \
      bcftools norm --no-version --output-type ~{if uncompressed then "u" else "b"} --output "~{filebase}.bcf" --check-ref x --fasta-ref "~{basename(fasta_ref)}"~{if cpu > 1 then " --threads " + (cpu - 1) else ""}
    bcftools index --force "~{filebase}.bcf"
    bcftools query --list-samples "~{filebase}.bcf" > "~{filebase}.sample_id.lines"
    echo "~{sep="\n" select_all([csv_file, fasta_ref, fasta_ref_fai, calls_file, confidences_file, summary_file, snp_file, sam_file])}" | \
      ~{if defined(chp_files) then "cat - $chp_files | " else ""}awk -F/ '{print $0"\n"$NF}' | tr '\n' '\0' | xargs -0 rm
    ~{if defined(chp_files) then "rm chp_file.lines" else ""}
    touch "~{filebase}.chp.tsv"
  >>>

  output {
    File? chp_tsv = filebase + ".chp.tsv"
    File vcf_file = filebase + ".bcf"
    File vcf_idx = filebase + ".bcf.csi"
    File sample_id_lines = filebase + ".sample_id.lines"
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
    Int n_chrs
    File fasta_ref_fai
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

  String filebase = basename(fasta_ref_fai, ".fai")

  command <<<
    ln -s "~{fasta_ref_fai}"
    ln -s "~{genetic_map_file}"
    head -n~{n_chrs} "~{basename(fasta_ref_fai)}" | cut -f1,2 > chr2len.tsv
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
      chr_cm_len = max(df_group['CM'])
      n_win = np.ceil((chr_cm_len - ~{overlap_size_cm})/(~{max_win_size_cm} - ~{overlap_size_cm}))
      win_size = (chr_cm_len - ~{overlap_size_cm}) / n_win + ~{overlap_size_cm}
      cm_begs = (win_size - ~{overlap_size_cm}) * np.arange(1, n_win)
      cm_ends = (win_size - ~{overlap_size_cm}) * np.arange(1, n_win) + ~{overlap_size_cm}
      pos_begs = np.concatenate(([1], np.interp(cm_begs, df_group['CM'], df_group['POS']).astype(int)))
      pos_ends = np.concatenate((np.interp(cm_ends, df_group['CM'], df_group['POS']).astype(int), [chr2len[fai_chr]]))
      df_out[fai_chr] = pd.DataFrame.from_dict({'CHR': fai_chr, 'BEG': pos_begs, 'END': pos_ends})
    df = pd.concat([df_out[fai_chr] for fai_chr in chr2len.keys()])
    df[['CHR', 'BEG', 'END']].to_csv('~{filebase}.bed', sep='\t', header = False, index = False)
    CODE
    rm chr2len.tsv
    rm "~{fasta_ref_fai}" \
      "~{basename(fasta_ref_fai)}"
    rm "~{genetic_map_file}" \
      "~{basename(genetic_map_file)}"
  >>>

  output {
    File intervals_bed = filebase + ".bed"
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

# https://support.terra.bio/hc/en-us/community/posts/360071476431-Terra-fails-to-delocalize-files-listed-through-read-lines-
task vcf_scatter {
  input {
    File vcf_file
    File intervals_bed
    String? filebase
    String other = "other"
    Boolean uncompressed = false

    String docker
    Int cpu = 1
    Int? disk_size_override
    Float? memory_override
    Int preemptible = 1
    Int maxRetries = 0
  }

  Float vcf_size = size(vcf_file, "GiB")
  Int disk_size = select_first([disk_size_override, ceil(10.0 + 3.0 * vcf_size)])
  Float memory = select_first([memory_override, 3.5])

  command <<<
    set -euo pipefail
    echo "~{sep="\n" flatten([[vcf_file, intervals_bed]])}" | \
      tr '\n' '\0' | xargs -0 -n 1 ln -s
    awk -F"\t" '{print $1":"$2"-"$3"\t"NR-1}' "~{basename(intervals_bed)}" > regions.lines
    bcftools query --list-samples "~{basename(vcf_file)}" > "~{basename(vcf_file, ".bcf")}.sample_id.lines"
    bcftools annotate \
      --no-version \
      --output-type u \
      --remove ID,QUAL,INFO,^FMT/GT \
      ~{if cpu > 1 then "--threads " + (cpu - 1) else ""} \
      "~{basename(vcf_file)}" | \
    bcftools +scatter \
      --no-version \
      --output-type ~{if uncompressed then "u" else "b"} \
      --output vcfs \
      ~{if cpu > 1 then "--threads " + (cpu - 1) else ""} \
      --scatter-file regions.lines \
      ~{"--prefix " + filebase} \
      --extra "~{other}"
    awk -F"\t" '{print "vcfs/~{filebase}"$2".bcf"}' regions.lines
    echo "~{sep="\n" flatten([[vcf_file, intervals_bed]])}" | \
      awk -F/ '{print $0"\n"$NF}' | tr '\n' '\0' | xargs -0 rm
    rm regions.lines
  >>>

  output {
    File sample_id_lines = basename(vcf_file, ".bcf") + ".sample_id.lines"
    Directory vcfs = "vcfs"
    Array[File] vcf_files = read_lines(stdout())
    File other_vcf_file = "vcfs/" + filebase + other + ".bcf"
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
    Boolean uncompressed = false

    String docker
    Int cpu = 1
    Int? disk_size_override
    Float? memory_override
    Int preemptible = 1
    Int maxRetries = 0
  }

  Float vcf_size = size(vcf_files, "GiB")
  Int disk_size = select_first([disk_size_override, ceil(10 + 2.0 * vcf_size)])
  Float memory = select_first([memory_override, 3.5])

  command <<<
    set -euo pipefail
    ~{if length(vcf_files) > 1 then "vcf_files=" else ""}~{if length(vcf_files) > 1 then write_lines(vcf_files) else ""}
    ~{if length(vcf_files) > 1 then "bcftools merge \\\n" +
      "  --no-version \\\n" +
      "  --output-type " + (if uncompressed then "u" else "b") + " \\\n" +
      "  --output \"" + filebase + ".bcf\" \\\n" +
      "  --file-list $vcf_files \\\n" +
      "  --merge none \\\n" +
      "  --no-index \\\n" +
      (if cpu > 1 then "  --threads " + (cpu - 1) else "")
    else "mv \"" + vcf_files[0] + "\" \"" + filebase + ".bcf\""}
    bcftools index --force "~{filebase}.bcf"
    ~{if length(vcf_files) > 1 then "cat $vcf_files | tr '\\n' '\\0' | xargs -0 rm" else ""}
  >>>

  output {
    File vcf_file = filebase + ".bcf"
    File vcf_idx = filebase + ".bcf.csi"
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

task vcf_qc {
  input {
    File vcf_file
    File vcf_idx
    String mhc_reg
    File dup_file
    File sample_id_file
    File computed_gender_file
    File call_rate_file
    Float sample_call_rate_thr = 0.97
    Float variant_call_rate_thr = 0.97
    File? extra_xcl_vcf_file
    Boolean uncompressed = false

    String docker
    Int cpu = 1
    Int? disk_size_override
    Float? memory_override
    Int preemptible = 1
    Int maxRetries = 0
  }

  Float vcf_size = size(vcf_file, "GiB")
  Int disk_size = select_first([disk_size_override, ceil(10.0 + vcf_size)])
  Float memory = select_first([memory_override, 3.5])

  String filebase = basename(vcf_file, ".bcf")

  command <<<
    set -euo pipefail
    echo "~{sep="\n" select_all([vcf_file, vcf_idx, dup_file, sample_id_file, computed_gender_file, call_rate_file, extra_xcl_vcf_file])}" | \
      tr '\n' '\0' | xargs -0 -n 1 ln -s
    mhc_chr=$(echo ~{mhc_reg} | cut -d: -f1)
    mhc_beg=$(echo ~{mhc_reg} | cut -d: -f2 | cut -d- -f1)
    mhc_end=$(echo ~{mhc_reg} | cut -d- -f2)
    paste -d $'\t' "~{basename(sample_id_file)}" "~{basename(computed_gender_file)}" > computed_gender.map
    paste -d $'\t' "~{basename(sample_id_file)}" "~{basename(call_rate_file)}" | \
      awk -F"\t" '$2<~{sample_call_rate_thr} {print $1}' > samples_xcl.lines
    bcftools query --format "\n" "~{basename(vcf_file)}" | wc -l
    echo '##INFO=<ID=JK,Number=1,Type=Float,Description="Jukes Cantor">' | \
      bcftools annotate --no-version --output-type u --annotations "~{basename(dup_file)}" --columns CHROM,FROM,TO,JK --header-lines /dev/stdin --samples-file ^samples_xcl.lines~{if cpu > 1 then " --threads " + (cpu - 1) else ""} "~{basename(vcf_file)}" | \
      bcftools filter --no-version --output-type u --exclude "CHROM==\"$mhc_chr\" && POS>$mhc_beg && POS<$mhc_end" --soft-filter MHC | \
      bcftools +fill-tags --no-version --output-type u --targets ^Y,MT,chrY,chrM~{if cpu > 1 then " --threads " + (cpu - 1) else ""} -- --tags ExcHet,F_MISSING | \
      bcftools +mochatools --no-version --output-type u~{if cpu > 1 then " --threads " + (cpu - 1) else ""} -- --sex computed_gender.map --drop-genotypes | \
      bcftools annotate --no-version --output-type ~{if uncompressed then "u" else "b"} --output "~{filebase}.xcl.bcf" \
      --include 'FILTER!="." && FILTER!="PASS" || INFO/JK<.02 || INFO/ExcHet<1e-6 || INFO/F_MISSING>1-~{variant_call_rate_thr} || INFO/AC_Sex_Test>6' \
      --remove ^INFO/JK,^INFO/ExcHet,^INFO/F_MISSING,^INFO/AC_Sex_Test~{if cpu > 1 then " --threads " + (cpu - 1) else ""}
    bcftools index --force "~{filebase}.xcl.bcf"
    ~{if defined(extra_xcl_vcf_file) then "mv \"" + filebase + ".xcl.bcf\" \"" + filebase + ".tmp.bcf\"\n" +
      "mv \"" + filebase + ".xcl.bcf.csi\" \"" + filebase + ".tmp.bcf.csi\"\n" +
      "bcftools index --force \"" + basename(select_first([extra_xcl_vcf_file])) + "\"\n" +
      "bcftools merge --no-version --output-type " + (if uncompressed then "u" else "b") + " --output \"" + filebase + ".xcl.bcf\" --merge none \"" + filebase + ".tmp.bcf\" \"" + basename(select_first([extra_xcl_vcf_file])) + "\"\n" +
      "bcftools index --force \"" + filebase + ".xcl.bcf\"\n" +
      "rm \"" + filebase + ".tmp.bcf\" \"" + filebase + ".tmp.bcf.csi\"" else ""}
    echo "~{sep="\n" select_all([vcf_file, vcf_idx, dup_file, sample_id_file, computed_gender_file, call_rate_file, extra_xcl_vcf_file])}" | \
      awk -F/ '{print $0"\n"$NF}' | tr '\n' '\0' | xargs -0 rm
    rm computed_gender.map samples_xcl.lines
  >>>

  output {
    Int n_vars = read_int(stdout())
    File xcl_vcf_file = filebase + ".xcl.bcf"
    File xcl_vcf_idx = filebase + ".xcl.bcf.csi"
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

# bcftools query --list-samples / --samples-file /dev/stdin should be removed once https://github.com/samtools/bcftools/issues/1259 is fixed
task vcf_phase {
  input {
    Int n_vars
    Int n_smpls
    File unphased_vcf_file
    File unphased_vcf_idx
    File genetic_map_file
    Int? n_ref_smpls
    File? ref_vcf_file
    File? ref_vcf_idx
    File? xcl_vcf_file
    File? xcl_vcf_idx
    String? chr
    Boolean sequencing = false
    Boolean eagle = false
    Int? pbwtIters
    Boolean uncompressed = false

    String docker
    Int cpu = 1
    Int? disk_size_override
    Float? memory_override
    Int preemptible = 1
    Int maxRetries = 0

    Int mult = 4 # how much more memory does SHAPEIT4 consume compared to Eagle
  }

  Float vcf_size = size(unphased_vcf_file, "GiB")
  Float ref_size = size(ref_vcf_file, "GiB")
  Float xcl_size = size(xcl_vcf_file, "GiB")
  Int disk_size = select_first([disk_size_override, ceil(10.0 + 3.0 * vcf_size + ref_size + xcl_size)])
  Float memory = select_first([memory_override, 3.5 + (if eagle then 1.5 else mult * 1.5) *
    n_vars * (n_smpls + (if defined(ref_vcf_file) then select_first([n_ref_smpls]) else 0)) / 1024 / 1024 / 1024])

  String filebase = basename(unphased_vcf_file, ".bcf")
  String dollar = "$"

  command <<<
    set -euo pipefail
    echo "~{sep="\n" select_all([unphased_vcf_file, unphased_vcf_idx, genetic_map_file, ref_vcf_file, ref_vcf_idx, xcl_vcf_file, xcl_vcf_idx])}" | \
      tr '\n' '\0' | xargs -0 -n 1 ln -s
    ~{if eagle then
      "bio-eagle \\\n" +
      "  --geneticMapFile \"" + basename(genetic_map_file) + "\" \\\n" +
      "  --outPrefix \"" + filebase + ".phased\" \\\n" +
      (if cpu > 1 then "  --numThreads " + cpu + " \\\n" else "") +
      (if defined(ref_vcf_file) then "  --vcfRef \"" + basename(select_first([ref_vcf_file])) + "\" \\\n" else "") +
      "  --vcfTarget \"" + basename(unphased_vcf_file) + "\" \\\n" +
      "  --vcfOutFormat " + (if uncompressed then "u" else "b") + " \\\n" +
      "  --noImpMissing \\\n" +
      "  --outputUnphased \\\n" +
      (if defined(xcl_vcf_file) then "  --vcfExclude \"" + basename(select_first([xcl_vcf_file])) + "\" \\\n" else "") +
      (if defined(chr) then " --chrom " + chr + " \\\n" else "") +
      (if defined(pbwtIters) then "  --pbwtIters " + pbwtIters + " \\\n" else "") +
      "  1>&2"
    else
      (if defined(xcl_vcf_file) then
        "bcftools isec --no-version --output-type " +(if uncompressed then "u" else "b") + " --output \"" + filebase + ".unphased.bcf\" --complement --write 1 \"" + basename(unphased_vcf_file) + "\" \"" + basename(select_first([xcl_vcf_file])) + "\"\n" +
        "bcftools index --force \"" + filebase + ".unphased.bcf\"\n" else "xcl_vcf_file") +
      "chr=" + chr + "; zcat \"" + basename(genetic_map_file) + "\" | \\\n" +
      "  sed 's/^23/X/' | awk -v chr=" + dollar + "{chr#chr} '$1==chr {print $2,$3,$4}' > genetic_map.txt\n" +
      "shapeit4 \\\n" +
      (if cpu > 1 then "  --thread " + cpu + " \\\n" else "") +
      "  --input \"" + (if defined(xcl_vcf_file) then filebase + ".unphased.bcf" else basename(unphased_vcf_file)) + "\" \\\n" +
      (if defined(ref_vcf_file) then "  --reference \"" + basename(select_first([ref_vcf_file])) + "\" \\\n" else "") +
      "  --map genetic_map.txt \\\n" +
      "  --region " + chr + " \\\n" +
      (if sequencing then "  --sequencing \\\n" else "") +
      "  --output \"" + filebase + ".shapeit4.bcf\" \\\n" +
      "  1>&2\n" +
      "bcftools index --force \"" + filebase + ".shapeit4.bcf\"\n" +
      "bcftools query --list-samples \"" + basename(unphased_vcf_file) + "\" | \\\n" +
      "  bcftools annotate \\\n" +
      "  --no-version \\\n" +
      "  --output-type " + (if uncompressed then "u" else "b") + " \\\n" +
      "  --output \"" + filebase + ".phased.bcf\" \\\n" +
      "  --annotations \"" + filebase + ".shapeit4.bcf\" \\\n" +
      "  --columns -FMT/GT \\\n" +
      "  --samples-file /dev/stdin \\\n" +
      "  \"" + basename(unphased_vcf_file) + "\" \\\n" +
      (if cpu > 1 then "  --threads " + (cpu - 1) + "\n" else "") +
      "rm genetic_map.txt \"" + filebase + ".shapeit4.bcf\" \"" + filebase + ".shapeit4.bcf.csi\"" +
      (if defined(xcl_vcf_file) then " \"" + filebase + ".unphased.bcf\" \"" + filebase + ".unphased.bcf.csi\"" else "")}
    echo "~{sep="\n" select_all([unphased_vcf_file, unphased_vcf_idx, genetic_map_file, ref_vcf_file, ref_vcf_idx, xcl_vcf_file, xcl_vcf_idx])}" | \
      awk -F/ '{print $0"\n"$NF}' | tr '\n' '\0' | xargs -0 rm
  >>>

  output {
    File phased_vcf_file = filebase + ".phased.bcf"
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

# once suffix() is implemented the code could be simplified
# https://support.terra.bio/hc/en-us/community/posts/360071476431-Terra-fails-to-delocalize-files-listed-through-read-lines-
task vcf_split {
  input {
    File vcf_file
    Array[String] batches
    File sample_id_file
    File? ped_file
    String? rule
    Boolean uncompressed = false

    String docker
    Int? disk_size_override
    Int cpu = 1
    Float? memory_override
    Int preemptible = 1
    Int maxRetries = 0
  }

  Float vcf_size = size(vcf_file, "GiB")
  Int disk_size = select_first([disk_size_override, ceil(10.0 + 2.0 * vcf_size)])
  Float memory = select_first([memory_override, 3.5])

  String filebase = basename(vcf_file, ".bcf")

  command <<<
    set -euo pipefail
    filebases=~{write_lines(prefix(filebase + '.', batches))}
    echo "~{sep="\n" select_all([vcf_file, sample_id_file, ped_file])}" | \
      tr '\n' '\0' | xargs -0 -n 1 ln -s
    ~{if defined(ped_file) && defined(rule) then "bcftools +mendelian --output " + filebase + ".mendel.tsv --rules \"" + select_first([rule]) + "\" --ped \"" + basename(select_first([ped_file])) + "\" \"" + basename(vcf_file) + "\"" else ""}
    sed 's/\t/,/g;s/$/\t-/' "~{basename(sample_id_file)}" | paste -d $'\t' - $filebases > samples_file.txt
    ~{if defined(ped_file) then "bcftools +trio-phase --no-version --output-type u" + (if cpu > 1 then " --threads " + (cpu - 1) else "") + " \"" + basename(vcf_file) + "\" -- --ped \"" + basename(select_first([ped_file])) + "\" | \\\n" +
      "  bcftools +split --output-type " + (if uncompressed then "u" else "b") + " --output vcfs --samples-file samples_file.txt"
      else "bcftools +split --output-type " + (if uncompressed then "u" else "b") + " --output vcfs --samples-file samples_file.txt \"" + basename(vcf_file) + "\""}
    awk '{print "vcfs/"$1".bcf"}' $filebases
    echo "~{sep="\n" select_all([vcf_file, sample_id_file, ped_file])}" | \
      awk -F/ '{print $0"\n"$NF}' | tr '\n' '\0' | xargs -0 rm
    rm samples_file.txt
    touch "~{filebase}.mendel.tsv"
  >>>

  output {
    File? mendel_tsv = filebase + ".mendel.tsv"
    Directory vcfs = "vcfs"
    Array[File] vcf_files = read_lines(stdout())
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

task vcf_concat {
  input {
    Array[File]+ vcf_files
    File? other_vcf_file
    Boolean ligate = true
    String filebase
    Boolean uncompressed = false

    String docker
    Int cpu = 1
    Int? disk_size_override
    Float? memory_override
    Int preemptible = 1
    Int maxRetries = 0
  }

  Float vcf_size = size(vcf_files, "GiB")
  Float other_vcf_size = size(other_vcf_file, "GiB")
  Int disk_size = select_first([disk_size_override, ceil(10.0 + 2.0 * vcf_size + 2.0 * other_vcf_size)])
  Float memory = select_first([memory_override, 3.5])

  command <<<
    set -euo pipefail
    vcf_files=~{write_lines(vcf_files)}
    ~{if defined(other_vcf_file) then "echo \"" + other_vcf_file + "\" >> $vcf_files" else ""}
    ~{if ligate then "cat $vcf_files | tr '\\n' '\\0' | xargs -0 -n 1 bcftools index --force" else ""}
    bcftools concat \
      --no-version \
      --output-type ~{if uncompressed then "u" else "b"} \
      --output "~{filebase}.bcf" \
      --file-list $vcf_files \
      ~{if ligate then "--ligate" else ""} \
      ~{if cpu > 1 then "--threads " + (cpu - 1) else ""}
    bcftools index --force "~{filebase}.bcf"
    ~{if ligate then "cat $vcf_files | tr '\\n' '\\0' | xargs -0 -i rm \"{}.csi\"" else ""}
    cat $vcf_files | tr '\n' '\0' | xargs -0 rm
  >>>

  output {
    File vcf_file = filebase + ".bcf"
    File vcf_idx = filebase + ".bcf.csi"
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

# bcftools query --list-samples / --samples-file /dev/stdin should be removed once https://github.com/samtools/bcftools/issues/1259 is fixed
task vcf_import {
  input {
    File phased_vcf_file
    File phased_vcf_idx
    File unphased_vcf_file
    File unphased_vcf_idx
    Boolean uncompressed = false

    String docker
    Int cpu = 1
    Int? disk_size_override
    Float? memory_override
    Int preemptible = 1
    Int maxRetries = 0
  }

  String filebase = basename(unphased_vcf_file, ".bcf")

  Float gt_size = size(phased_vcf_file, "GiB")
  Float vcf_size = size(unphased_vcf_file, "GiB")
  Int disk_size = select_first([disk_size_override, ceil(10.0 + 2.0 * gt_size + 2.0 * vcf_size)])
  Float memory = select_first([memory_override, 3.5])

  command <<<
    set -euo pipefail
    echo "~{sep="\n" flatten([[phased_vcf_file, phased_vcf_idx, unphased_vcf_file, unphased_vcf_idx]])}" | \
      tr '\n' '\0' | xargs -0 -n 1 ln -s
    bcftools query --list-samples "~{basename(unphased_vcf_file)}" | \
      bcftools annotate \
      --no-version \
      --output-type ~{if uncompressed then "u" else "b"} \
      --output "~{filebase}.phased.bcf" \
      --annotations "~{basename(phased_vcf_file)}" \
      --columns -FMT/GT \
      --samples-file /dev/stdin \
      ~{if cpu > 1 then "--threads " + (cpu - 1) else ""} \
      "~{basename(unphased_vcf_file)}"
    bcftools index --force "~{filebase}.phased.bcf"
    echo "~{sep="\n" flatten([[phased_vcf_file, phased_vcf_idx, unphased_vcf_file, unphased_vcf_idx]])}" | \
      awk -F/ '{print $0"\n"$NF}' | tr '\n' '\0' | xargs -0 rm
  >>>

  output {
    File vcf_file = filebase + ".phased.bcf"
    File vcf_idx = filebase + ".phased.bcf.csi"
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

task vcf_mocha {
  input {
    String filebase
    String rule
    File vcf_file
    File vcf_idx
    File? sample_id_file
    File? computed_gender_file
    File? xcl_vcf_file
    File? xcl_vcf_idx
    File? cnp_file
    String? bdev_lrr_baf
    String? bdev_baf_phase
    Int? min_dist
    Int? adjust_baf_lrr
    Int? regress_baf_lrr
    Int? lrr_gc_order
    Float? xy_prob
    Float? err_prob
    Float? flip_prob
    Float? telomere_advantage
    Float? centromere_penalty
    String? short_arm_chrs
    Boolean use_short_arms = false
    Boolean use_centromeres = false
    Boolean use_no_rules_chrs = false
    Float? lrr_weight
    Float? lrr_hap2dip
    Float? lrr_cutoff
    Boolean uncompressed = false

    String docker
    Int cpu = 1
    Int? disk_size_override
    Float? memory_override
    Int preemptible = 1
    Int maxRetries = 0
  }

  Float vcf_size = size(vcf_file, "GiB")
  Float xcl_size = size(xcl_vcf_file, "GiB")
  Int disk_size = select_first([disk_size_override, ceil(10.0 + 2.0 * vcf_size + xcl_size)])
  Float memory = select_first([memory_override, ceil(4.0 + 0.1 * vcf_size)])

  command <<<
    set -euo pipefail
    echo "~{sep="\n" select_all([vcf_file, vcf_idx, sample_id_file, computed_gender_file, xcl_vcf_file, xcl_vcf_idx, cnp_file])}" | \
      tr '\n' '\0' | xargs -0 -n 1 ln -s
    ~{if defined(sample_id_file) && defined(computed_gender_file) then "paste -d $'\\t' \"" + basename(select_first([sample_id_file])) + "\" \"" +
      basename(select_first([computed_gender_file])) + "\" > computed_gender.map" else ""}
    bcftools +mocha \
      --rules "~{rule}" \
      --no-version \
      ~{if defined(sample_id_file) && defined(computed_gender_file) then "--sex computed_gender.map" else ""} \
      ~{if defined(xcl_vcf_file) then "--variants \"^" + basename(select_first([xcl_vcf_file])) + "\"" else ""} \
      ~{if defined(cnp_file) then "--cnp \"" + basename(select_first([cnp_file])) + "\"" else ""} \
      ~{if cpu > 1 then "--threads " + (cpu - 1) else ""} \
      --output-type ~{if uncompressed then "u" else "b"} \
      --output "~{filebase}.mocha.bcf" \
      --mosaic-calls "~{filebase}.mocha.tsv" \
      --genome-stats "~{filebase}.stats.tsv" \
      --ucsc-bed "~{filebase}.ucsc.bed" \
      "~{basename(vcf_file)}" \
      ~{"--bdev-LRR-BAF " + bdev_lrr_baf} \
      ~{"--bdev-BAF-phase " + bdev_baf_phase} \
      ~{"--min-dist " + min_dist} \
      ~{"--adjust-BAF-LRR " + adjust_baf_lrr} \
      ~{"--regress-BAF-LRR " + regress_baf_lrr} \
      ~{"--LRR-GC-order " + lrr_gc_order} \
      ~{"--xy-prob " + xy_prob} \
      ~{"--err-prob " + err_prob} \
      ~{"--flip-prob " + flip_prob} \
      ~{"--telomere-advantage " + telomere_advantage} \
      ~{"--centromere-penalty " + centromere_penalty} \
      ~{"--short-arm-chrs " + short_arm_chrs} \
      ~{if use_short_arms then "--use-short-arms" else ""} \
      ~{if use_centromeres then "--use-centromeres" else ""} \
      ~{if use_no_rules_chrs then "--use-no-rules-chrs" else ""} \
      ~{"--LRR-weight " + lrr_weight} \
      ~{"--LRR-hap2dip " + lrr_hap2dip} \
      ~{"--LRR-cutoff " + lrr_cutoff}
    bcftools index --force "~{filebase}.mocha.bcf"
    echo "~{sep="\n" select_all([vcf_file, vcf_idx, sample_id_file, computed_gender_file, xcl_vcf_file, xcl_vcf_idx, cnp_file])}" | \
      awk -F/ '{print $0"\n"$NF}' | tr '\n' '\0' | xargs -0 rm
    ~{if defined(sample_id_file) && defined(computed_gender_file) then "rm computed_gender.map" else ""}
  >>>

  output {
    File mocha_vcf_file = filebase + ".mocha.bcf"
    File mocha_vcf_idx = filebase + ".mocha.bcf.csi"
    File calls_tsv = filebase + ".mocha.tsv"
    File stats_tsv = filebase + ".stats.tsv"
    File ucsc_bed = filebase + ".ucsc.bed"
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

task mocha_summary {
  input {
    File calls_tsv
    File stats_tsv
    Array[File]+ ucsc_beds
    String filebase

    String docker
    Int cpu = 1
    Int disk_size = 10
    Float? memory = 3.5
    Int preemptible = 1
    Int maxRetries = 0
  }

  command <<<
    set -euo pipefail
    ucsc_files=~{write_lines(ucsc_beds)}
    echo "~{sep="\n" flatten([[calls_tsv, stats_tsv]])}" | \
      tr '\n' '\0' | xargs -0 -n 1 ln -s
    summary_plot.R \
      --pdf "~{filebase}.pdf" \
      --stats "~{basename(stats_tsv)}" \
      --calls "~{basename(calls_tsv)}"
    cat $ucsc_files | tr '\n' '\0' | xargs -0 cat | \
      awk '{if ($0~"^track") track=$0; else bed[track]=bed[track]$0"\n"}
      END {for (track in bed) printf track"\n"bed[track]}' > "~{filebase}.ucsc.bed"
    echo "~{sep="\n" flatten([[calls_tsv, stats_tsv]])}" | \
      awk -F/ '{print $0"\n"$NF}' | tr '\n' '\0' | xargs -0 rm
  >>>

  output {
    File ucsc_bed = filebase + ".ucsc.bed"
    File pdf = filebase + ".pdf"
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
