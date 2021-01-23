version development

## Copyright (c) 2021 Giulio Genovese
##
## Version 2021-01-20
##
## Contact Giulio Genovese <giulio.genovese@gmail.com>
##
## This WDL workflow runs impute5 or beagle5 on a set of VCFs
##
## Cromwell version support
## - Successfully tested on v55
##
## Distributed under terms of the MIT License

struct Reference {
  String name
  File fasta_fai
  Int n_chrs
  File genetic_map_file
  String panel_pfx
  String panel_sfx
  String panel_idx
  Int n_panel_smpls
  Array[Int] n_panel_vars
}

workflow impute {
  input {
    String sample_set_id
    String mode = "pvcf" # pvcf imp
    String target = "imp" # imp ext
    String data_path = ""
    File batch_tsv_file # batch_id path vcf vcf_index xcl_vcf xcl_vcf_index
    Float max_win_size_cm = 10.0
    Float overlap_size_cm = 2.0
    Boolean? convert_panel_override
    String format_id = "Bdev_Phase"
    String ext_string = "bdev"

    String ref_name = "GRCh38"
    String ref_path = ""
    String? ref_fasta_fai
    Int? n_chrs
    String? genetic_map_file
    String? panel_pfx
    String? panel_sfx
    String? panel_idx
    Int? n_panel_smpls
    Array[Int]? n_panel_vars

    Boolean beagle = true
    Boolean out_gp = true
    Boolean out_ap = false
    String? impute_extra_args
    String basic_bash_docker = "ubuntu:latest"
    String pandas_docker = "amancevice/pandas:slim"
    String bcftools_docker = "us.gcr.io/mccarroll-mocha/mocha:1.11-20210120"
    String beagle5_docker = "us.gcr.io/mccarroll-mocha/beagle5:1.11-20210120"
    String impute5_docker = "us.gcr.io/mccarroll-mocha/impute5:1.11-20210120"
  }

  String ref_path_with_sep = ref_path + (if ref_path == "" || sub(ref_path, "/$", "") != ref_path then "" else "/")
  Reference ref = object {
    name: ref_name,
    fasta_fai: ref_path_with_sep + select_first([ref_fasta_fai, if ref_name == "GRCh38" then "GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai" else if ref_name == "GRCh37" then "human_g1k_v37.fasta.fai" else None]),
    n_chrs: select_first([n_chrs, 23]),
    genetic_map_file: ref_path_with_sep + select_first([genetic_map_file, if ref_name == "GRCh38" then "genetic_map_hg38_withX.txt.gz" else if ref_name == "GRCh37" then "genetic_map_hg19_withX.txt.gz" else None]),
    panel_pfx: ref_path_with_sep + select_first([panel_pfx, if ref_name == "GRCh38" then "ALL." else if ref_name == "GRCh37" then "ALL.chr" else None]),
    panel_sfx: select_first([panel_sfx, if ref_name == "GRCh38" then "_GRCh38.genotypes.20170504.bcf" else if ref_name == "GRCh37" then ".phase3_integrated.20130502.genotypes.bcf" else None]),
    panel_idx: select_first([panel_idx, ".csi"]),
    n_panel_smpls: select_first([n_panel_smpls, 2504]),
    n_panel_vars: select_first([n_panel_vars, if ref_name == "GRCh38" then [3757099, 4053723, 3353433, 3334752, 3029405, 2951809, 2749818, 2649365, 2053531, 2352560, 2330076, 2183796, 1660088, 1532468, 1402526, 1547663, 1342968, 1318678, 1081617, 1046636, 654342, 650718, 1140847]
      else if ref_name == "GRCh37" then [3739969, 4059474, 3357899, 3339176, 3033495, 2955302, 2755302, 2653225, 2063366, 2335239, 2335653, 2243891, 1662376, 1537461, 1405068, 1549896, 1346874, 1320491, 1085279, 1048351, 656422, 653034, 1911917] else None])
  }

  Array[Array[String]] ref_fasta_fai_tbl = transpose(read_tsv(ref.fasta_fai))

  # read table with batches information (scatter could be avoided if there was a tail() function)
  Array[Array[String]] batch_tsv = read_tsv(batch_tsv_file)
  Int n_batches = length(batch_tsv)-1
  scatter (idx in range(n_batches)) { Array[String] batch_tsv_rows = batch_tsv[(idx+1)] }
  Map[String, Array[String]] batch_tbl = as_map(zip(batch_tsv[0], transpose(batch_tsv_rows)))

  Boolean convert_panel = beagle || select_first([convert_panel_override, length(batch_tbl["vcf"]) > 2])

  # compute data paths for each batch, if available (scatter could be avoided if there was a contains_key() function)
  scatter (key in keys(batch_tbl)) { Boolean? is_key_equal_path = if key == "path" then true else None }
  scatter (idx in range(n_batches)) {
    String data_paths = if length(select_all(is_key_equal_path))>0 then batch_tbl["path"][idx] else data_path
    String data_paths_with_sep = data_paths + (if data_paths == "" || sub(data_paths, "/$", "") != data_paths then "" else "/")
  }

  if (mode == "pvcf") {
    call ref_scatter {
      input:
        n_chrs = ref.n_chrs,
        ref_fasta_fai = ref.fasta_fai,
        genetic_map_file = ref.genetic_map_file,
        max_win_size_cm = max_win_size_cm,
        overlap_size_cm = overlap_size_cm,
        docker = pandas_docker
    }
    Array[Array[String]] intervals_tbl = transpose(read_tsv(ref_scatter.intervals_bed))
    Array[String] chrs = intervals_tbl[0]
    Array[String] begs = intervals_tbl[1]
    Array[String] ends = intervals_tbl[2]
    # this is a trick to table how many intervals you will use for each chromosome
    Map[String, Array[Int]] chr_map = collect_by_key(zip(chrs, range(length(chrs))))

    # scatter reference panel
    scatter (idx in range(ref.n_chrs)) {
      String chr = ref_fasta_fai_tbl[0][idx]
      if (length(chr_map[chr]) > 1) {
        call vcf_scatter as panel_scatter {
          input:
            vcf_file = ref.panel_pfx + chr + ref.panel_sfx,
            intervals_bed = ref_scatter.intervals_bed,
            chr = chr,
            docker = bcftools_docker
        }
      }
      Array[File] panel_scatter_vcf_files = select_first([panel_scatter.vcf_files, [ref.panel_pfx + chr + ref.panel_sfx]])
      Array[File] panel_scatter_vcf_idxs = select_first([panel_scatter.vcf_idxs, [ref.panel_pfx + chr + ref.panel_sfx + ref.panel_idx]])
      Array[Int] n_panel_scatter_vars = select_first([panel_scatter.n_vars, [ref.n_panel_vars[idx]]])
    }
    Array[File] panel_flatten_vcf_files = flatten(panel_scatter_vcf_files)
    Array[File] panel_flatten_vcf_idxs = flatten(panel_scatter_vcf_idxs)

    # convert reference panel
    scatter (idx in range(length(chrs))) {
      call init_panel {
        input:
          vcf_file = panel_flatten_vcf_files[idx],
          vcf_idx = panel_flatten_vcf_idxs[idx],
          chr = chrs[idx],
          convert_panel = convert_panel,
          beagle = beagle,
          docker = if beagle then beagle5_docker else impute5_docker
      }
      File panel_genotypes_files = select_first([init_panel.panel_genotypes_file, panel_flatten_vcf_files[idx]])
      File? panel_genotypes_idxs = if convert_panel && beagle then None else select_first([init_panel.panel_genotypes_idx, panel_flatten_vcf_idxs[idx]])
    }

     # scatter target genotypes
    scatter (idx in range(n_batches)) {
      call vcf_scatter {
        input:
          vcf_file = data_paths_with_sep[idx] + batch_tbl["vcf"][idx],
          intervals_bed = ref_scatter.intervals_bed,
          docker = bcftools_docker
      }
    }

    # impute genotypes
    Array[Int] n_panel_flatten_vars = flatten(n_panel_scatter_vars)
    Map[String, Int] chr2int = as_map(zip(ref_fasta_fai_tbl[0], range(length(ref_fasta_fai_tbl[0]))))
    scatter (p in cross(range(n_batches), range(length(chrs)))) {
      if (vcf_scatter.n_vars[p.left][p.right] > 0) {
        Int cross_idx = p.left * ref.n_chrs + chr2int[(chrs[p.right])]
        call vcf_impute {
          input:
            n_vars = vcf_scatter.n_vars[p.left][p.right],
            n_smpls = vcf_scatter.n_smpls[p.left],
            pvcf_file = vcf_scatter.vcf_files[p.left][p.right],
            genetic_map_file = ref.genetic_map_file,
            n_panel_vars = n_panel_flatten_vars[p.right],
            n_panel_smpls = ref.n_panel_smpls,
            panel_genotypes_file = panel_genotypes_files[p.right],
            panel_genotypes_idx = panel_genotypes_idxs[p.right],
            panel_sites_file = init_panel.panel_sites_file[p.right],
            panel_sites_idx = init_panel.panel_sites_idx[p.right],
            xcl_vcf_file = data_paths_with_sep[p.left] + batch_tbl["xcl_vcf"][p.left],
            xcl_vcf_idx = data_paths_with_sep[p.left] + batch_tbl["xcl_vcf_index"][p.left],
            chr = chrs[p.right],
            region = chrs[p.right] + ":" + begs[p.right] + "-" + ends[p.right],
            beagle = beagle,
            out_gp= out_gp,
            out_ap= out_ap,
            impute_extra_args = impute_extra_args,
            docker = if beagle then beagle5_docker else impute5_docker
        }
      }
    }

    scatter (p in cross(range(n_batches), range(length(chrs)))) {  }
    Map[Int, Array[File]] idx2vcf_files = collect_by_key(zip(select_all(cross_idx), select_all(vcf_impute.imp_vcf_file)))
    Map[Int, Array[File]] idx2vcf_idxs = collect_by_key(zip(select_all(cross_idx), select_all(vcf_impute.imp_vcf_idx)))
    scatter (p in cross(range(n_batches), range(ref.n_chrs))) {
      Int batch_idx = p.left
      if (length(idx2vcf_files[(batch_idx * ref.n_chrs + p.right)]) > 1) {
        call vcf_concat as chr_concat {
          input:
          vcf_files = idx2vcf_files[(batch_idx * ref.n_chrs + p.right)],
          ref_fasta_fai = ref.fasta_fai,
          filebase = basename(basename(batch_tbl["vcf"][batch_idx], ".bcf"), ".vcf.gz") + ".imp" + "." + ref_fasta_fai_tbl[0][p.right],
          docker = bcftools_docker
        }
      }

      if (target == "ext") {
        call vcf_extend as chr_extend {
          input:
            vcf_file = select_first([chr_concat.vcf_file, idx2vcf_files[(batch_idx * ref.n_chrs + p.right)][0]]),
            vcf_idx = select_first([chr_concat.vcf_idx, idx2vcf_idxs[(batch_idx * ref.n_chrs + p.right)][0]]),
            annot_vcf_file = data_paths_with_sep[batch_idx] + batch_tbl["vcf"][batch_idx],
            annot_vcf_idx = data_paths_with_sep[batch_idx] + batch_tbl["vcf_index"][batch_idx],
            xcl_vcf_file = data_paths_with_sep[batch_idx] + batch_tbl["xcl_vcf"][batch_idx],
            xcl_vcf_idx = data_paths_with_sep[batch_idx] + batch_tbl["xcl_vcf_index"][batch_idx],
            format_id = format_id,
            ext_string = ext_string,
            docker = bcftools_docker
        }
      }

      File chr_concat_vcf_files = select_first([chr_extend.ext_vcf_file, chr_concat.vcf_file, idx2vcf_files[(batch_idx * ref.n_chrs + p.right)][0]])
    }

    Map[Int, Array[File]] batch2dose_vcf_files = collect_by_key(zip(batch_idx, chr_concat_vcf_files))
    scatter (idx in range(n_batches)) {
      call vcf_concat {
        input:
          vcf_files = batch2dose_vcf_files[idx],
          ligate = false,
          filebase = basename(basename(batch_tbl["vcf"][idx], ".bcf"), ".vcf.gz") + ".imp" + if target == "ext" then ext_string else "",
          docker = bcftools_docker
      }
    }
  }

  if (mode != "pvcf" && target == "ext") {
    scatter (idx in range(n_batches)) {
      call vcf_extend {
        input:
          vcf_file = data_paths_with_sep[idx] + batch_tbl["imp_vcf"][idx],
          vcf_idx = data_paths_with_sep[idx] + batch_tbl["imp_vcf_idx"][idx],
          annot_vcf_file = data_paths_with_sep[idx] + batch_tbl["vcf"][idx],
          annot_vcf_idx = data_paths_with_sep[idx] + batch_tbl["vcf_index"][idx],
          xcl_vcf_file = data_paths_with_sep[idx] + batch_tbl["xcl_vcf"][idx],
          xcl_vcf_idx = data_paths_with_sep[idx] + batch_tbl["xcl_vcf_index"][idx],
          format_id = format_id,
          ext_string = ext_string,
          docker = bcftools_docker
      }
    }
  }

  # generate a table summarizing the main output files and serialize the table to disk
  # vcf_files and vcf_idxs are defined in the output section
  scatter (idx in range(n_batches)) {
    String basename_vcf_files = basename(vcf_files[idx])
    String basename_vcf_idxs = basename(vcf_idxs[idx])
    String basename_xcl_vcf_files = basename(batch_tbl["xcl_vcf"][idx])
    String basename_xcl_vcf_idxs = basename(batch_tbl["xcl_vcf_index"][idx])
  }
  Map[String, Array[String]] output_map = {
    "batch_id": batch_tbl["batch_id"],
    "imp_vcf": basename_vcf_files,
    "imp_vcf_index": basename_vcf_idxs,
    "xcl_vcf": basename_xcl_vcf_files,
    "xcl_vcf_index": basename_xcl_vcf_idxs,
  }
  scatter (key in ["batch_id", "imp_vcf", "imp_vcf_index", "xcl_vcf", "xcl_vcf_index"]) { Array[String] output_tsv_cols = output_map[key] }
  call write_tsv {
    input:
      tsv = flatten([[["batch_id", "imp_vcf", "imp_vcf_index", "xcl_vcf", "xcl_vcf_index"]], transpose(output_tsv_cols)]),
      filebase = sample_set_id + ".batch",
      docker = basic_bash_docker
  }

  output {
    Array[File] vcf_files = select_first([vcf_extend.ext_vcf_file, vcf_concat.vcf_file])
    Array[File] vcf_idxs = select_first([vcf_extend.ext_vcf_idx, vcf_concat.vcf_idx])
    Array[File]? logs = if mode == "pvcf" then select_all(select_first([vcf_impute.log])) else None
    File output_tsv_file = write_tsv.file
  }
}

task write_tsv {
  input {
    Array[Array[String]] tsv
    String filebase

    String docker
    Int cpu = 1
    Int disk_size = 10
    Float memory = 3.5
    Int preemptible = 1
    Int maxRetries = 0
  }

  command <<<
    mv ~{write_tsv(tsv)} "~{filebase}.tsv"
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

task ref_scatter {
  input {
    Int n_chrs
    File ref_fasta_fai
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

  String filebase = basename(ref_fasta_fai, ".fai")

  command <<<
    set -euo pipefail
    mv "~{ref_fasta_fai}" .
    mv "~{genetic_map_file}" .
    head -n~{n_chrs} "~{basename(ref_fasta_fai)}" | cut -f1,2 > chr2len.tsv
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
      pos_begs = np.concatenate(([1], np.interp(cm_begs, df_group['CM'], df_group['POS'], period = np.inf).astype(int)))
      pos_ends = np.concatenate((np.interp(cm_ends, df_group['CM'], df_group['POS'], period = np.inf).astype(int), [chr2len[fai_chr]]))
      df_out[fai_chr] = pd.DataFrame.from_dict({'CHR': fai_chr, 'BEG': pos_begs, 'END': pos_ends})
    df = pd.concat([df_out[fai_chr] for fai_chr in chr2len.keys()])
    df[['CHR', 'BEG', 'END']].to_csv('~{filebase}.bed', sep='\t', header = False, index = False)
    CODE
    rm chr2len.tsv
    rm "~{basename(ref_fasta_fai)}"
    rm "~{basename(genetic_map_file)}"
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

task vcf_scatter {
  input {
    File vcf_file
    File intervals_bed
    String? chr
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
    --remove ID,QUAL,INFO,^FMT/GT \
    ~{if cpu > 1 then "--threads " + (cpu - 1) else ""} \
    "~{basename(vcf_file)}" | \
    bcftools +scatter \
    --no-version \
    --output-type ~{if uncompressed then "u" else "b"} \
    --output vcfs \
    ~{if cpu > 1 then "--threads " + (cpu - 1) else ""} \
    --scatter-file regions.lines \
    --prefix "~{filebase}."
    while read reg i; do
      bcftools query --format "\n" "vcfs/~{filebase}.$i.bcf" | wc -l
      bcftools index --force "vcfs/~{filebase}.$i.bcf"
    done < regions.lines > n_vars.lines
    cut -f2 regions.lines | sed 's/^/vcfs\/~{filebase}./;s/$/.bcf/'
    cut -f2 regions.lines | sed 's/^/vcfs\/~{filebase}./;s/$/.bcf.csi/' > vcf_idxs.lines
    rm "~{basename(vcf_file)}"
    rm "~{basename(intervals_bed)}"
    rm regions.lines
  >>>

  output {
    Int n_smpls = read_int("n_smpls.int")
    Array[Int] n_vars = read_lines("n_vars.lines")
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

task init_panel {
  input {
    File vcf_file
    File? vcf_idx
    String chr
    Boolean convert_panel = false
    Boolean beagle = false

    String docker
    Int cpu = 1
    Int? disk_size_override
    Float memory = 3.5
    Int preemptible = 1
    Int maxRetries = 0
  }

  Float vcf_size = size(vcf_file, "GiB")
  Int disk_size = select_first([disk_size_override, ceil(10.0 + 3.0 * vcf_size)])
  String filebase = basename(basename(vcf_file, ".bcf"), ".vcf.gz")

  command <<<
    set -euo pipefail
    mv "~{vcf_file}" .
    ~{if defined(vcf_idx) then "mv \"" + select_first([vcf_idx]) + "\" ." else ""}
    bcftools view \
      --no-version \
      --drop-genotypes \
      --output-type u \
      "~{basename(vcf_file)}" | \
      bcftools annotate \
        --no-version \
        --output "~{filebase}.sites.bcf" \
        --output-type b \
        --remove ID,QUAL,FILTER,INFO
    bcftools index --force "~{filebase}.sites.bcf"
    ~{if convert_panel then
        (if beagle then
          (if chr != "X" && chr != "chrX" then
            "bcftools view --no-version \"" + basename(vcf_file) + "\" | \\\n"
          else
            "bcftools +fixploidy --no-version \"" + basename(vcf_file) + "\" | \\\n" +
            "  sed 's/0\\/0/0|0/g;s/1\\/1/1|1/g' | \\\n") +
          "  java -jar /usr/bin/bref3.jar > \"" + filebase +  ".bref3\""
        else
          "imp5Converter \\\n" +
          "  --h \"" + basename(vcf_file) + "\" \\\n" +
          "  --r " + chr + " \\\n" +
          "  --o \"" + filebase + ".imp5\"")
      else ""}
    rm "~{basename(vcf_file)}"
    ~{if defined(vcf_idx) then "rm \"" + basename(select_first([vcf_idx])) + "\"" else ""}
  >>>

  output {
    File? panel_genotypes_file = if convert_panel then filebase + (if beagle then ".bref3" else ".imp5") else None
    File? panel_genotypes_idx = if convert_panel && !beagle then filebase + ".imp5.idx" else None
    File panel_sites_file = filebase + ".sites.bcf"
    File panel_sites_idx = filebase + ".sites.bcf.csi"
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

task vcf_impute {
  input {
    Int n_vars
    Int n_smpls
    File pvcf_file
    File genetic_map_file
    Int n_panel_vars
    Int n_panel_smpls
    File panel_genotypes_file
    File? panel_genotypes_idx
    File panel_sites_file
    File panel_sites_idx
    File? xcl_vcf_file
    File? xcl_vcf_idx
    String region
    String chr
    Boolean beagle = true
    Boolean out_gp = true
    Boolean out_ap = false
    String? impute_extra_args

    String docker
    Int cpu = 1
    Int? disk_size_override
    Float? memory_override
    Int preemptible = 1
    Int maxRetries = 0

    Int max_ram_percentage = 90
    Float? mult_override
  }

  Float pvcf_size = size(pvcf_file, "GiB")
  Float panel_size = size(panel_genotypes_file, "GiB")
  Float xcl_size = size(xcl_vcf_file, "GiB")
  Int disk_size = select_first([disk_size_override, ceil(10.0 + 3.0 * pvcf_size + (1.0 + 2.0 * n_smpls / n_panel_smpls) * panel_size + xcl_size)])
  Float mult = select_first([mult_override, if beagle then 30.0 else 10.0])
  Float memory = select_first([memory_override, 3.5 + mult * n_panel_vars * (n_smpls + n_panel_smpls) / 1024 / 1024 / 1024])
  String filebase = basename(basename(pvcf_file, ".bcf"), ".vcf.gz")

  command <<<
    set -euo pipefail
    echo "~{sep="\n" select_all([pvcf_file, genetic_map_file, panel_genotypes_file, panel_genotypes_idx, panel_sites_file, panel_sites_idx, xcl_vcf_file, xcl_vcf_idx])}" | \
      tr '\n' '\0' | xargs -0 mv -t .
    bcftools index --force "~{basename(pvcf_file)}"
    ~{if defined(xcl_vcf_file) then
      "bcftools isec \\\n" +
      "  --no-version \\\n" +
      "  --output-type b \\\n" +
      "  --output \"" + filebase + ".isec.bcf\" \\\n" +
      "  --complement \\\n" +
      "  --exclude \"N_ALT>1\" \\\n" +
      "  --write 1 \\\n" +
      "  \"" + basename(pvcf_file) + "\" \\\n" +
      "  \"" + basename(select_first([xcl_vcf_file])) + "\"\n" +
      "mv \"" + filebase + ".isec.bcf\" \"" + basename(pvcf_file) + "\"\n" +
      "bcftools index --force \"" + basename(pvcf_file) + "\"\n" else ""}
    n_vars=$(bcftools isec \
      --no-version \
      --output-type u \
      --nfiles 2 \
      --write 1 \
      "~{basename(pvcf_file)}" \
      "~{basename(panel_sites_file)}" | \
      bcftools query --format "\n" | wc -l)
    if [ $n_vars == 0 ]; then
      cp "~{basename(pvcf_file)}" "~{filebase}.imp.bcf"
    else
      mkdir logs
    ~{if beagle then
      "  bcftools norm \\\n" +
      "    --no-version \\\n" +
      "    --rm-dup exact \\\n" +
      "    --output-type z \\\n" +
      "    --output \"" + filebase + ".vcf.gz\" \\\n" +
      "    \"" + basename(pvcf_file) + "\"\n" +
      "  chr=" + chr +"; zcat \"" + basename(genetic_map_file) + "\" | \\\n" +
      "    sed 's/^23/X/' | awk -v chr=$chr '$1==chr || \"chr\"$1==chr {print chr,\".\",$4,$2}' > genetic_map.txt\n" +
      "  java -XX:MaxRAMPercentage=" + max_ram_percentage + " \\\n" +
      "    -jar /usr/bin/beagle.jar \\\n" +
      "    gt=\"" + filebase + ".vcf.gz\" \\\n" +
      "    ref=\"" + basename(panel_genotypes_file) + "\" \\\n" +
      "    out=\"" + filebase + ".imp\" \\\n" +
      "    map=genetic_map.txt \\\n" +
      "    chrom=" + region + " \\\n" +
      (if cpu > 1 then "    nthread=" + cpu + " \\\n" else "") +
      (if out_ap then "    ap=true \\\n" else "") +
      (if out_gp then "    gp=true \\\n" else "") +
      (if defined(impute_extra_args) then "    " + impute_extra_args + " \\\n" else "") +
      "    1>&2\n" +
      "  bcftools index --force --tbi \"" + filebase + ".imp.vcf.gz\"\n" +
      "  bcftools view \\\n" +
      "    --no-version \\\n" +
      "    --output-type b \\\n" +
      "    --output \"" + filebase + ".imp.bcf\" \\\n" +
      "    \"" + filebase + ".imp.vcf.gz\"\n" +
      "  rm \"" + filebase + ".vcf.gz\" \"" + filebase + ".imp.vcf.gz\" \"" + filebase + ".imp.vcf.gz.tbi\"\n" +
      "  mv \"" + filebase + ".imp.log\" \"logs/" + filebase + ".imp.log\""
    else
      "  chr=" + chr +"; zcat \"" + basename(genetic_map_file) + "\" | \\\n" +
      "  sed 's/^23/X/' | awk -v chr=$chr -v OFS=\"\\t\" 'BEGIN {print \"pos\",\"chr\",\"cM\"}\n" +
      "    $1==chr || \"chr\"$1==chr {print $2,chr,$4}' > genetic_map.txt\n" +
      "  impute5 \\\n" +
      "    --h \"" + basename(panel_genotypes_file) + "\" \\\n" +
      "    --m genetic_map.txt \\\n" +
      "    --g \"" + basename(pvcf_file) + "\" \\\n" +
      "    --r " + region + " \\\n" +
      "    --l \"logs/" + filebase + ".imp.log\" \\\n" +
      (if out_gp then "    --out-gp-field \\\n" else "") +
      (if out_ap then "    --out-ap-field \\\n" else "") +
      "    --o \"" + filebase + ".imp.bcf\" \\\n" +
      (if cpu > 1 then "  --threads " + cpu + " \\\n" else "") +
      (if defined(impute_extra_args) then "  " + impute_extra_args + " \\\n" else "") +
      "    1>&2"}
    fi
    rm "~{basename(pvcf_file)}.csi" genetic_map.txt
    bcftools index --force "~{filebase}.imp.bcf"
    echo "~{sep="\n" select_all([pvcf_file, genetic_map_file, panel_genotypes_file, panel_genotypes_idx, panel_sites_file, panel_sites_idx, xcl_vcf_file, xcl_vcf_idx])}" | \
      sed 's/^.*\///' | tr '\n' '\0' | xargs -0 rm
  >>>

  output {
    File imp_vcf_file = filebase + ".imp.bcf"
    File imp_vcf_idx = filebase + ".imp.bcf.csi"
    File log = "logs/" + filebase + ".imp.log"
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
    File? ref_fasta_fai
    Boolean ligate = true
    String filebase
    Boolean uncompressed = false

    String docker
    Int cpu = 1
    Int? disk_size_override
    Float? memory_override
    Int preemptible = 1
    Int maxRetries = 0

    Float mult = 4.0
  }

  Float vcf_size = size(vcf_files, "GiB")
  Int disk_size = select_first([disk_size_override, ceil(10.0 + 2.0 * vcf_size)])
  Float memory = select_first([memory_override, 3.5 + if ligate then mult * 2.0 * vcf_size else 0.0])

  command <<<
    set -euo pipefail
    vcf_files=~{write_lines(vcf_files)}
    ~{if defined(ref_fasta_fai) then "mv \"" + select_first([ref_fasta_fai]) + "\" ." else ""}
    cat $vcf_files | tr '\n' '\0' | xargs -0 mv -t .
    sed -i 's/^.*\///' $vcf_files
    ~{if ligate then "cat $vcf_files | tr '\\n' '\\0' | xargs -0 -n 1 bcftools index --force" else ""}
    bcftools concat \
      --no-version \
      --output-type ~{if uncompressed then "u" else "b"} \
      --file-list $vcf_files \
      ~{if ligate then "--ligate" else ""} \
      ~{if cpu > 1 then "--threads " + (cpu - 1) else ""} \
      --output "~{filebase}.bcf"
    ~{if defined(ref_fasta_fai) then
      "mv  \"" + filebase + ".bcf\" \"" + filebase + ".tmp.bcf\"\n" +
      "bcftools reheader \\\n" +
      "  --fai \"" + basename(select_first([ref_fasta_fai])) + "\" \\\n" +
      "  --output \"" + filebase + ".bcf\" \\\n" +
      "  \"" + filebase + ".tmp.bcf\"\n" +
      "rm \"" + filebase + ".tmp.bcf\""
      else ""}
    bcftools index --force "~{filebase}.bcf"
    ~{if ligate then "cat $vcf_files | sed 's/$/.csi/' | tr '\\n' '\\0' | xargs -0 rm" else ""}
    ~{if defined(ref_fasta_fai) then "rm \"" + basename(select_first([ref_fasta_fai])) + "\"" else ""}
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

task vcf_extend {
  input {
    File vcf_file
    File vcf_idx
    File annot_vcf_file
    File annot_vcf_idx
    File? xcl_vcf_file
    File? xcl_vcf_idx
    String format_id
    String ext_string
    Int dist = 500000
    Boolean uncompressed = false

    String docker
    Int cpu = 1
    Int? disk_size_override
    Float? memory_override
    Int preemptible = 1
    Int maxRetries = 0
  }

  Float vcf_size = size(vcf_file, "GiB")
  Float annot_vcf_size = size(annot_vcf_file, "GiB")
  Int disk_size = select_first([disk_size_override, ceil(10.0 + 2.0 * (vcf_size + annot_vcf_size))])
  Float memory = select_first([memory_override, 3.5])
  String filebase = basename(basename(vcf_file, ".bcf"), ".vcf.gz")
  String annot_filebase = basename(basename(annot_vcf_file, ".bcf"), ".vcf.gz")

  command <<<
    set -euo pipefail
    echo "~{sep="\n" select_all([vcf_file, vcf_idx, annot_vcf_file, annot_vcf_idx, xcl_vcf_file, xcl_vcf_idx])}" | \
      tr '\n' '\0' | xargs -0 mv -t .
    ~{if defined(xcl_vcf_file) then
      "bcftools isec \\\n" +
      "  --no-version \\\n" +
      "  --output-type b \\\n" +
      "  --output \"" + annot_filebase + ".isec.bcf\" \\\n" +
      "  --complement \\\n" +
      "  --exclude \"N_ALT>1\" \\\n" +
      "  --write 1 \\\n" +
      "  \"" + basename(annot_vcf_file) + "\" \\\n" +
      "  \"" + basename(select_first([xcl_vcf_file])) + "\"\n" +
      "bcftools index --force \"" + annot_filebase + ".isec.bcf\"\n" else ""}
    bcftools annotate \
      --no-version \
      --output-type u \
      --columns FMT/~{format_id} \
      ~{basename(vcf_file)} \
      --annotations "~{if defined(xcl_vcf_file) then annot_filebase + ".isec.bcf" else basename(annot_vcf_file)}" |
    bcftools +extendFMT \
      --no-version \
      --output-type ~{if uncompressed then "u" else "b"} \
      --output "~{filebase}.~{ext_string}.bcf" \
      --format ~{format_id} \
      --phase \
      --dist ~{dist}
    bcftools index --force "~{filebase}.~{ext_string}.bcf"
    ~{if defined(xcl_vcf_file) then "rm \"" + annot_filebase + ".isec.bcf\" \"" + annot_filebase + ".isec.bcf.csi\"" else ""}
    echo "~{sep="\n" select_all([vcf_file, vcf_idx, annot_vcf_file, annot_vcf_idx, xcl_vcf_file, xcl_vcf_idx])}" | \
      sed 's/^.*\///' | tr '\n' '\0' | xargs -0 rm
  >>>

  output {
    File ext_vcf_file = filebase + "." + ext_string + ".bcf"
    File ext_vcf_idx = filebase + "." + ext_string + ".bcf.csi"
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
