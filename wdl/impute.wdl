version development

## Copyright (c) 2021-2022 Giulio Genovese
##
## Version 2022-05-18
##
## Contact Giulio Genovese <giulio.genovese@gmail.com>
##
## This WDL workflow runs impute5 or beagle5 on a set of VCFs
##
## Cromwell version support
## - Successfully tested on v79
##
## Distributed under terms of the MIT License

struct Reference {
  File fasta_fai
  Int n_chrs
  String mhc_reg
  File genetic_map_file
  String panel_pfx
  String panel_sfx
  String panel_idx
  Int n_panel_smpls
}

workflow impute {
  input {
    String sample_set_id
    String mode = "pgt" # pgt imp
    String target = "imp" # imp ext
    Float max_win_size_cm = 10.0
    Float overlap_size_cm = 2.0
    String format_id = "AS"
    String ext_string = "as"
    Array[String]? target_chrs

    String ref_name = "GRCh38"
    String? ref_path
    String? ref_fasta_fai
    Int? ref_n_chrs
    String? mhc_reg
    String? genetic_map_file
    String? panel_pfx
    String? panel_sfx
    String? panel_idx
    Int? n_panel_smpls

    File mocha_tsv_file # batch_id n_smpls path vcf vcf_index pgt_vcf pgt_vcf_index chr1_imp_vcf_index chr2_imp_vcf chr2_imp_vcf_index ...
    String? mocha_data_path
    String? impute_data_path
    File? remove_samples_file
    Boolean convert_panel = true
    Boolean beagle = false
    Boolean out_ds = true
    Boolean out_gp = false
    Boolean out_ap = false
    String? impute_extra_args
    String basic_bash_docker = "debian:stable-slim"
    String pandas_docker = "amancevice/pandas:slim"
    String docker_repository = "us.gcr.io/mccarroll-mocha"
    String bcftools_docker = "bcftools:1.15.1-20220518"
    String impute5_docker = "impute5:1.15.1-20220518"
    String beagle5_docker = "beagle5:1.15.1-20220518"
  }

  String docker_repository_with_sep = docker_repository + if docker_repository != "" && docker_repository == sub(docker_repository, "/$", "") then "/" else ""

  String ref_path_with_sep = select_first([ref_path, ""]) + if defined(ref_path) && select_first([ref_path]) == sub(select_first([ref_path]), "/$", "") then "/" else ""
  Reference ref = object {
    fasta_fai: ref_path_with_sep + select_first([ref_fasta_fai, if ref_name == "GRCh38" then "GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai" else if ref_name == "GRCh37" then "human_g1k_v37.fasta.fai" else None]),
    n_chrs: select_first([ref_n_chrs, 23]),
    mhc_reg: select_first([mhc_reg, if ref_name == "GRCh38" then "chr6:27518932-33480487" else if ref_name == "GRCh37" then "6:27486711-33448264" else None]),
    genetic_map_file: ref_path_with_sep + select_first([genetic_map_file, if ref_name == "GRCh38" then "genetic_map_hg38_withX.txt.gz" else if ref_name == "GRCh37" then "genetic_map_hg19_withX.txt.gz" else None]),
    panel_pfx: ref_path_with_sep + select_first([panel_pfx, if ref_name == "GRCh38" then "CCDG_14151_B01_GRM_WGS_2020-08-05_" else if ref_name == "GRCh37" then "ALL.chr" else None]),
    panel_sfx: select_first([panel_sfx, if ref_name == "GRCh38" then ".filtered.phased.bcf" else if ref_name == "GRCh37" then ".phase3_integrated.20130502.genotypes.bcf" else None]),
    panel_idx: select_first([panel_idx, ".csi"]),
    n_panel_smpls: select_first([n_panel_smpls, if ref_name == "GRCh38" then 3202 else if ref_name == "GRCh37" then 2504 else None]),
  }

  # read table with batches information (scatter could be avoided if there was a tail() function)
  Array[Array[String]] mocha_tsv = read_tsv(mocha_tsv_file)
  Int n_mocha_batches = length(mocha_tsv)-1
  scatter (idx in range(n_mocha_batches)) { Array[String] mocha_tsv_rows = mocha_tsv[(idx+1)] }
  Map[String, Array[String]] mocha_tbl = as_map(zip(mocha_tsv[0], transpose(mocha_tsv_rows)))

  # compute data paths for each batch, if available (scatter could be avoided if there was a contains_key() function)
  scatter (key in keys(mocha_tbl)) { Boolean? is_key_equal_path = if key == "path" then true else None }
  scatter (idx in range(n_mocha_batches)) {
    String mocha_data_paths = select_first([mocha_data_path, if length(select_all(is_key_equal_path))>0 then mocha_tbl["path"][idx] else ""])
    String mocha_data_paths_with_sep = mocha_data_paths + (if mocha_data_paths == "" || sub(mocha_data_paths, "/$", "") != mocha_data_paths then "" else "/")
    String impute_data_paths = select_first([impute_data_path, if length(select_all(is_key_equal_path))>0 then mocha_tbl["path"][idx] else "", ""])
    String impute_data_paths_with_sep = impute_data_paths + (if impute_data_paths == "" || sub(impute_data_paths, "/$", "") != impute_data_paths then "" else "/")
  }

  Array[Array[String]] ref_fasta_fai_tbl = transpose(read_tsv(ref.fasta_fai))
  scatter (idx in range(ref.n_chrs)) { String ref_chrs = ref_fasta_fai_tbl[0][idx] }
  Array[String] chrs = select_first([target_chrs, ref_chrs])
  Int n_chrs = length(chrs)
  scatter (idx in range(n_chrs)) { String chr_strings = sub(chrs[idx], "^chr", "") }

  if (mode == "pgt") {
    Map[String, Int] chr2len = as_map(zip(ref_fasta_fai_tbl[0], ref_fasta_fai_tbl[1]))
    scatter (chr in chrs) { Int lens = chr2len[chr] }
    call ref_scatter {
      input:
        chrs = chrs,
        lens = lens,
        genetic_map_file = ref.genetic_map_file,
        max_win_size_cm = max_win_size_cm,
        overlap_size_cm = overlap_size_cm,
        docker = pandas_docker
    }
    Array[Array[String]] intervals_tbl = transpose(read_tsv(ref_scatter.intervals_tsv))
    # this is a trick to table how many intervals you will use for each chromosome
    Map[String, Array[Int]] chr_map = collect_by_key(zip(intervals_tbl[0], range(length(intervals_tbl[0]))))

    # scatter reference panel
    scatter (idx in range(n_chrs)) {
      call get_nrecords as get_panel_markers { input: vcf_idx = ref.panel_pfx + chrs[idx] + ref.panel_sfx + ref.panel_idx, docker = docker_repository_with_sep + bcftools_docker }
      if (length(chr_map[(chrs[idx])]) > 1) {
        call vcf_scatter as panel_scatter {
          input:
            vcf_file = ref.panel_pfx + chrs[idx] + ref.panel_sfx,
            intervals_tsv = ref_scatter.intervals_tsv,
            chr = chrs[idx],
            docker = docker_repository_with_sep + bcftools_docker
        }
      }
      Array[File] panel_scatter_vcf_files = select_first([panel_scatter.vcf_files, [ref.panel_pfx + chrs[idx] + ref.panel_sfx]])
      Array[File] panel_scatter_vcf_idxs = select_first([panel_scatter.vcf_idxs, [ref.panel_pfx + chrs[idx] + ref.panel_sfx + ref.panel_idx]])
      Array[Int] n_panel_scatter_markers = select_first([panel_scatter.n_markers, [get_panel_markers.n]])
    }
    Array[File] panel_flatten_vcf_files = flatten(panel_scatter_vcf_files)
    Array[File] panel_flatten_vcf_idxs = flatten(panel_scatter_vcf_idxs)
    Array[Int] n_panel_flatten_markers = flatten(n_panel_scatter_markers)

    # convert reference panel
    scatter (idx in range(length(intervals_tbl[0]))) {
      call init_panel {
        input:
          vcf_file = panel_flatten_vcf_files[idx],
          vcf_idx = panel_flatten_vcf_idxs[idx],
          chr = intervals_tbl[0][idx],
          convert = convert_panel,
          beagle = beagle,
          docker = if beagle then docker_repository_with_sep + beagle5_docker else docker_repository_with_sep + impute5_docker
      }
    }

    # scatter target genotypes
    scatter (idx in range(n_mocha_batches)) {
      call vcf_scatter {
        input:
          vcf_file = mocha_data_paths_with_sep[idx] + mocha_tbl["pgt_vcf"][idx],
          intervals_tsv = ref_scatter.intervals_tsv,
          remove_samples_file = remove_samples_file,
          docker = docker_repository_with_sep + bcftools_docker
      }
    }

    # impute genotypes
    Map[String, Int] chr2int = as_map(zip(chrs, range(length(chrs))))
    String mhc_chr = sub(ref.mhc_reg, ":.*$", "")
    Int mhc_beg = sub(sub(ref.mhc_reg, "-.*$", ""), "^.*:", "")
    Int mhc_end = sub(ref.mhc_reg, "^.*-", "")
    scatter (p in cross(range(n_mocha_batches), range(length(intervals_tbl[0])))) {
      if (vcf_scatter.n_markers[p.left][p.right] > 0) {
        Int cross_idx = p.left * n_chrs + chr2int[(intervals_tbl[0][p.right])]
        Int buffer_beg = intervals_tbl[1][p.right] # cast string to integer
        Int buffer_end = intervals_tbl[2][p.right] # cast string to integer
        Int beg = intervals_tbl[3][p.right] # cast string to integer
        Int end = intervals_tbl[4][p.right] # cast string to integer
        call vcf_impute {
          input:
            n_smpls = vcf_scatter.n_smpls[p.left],
            n_markers = vcf_scatter.n_markers[p.left][p.right],
            pgt_file = vcf_scatter.vcf_files[p.left][p.right],
            genetic_map_file = ref.genetic_map_file,
            n_panel_smpls = ref.n_panel_smpls,
            n_panel_markers = n_panel_flatten_markers[p.right],
            panel_genotypes_file = select_first([init_panel.panel_genotypes_file[p.right], panel_flatten_vcf_files[p.right]]),
            panel_genotypes_idx = select_first([init_panel.panel_genotypes_idx[p.right], panel_flatten_vcf_idxs[p.right]]),
            panel_sites_file = init_panel.panel_sites_file[p.right],
            panel_sites_idx = init_panel.panel_sites_idx[p.right],
            ref_fasta_fai = ref.fasta_fai,
            chr = intervals_tbl[0][p.right],
            mhc = intervals_tbl[0][p.right] == mhc_chr && buffer_beg <= mhc_end && buffer_end > mhc_beg,
            region = intervals_tbl[0][p.right] + ":" + (1 + beg) + "-" + end,
            buffer_region = intervals_tbl[0][p.right] + ":" + (1 + buffer_beg) + "-" + buffer_end,
            beagle = beagle,
            out_ds= out_ds,
            out_gp= out_gp,
            out_ap= out_ap,
            impute_extra_args = impute_extra_args,
            docker = if beagle then docker_repository_with_sep + beagle5_docker else docker_repository_with_sep + impute5_docker
        }
      }
    }

    Map[Int, Array[File]] idx2vcf_files = collect_by_key(zip(select_all(cross_idx), select_all(vcf_impute.imp_vcf_file)))
    Map[Int, Array[File]] idx2vcf_idxs = collect_by_key(zip(select_all(cross_idx), select_all(vcf_impute.imp_vcf_idx)))
    scatter (p in cross(range(n_mocha_batches), range(n_chrs))) {
      if (length(idx2vcf_files[(p.left * n_chrs + p.right)]) > 1) {
        call vcf_concat {
          input:
            vcf_files = idx2vcf_files[(p.left * n_chrs + p.right)],
            filebase = basename(basename(basename(mocha_tbl["pgt_vcf"][p.left], ".bcf"), ".vcf.gz"), ".pgt") + ".chr" + chr_strings[p.right] + ".imp",
            docker = docker_repository_with_sep + bcftools_docker
        }
      }
      File chr_imp_vcf_files = select_first([vcf_concat.vcf_file, idx2vcf_files[(p.left * n_chrs + p.right)][0]])
      File chr_imp_vcf_idxs = select_first([vcf_concat.vcf_idx, idx2vcf_idxs[(p.left * n_chrs + p.right)][0]])
    }
  }

  Array[Int] n_smpls = if mode == "pgt" then select_first([vcf_scatter.n_smpls]) else mocha_tbl["n_smpls"]

  if (target == "ext") {
    scatter (p in cross(range(n_mocha_batches), range(n_chrs))) {
      File vcf_file = if mode == "pgt" then select_first([chr_imp_vcf_files])[(p.left * n_chrs + p.right)]
                      else impute_data_paths_with_sep[p.left] + mocha_tbl[("chr" +  chr_strings[p.right] + "_imp_vcf")][p.left]
      File vcf_idx = if mode == "pgt" then select_first([chr_imp_vcf_idxs])[(p.left * n_chrs + p.right)]
                     else impute_data_paths_with_sep[p.left] + mocha_tbl[("chr" +  chr_strings[p.right] + "_imp_vcf_index")][p.left]
      call get_nrecords { input: vcf_idx = vcf_idx, docker = docker_repository_with_sep + bcftools_docker }
      call vcf_extend {
        input:
          n_smpls = n_smpls[p.left],
          n_markers = get_nrecords.n,
          vcf_file = vcf_file,
          vcf_idx = vcf_idx,
          annot_vcf_file = mocha_data_paths_with_sep[p.left] + mocha_tbl["vcf"][p.left],
          annot_vcf_idx = mocha_data_paths_with_sep[p.left] + mocha_tbl["vcf_index"][p.left],
          format_id = format_id,
          ext_string = ext_string,
          chr_string = chr_strings[p.right],
          docker = docker_repository_with_sep + bcftools_docker
      }
    }
  }

  # generate a table summarizing the main output files and serialize the table to disk
  # vcf_files and vcf_idxs are defined in the output section
  scatter (p in cross(range(n_mocha_batches), range(n_chrs))) {
    Pair[String, String] output_vcfs = ("chr" + chr_strings[p.right] + "_imp_vcf", basename(vcf_files[(p.left * n_chrs + p.right)]))
    Pair[String, String] output_idxs = ("chr" + chr_strings[p.right] + "_imp_vcf_index", basename(vcf_idxs[(p.left * n_chrs + p.right)]))
  }
  Map[String, Array[String]] output_map = as_map(flatten([[("batch_id", mocha_tbl["batch_id"]),
    ("n_smpls", n_smpls)], as_pairs(collect_by_key(output_vcfs)), as_pairs(collect_by_key(output_idxs))]))
  # cannot use output_keys = keys(output_map) because of unresolved Cromwell bug
  # https://github.com/broadinstitute/cromwell/issues/5559
  scatter (idx in range(n_chrs)) { Array[String] keys = ["chr" + chr_strings[idx] + "_imp_vcf", "chr" + chr_strings[idx] + "_imp_vcf_index"] }
  Array[String] output_keys = flatten([["batch_id", "n_smpls"], flatten(keys)])
  scatter (key in output_keys) { Array[String] output_tsv_cols = output_map[key] }
  # this is run as a separate task rather than using write_tsv() as Cromwell can break the WDL specification
  # https://support.terra.bio/hc/en-us/community/posts/360071465631-write-lines-write-map-write-tsv-write-json-fail-when-run-in-a-workflow-rather-than-in-a-task
  call write_tsv {
    input:
      tsv = flatten([[output_keys], transpose(output_tsv_cols)]),
      filebase = sample_set_id + ".impute",
      docker = basic_bash_docker
  }

  output {
    Array[File] vcf_files = select_first([vcf_extend.ext_vcf_file, chr_imp_vcf_files])
    Array[File] vcf_idxs = select_first([vcf_extend.ext_vcf_idx, chr_imp_vcf_idxs])
    Array[File]? logs = if mode == "pgt" then select_all(select_first([vcf_impute.log])) else None
    File impute_tsv_file = write_tsv.file
  }

  meta {
    author: "Giulio Genovese"
    email: "giulio.genovese@gmail.com"
    description: "See the [MoChA](https://github.com/freeseek/mocha) website for more information"
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
    set -euo pipefail
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

task get_nrecords {
  input {
    File vcf_idx
    Boolean binary_vcf = true

    String docker
    Int cpu = 1
    Int disk_size = 10
    Float memory = 3.5
    Int preemptible = 1
    Int maxRetries = 0
  }

  String ext = if binary_vcf then "bcf" else "vcf.gz"

  command <<<
    set -euo pipefail
    mv "~{vcf_idx}" .
    bcftools index --nrecords "~{basename(vcf_idx)}"
    rm "~{basename(vcf_idx)}"
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
        cm_mids = (win_size - 0.5 * ~{overlap_size_cm}) * np.arange(1, n_win - 1)
        pos_begs = np.concatenate(([0], 0 + np.interp(cm_begs, df_group['CM'], df_group['POS'], period = np.inf).astype(int)))
        pos_ends = np.concatenate((np.interp(cm_ends, df_group['CM'], df_group['POS'], period = np.inf).astype(int), [chr2len[fai_chr]]))
        pos_begs2 = np.concatenate(([0], 0 + np.interp(cm_begs + 0.5 * ~{overlap_size_cm}, df_group['CM'], df_group['POS'], period = np.inf).astype(int)))
        pos_ends2 = np.concatenate((np.interp(cm_ends - 0.5 * ~{overlap_size_cm}, df_group['CM'], df_group['POS'], period = np.inf).astype(int), [chr2len[fai_chr]]))
        df_out[fai_chr] = pd.DataFrame.from_dict({'CHR': fai_chr, 'BEG': pos_begs, 'END': pos_ends, 'BEG2': pos_begs2, 'END2': pos_ends2})
    df = pd.concat([df_out[fai_chr] for fai_chr in chr2len.keys()])
    df[['CHR', 'BEG', 'END', 'BEG2', 'END2']].to_csv('ref_scatter.tsv', sep='\t', header = False, index = False)
    CODE
    rm chr2len.tsv
    rm "~{basename(genetic_map_file)}"
  >>>

  output {
    File intervals_tsv = "ref_scatter.tsv"
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
    File intervals_tsv # zero-based intervals
    Int clevel = 2
    File? remove_samples_file
    String? chr

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
    echo "~{sep="\n" select_all([vcf_file, intervals_tsv, remove_samples_file])}" | \
      tr '\n' '\0' | xargs -0 mv -t .
    ~{if defined(chr) then
      "mv \"" + basename(intervals_tsv) + "\" \"" + basename(intervals_tsv, ".tsv") + ".all.tsv\"\n" +
      "awk -v chr=\"" + chr + "\" '$1==chr' \"" + basename(intervals_tsv, ".tsv") + ".all.tsv\" > \"" + basename(intervals_tsv) + "\"\n" +
      "rm \"" + basename(intervals_tsv, ".tsv") + ".all.tsv\""
      else ""}
    bcftools query --list-samples "~{basename(vcf_file)}" | wc -l > n_smpls.int
    awk -F"\t" '{print $1":"1+$2"-"$3"\t"NR-1}' "~{basename(intervals_tsv)}" > regions.lines
    bcftools annotate \
      --no-version \
      --output-type u \
      --remove ID,QUAL,FILTER,INFO,^FMT/GT \
      ~{if cpu > 1 then "--threads " + (cpu - 1) else ""} \
      "~{basename(vcf_file)}" | \
    ~{if defined(remove_samples_file) then
      "bcftools view \\\n" +
      "  --no-version \\\n" +
      "  --output-type u \\\n" +
      "  --samples-file ^\"" + basename(select_first([remove_samples_file])) + "\" \\\n" +
      "  --force-samples | \\\n"
      else ""}bcftools norm \
      --no-version \
      --output-type u \
      --rm-dup exact \
      ~{if cpu > 1 then "--threads " + (cpu - 1) else ""} | \
    bcftools +scatter \
      --no-version \
      --output-type b~{clevel} \
      --output vcfs \
      ~{if cpu > 1 then "--threads " + (cpu - 1) else ""} \
      --scatter-file regions.lines \
      --prefix "~{filebase}."
    while read reg i; do
      bcftools index --force "vcfs/~{filebase}.$i.bcf"
      bcftools index --nrecords "vcfs/~{filebase}.$i.bcf.csi"
    done < regions.lines > n_markers.lines
    cut -f2 regions.lines | sed 's/^/vcfs\/~{filebase}./;s/$/.bcf/'
    cut -f2 regions.lines | sed 's/^/vcfs\/~{filebase}./;s/$/.bcf.csi/' > vcf_idxs.lines
    echo "~{sep="\n" select_all([vcf_file, intervals_tsv, remove_samples_file])}" | \
      sed 's/^.*\///' | tr '\n' '\0' | xargs -0 rm
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

task init_panel {
  input {
    File vcf_file
    File? vcf_idx
    String chr
    Boolean convert = true
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
        --output-type b \
        --remove ID,QUAL,FILTER,INFO,^FMT/GT | \
      tee "~{filebase}.sites.bcf" | \
      bcftools index --force --output "~{filebase}.sites.bcf.csi"
    ~{if beagle then
        (if chr != "X" && chr != "chrX" then
          "bcftools view --no-version \"" + basename(vcf_file) + "\" | \\\n"
        else
          "bcftools +fixploidy --no-version \"" + basename(vcf_file) + "\" | \\\n" +
          "  sed 's/0\\/0/0|0/g;s/1\\/1/1|1/g' | \\\n") +
        "  java -jar /usr/share/beagle/bref3.jar > \"" + filebase +  ".bref3\""
      else if convert then
        "imp5Converter \\\n" +
        "  --h \"" + basename(vcf_file) + "\" \\\n" +
        "  --r " + chr + " \\\n" +
        "  --o \"" + filebase + ".imp5\""
      else ""}
    rm "~{basename(vcf_file)}"
    ~{if defined(vcf_idx) then "rm \"" + basename(select_first([vcf_idx])) + "\"" else ""}
  >>>

  output {
    File? panel_genotypes_file = if beagle then filebase + ".bref3" else if convert then filebase + ".imp5" else None
    File? panel_genotypes_idx = if !beagle && convert then filebase + ".imp5.idx" else None
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

# hack https://github.com/samtools/bcftools/issues/1425 is employed in the end to fix the header
# the command requires BCFtools 1.14 due to bug https://github.com/samtools/bcftools/issues/1497
task vcf_impute {
  input {
    Int n_smpls
    Int n_markers
    File pgt_file
    File genetic_map_file
    Int n_panel_smpls
    Int n_panel_markers
    File panel_genotypes_file
    File? panel_genotypes_idx
    File panel_sites_file
    File panel_sites_idx
    File? ref_fasta_fai
    String region
    String buffer_region
    String chr
    Boolean mhc = false # requires additional memory if imputing the MHC region
    Boolean beagle = false
    Boolean out_ds = true
    Boolean out_gp = false
    Boolean out_ap = false
    Int max_ram_percentage = 90
    String? impute_extra_args

    String docker
    Int? cpu_override
    Int? disk_size_override
    Float? mult_override
    Float? memory_override
    Int preemptible = 1
    Int maxRetries = 0
  }

  Float pgt_size = size(pgt_file, "GiB")
  Float panel_size = size(panel_genotypes_file, "GiB")
  Int disk_size = select_first([disk_size_override, ceil(10.0 + 3.0 * pgt_size + (1.0 + 2.0 * n_smpls / n_panel_smpls) * panel_size)])
  Float mult = select_first([mult_override, (if beagle then 15.0 else 10.0) * (if mhc then 2.0 else if chr == "X" || chr == "chrX" then 1.5 else 1.0)])
  Float memory = select_first([memory_override, 3.5 + mult * n_panel_markers * (n_smpls + n_panel_smpls) / 1024 / 1024 / 1024])
  Int cpu = select_first([cpu_override, if memory > 6.5 then 2 * ceil(memory / 13) else 1])
  String filebase = basename(basename(pgt_file, ".bcf"), ".vcf.gz")

  command <<<
    set -euo pipefail
    echo "~{sep="\n" select_all([pgt_file, genetic_map_file, panel_genotypes_file, panel_genotypes_idx, panel_sites_file, panel_sites_idx])}" | \
      tr '\n' '\0' | xargs -0 mv -t .
    ~{if defined(ref_fasta_fai) then "mv \"" + select_first([ref_fasta_fai]) + "\" ." else ""}
    bcftools index --force "~{basename(pgt_file)}"
    n_markers=$(bcftools isec \
      --no-version \
      --output-type u \
      --nfiles 2 \
      --write 1 \
      "~{basename(pgt_file)}" \
      "~{basename(panel_sites_file)}" | \
      bcftools query --format "\n" | wc -l)
    if [ $n_markers == 0 ]; then
      cp "~{basename(pgt_file)}" "~{filebase}.imp.bcf"
    else
      mkdir logs
    ~{if beagle then
      "  bcftools norm \\\n" +
      "    --no-version \\\n" +
      "    --rm-dup exact \\\n" +
      "    --output-type z \\\n" +
      "    --output \"" + filebase + ".vcf.gz\" \\\n" +
      "    \"" + basename(pgt_file) + "\"\n" +
      "  chr=" + chr +"; zcat \"" + basename(genetic_map_file) + "\" | \\\n" +
      "    sed 's/^23/X/' | awk -v chr=$chr '$1==chr || \"chr\"$1==chr {print chr,\".\",$4,$2}' > genetic_map.txt\n" +
      "  java -XX:MaxRAMPercentage=" + max_ram_percentage + " \\\n" +
      "    -jar /usr/share/beagle/beagle.jar \\\n" +
      "    gt=\"" + filebase + ".vcf.gz\" \\\n" +
      "    ref=\"" + basename(panel_genotypes_file) + "\" \\\n" +
      "    out=\"" + filebase + ".imp\" \\\n" +
      "    map=genetic_map.txt \\\n" +
      "    chrom=" + buffer_region + " \\\n" +
      (if cpu > 1 then "    nthreads=" + cpu + " \\\n" else "") +
      (if out_ap then "    ap=true \\\n" else "") +
      (if out_gp then "    gp=true \\\n" else "") +
      (if defined(impute_extra_args) then "    " + impute_extra_args + " \\\n" else "") +
      "    1>&2\n" +
      "  bcftools index --force --tbi \"" + filebase + ".imp.vcf.gz\"\n" +
      "  bcftools view \\\n" +
      "    --no-version \\\n" +
      "    --output-type b \\\n" +
      "    --regions " + region + " \\\n" +
      "    \"" + filebase + ".imp.vcf.gz\" | \\\n" +
      "  tee \"" + filebase + ".imp.bcf\" | \\\n" +
      "  bcftools index --force --output \"" + filebase + ".imp.bcf.csi\"\n" +
      "  rm \"" + filebase + ".vcf.gz\" \"" + filebase + ".imp.vcf.gz\" \"" + filebase + ".imp.vcf.gz.tbi\"\n" +
      "  mv \"" + filebase + ".imp.log\" \"logs/" + filebase + ".imp.log\""
    else
      "  chr=" + chr +"; zcat \"" + basename(genetic_map_file) + "\" | \\\n" +
      "  sed 's/^23/X/' | awk -v chr=$chr -v OFS=\"\\t\" 'BEGIN {print \"pos\",\"chr\",\"cM\"}\n" +
      "    $1==chr || \"chr\"$1==chr {print $2,chr,$4}' > genetic_map.txt\n" +
      "  impute5 \\\n" +
      "    --h \"" + basename(panel_genotypes_file) + "\" \\\n" +
      "    --m genetic_map.txt \\\n" +
      "    --g \"" + basename(pgt_file) + "\" \\\n" +
      "    --r " + region + " \\\n" +
      "    --buffer-region " + buffer_region + " \\\n" +
      "    --l \"logs/" + filebase + ".imp.log\" \\\n" +
      (if out_gp then "    --out-gp-field \\\n" else "") +
      (if out_ap then "    --out-ap-field \\\n" else "") +
      "    --o \"" + filebase + ".imp.bcf\" \\\n" +
      (if cpu > 1 then "    --threads " + cpu + " \\\n" else "") +
      (if defined(impute_extra_args) then "  " + impute_extra_args + " \\\n" else "") +
      "    1>&2"}
    fi
    ~{if (!out_ds) then
      "mv \"" + filebase + ".imp.bcf\" \"" + filebase + ".tmp.bcf\"\n" +
      "bcftools annotate \\\n" +
      "  --no-version \\\n" +
      "  --output-type b \\\n" +
      "  --remove FMT/DS \\\n" +
      "  \"" + filebase + ".tmp.bcf\" | \\\n" +
      "tee \"" + filebase + ".imp.bcf\" | \\\n" +
      "bcftools index --force --output \"" + filebase + ".imp.bcf.csi\"\n" +
      "rm \"" + filebase + ".tmp.bcf\""
      else ""}
    ~{if defined(ref_fasta_fai) then
      "(echo -en \"##fileformat=VCFv4.2\\n#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\tFORMAT\\t\"\n" +
      "bcftools query -l \"" + filebase + ".imp.bcf\" | tr '\\n' '\\t' | sed 's/\\t$/\\n/') > tmp.vcf\n" +
      "bcftools reheader --fai \"" + basename(select_first([ref_fasta_fai])) + "\" --output fai.vcf --temp-prefix ./bcftools.XXXXXX tmp.vcf\n" +
      "mv \"" + filebase + ".imp.bcf\" \"" + filebase + ".tmp.bcf\"\n" +
      "bcftools concat \\\n" +
      "  --no-version \\\n" +
      "  --output-type b \\\n" +
      "  fai.vcf \"" + filebase + ".tmp.bcf\" | \\\n" +
      "tee \"" + filebase + ".imp.bcf\" | \\\n" +
      "bcftools index --force --output \"" + filebase + ".imp.bcf.csi\"\n" +
      "rm \"" + filebase + ".tmp.bcf\" fai.vcf tmp.vcf \"" + basename(select_first([ref_fasta_fai])) + "\""
      else ""}
    ~{if (out_ds && !defined(ref_fasta_fai)) then
      "bcftools index --force \"" + filebase + ".imp.bcf\""
      else ""}
    rm "~{basename(pgt_file)}.csi" genetic_map.txt
    echo "~{sep="\n" select_all([pgt_file, genetic_map_file, panel_genotypes_file, panel_genotypes_idx, panel_sites_file, panel_sites_idx])}" | \
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
    String filebase

    String docker
    Int cpu = 1
    Int? disk_size_override
    Float memory = 3.5
    Int preemptible = 1
    Int maxRetries = 0
  }

  Float vcf_size = size(vcf_files, "GiB")
  Int disk_size = select_first([disk_size_override, ceil(10.0 + 2.0 * vcf_size)])

  command <<<
    set -euo pipefail
    vcf_files=~{write_lines(vcf_files)}
    cat $vcf_files | tr '\n' '\0' | xargs -0 mv -t .
    sed -i 's/^.*\///' $vcf_files
    bcftools concat \
      --no-version \
      --output-type b \
      --file-list $vcf_files \
      ~{if cpu > 1 then "--threads " + (cpu - 1) else ""} | \
    tee "~{filebase}.bcf" | \
    bcftools index --force --output "~{filebase}.bcf.csi"
    cat $vcf_files | tr '\n' '\0' | xargs -0 rm
    bcftools query --list-samples "~{filebase}.bcf" | wc -l
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

task vcf_extend {
  input {
    Int n_smpls
    Int n_markers
    File vcf_file
    File vcf_idx
    File annot_vcf_file
    File annot_vcf_idx
    String format_id
    String ext_string
    String chr_string
    Int dist = 500000

    String docker
    Int? cpu_override
    Int? disk_size_override
    Float? memory_override
    Int preemptible = 1
    Int maxRetries = 0

    Float mult = 0.2
  }

  Float vcf_size = size(vcf_file, "GiB")
  Float annot_vcf_size = size(annot_vcf_file, "GiB")
  Int disk_size = select_first([disk_size_override, ceil(10.0 + 2.0 * (vcf_size + annot_vcf_size))])
  Float memory = select_first([memory_override, 3.5 + mult * n_smpls * n_markers / 1024 / 1024 / 1024])
  Int cpu = select_first([cpu_override, if memory > 6.5 then 2 * ceil(memory / 13) else 1])
  String filebase = basename(basename(basename(basename(basename(vcf_file, ".bcf"), ".vcf.gz"), ".imp"), "." + ext_string), ".chr" + chr_string) + ".chr" + chr_string
  String input_vcf_file = sub(sub(basename(vcf_file), "\\." + ext_string + "\\.bcf$", ".bcf"), "\\." + ext_string + "\\.vcf\\.gz$", ".vcf.gz")
  String input_vcf_idx = sub(sub(sub(basename(vcf_idx), "\\." + ext_string + "\\.bcf\\.csi$", ".bcf.csi"), "\\." + ext_string + "\\.vcf\\.gz\\.tbi$", ".vcf.gz.tbi"), "\\." + ext_string + "\\.vcf\\.gz\\.csi$", ".vcf.gz.csi")

  command <<<
    set -euo pipefail
    mv "~{vcf_file}" ~{if basename(vcf_file) != input_vcf_file then "\"" + input_vcf_file + "\"" else "."}
    mv "~{vcf_idx}" ~{if basename(vcf_file) != input_vcf_file then "\"" + input_vcf_idx + "\"" else "."}
    mv "~{annot_vcf_file}" .
    mv "~{annot_vcf_idx}" .
    bcftools annotate \
      --no-version \
      --output-type u \
      --annotations "~{basename(annot_vcf_file)}" \
      --columns "FMT/~{format_id}" \
      --remove ID,QUAL,FILTER,INFO,^FMT/GT \
      "~{input_vcf_file}" | \
    bcftools +extendFMT \
      --no-version \
      --output-type b \
      --format "~{format_id}" \
      --phase \
      --dist ~{dist} | \
    tee "~{filebase}.~{ext_string}.pgt.bcf" | \
    bcftools index --force --output "~{filebase}.~{ext_string}.pgt.bcf.csi"
    bcftools annotate \
      --no-version \
      --output-type b \
      --annotations "~{filebase}.~{ext_string}.pgt.bcf" \
      --columns "FMT/~{format_id}" \
      ~{if cpu > 1 then "--threads " + (cpu - 1) else ""} \
      "~{input_vcf_file}" | \
    tee "~{filebase}.~{ext_string}.bcf" | \
    bcftools index --force --output "~{filebase}.~{ext_string}.bcf.csi"
    rm "~{filebase}.~{ext_string}.pgt.bcf" "~{filebase}.~{ext_string}.pgt.bcf.csi"
    rm "~{input_vcf_file}"
    rm "~{input_vcf_idx}"
    rm "~{basename(annot_vcf_file)}"
    rm "~{basename(annot_vcf_idx)}"
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
