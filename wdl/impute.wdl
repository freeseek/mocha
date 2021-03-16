version development

## Copyright (c) 2021 Giulio Genovese
##
## Version 2021-03-12
##
## Contact Giulio Genovese <giulio.genovese@gmail.com>
##
## This WDL workflow runs impute5 or beagle5 on a set of VCFs
##
## Cromwell version support
## - Successfully tested on v58
##
## Distributed under terms of the MIT License

struct Reference {
  File fasta_fai
  Int n_chrs
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
    String? genetic_map_file
    String? panel_pfx
    String? panel_sfx
    String? panel_idx
    Int? n_panel_smpls

    String? data_path
    File batch_tsv_file # batch_id n_smpls path vcf vcf_index pgt_vcf pgt_vcf_index chr1_imp_vcf chr1_imp_vcf_index ...
    Boolean beagle = false
    Boolean out_ds = true
    Boolean out_gp = false
    Boolean out_ap = false
    String? impute_extra_args
    String basic_bash_docker = "ubuntu:latest"
    String pandas_docker = "amancevice/pandas:slim"
    String bcftools_docker = "us.gcr.io/mccarroll-mocha/bcftools:1.11-20210315"
    String impute5_docker = "us.gcr.io/mccarroll-mocha/impute5:1.11-20210315"
    String beagle5_docker = "us.gcr.io/mccarroll-mocha/beagle5:1.11-20210315"
  }

  String ref_path_with_sep = if defined(ref_path) then select_first([ref_path]) + (if sub(select_first([ref_path]), "/$", "") != select_first([ref_path]) then "" else "/") else ""
  Reference ref = object {
    fasta_fai: ref_path_with_sep + select_first([ref_fasta_fai, if ref_name == "GRCh38" then "GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai" else if ref_name == "GRCh37" then "human_g1k_v37.fasta.fai" else None]),
    n_chrs: select_first([ref_n_chrs, 23]),
    genetic_map_file: ref_path_with_sep + select_first([genetic_map_file, if ref_name == "GRCh38" then "genetic_map_hg38_withX.txt.gz" else if ref_name == "GRCh37" then "genetic_map_hg19_withX.txt.gz" else None]),
    panel_pfx: ref_path_with_sep + select_first([panel_pfx, if ref_name == "GRCh38" then "CCDG_14151_B01_GRM_WGS_2020-08-05_" else if ref_name == "GRCh37" then "ALL.chr" else None]),
    panel_sfx: select_first([panel_sfx, if ref_name == "GRCh38" then ".filtered.phased.bcf" else if ref_name == "GRCh37" then ".phase3_integrated.20130502.genotypes.bcf" else None]),
    panel_idx: select_first([panel_idx, ".csi"]),
    n_panel_smpls: select_first([n_panel_smpls, if ref_name == "GRCh38" then 3202 else if ref_name == "GRCh37" then 2504 else None]),
  }

  # read table with batches information (scatter could be avoided if there was a tail() function)
  Array[Array[String]] batch_tsv = read_tsv(batch_tsv_file)
  Int n_batches = length(batch_tsv)-1
  scatter (idx in range(n_batches)) { Array[String] batch_tsv_rows = batch_tsv[(idx+1)] }
  Map[String, Array[String]] batch_tbl = as_map(zip(batch_tsv[0], transpose(batch_tsv_rows)))

  # compute data paths for each batch, if available (scatter could be avoided if there was a contains_key() function)
  scatter (key in keys(batch_tbl)) { Boolean? is_key_equal_path = if key == "path" then true else None }
  scatter (idx in range(n_batches)) {
    String data_paths = select_first([data_path, if length(select_all(is_key_equal_path))>0 then batch_tbl["path"][idx] else ""])
    String data_paths_with_sep = data_paths + (if data_paths == "" || sub(data_paths, "/$", "") != data_paths then "" else "/")
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
    Array[Array[String]] intervals_tbl = transpose(read_tsv(ref_scatter.intervals_bed))
    # this is a trick to table how many intervals you will use for each chromosome
    Map[String, Array[Int]] chr_map = collect_by_key(zip(intervals_tbl[0], range(length(intervals_tbl[0]))))

    # scatter reference panel
    scatter (idx in range(n_chrs)) {
      call get_markers as get_panel_markers { input: vcf_idx = ref.panel_pfx + chrs[idx] + ref.panel_sfx + ref.panel_idx, docker = bcftools_docker }
      if (length(chr_map[(chrs[idx])]) > 1) {
        call vcf_scatter as panel_scatter {
          input:
            vcf_file = ref.panel_pfx + chrs[idx] + ref.panel_sfx,
            intervals_bed = ref_scatter.intervals_bed,
            chr = chrs[idx],
            docker = bcftools_docker
        }
      }
      Array[File] panel_scatter_vcf_files = select_first([panel_scatter.vcf_files, [ref.panel_pfx + chrs[idx] + ref.panel_sfx]])
      Array[File] panel_scatter_vcf_idxs = select_first([panel_scatter.vcf_idxs, [ref.panel_pfx + chrs[idx] + ref.panel_sfx + ref.panel_idx]])
      Array[Int] n_panel_scatter_markers = select_first([panel_scatter.n_markers, [get_panel_markers.n]])
    }
    Array[File] panel_flatten_vcf_files = flatten(panel_scatter_vcf_files)
    Array[File] panel_flatten_vcf_idxs = flatten(panel_scatter_vcf_idxs)

    # convert reference panel
    scatter (idx in range(length(intervals_tbl[0]))) {
      call init_panel {
        input:
          vcf_file = panel_flatten_vcf_files[idx],
          vcf_idx = panel_flatten_vcf_idxs[idx],
          chr = intervals_tbl[0][idx],
          beagle = beagle,
          docker = if beagle then beagle5_docker else impute5_docker
      }
    }

    # scatter target genotypes
    scatter (idx in range(n_batches)) {
      call vcf_scatter {
        input:
          vcf_file = data_paths_with_sep[idx] + batch_tbl["pgt_vcf"][idx],
          intervals_bed = ref_scatter.intervals_bed,
          docker = bcftools_docker
      }
    }

    # impute genotypes
    Array[Int] n_panel_flatten_markers = flatten(n_panel_scatter_markers)
    Map[String, Int] chr2int = as_map(zip(chrs, range(length(chrs))))
    scatter (p in cross(range(n_batches), range(length(intervals_tbl[0])))) {
      if (vcf_scatter.n_markers[p.left][p.right] > 0) {
        Int cross_idx = p.left * n_chrs + chr2int[(intervals_tbl[0][p.right])]
        call vcf_impute {
          input:
            n_smpls = vcf_scatter.n_smpls[p.left],
            n_markers = vcf_scatter.n_markers[p.left][p.right],
            pgt_file = vcf_scatter.vcf_files[p.left][p.right],
            genetic_map_file = ref.genetic_map_file,
            n_panel_smpls = ref.n_panel_smpls,
            n_panel_markers = n_panel_flatten_markers[p.right],
            panel_genotypes_file = init_panel.panel_genotypes_file[p.right],
            panel_genotypes_idx = init_panel.panel_genotypes_idx[p.right],
            panel_sites_file = init_panel.panel_sites_file[p.right],
            panel_sites_idx = init_panel.panel_sites_idx[p.right],
            ref_fasta_fai = ref.fasta_fai,
            chr = intervals_tbl[0][p.right],
            region = intervals_tbl[0][p.right] + ":" + intervals_tbl[1][p.right] + "-" + intervals_tbl[2][p.right],
            beagle = beagle,
            out_ds= out_ds,
            out_gp= out_gp,
            out_ap= out_ap,
            impute_extra_args = impute_extra_args,
            docker = if beagle then beagle5_docker else impute5_docker
        }
      }
    }

    Map[Int, Array[File]] idx2vcf_files = collect_by_key(zip(select_all(cross_idx), select_all(vcf_impute.imp_vcf_file)))
    Map[Int, Array[File]] idx2vcf_idxs = collect_by_key(zip(select_all(cross_idx), select_all(vcf_impute.imp_vcf_idx)))
    scatter (p in cross(range(n_batches), range(n_chrs))) {
      if (length(idx2vcf_files[(p.left * n_chrs + p.right)]) > 1) {
        call vcf_ligate {
          input:
            n_smpls = vcf_scatter.n_smpls[p.left],
            n_markers = get_panel_markers.n[p.right],
            vcf_files = idx2vcf_files[(p.left * n_chrs + p.right)],
            filebase = basename(basename(basename(batch_tbl["pgt_vcf"][p.left], ".bcf"), ".vcf.gz"), ".pgt") + ".chr" + chr_strings[p.right] + ".imp",
            docker = bcftools_docker
        }
      }
      File chr_imp_vcf_files = select_first([vcf_ligate.vcf_file, idx2vcf_files[(p.left * n_chrs + p.right)][0]])
      File chr_imp_vcf_idxs = select_first([vcf_ligate.vcf_idx, idx2vcf_idxs[(p.left * n_chrs + p.right)][0]])
    }
  }

  Array[Int] n_smpls = if mode == "pgt" then select_first([vcf_scatter.n_smpls]) else batch_tbl["n_smpls"]

  if (target == "ext") {
    scatter (p in cross(range(n_batches), range(n_chrs))) {
      File vcf_file = if mode == "pgt" then select_first([chr_imp_vcf_files])[(p.left * n_chrs + p.right)]
                      else data_paths_with_sep[p.left] + batch_tbl[("chr" +  chr_strings[p.right] + "_imp_vcf")][p.left]
      File vcf_idx = if mode == "pgt" then select_first([chr_imp_vcf_idxs])[(p.left * n_chrs + p.right)]
                     else data_paths_with_sep[p.left] + batch_tbl[("chr" +  chr_strings[p.right] + "_imp_vcf_index")][p.left]
      call get_markers { input: vcf_idx = vcf_idx, docker = bcftools_docker }
      call vcf_extend {
        input:
          n_smpls = n_smpls[p.left],
          n_markers = get_markers.n,
          vcf_file = vcf_file,
          vcf_idx = vcf_idx,
          annot_vcf_file = data_paths_with_sep[p.left] + batch_tbl["vcf"][p.left],
          annot_vcf_idx = data_paths_with_sep[p.left] + batch_tbl["vcf_index"][p.left],
          format_id = format_id,
          ext_string = ext_string,
          docker = bcftools_docker
      }
    }
  }

  # generate a table summarizing the main output files and serialize the table to disk
  # vcf_files and vcf_idxs are defined in the output section
  scatter (p in cross(range(n_batches), range(n_chrs))) {
    Pair[String, String] output_vcfs = ("chr" + chr_strings[p.right] + "_imp_vcf", basename(vcf_files[(p.left * n_chrs + p.right)]))
    Pair[String, String] output_idxs = ("chr" + chr_strings[p.right] + "_imp_vcf_index", basename(vcf_idxs[(p.left * n_chrs + p.right)]))
  }
  Map[String, Array[String]] output_map = as_map(flatten([[("batch_id", batch_tbl["batch_id"]),
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
      filebase = sample_set_id + ".output",
      docker = basic_bash_docker
  }

  output {
    Array[File] vcf_files = select_first([vcf_extend.ext_vcf_file, chr_imp_vcf_files])
    Array[File] vcf_idxs = select_first([vcf_extend.ext_vcf_idx, chr_imp_vcf_idxs])
    Array[File]? logs = if mode == "pgt" then select_all(select_first([vcf_impute.log])) else None
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

# uses hack from https://github.com/samtools/bcftools/issues/1418
task get_markers {
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
    echo -e "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" | \
      bcftools view --no-version -O~{if binary_vcf then "b" else "z"} -o in.~{ext}
    bcftools index --nrecords "in.~{ext}##idx##~{basename(vcf_idx)}"
    rm in.~{ext} "~{basename(vcf_idx)}"
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
    Boolean uncompressed = false

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

task init_panel {
  input {
    File vcf_file
    File? vcf_idx
    String chr
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
    ~{if beagle then
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
        "  --o \"" + filebase + ".imp5\""}
    rm "~{basename(vcf_file)}"
    ~{if defined(vcf_idx) then "rm \"" + basename(select_first([vcf_idx])) + "\"" else ""}
  >>>

  output {
    File panel_genotypes_file = filebase + (if beagle then ".bref3" else ".imp5")
    File? panel_genotypes_idx = if beagle then None else filebase + ".imp5.idx"
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
    String chr
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
  Float mult = select_first([mult_override, (if beagle then 15.0 else 8.0) * (if chr == "X" || chr == "chrX" then 1.5 else 1.0)])
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
      "    -jar /usr/bin/beagle.jar \\\n" +
      "    gt=\"" + filebase + ".vcf.gz\" \\\n" +
      "    ref=\"" + basename(panel_genotypes_file) + "\" \\\n" +
      "    out=\"" + filebase + ".imp\" \\\n" +
      "    map=genetic_map.txt \\\n" +
      "    chrom=" + region + " \\\n" +
      (if cpu > 1 then "    nthreads=" + cpu + " \\\n" else "") +
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
      "    --g \"" + basename(pgt_file) + "\" \\\n" +
      "    --r " + region + " \\\n" +
      "    --l \"logs/" + filebase + ".imp.log\" \\\n" +
      (if out_gp then "    --out-gp-field \\\n" else "") +
      (if out_ap then "    --out-ap-field \\\n" else "") +
      "    --o \"" + filebase + ".imp.bcf\" \\\n" +
      (if cpu > 1 then "  --threads " + cpu + " \\\n" else "") +
      (if defined(impute_extra_args) then "  " + impute_extra_args + " \\\n" else "") +
      "    1>&2"}
    fi
    ~{if (!out_ds) then
      "mv \"" + filebase + ".imp.bcf\" \"" + filebase + ".tmp.bcf\"\n" +
      "bcftools annotate \\\n" +
      "  --no-version \\\n" +
      "  --output-type b \\\n" +
      "  --output \"" + filebase + ".imp.bcf\" \\\n" +
      "  --remove FMT/DS \\\n" +
      "  \"" + filebase + ".tmp.bcf\"\n" +
      "rm \"" + filebase + ".tmp.bcf\""
      else ""}
    rm "~{basename(pgt_file)}.csi" genetic_map.txt
    ~{if defined(ref_fasta_fai) then
      "(echo -en \"##fileformat=VCFv4.2\\n#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\tFORMAT\\t\"\n" +
      "bcftools query -l \"" + filebase + ".imp.bcf\" | tr '\\n' '\\t' | sed 's/\\t$/\\n/') > tmp.vcf\n" +
      "bcftools reheader --fai \"" + basename(select_first([ref_fasta_fai])) + "\" --output fai.vcf tmp.vcf\n" +
      "mv \"" + filebase + ".imp.bcf\" \"" + filebase + ".tmp.bcf\"\n" +
      "bcftools concat \\\n" +
      "  --no-version \\\n" +
      "  --output-type b \\\n" +
      "  --output \"" + filebase + ".imp.bcf\" \\\n" +
      "  fai.vcf \"" + filebase + ".tmp.bcf\"\n" +
      "rm \"" + filebase + ".tmp.bcf\" fai.vcf tmp.vcf \"" + basename(select_first([ref_fasta_fai])) + "\""
      else ""}
    bcftools index --force "~{filebase}.imp.bcf"
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

task vcf_ligate {
  input {
    Int n_smpls
    Int n_markers
    Array[File]+ vcf_files
    String filebase
    Boolean uncompressed = false

    String docker
    Int? cpu_override
    Int? disk_size_override
    Float? memory_override
    Int preemptible = 1
    Int maxRetries = 0

    Float mult = 4.0 # how many bytes per genotype are required
  }

  Float vcf_size = size(vcf_files, "GiB")
  Int disk_size = select_first([disk_size_override, ceil(10.0 + 2.0 * vcf_size)])
  Float memory = select_first([memory_override, 3.5 + mult * n_smpls * n_markers / length(vcf_files) / 1024 / 1024 / 1024])
  Int cpu = select_first([cpu_override, if memory > 6.5 then 2 * ceil(memory / 13) else 1])

  command <<<
    set -euo pipefail
    vcf_files=~{write_lines(vcf_files)}
    cat $vcf_files | tr '\n' '\0' | xargs -0 mv -t .
    sed -i 's/^.*\///' $vcf_files
    cat $vcf_files | tr '\n' '\0' | xargs -0 -n 1 bcftools index --force
    bcftools concat \
      --no-version \
      --output-type ~{if uncompressed then "u" else "b"} \
      --output "~{filebase}.unligated.bcf" \
      --allow-overlaps \
      --rm-dups none \
      --file-list $vcf_files \
      ~{if cpu > 1 then "--threads " + (cpu - 1) else ""}
    bcftools index --force "~{filebase}.unligated.bcf"
    cat $vcf_files | tr '\n' '\0' | \
    xargs -0 -i bcftools annotate \
      --no-version \
      --output-type ~{if uncompressed then "u" else "b"} \
      --output "{}.pgt.bcf" \
      --remove ID,QUAL,FILTER,INFO,^FMT/GT \
      ~{if cpu > 1 then "--threads " + (cpu - 1) else ""} \
      "{}"
    cat $vcf_files | tr '\n' '\0' | xargs -0 rm
    cat $vcf_files | sed 's/$/.csi/' | tr '\n' '\0' | xargs -0 rm
    sed -i 's/$/.pgt.bcf/' $vcf_files
    cat $vcf_files | tr '\n' '\0' | xargs -0 -n 1 bcftools index --force
    bcftools concat \
      --no-version \
      --output-type ~{if uncompressed then "u" else "b"} \
      --compact-PS \
      --file-list $vcf_files \
      --ligate \
      ~{if cpu > 1 then "--threads " + (cpu - 1) else ""} \
      --output "~{filebase}.pgt.bcf"
    bcftools index --force "~{filebase}.pgt.bcf"
    cat $vcf_files | tr '\n' '\0' | xargs -0 rm
    cat $vcf_files | sed 's/$/.csi/' | tr '\n' '\0' | xargs -0 rm
    bcftools annotate \
      --no-version \
      --output-type ~{if uncompressed then "u" else "b"} \
      --output "~{filebase}.bcf" \
      --annotations "~{filebase}.pgt.bcf" \
      --columns FMT/GT \
      ~{if cpu > 1 then "--threads " + (cpu - 1) else ""} \
      "~{filebase}.unligated.bcf"
    rm "~{filebase}.unligated.bcf"
    rm "~{filebase}.unligated.bcf.csi"
    rm "~{filebase}.pgt.bcf"
    rm "~{filebase}.pgt.bcf.csi"
    bcftools index --force "~{filebase}.bcf"
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

# due to bug https://github.com/samtools/bcftools/issues/1414 I need two bcftools annotate invocations
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
    Int dist = 500000
    Boolean uncompressed = false

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
  String filebase = basename(basename(basename(vcf_file, ".bcf"), ".vcf.gz"), ".imp")

  command <<<
    set -euo pipefail
    mv "~{vcf_file}" .
    mv "~{vcf_idx}" .
    mv "~{annot_vcf_file}" .
    mv "~{annot_vcf_idx}" .
    bcftools annotate \
      --no-version \
      --output-type b \
      --output "~{filebase}.pgt.bcf" \
      --remove ID,QUAL,FILTER,INFO,^FMT/GT \
      "~{basename(vcf_file)}"
    bcftools index --force "~{filebase}.pgt.bcf"
    bcftools annotate \
      --no-version \
      --output-type u \
      --columns "FMT/~{format_id}" \
      --annotations "~{basename(annot_vcf_file)}" \
      "~{filebase}.pgt.bcf" | \
    bcftools +extendFMT \
      --no-version \
      --output-type ~{if uncompressed then "u" else "b"} \
      --output "~{filebase}.pgt.~{ext_string}.bcf" \
      --format "~{format_id}" \
      --phase \
      --dist ~{dist}
    bcftools index --force "~{filebase}.pgt.~{ext_string}.bcf"
    rm "~{filebase}.pgt.bcf" "~{filebase}.pgt.bcf.csi"
    bcftools annotate \
      --no-version \
      --output-type ~{if uncompressed then "u" else "b"} \
      --output "~{filebase}.~{ext_string}.bcf" \
      --annotations "~{filebase}.pgt.~{ext_string}.bcf" \
      --columns "FMT/~{format_id}" \
      ~{if cpu > 1 then "--threads " + (cpu - 1) else ""} \
      "~{basename(vcf_file)}"
    bcftools index --force "~{filebase}.~{ext_string}.bcf"
    rm "~{filebase}.pgt.~{ext_string}.bcf" "~{filebase}.pgt.~{ext_string}.bcf.csi"
    rm "~{basename(vcf_file)}"
    rm "~{basename(vcf_idx)}"
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
