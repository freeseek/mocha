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
  File? cyto_file
}

workflow shift {
  input {
    String sample_set_id
    String region
    String as_id = "AS"
    String ext_string = "as"

    String? ref_path
    File? cyto_file

    String? data_path
    File batch_tsv_file # batch_id path chr1_imp_vcf chr1_imp_vcf_index chr2_imp_vcf chr2_imp_vcf_index ...
    File? samples_file
    Boolean drop_genotypes = true
    Boolean phred_score = false
    Boolean plot = true
    String basic_bash_docker = "ubuntu:latest"
    String bcftools_docker = "us.gcr.io/mccarroll-mocha/bcftools:1.11-20210315"
    String mocha_plot_docker = "us.gcr.io/mccarroll-mocha/mocha_plot:1.11-20210315"
  }

  String ref_path_with_sep = if defined(ref_path) then select_first([ref_path]) + (if sub(select_first([ref_path]), "/$", "") != select_first([ref_path]) then "" else "/") else ""
  Reference ref = object {
    cyto_file: if defined(ref_path) || defined(cyto_file) then ref_path_with_sep + select_first([cyto_file, "cytoBand.txt.gz"]) else None
  }
  # call the relevant chromosome
  String chr_string = sub(sub(region, ":.*$", ""), "^chr", "")

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

  # scatter target genotypes
  scatter (idx in range(n_batches)) {
    call vcf_summary {
      input:
        vcf_file = data_paths_with_sep[idx] + batch_tbl[("chr" + chr_string + "_imp_vcf")][idx],
        vcf_idx = data_paths_with_sep[idx] + batch_tbl[("chr" + chr_string + "_imp_vcf_index")][idx],
        region = region,
        samples_file = samples_file,
        as_id = as_id,
        drop_genotypes = drop_genotypes,
        docker = bcftools_docker
    }
  }

  call vcf_merge {
    input:
      vcf_files = vcf_summary.as_vcf_file,
      vcf_idxs = vcf_summary.as_vcf_idx,
      as_id = as_id,
      ext_string = ext_string,
      filebase = sample_set_id + "." + sub(region, "[:-]", "_"),
      phred_score = phred_score,
      docker = bcftools_docker
  }

  if (plot) {
    call shift_plot {
      input:
        vcf_file = vcf_merge.as_vcf_file,
        vcf_idx = vcf_merge.as_vcf_idx,
        region = region,
        cyto_file = ref.cyto_file,
        filebase = sample_set_id + "." + sub(region, "[:-]", "_"),
        docker = mocha_plot_docker
    }
  }

  output {
    File vcf_file = vcf_merge.as_vcf_file
    File vcf_idx = vcf_merge.as_vcf_idx
    File? png_file = shift_plot.png_file
  }
}

task vcf_summary {
  input {
    File vcf_file
    File vcf_idx
    String region
    File? samples_file
    String as_id
    String ext_string = "sites"
    Boolean drop_genotypes = true

    String docker
    Int cpu = 1
    Int? disk_size_override
    Float memory = 3.5
    Int preemptible = 1
    Int maxRetries = 0
  }

  Float vcf_size = size(vcf_file, "GiB")
  Int disk_size = select_first([disk_size_override, ceil(10.0 + 2.0 * vcf_size)])
  String filebase = basename(basename(vcf_file, ".bcf"), ".vcf.gz")

  command <<<
    set -euo pipefail
    mv "~{vcf_file}" .
    mv "~{vcf_idx}" .
    ~{if defined(samples_file) then "mv \"" + select_first([samples_file]) + "\" ." else ""}
    bcftools +mochatools \
      --no-version \
      --output-type b \
      --output "~{filebase}.~{ext_string}.bcf" \
      --regions "~{region}" \
      "~{basename(vcf_file)}" \
      -- --summary ~{as_id} \
      ~{if defined(samples_file) then "--samples-file \"" + basename(select_first([samples_file])) + "\" \\\n" +
      "  --force-samples" else ""} \
      ~{if drop_genotypes then "--drop-genotypes" else ""}
    bcftools index --force "~{filebase}.~{ext_string}.bcf"
    rm "~{basename(vcf_file)}"
    rm "~{basename(vcf_idx)}"
    ~{if defined(samples_file) then "rm \"" + basename(select_first([samples_file])) + "\"" else ""}
  >>>

  output {
    File as_vcf_file = filebase + "." + ext_string + ".bcf"
    File as_vcf_idx = filebase + "." + ext_string + ".bcf.csi"
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

task vcf_merge {
  input {
    Array[File] vcf_files
    Array[File] vcf_idxs
    String as_id
    String ext_string
    String filebase
    Boolean phred_score = false

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
    vcf_idxs=~{write_lines(vcf_idxs)}
    cat $vcf_files $vcf_idxs | tr '\n' '\0' | xargs -0 mv -t .
    sed -i 's/^.*\///' $vcf_files $vcf_idxs
    bcftools merge \
      --no-version \
      --output-type u \
      --info-rules ~{as_id}:sum \
      --file-list $vcf_files \
      --merge none | \
    bcftools +mochatools \
    --no-version \
    --output-type b \
    --output "~{filebase}.~{ext_string}.bcf" \
    -- --test ~{as_id} \
    ~{if phred_score then "--phred" else ""}
    bcftools index --force "~{filebase}.~{ext_string}.bcf"
    cat $vcf_files $vcf_idxs | tr '\n' '\0' | xargs -0 rm
  >>>

  output {
    File as_vcf_file = filebase + "." + ext_string + ".bcf"
    File as_vcf_idx = filebase + "." + ext_string + ".bcf.csi"
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

task shift_plot {
  input {
    File vcf_file
    File vcf_idx
    String region
    File? cyto_file
    String filebase

    String docker
    Int cpu = 1
    Int? disk_size_override
    Float memory = 3.5
    Int preemptible = 1
    Int maxRetries = 0
  }

  Float vcf_size = size(vcf_file, "GiB")
  Int disk_size = select_first([disk_size_override, ceil(10.0 + vcf_size)])

  command <<<
    set -euo pipefail
    mv "~{vcf_file}" .
    mv "~{vcf_idx}" .
    ~{if defined(cyto_file) then "mv \"" + select_first([cyto_file]) + "\" ." else ""}
    shift_plot.R \
      ~{if defined(cyto_file) then "--cytoband \"" + basename(select_first([cyto_file])) + "\"" else ""} \
      --vcf "~{basename(vcf_file)}" \
      --region ~{region} \
      --png "~{filebase}.png"
    rm "~{basename(vcf_file)}"
    rm "~{basename(vcf_idx)}"
    ~{if defined(cyto_file) then "rm \"" + basename(select_first([cyto_file])) + "\"" else ""}
  >>>

  output {
    File png_file = filebase + '.png'
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
