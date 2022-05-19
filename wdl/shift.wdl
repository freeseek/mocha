version development

## Copyright (c) 2021-2022 Giulio Genovese
##
## Version 2022-05-18
##
## Contact Giulio Genovese <giulio.genovese@gmail.com>
##
## This WDL workflow runs allelic shift imbalance analysis in a given region
##
## Cromwell version support
## - Successfully tested on v79
##
## Distributed under terms of the MIT License

struct Reference {
  File? cyto_file
  String? chr_prefix
  Array[Int] len
  Array[Int] pcen
  Array[Int] qcen
}

workflow shift {
  input {
    String sample_set_id
    File? keep_samples_file
    File? remove_samples_file
    File pheno_tsv_file
    String as_id = "AS"
    String ext_string = "as"

    String ref_name = "GRCh38"
    String? ref_path
    String? chr_prefix
    File? cyto_file

    File impute_tsv_file # batch_id path chr1_imp_vcf chr1_imp_vcf_index chr2_imp_vcf chr2_imp_vcf_index ...
    String? impute_data_path
    Boolean fisher_exact = true
    Boolean drop_genotypes = true
    Boolean phred_score = true
    Boolean plot = true
    String basic_bash_docker = "debian:stable-slim"
    String docker_repository = "us.gcr.io/mccarroll-mocha"
    String bcftools_docker = "bcftools:1.15.1-20220518"
    String r_mocha_docker = "r_mocha:1.15.1-20220518"
  }

  String docker_repository_with_sep = docker_repository + if docker_repository != "" && docker_repository == sub(docker_repository, "/$", "") then "/" else ""

  String ref_path_with_sep = select_first([ref_path, ""]) + if defined(ref_path) && select_first([ref_path]) == sub(select_first([ref_path]), "/$", "") then "/" else ""
  Reference ref = object {
    cyto_file: if defined(ref_path) || defined(cyto_file) then ref_path_with_sep + select_first([cyto_file, "cytoBand.txt.gz"]) else None,
    chr_prefix: if ref_name == "GRCh38" then "chr" else if ref_name == "GRCh37" then "" else chr_prefix,
    len: if ref_name == "GRCh38" then
      [248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973, 145138636, 138394717, 133797422, 135086622, 133275309, 114364328, 107043718, 101991189, 90338345, 83257441, 80373285, 58617616, 64444167, 46709983, 50818468, 156040895]
    else if ref_name == "GRCh37" then
      [249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 141213431, 135534747, 135006516, 133851895, 115169878, 107349540, 102531392, 90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 51304566, 155270560]
    else None,
    pcen: if ref_name == "GRCh38" then
      [122026459, 92188145, 90772458, 49712061, 46485900, 58553888, 58169653, 44033744, 43389635, 39686682, 51078348, 34769407, 16000000, 16000000, 17083673, 36311158, 22813679, 15460899, 24498980, 26436232, 10864560, 12954788, 58605579]
    else if ref_name == "GRCh37" then
      [121535434, 92326171, 90504854, 49660117, 46405641, 58830166, 58054331, 43838887, 47367679, 39254935, 51644205, 34856694, 16000000, 16000000, 17000000, 35335801, 22263006, 15460898, 24681782, 26369569, 11288129, 13000000, 58632012]
    else None,
    qcen: if ref_name == "GRCh38" then
      [124932724, 94090557, 93655574, 51743951, 50059807, 59829934, 61528020, 45877265, 45518558, 41593521, 54425074, 37185252, 18051248, 18173523, 19725254, 38265669, 26616164, 20861206, 27190874, 30038348, 12915808, 15054318, 62412542]
    else if ref_name == "GRCh37" then
      [124535434, 95326171, 93504854, 52660117, 49405641, 61830166, 61054331, 46838887, 50367679, 42254935, 54644205, 37856694, 19000000, 19000000, 20000000, 38335801, 25263006, 18460898, 27681782, 29369569, 14288129, 16000000, 61632012]
    else None
  }

  # read table with batches information (scatter could be avoided if there was a tail() function)
  Array[Array[String]] impute_tsv = read_tsv(impute_tsv_file)
  Int n_batches = length(impute_tsv)-1
  scatter (idx in range(n_batches)) { Array[String] impute_tsv_rows = impute_tsv[(idx+1)] }
  Map[String, Array[String]] impute_tbl = as_map(zip(impute_tsv[0], transpose(impute_tsv_rows)))

  # compute data paths for each batch, if available (scatter could be avoided if there was a contains_key() function)
  scatter (key in keys(impute_tbl)) { Boolean? is_key_equal_path = if key == "path" then true else None }
  scatter (idx in range(n_batches)) {
    String impute_data_paths = select_first([impute_data_path, if length(select_all(is_key_equal_path))>0 then impute_tbl["path"][idx] else ""])
    String impute_data_paths_with_sep = impute_data_paths + (if impute_data_paths == "" || sub(impute_data_paths, "/$", "") != impute_data_paths then "" else "/")
  }

  call lst_header { input: pheno_tsv_file = pheno_tsv_file, docker = basic_bash_docker }
  # compute phenotype regions to test
  scatter (idx in range(length(lst_header.phenos))) {
    String region_name = sub(lst_header.phenos[idx], "_.*$", "")
    String chr_string = sub(sub(region_name, "[pq]*$", ""), "Y", "X")
    if (sub(chr_string, "[0-9X]+", "") == "") {
      Int pheno_idx = idx
      Int chr_idx = sub(chr_string, "^X$", "23")
      String arm = sub(region_name, "^[0-9XY]*", "")
      String pheno_regions = ref.chr_prefix + chr_string + if arm == "p" then ":1-" + ref.pcen[(chr_idx - 1)] else if arm == "q" then ":" + ref.qcen[(chr_idx - 1)] + "-" + ref.len[(chr_idx - 1)] else ""
    }
  }

  # scatter target genotypes
  scatter (p in cross(range(n_batches), select_all(pheno_idx))) {
    String cross_idx = p.right
    call vcf_summary {
      input:
        vcf_file = impute_data_paths_with_sep[p.left] + impute_tbl[("chr" + chr_string[p.right] + "_imp_vcf")][p.left],
        vcf_idx = impute_data_paths_with_sep[p.left] + impute_tbl[("chr" + chr_string[p.right] + "_imp_vcf_index")][p.left],
        pheno_name = lst_header.phenos[p.right],
        region = select_first([pheno_regions[p.right]]),
        keep_samples_file = keep_samples_file,
        remove_samples_file = remove_samples_file,
        pheno_tsv_file = pheno_tsv_file,
        fisher_exact = fisher_exact,
        as_id = as_id,
        ext_string = ext_string,
        drop_genotypes = drop_genotypes,
        docker = docker_repository_with_sep + bcftools_docker
    }
  }

  scatter (idx in select_all(pheno_idx)) {
    Map[Int, Array[File]] idx2vcf_files = collect_by_key(zip(cross_idx, vcf_summary.as_vcf_file))
    Map[Int, Array[File]] idx2vcf_idxs = collect_by_key(zip(cross_idx, vcf_summary.as_vcf_idx))
    call vcf_merge {
      input:
        vcf_files = idx2vcf_files[idx],
        vcf_idxs = idx2vcf_idxs[idx],
        as_id = as_id,
        fisher_exact = fisher_exact,
        ext_string = ext_string,
        filebase = sample_set_id + "." + lst_header.phenos[idx],
        phred_score = phred_score,
        docker = docker_repository_with_sep + bcftools_docker
    }

    if (plot) {
      call assoc_plot {
        input:
          vcf_file = vcf_merge.as_vcf_file,
          vcf_idx = vcf_merge.as_vcf_idx,
          region = select_first([pheno_regions[idx]]),
          cyto_file = ref.cyto_file,
          filebase = sample_set_id + "." + lst_header.phenos[idx],
          docker = docker_repository_with_sep + r_mocha_docker
      }
    }
  }

  output {
    Array[File] vcf_files = vcf_merge.as_vcf_file
    Array[File] vcf_idxs = vcf_merge.as_vcf_idx
    Array[File]? png_files = if plot then select_all(assoc_plot.png_file) else None
  }

  meta {
    author: "Giulio Genovese"
    email: "giulio.genovese@gmail.com"
    description: "See the [MoChA](https://github.com/freeseek/mocha) website for more information"
  }
}

task lst_header {
  input {
    File pheno_tsv_file
    String space_character = '_'

    String docker
    Int cpu = 1
    Int disk_size = 10
    Float memory = 3.5
    Int preemptible = 1
    Int maxRetries = 0
  }

  command <<<
    set -euo pipefail
    mv "~{pheno_tsv_file}" .
    head -n1 "~{basename(pheno_tsv_file)}" | cut -f2- | tr '\t ' '\n~{space_character}'
    rm "~{basename(pheno_tsv_file)}"
  >>>

  output {
    Array[String] phenos = read_lines(stdout())
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

# the command requires BCFtools 1.14 due to bug https://github.com/samtools/bcftools/issues/1566
task vcf_summary {
  input {
    File vcf_file
    File vcf_idx
    String pheno_name
    String region
    File? keep_samples_file
    File? remove_samples_file
    File pheno_tsv_file
    Boolean fisher_exact
    String as_id
    String ext_string
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
  String filebase = basename(basename(basename(basename(vcf_file, ".bcf"), ".vcf.gz"), ".as"), ".chr" + sub(sub(region, ":.*$", ""), "^chr", "")) + '.' + pheno_name

  command <<<
    set -euo pipefail
    echo "~{sep="\n" select_all([vcf_file, vcf_idx, keep_samples_file, remove_samples_file, pheno_tsv_file])}" | \
      tr '\n' '\0' | xargs -0 mv -t .
    ~{if fisher_exact then
      "awk -F\"\\t\" 'NR==1 {for (i=1; i<=NF; i++) f[$i] = i}\n" +
      "  NR>1 && $(f[\"" + pheno_name + "\"])==0 {print $(f[\"sample_id\"])}' \"" + basename(pheno_tsv_file) + "\" > \"" + filebase + ".controls.lines\""
    else ""}
    awk -F"\t" 'NR==1 {for (i=1; i<=NF; i++) f[$i] = i}
      NR>1 && $(f["~{pheno_name}"])==1 {print $(f["sample_id"])}' "~{basename(pheno_tsv_file)}" > "~{filebase}.cases.lines"
    bcftools annotate \
      --no-version \
      --output-type u \
      --regions "~{region}" \
      --remove ID,QUAL,FILTER,INFO,^FMT/GT,FMT/~{as_id} \
      "~{basename(vcf_file)}" | \
    ~{if defined(keep_samples_file) then
      "bcftools view \\\n" +
      "  --no-version \\\n" +
      "  --output-type u \\\n" +
      "  --samples-file \"" + basename(select_first([keep_samples_file])) + "\" \\\n" +
      "  --force-samples |"
      else ""} \
    ~{if defined(remove_samples_file) then
      "bcftools view \\\n" +
      "  --no-version \\\n" +
      "  --output-type u \\\n" +
      "  --samples-file \"^" + basename(select_first([remove_samples_file])) + "\" \\\n" +
      "  --force-samples |"
    else ""} \
    ~{if fisher_exact then
      "bcftools +contrast \\\n" +
      "  --output-type u \\\n" +
      "  --annots NASSOC \\\n" +
      "  --control-samples \"" + filebase + ".controls.lines\" \\\n" +
      "  --case-samples \"" + filebase + ".cases.lines\" \\\n" +
      "  --force-samples |"
    else ""} \
    bcftools +mochatools \
      --no-version \
      --output-type b \
      -- --summary ~{as_id} \
      --samples-file "~{filebase}.cases.lines" \
      --force-samples \
      ~{if drop_genotypes then "--drop-genotypes" else ""} | \
    tee "~{filebase}.~{ext_string}.bcf" | \
    bcftools index --force --output "~{filebase}.~{ext_string}.bcf.csi"
    rm~{if fisher_exact then " \"" + filebase + ".controls.lines\"" else ""} "~{filebase}.cases.lines"
    echo "~{sep="\n" select_all([vcf_file, vcf_idx, keep_samples_file, remove_samples_file, pheno_tsv_file])}" | \
      sed 's/^.*\///' | tr '\n' '\0' | xargs -0 rm
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
    Boolean fisher_exact
    String as_id
    String ext_string
    String filebase
    Boolean phred_score = true

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
      --info-rules ~{if fisher_exact then "NASSOC:sum," else ""}~{as_id}:sum \
      --file-list $vcf_files \
      --merge none | \
    ~{if fisher_exact then
      "bcftools +mochatools \\\n" +
      "  --no-version \\\n" +
      "  --output-type u \\\n" +
      "  -- --test NASSOC \\\n" +
      (if phred_score then "  --phred" else "") + " |"
    else ""} \
    bcftools +mochatools \
      --no-version \
      --output-type b \
      -- --test ~{as_id} \
      ~{if phred_score then "--phred" else ""} | \
    tee "~{filebase}.~{ext_string}.bcf" | \
    bcftools index --force --output "~{filebase}.~{ext_string}.bcf.csi"
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

task assoc_plot {
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
    assoc_plot.R \
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
