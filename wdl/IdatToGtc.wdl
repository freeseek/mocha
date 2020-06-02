version 1.0

workflow IdatToGtc {

  input {
    String sample
    File red_idat
    File green_idat
    File bpm_manifest
    File egt_cluster_file
  
    String? mocha_docker 
    Int? disk_size 
    Int? memory 
    Int? preemptible_attempts
  }

  call GenCall {
    input:
      red_idat = red_idat,
      green_idat = green_idat,
      bpm_manifest = bpm_manifest,
      egt_cluster_file = egt_cluster_file,
      out_prefix = sample,

      disk_size = disk_size,
      docker = mocha_docker,
      memory = memory,
      preemptible_attempts = preemptible_attempts

  }


  output {
    File gtc = GenCall.gtc
  }
}

task GenCall {
  input {
    File red_idat
    File green_idat 
    File egt_cluster_file
    File bpm_manifest

    String out_prefix

    Int disk_size = 10
    String docker = "cwhelan/mocha:v1.0"
    Int memory = 6
    Int preemptible_attempts = 3
  }

  String idats_dir = "~{out_prefix}_idats"
  String filebase = basename(red_idat, "_Red.idat")
  String outfile = "~{filebase}.gtc"

  output {
    File gtc = "~{outfile}"
  }

  command <<<

    set -euo pipefail

    mkdir idats
    mv ~{green_idat} idats
    mv ~{red_idat} idats
    

    iaap-cli gencall \
      ~{bpm_manifest} \
      ~{egt_cluster_file} \
      . \
      -f idats -g  

  >>>
  runtime {
    memory: memory + " GiB"
    disks: "local-disk " + disk_size + " HDD"
    docker: docker
    preemptible_attempts: preemptible_attempts
    maxRetries: 1
  }
}

