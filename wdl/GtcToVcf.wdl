version 1.0

workflow Gtc2Vcf {

  input {
    String out_prefix

    Array[File] gtcs

    File bpm_manifest
    File egt_cluster_file
    File csv_manifest
    
    File ref_fasta
    File ref_fasta_fai

    # bwa indices
    File ref_fasta_sa
    File ref_fasta_amb
    File ref_fasta_bwt
    File ref_fasta_ann
    File ref_fasta_pac
  
    String? gtc2vcf_mocha_docker 
    Int? gtc2vcf_disk_size 
    Int? gtc2vcf_preemptible_attempts
    Int? gtc2vcf_memory 

    String? align_flank_mocha_docker 
    Int? align_flank_disk_size 
    Int? align_flank_preemptible_attempts
    Int? align_flank_memory 

  }

  call AlignFlankingSequences {
    input:
      csv_manifest = csv_manifest,
      ref_fasta = ref_fasta,
      ref_fasta_fai = ref_fasta_fai,
      ref_fasta_sa = ref_fasta_sa,
      ref_fasta_amb = ref_fasta_amb,
      ref_fasta_bwt = ref_fasta_bwt,
      ref_fasta_ann = ref_fasta_ann,
      ref_fasta_pac = ref_fasta_pac,
      
      disk_size = align_flank_disk_size,
      preemptible_attempts = align_flank_preemptible_attempts,
      docker = align_flank_mocha_docker,
      memory = align_flank_memory
  }

  call Gtc2Vcf {
    input:
      gtcs = gtcs,
      bpm_manifest = bpm_manifest,
      egt_cluster_file = egt_cluster_file,
      csv_manifest = csv_manifest,
      flanking_sequence_alignment_bam = AlignFlankingSequences.flanking_sequence_alignment_bam,

      ref_fasta = ref_fasta,
      ref_fasta_fai = ref_fasta_fai,

      out_prefix = out_prefix,      

      disk_size = gtc2vcf_disk_size,
      preemptible_attempts = gtc2vcf_preemptible_attempts,
      docker = gtc2vcf_mocha_docker,
      memory = gtc2vcf_memory
  }

  output {
    File flanking_sequence_alignment_bam = AlignFlankingSequences.flanking_sequence_alignment_bam
    File vcf = Gtc2Vcf.bcf
    File vcf_idx = Gtc2Vcf.bcf_idx
    File sex_file = Gtc2Vcf.sex_file
  }
}

task AlignFlankingSequences {
  input {
    File csv_manifest

    File ref_fasta
    File ref_fasta_fai
    File ref_fasta_sa
    File ref_fasta_amb
    File ref_fasta_bwt
    File ref_fasta_ann
    File ref_fasta_pac

    Int disk_size = 50
    String docker = "cwhelan/mocha:v1.0"
    Int threads = 1
    Int memory = 16
    Int preemptible_attempts = 3
  }

  String filebase = basename(ref_fasta)

  String outfile = "~{filebase}.bam"

  output {
    File flanking_sequence_alignment_bam = "~{outfile}"
  }

  command <<<

    set -euo pipefail

    bcftools +gtc2vcf --csv ~{csv_manifest} --fasta-flank | \
      bwa mem -M ~{ref_fasta} - | \
      samtools view -bS -o ~{outfile}

  >>>

  runtime {
    memory: memory + " GiB"
    disks: "local-disk " + disk_size + " HDD"
    cpu: threads
    docker: docker
    preemptible_attempts: preemptible_attempts
  }
}

task Gtc2Vcf {
  input {
    Array[File] gtcs

    File bpm_manifest
    File csv_manifest
    File egt_cluster_file
    File flanking_sequence_alignment_bam

    File ref_fasta
    File ref_fasta_fai

    String out_prefix

    Int disk_size = 50
    String docker = "cwhelan/mocha:v1.0"
    Int threads = 1
    Int memory = 16
    Int preemptible_attempts = 3
  }

  String gtc_dir = "~{out_prefix}_gtcs"

  String bcf_filename = "${out_prefix}.bcf"

  output {
    File bcf = "~{bcf_filename}"
    File bcf_idx = "${out_prefix}.bcf.csi"
    File sex_file = "${out_prefix}.sex"
  }

  command <<<

    set -euo pipefail

    while read gtc_file
    do
      ln -s $gtc_file 
    done < ~{write_lines(gtcs)}

    bcftools +gtc2vcf \
      --no-version -Ou \
      -b ~{bpm_manifest} \
      -c ~{csv_manifest} \
      -e ~{egt_cluster_file} \
      --sam ~{flanking_sequence_alignment_bam} \
      -g . \
      -f ~{ref_fasta} \
      -x ~{out_prefix}.sex | \
      bcftools sort -Ou -T ./bcftools-sort.XXXXXX | \
      bcftools norm --no-version -Ob -o ~{out_prefix}.bcf -c x -f ~{ref_fasta} && \
      bcftools index -f ~{bcf_filename}

  >>>
  runtime {
    memory: memory + " GiB"
    disks: "local-disk " + disk_size + " HDD"
    cpu: threads
    docker: docker
    preemptible_attempts: preemptible_attempts
  }
}
