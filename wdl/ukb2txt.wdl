version development

## Copyright (c) 2021 Giulio Genovese
##
## Version 2021-05-14
##
## Contact Giulio Genovese <giulio.genovese@gmail.com>
##
## Converts UK biobank genotype and intensity files back to Affymetrix format suitable for the MoChA WDL pipeline
##
## Cromwell version support
## - Successfully tested on v61
##
## Distributed under terms of the MIT License

workflow ukb2txt {
  input {
    String sqc_path
    String sqc_filename = "ukb_sqc_v2.txt"
    String cal_path
    String cal_filename_pfx = "ukb_cal_chr"
    String cal_filename_sfx = "_v2.bed"
    String int_path
    String int_filename_pfx = "ukb_int_chr"
    String int_filename_sfx = "_v2.bin"
    Array[String] chrs = [ "1",  "2",  "3",  "4",  "5",  "6",  "7",  "8",  "9", "10", "11", "12", "13",
                          "14", "15", "16", "17", "18", "19", "20", "21", "22",  "X",  "Y", "XY", "MT"]
    File? snp_posterior_tar
    File? snp_qc_file
    File? UKBL_zip
    File? UKBB_zip
    String docker = "ubuntu:latest"
  }

  File sqc_file = sqc_path + (if sub(sqc_path, "/$", "") != sqc_path then "" else "/") + sqc_filename

  call init {
    input:
      sqc_file = sqc_file,
      docker = docker
  }

  if (!defined(snp_posterior_tar) || !defined(snp_qc_file) || !defined(UKBL_zip) || !defined(UKBB_zip)) {
    call wget {
      input:
        docker = docker
    }
  }

  call csv {
    input:
      UKBL_zip = select_first([UKBL_zip, wget.UKBL_zip]),
      UKBB_zip = select_first([UKBB_zip, wget.UKBB_zip]),
      docker = docker
  }

  call split_report {
    input:
      sqc_file = sqc_file,
      docker = docker
  }

  call split_snp_posterior {
    input:
      snp_posterior_tar = select_first([snp_posterior_tar, wget.snp_posterior_tar]),
      sqc_file = sqc_file,
      snp_qc_file = select_first([snp_qc_file, wget.snp_qc_file]),
      split_elf = init.split_elf,
      dump_elf = init.dump_elf,
      docker = docker
  }

  scatter(idx in range(length(chrs))) {
    call split_cal {
      input:
        chr = chrs[idx],
        UKBL_zip = select_first([UKBL_zip, wget.UKBL_zip]),
        UKBB_zip = select_first([UKBB_zip, wget.UKBB_zip]),
        snp_qc_file = select_first([snp_qc_file, wget.snp_qc_file]),
        sqc_file = sqc_file,
        cal_file = cal_path + (if sub(cal_path, "/$", "") != cal_path then "" else "/") + cal_filename_pfx + chrs[idx] + cal_filename_sfx,
        unpack_elf = init.unpack_elf,
        split_elf = init.split_elf,
        docker = docker
    }

    call split_int {
      input:
        sqc_file = sqc_file,
        chr = chrs[idx],
        int_file = int_path + (if sub(int_path, "/$", "") != int_path then "" else "/") + int_filename_pfx + chrs[idx] + int_filename_sfx,
        split_elf = init.split_elf,
        docker = docker
    }
  }

  Array[Array[String]] tsv = transpose(init.tsv)
  Array[Array[File]] split_cal_files = transpose(split_cal.files)
  Array[Array[File]] split_int_files = transpose(split_int.files)

  scatter(idx in range(length(init.tsv))) {
    call revert_calls {
      input:
        n = init.tsv[idx][0],
        array = init.tsv[idx][1],
        batch = init.tsv[idx][2],
        snp_qc_file = select_first([snp_qc_file, wget.snp_qc_file]),
        sqc_file = sqc_file,
        cal_files = split_cal_files[idx],
        docker = docker
    }

    call revert_summary {
      input:
        n = init.tsv[idx][0],
        array = init.tsv[idx][1],
        batch = init.tsv[idx][2],
        snp_qc_file = select_first([snp_qc_file, wget.snp_qc_file]),
        sqc_file = sqc_file,
        int_files = split_int_files[idx],
        dump_elf = init.dump_elf,
        docker = docker
    }
  }

  output {
    File UKBL_csv = csv.UKBL_csv
    File UKBB_csv = csv.UKBB_csv
    Array[File] report_files = split_report.files
    Array[File] snp_posteriors_files = split_snp_posterior.files
    Array[File] calls_files = revert_calls.file
    Array[File] summary_files = revert_summary.file
    File sample_tsv = init.sample_tsv
    File batch_tsv = init.batch_tsv
  }
}

task init {
  input {
    File sqc_file
    String docker
    Int cpu = 1
    Int disk_size = 10
    Float memory = 3.5
    Int preemptible = 1
    Int maxRetries = 0
  }

  command <<<
    set -euo pipefail
    apt-get -qqy update --fix-missing
    apt-get -qqy install --no-install-recommends gcc libc6-dev 1>&2

    echo '#include <stdio.h>
    #include <stdlib.h>
    #include <string.h>
    int main(int argc, char **argv) {
      if (argc != 3)
        return 1;
      int i, n, nread, swap;
      if (sscanf(argv[1], "%d", &n) != 1)
        return 1;
      char *in_buffer = (char *)malloc((n + 3) / 4 * sizeof(char));
      char *out_buffer = (char *)malloc((n + 3) / 4 * 4 * sizeof(char));
      FILE *in = fopen(argv[2], "r");
      while ((nread = fread(in_buffer, 1, (n + 3) / 4, stdin))) {
        if (nread != (n + 3) / 4)
          return 1;
        if (fscanf(in, "%d\n", &swap) != 1)
          return 1;
        for (i = 0; i < (n + 3) / 4; i++) {
          char x = in_buffer[i];
          char y = (swap ? (~x & 0xAA) : (x & 0x55) << 1) ^
                   ((x & 0x55) ^ ((x & 0xAA) >> 1));
          int z = (int)(y & 0x03) ^ ((int)(y & 0x0C) << 6) ^
                  ((int)(y & 0x30) << 12) ^ ((int)(y & 0xC0) << 18) ^ 0x30303030;
          memcpy(&out_buffer[4 * i], &z, 4);
        }
        fwrite(out_buffer, 1, n, stdout);
      }
      free(in_buffer);
      free(out_buffer);
      return 0;
    }' > unpack.c
    gcc -O3 -o unpack unpack.c

    echo '#include <stdio.h>
    #include <stdlib.h>
    #include <string.h>
    int main(int argc, char **argv) {
      if (argc != 3)
        return 1;
      int i, nread, n = 0, m = 4;
      int *offsets = (int *)malloc(m * sizeof(int));
      FILE **batches = (FILE **)malloc(m * sizeof(FILE *));
      FILE *in = fopen(argv[2], "r");
      char buffer[BUFSIZ];
      offsets[0] = 0;
      while (fscanf(in, "%d %s\n", &i, buffer) == 2) {
        sprintf(buffer + strlen(buffer), "%s", argv[1]);
        batches[n] = fopen(buffer, "w");
        n++;
        if (n == m) {
          m <<= 1;
          offsets = (int *)realloc(offsets, m * sizeof(int));
          batches = (FILE **)realloc(batches, m * sizeof(FILE *));
        }
        offsets[n] = offsets[n - 1] + i;
      }
      fclose(in);
      char *line = (char *)malloc(offsets[n] * sizeof(char *));
      while ((nread = fread(line, 1, offsets[n], stdin))) {
        if (nread != offsets[n])
          return 1;
        for (i = 0; i < n; i++)
          fwrite(line + offsets[i], 1, offsets[i + 1] - offsets[i], batches[i]);
      }
      for (i = 0; i < n; i++)
        fclose(batches[i]);
      free(offsets);
      free(batches);
      free(line);
      return 0;
    }' > split.c
    gcc -O3 -o split split.c

    echo '#include <stdio.h>
    #include <stdlib.h>
    #include <string.h>
    int main(int argc, char **argv) {
      if (argc != 3)
        return 1;
      int i, j, n, m, nread;
      if (sscanf(argv[1], "%d", &n) != 1)
        return 1;
      if (sscanf(argv[2], "%d", &m) != 1)
        return 1;
      char *buffer = (char *)malloc(n * m * sizeof(float));
      while ((nread = fread(buffer, sizeof(float), n * m, stdin))) {
        if (nread != n * m)
          return 1;
        for (i = 0; i < n; i++)
          for (j = 0; j < m; j++)
            fprintf(stdout, "%g%c", *(float *)&buffer[4 * (n * j + i)],
                    j + 1 == m ? 0x0A : 0x09);
      }
      free(buffer);
      return 0;
    }' > dump.c
    gcc -O3 -o dump dump.c

    mv "~{sqc_file}" .
    mkdir -p out
    awk 'BEGIN {print "sample_id\tbatch_id\tcel"} {printf "%s\t%s\t%s\n",$1,$4,$1}' "~{basename(sqc_file)}" > out/ukb.sample.tsv
    cut -d" " -f3,4 "~{basename(sqc_file)}" | uniq -c | \
      awk 'BEGIN {csv["UKBB"]="Axiom_UKB_WCSG.na34.annot.csv.gz"; csv["UKBL"]="Axiom_UKBiLEVE.na34.annot.csv.gz"
      print "batch_id\tn_smpls\tcsv\tsnp\treport\tcalls\tsummary"}
      {printf "%s\t%s\t%s\t%s.AxiomGT1.snp-posteriors.txt.gz\t%s.AxiomGT1.report.txt.gz\t%s.AxiomGT1.calls.txt.gz\t%s.AxiomGT1.summary.txt.gz\n",
      $3,$1,csv[$2],$3,$3,$3,$3}' > out/ukb.batch.tsv
    cut -d" " -f3,4 "~{basename(sqc_file)}" | uniq -c | awk '{print $1"\t"$2"\t"$3}'
  >>>

  output {
    File unpack_elf = "unpack"
    File split_elf = "split"
    File dump_elf = "dump"
    File sample_tsv = "out/ukb.sample.tsv"
    File batch_tsv = "out/ukb.batch.tsv"
    Array[Array[String]] tsv = read_tsv(stdout())
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

task wget {
  input {
    String docker
    Int cpu = 1
    Int disk_size = 20
    Float memory = 3.5
    Int preemptible = 1
    Int maxRetries = 0
  }

  command <<<
    set -euo pipefail
    apt-get -qqy update --fix-missing
    apt-get -qqy install --no-install-recommends wget 1>&2
    wget -nd --no-check-certificate biobank.ctsu.ox.ac.uk/crystal/crystal/auxdata/ukb_snp_posterior.tar
    wget -nd --no-check-certificate biobank.ctsu.ox.ac.uk/crystal/crystal/auxdata/ukb_snp_qc.txt
    wget -nd --no-check-certificate biobank.ctsu.ox.ac.uk/crystal/crystal/docs/Array_BIL_34.zip
    wget -nd --no-check-certificate biobank.ctsu.ox.ac.uk/crystal/crystal/docs/Array_UKB_34.zip
  >>>

  output {
    File UKBL_zip = "Array_BIL_34.zip"
    File UKBB_zip = "Array_UKB_34.zip"
    File snp_posterior_tar = "ukb_snp_posterior.tar"
    File snp_qc_file = "ukb_snp_qc.txt"
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

task csv {
  input {
    File UKBL_zip
    File UKBB_zip
    String docker
    Int cpu = 1
    Int disk_size = 10
    Float memory = 3.5
    Int preemptible = 1
    Int maxRetries = 0
  }

  command <<<
    set -euo pipefail
    mv "~{UKBL_zip}" .
    mv "~{UKBB_zip}" .
    mkdir -p out
    zcat "~{basename(UKBL_zip)}" | \
      sed -e 's/ACTGAAGCATGACTCTTCTAAAACAAGTTAATAAT\[-\/GCATTCTACAGAAAACCTTTTAGACATTACAACCTCTAACCTAAAAAATGGAACTCTAAATACCTTTTCACCATCTCAGTAGTAATTTCTGACACATCACTTTCAGAAGGCTCATAAGCCAACTCAGATACAACCGTAATATTAATGAGGAAAGCAGTTAGTGTCAAAGGAAGAGGATTGAAAGTAGAAGTGTTATGAAACTCCTCTTTTAC/ACTGAAGCATGACTCTTCTAAAACAAGTTAATAAT\[-\/GCATTCTACAGAAAACCTTTTAGACATTACAACCTCTAACCTAAAAAATGGAACTCTAAATACCTTTTCACCATCTCAGTAGTAATTTCTGACACATCACTTTCAGAAGGCTCATAAGCCAACTCAGATACAACCGTAATATTAATGAGGAAAGCAGTTAGTGTCAAAGGAAGAGGATTGAAAGTAGAAGTGTTATGAAACTCCTCTTTTACTTTTTTGCATGCTATCAAAATCACTGCCAATAATATGAATGTCAGGTCAATTTCCATAGGTAAATCCGTTACCTTTTACCTCTTTAAAAGAAAAGTTTTGCAGAAGAGGGCTGAAAATTTCTTGAGCCATTTCAGCACAAAGAATGGAAGTTCATTTCTCACCATGATACAACTCTACCCTGCTGTCATCTTCATTGTGATGGTGGCAGAAGTTTAGCAGGGTGCAAGTGACCACTA\]AATGACATCTTTTCATGAACAATTGATAAATCTTT/' | \
        gzip > out/Axiom_UKBiLEVE.na34.annot.csv.gz
    zcat "~{basename(UKBB_zip)}" | \
      sed -e 's/ACTGAAGCATGACTCTTCTAAAACAAGTTAATAAT\[-\/GCATTCTACAGAAAACCTTTTAGACATTACAACCTCTAACCTAAAAAATGGAACTCTAAATACCTTTTCACCATCTCAGTAGTAATTTCTGACACATCACTTTCAGAAGGCTCATAAGCCAACTCAGATACAACCGTAATATTAATGAGGAAAGCAGTTAGTGTCAAAGGAAGAGGATTGAAAGTAGAAGTGTTATGAAACTCCTCTTTTAC/ACTGAAGCATGACTCTTCTAAAACAAGTTAATAAT\[-\/GCATTCTACAGAAAACCTTTTAGACATTACAACCTCTAACCTAAAAAATGGAACTCTAAATACCTTTTCACCATCTCAGTAGTAATTTCTGACACATCACTTTCAGAAGGCTCATAAGCCAACTCAGATACAACCGTAATATTAATGAGGAAAGCAGTTAGTGTCAAAGGAAGAGGATTGAAAGTAGAAGTGTTATGAAACTCCTCTTTTACTTTTTTGCATGCTATCAAAATCACTGCCAATAATATGAATGTCAGGTCAATTTCCATAGGTAAATCCGTTACCTTTTACCTCTTTAAAAGAAAAGTTTTGCAGAAGAGGGCTGAAAATTTCTTGAGCCATTTCAGCACAAAGAATGGAAGTTCATTTCTCACCATGATACAACTCTACCCTGCTGTCATCTTCATTGTGATGGTGGCAGAAGTTTAGCAGGGTGCAAGTGACCACTA\]AATGACATCTTTTCATGAACAATTGATAAATCTTT/' \
          -e 's/TATTGGCTCAATTTCTTGGGGAGGGGGTGCTGTCA\[-\/GAGATTGTTATCTGAGGATGTGACATAGATTTCTCAGGGCACAATTTCAACTACTTTTTCAGCTTTAGGGTTTTTAGATACGTTTGTACCACAATTGAGCATGGGAGGGAGAGGGGTGAGCCTAAGCAGTGATGGCTGATTTCTGTCATGTCTGTCATGTGTCCCCCAGTACCTCCAGAGGTAACTGTGCTCACAAACAGCCCTGTGGAACT/TATTGGCTCAATTTCTTGGGGAGGGGGTGCTGTCA\[-\/GAGATTGTTATCTGAGGATGTGACATAGATTTCTCAGGGCACAATTTCAACTACTTTTTCAGCTTTAGGGTTTTTAGATACGTTTGTACCACAATTGAGCATGGGAGGGAGAGGGGTGAGCCTAAGCAGTGATGGCTGATTTCTGTCATGTCTGTCATGTGTCCCCCAGTACCTCCAGAGGTAACTGTGCTCACAAACAGCCCTGTGGAACTGAGAGAGCCCAACGTCCTCATCTGTTTCATA\]GACAAGTTCACCCCACCAGTGGTCAATGTCACGTG/' \
          -e 's/GTAATAGGACTCCTTGAAGAGACTTTCGGCGGGGC\[-\/CTGGGGAAGACAGAGGGAGACACACTGAACAGTCTGGCTGTGACCTCCAGGGCCACTGTCACCCACACAGCTACCAGGTGGCCAGAGCCAGGATCTGAACCCAGGTCTGTGGGGGATCCACACCTGAATCCCATTCTTGGGGGAGTCTCATTGGCACCACAGCAGAGGAACCTCTAACCTAGGCCTTCGTTCAAGACTAGAACCTGCCCCCA/GTAATAGGACTCCTTGAAGAGACTTTCGGCGGGGC\[-\/CTGGGGAAGACAGAGGGAGACACACTGAACAGTCTGGCTGTGACCTCCAGGGCCACTGTCACCCACACAGCTACCAGGTGGCCAGAGCCAGGATCTGAACCCAGGTCTGTGGGGGATCCACACCTGAATCCCATTCTTGGGGGAGTCTCATTGGCACCACAGCAGAGGAACCTCTAACCTAGGCCTTCGTTCAAGACTAGAACCTGCCCCCACTCCAGTGCTGACCACTTCAGAGCAGAGGGGAGGCTGAAGAGGACACAGGGTCCTCAGTGTCCCAATGCCAGATCCCCACTCTCCTTGGTCACCAGCTTGTGAATCTGGGCAGTCGCCTGGCTCCTGCCTACTGTCCTGAGCCATGTTTCAGAGGGCAGGTAACAAATGAGAAGGGAAAAGTACAGCTCTAGTTCGGGGGGTGGGAGGCCGCTCTATCCTTTACTCTGAAGGCCTGGGGGAGGCTGACCTCCAGACCTGCAGCTGCCAGAAAACCCTGGGGCCCATCCACTGCTTAC\]CCATGGGGTCTGAGGAGTCAGTGATGATCACGTCG/' | \
      gzip > out/Axiom_UKB_WCSG.na34.annot.csv.gz
  >>>

  output {
    File UKBL_csv = "out/Axiom_UKBiLEVE.na34.annot.csv.gz"
    File UKBB_csv = "out/Axiom_UKB_WCSG.na34.annot.csv.gz"
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

task split_report {
  input {
    File sqc_file
    String docker
    Int cpu = 1
    Int disk_size = 10
    Float memory = 3.5
    Int preemptible = 1
    Int maxRetries = 0
  }

  command <<<
    set -euo pipefail
    mv "~{sqc_file}" .
    cut -d" " -f4 "~{basename(sqc_file)}" | uniq > ukb_snp_posterior.batch
    awk 'NR==FNR {print "cel_files\tcomputed_gender\tcall_rate\tcn-probe-chrXY-ratio_gender_meanX\tcn-probe-chrXY-ratio_gender_meanY" > $1".AxiomGT1.report.txt"}
      NR>FNR {printf "%s\t%s\t%s\t%s\t%s\n",$1,$11,$7,$12,$13 > $4".AxiomGT1.report.txt"}' ukb_snp_posterior.batch "~{basename(sqc_file)}"
    sed 's/$/.AxiomGT1.report.txt/' ukb_snp_posterior.batch | tr '\n' '\0' | xargs -0 gzip --force --no-name
    mkdir -p out
    sed 's/$/.AxiomGT1.report.txt.gz/' ukb_snp_posterior.batch | tr '\n' '\0' | xargs -0 mv -t out/
    sed 's/^/out\//;s/$/.AxiomGT1.report.txt.gz/' ukb_snp_posterior.batch
  >>>

  output {
    Directory out = "out"
    Array[File] files = read_lines(stdout())
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

task split_snp_posterior {
  input {
    File snp_posterior_tar
    File sqc_file
    File snp_qc_file
    File split_elf
    File dump_elf
    String docker
    Int cpu = 1
    Int? disk_size_override
    Float memory = 3.5
    Int preemptible = 1
    Int maxRetries = 0
  }

  Float snp_posterior_size = size(snp_posterior_tar, "GiB")
  Int disk_size = select_first([disk_size_override, ceil(10.0 + 3.0 * snp_posterior_size)])

  command <<<
    set -euo pipefail
    mv "~{snp_posterior_tar}" .
    mv "~{snp_qc_file}" .
    mv "~{sqc_file}" .
    mv "~{split_elf}" .
    mv "~{dump_elf}" .
    chmod a+x "~{basename(split_elf)}" "~{basename(dump_elf)}"
    tar xvf "~{basename(snp_posterior_tar)}" 1>&2
    cat ukb_snp_posterior_chr{{1..22},X,Y,XY,MT}.bin | "./~{basename(split_elf)}" .snp-posteriors.bin <(sed 's/^/132 /' ukb_snp_posterior.batch)
    awk 'NR==FNR {x[$2]++} NR>FNR && FNR>1 {array=$9==0 || $9==2; print array"\t"$3;
      if ($1 in x) print array"\t"$3":1"}' ukb_snp_posterior_chrX_haploid.bim "~{basename(snp_qc_file)}" > snp-posteriors.UKBL.tsv
    awk 'NR==FNR {x[$2]++} NR>FNR && FNR>1 {array=$9==1 || $9==2; print array"\t"$3;
      if ($1 in x) print array"\t"$3":1"}' ukb_snp_posterior_chrX_haploid.bim "~{basename(snp_qc_file)}" > snp-posteriors.UKBB.tsv
    cut -d" " -f3,4 "~{basename(sqc_file)}" | uniq | while read array batch; do
      (echo -e "#%SnpPosteriorFormatVer=1\n#%data-order=meanX,varX,nObsMean,nObsVar,meanY,varY,covarXY\nid\tBB\tAB\tAA\tCV";
      cat $batch.snp-posteriors.bin | "./~{basename(dump_elf)}" 1 33 | paste "snp-posteriors.$array.tsv" - | sed '/^0/d;s/^..//' | \
      awk '{fmt="%s\t%s,%s,%s,%s,%s,%s,%s\t%s,%s,%s,%s,%s,%s,%s\t%s,%s,%s,%s,%s,%s,%s\t%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n"
      if ($2<$16) printf fmt,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,$29,$30,$31,$32,$33,$34;
             else printf fmt,$1,$16,$17,$18,$19,$20,$21,$22,$9,$10,$11,$12,$13,$14,$15,$2,$3,$4,$5,$6,$7,$8,$23,$24,$25,$26,$27,$28,$29,$30,$31,$32,$33,$34}') | \
        gzip > "$batch.AxiomGT1.snp-posteriors.txt.gz"
    done
    mkdir -p out
    sed 's/$/.AxiomGT1.snp-posteriors.txt.gz/' ukb_snp_posterior.batch | tr '\n' '\0' | xargs -0 mv -t out/
    sed 's/^/out\//;s/$/.AxiomGT1.snp-posteriors.txt.gz/' ukb_snp_posterior.batch
  >>>

  output {
    Directory out = "out"
    Array[File] files = read_lines(stdout())
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

task split_cal {
  input {
    String chr
    File UKBL_zip
    File UKBB_zip
    File snp_qc_file
    File sqc_file
    File cal_file
    File unpack_elf
    File split_elf
    String docker
    Int cpu = 1
    Int? disk_size_override
    Float memory = 3.5
    Int preemptible = 1
    Int maxRetries = 0
  }

  Float cal_size = size(cal_file, "GiB")
  Int disk_size = select_first([disk_size_override, ceil(10.0 + 5.0 * cal_size)])
  Int plink_chr = sub(sub(sub(sub(chr, "MT", "26"), "XY", "25"), "Y", "24"), "X", "23")

  command <<<
    set -euo pipefail
    mv "~{UKBL_zip}" .
    mv "~{UKBB_zip}" .
    mv "~{snp_qc_file}" .
    mv "~{sqc_file}" .
    mv "~{cal_file}" .
    mv "~{unpack_elf}" .
    mv "~{split_elf}" .
    chmod a+x "~{basename(unpack_elf)}" "~{basename(split_elf)}"
    (zcat "~{basename(UKBL_zip)}" | grep -v ^# | tail -n+2;
     zcat "~{basename(UKBB_zip)}" | grep -v ^# | tail -n+2) | tr -d '"' | \
      awk -F, '{if ($12==$14 && $13==$15) x=0; else if ($12==$15 && $13==$14) x=1; else x=0; print $1,x}' | \
      awk 'NR==FNR {x[$1]=$2} NR>FNR && $4=="~{plink_chr}" {print x[$3]}' - "~{basename(snp_qc_file)}" > swap.lines
    tail -qc+4 "~{basename(cal_file)}" | "./~{basename(unpack_elf)}" $(cat "~{basename(sqc_file)}" | wc -l) swap.lines | \
      "./~{basename(split_elf)}" ".chr~{chr}.calls.bin" <(cut -d" " -f4 "~{basename(sqc_file)}" | uniq -c | sed 's/^ *//')
    mkdir -p out
    cut -d" " -f4 "~{basename(sqc_file)}" | uniq | sed 's/$/.chr~{chr}.calls.bin/' | tr '\n' '\0' | xargs -0 mv -t out/
    cut -d" " -f4 "~{basename(sqc_file)}" | uniq | sed 's/^/out\//;s/$/.chr~{chr}.calls.bin/'
  >>>

  output {
    Directory out = "out"
    Array[File] files = read_lines(stdout())
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

task split_int {
  input {
    String chr
    File sqc_file
    File int_file
    File split_elf
    String docker
    Int cpu = 1
    Int? disk_size_override
    Float memory = 3.5
    Int preemptible = 1
    Int maxRetries = 0
  }

  Float int_size = size(int_file, "GiB")
  Int disk_size = select_first([disk_size_override, ceil(10.0 + 2.0 * int_size)])

  command <<<
    set -euo pipefail
    mv "~{sqc_file}" .
    mv "~{int_file}" .
    mv "~{split_elf}" .
    chmod a+x "~{basename(split_elf)}"
    cat "~{basename(int_file)}" | "./~{basename(split_elf)}" ".chr~{chr}.summary.bin" <(cut -d" " -f4 "~{basename(sqc_file)}" | uniq -c | awk '{print 8*$1,$2}')
    mkdir -p out/
    cut -d" " -f4 "~{basename(sqc_file)}" | uniq | sed 's/$/.chr~{chr}.summary.bin/' | tr '\n' '\0' | xargs -0 mv -t out/
    cut -d" " -f4 "~{basename(sqc_file)}" | uniq | sed 's/^/out\//;s/$/.chr~{chr}.summary.bin/'
  >>>

  output {
    Directory out = "out"
    Array[File] files = read_lines(stdout())
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

task revert_calls {
  input {
    Int n
    String array
    String batch
    File snp_qc_file
    File sqc_file
    Array[File] cal_files
    String docker
    Int cpu = 1
    Int? disk_size_override
    Float memory = 3.5
    Int preemptible = 1
    Int maxRetries = 0
  }

  Float cal_size = size(cal_files, "GiB")
  Int disk_size = select_first([disk_size_override, ceil(10.0 + cal_size)])

  command <<<
    set -euo pipefail
    mv "~{snp_qc_file}" .
    mv "~{sqc_file}" .
    cal_files=~{write_lines(cal_files)}
    cat $cal_files | tr '\n' '\0' | xargs -0 mv -t .
    sed -i 's/^.*\///' $cal_files
    awk 'NR>1 {array=$9==0 || $9==2; print array"\t"$3}' "~{basename(snp_qc_file)}" > calls.UKBL.tsv
    awk 'NR>1 {array=$9==1 || $9==2; print array"\t"$3}' "~{basename(snp_qc_file)}" > calls.UKBB.tsv
    (echo -en "#Calls: -1=NN, 0=AA, 1=AB, 2=BB\nprobeset_id\t";
    awk -v batch="~{batch}" '$4==batch {print $1}' "~{basename(sqc_file)}" | tr '\n' '\t' | sed 's/\t$/\n/';
    cat $cal_files | tr '\n' '\0' | xargs -0 cat | fold -w~{n} | tr '\n' '\r' | fold -w1 | tr '\n\r' '\t\n' | sed 's/3/-1/g' | \
    paste "calls.~{array}.tsv" - | sed '/^0/d;s/^..//') | \
      gzip > "~{batch}.AxiomGT1.calls.txt.gz"
    mkdir -p out
    mv "~{batch}.AxiomGT1.calls.txt.gz" out/
  >>>

  output {
    File file = "out/" + batch + ".AxiomGT1.calls.txt.gz"
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

task revert_summary {
  input {
    Int n
    String array
    String batch
    File snp_qc_file
    File sqc_file
    Array[File] int_files
    File dump_elf
    String docker
    Int cpu = 1
    Int? disk_size_override
    Float memory = 3.5
    Int preemptible = 1
    Int maxRetries = 0
  }

  Float int_size = size(int_files, "GiB")
  Int disk_size = select_first([disk_size_override, ceil(10.0 + 2.0 * int_size)])

  command <<<
    set -euo pipefail
    mv "~{snp_qc_file}" .
    mv "~{sqc_file}" .
    mv "~{dump_elf}" .
    int_files=~{write_lines(int_files)}
    cat $int_files | tr '\n' '\0' | xargs -0 mv -t .
    sed -i 's/^.*\///' $int_files
    chmod a+x "~{basename(dump_elf)}"
    awk 'NR>1 {array=$9==0 || $9==2; print array"\t"$3"-A"; print array"\t"$3"-B"}' "~{basename(snp_qc_file)}" > summary.UKBL.tsv
    awk 'NR>1 {array=$9==1 || $9==2; print array"\t"$3"-A"; print array"\t"$3"-B"}' "~{basename(snp_qc_file)}" > summary.UKBB.tsv
    (echo -en "probeset_id\t";
    awk -v batch="~{batch}" '$4==batch {print $1}' "~{basename(sqc_file)}" | tr '\n' '\t' | sed 's/\t$/\n/';
    cat $int_files | tr '\n' '\0' | xargs -0 cat | "./~{basename(dump_elf)}" 2 ~{n} | paste "summary.~{array}.tsv" - | sed '/^0/d;s/^..//') | \
      gzip > "~{batch}.AxiomGT1.summary.txt.gz"
    mkdir -p out
    mv "~{batch}.AxiomGT1.summary.txt.gz" out/
  >>>

  output {
    File file = "out/" + batch + ".AxiomGT1.summary.txt.gz"
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
