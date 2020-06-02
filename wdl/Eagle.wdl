version 1.0

workflow Eagle {

	input {
		File unphased_bcf
		File unphased_bcf_idx

		File sex_file

		File ref_fasta
		File ref_fasta_fai

		File segdup_exclusion_annotations
		File genetic_map_file

		File primary_contigs_list

		Array[File] reference_vcfs
		Array[File] reference_vcfs_idx
	}

	Array[String] contigs = read_lines(primary_contigs_list)

	call MinimizeVcf {
		input:
			bcf = unphased_bcf,
			bcf_idx = unphased_bcf_idx			
	}

	call ExcludeVariants {
		input:
			bcf = MinimizeVcf.out,
			bcf_idx = MinimizeVcf.out_idx,
			sex_file = sex_file,
			segdup_exclusion_annotations = segdup_exclusion_annotations
	}

	scatter (idx in range(length(contigs))) {
		call EagleChromosome {
			input:
				chromosome = contigs[idx],
				unphased_bcf = MinimizeVcf.out,
				unphased_bcf_idx = MinimizeVcf.out_idx,
				genetic_map_file = genetic_map_file,
				vcf_ref = reference_vcfs[idx],
				vcf_ref_idx = reference_vcfs_idx[idx],
				excluded_variants = ExcludeVariants.excluded_variants,
				excluded_variants_idx = ExcludeVariants.excluded_variants_idx,

		}
	}

	call ExtractUnphasedChromosomeData {
		input:
			bcf = MinimizeVcf.out,
			bcf_idx = MinimizeVcf.out_idx,
			phased_chromosomes = contigs
	}

	call GatherFinalVcf {
		input:
			phased_chrom_files = EagleChromosome.phased_bcf,
			other_chrom_file = ExtractUnphasedChromosomeData.out,
			ref_fasta = ref_fasta
	}

}

task MinimizeVcf {
	input {
		File bcf
		File bcf_idx

		Int disk_size = 10
		String docker = "cwhelan/mocha:v1.0"
		Int threads = 1
		Int memory = 2
		Int preemptible_attempts = 3

	}

	String filebase = basename(bcf, ".bcf")

	output {
		File out = "~{filebase}.unphased.bcf"
		File out_idx = "~{filebase}.unphased.bcf.csi"
	}

	command <<<
		set -euo pipefail

		bcftools annotate --no-version -Ob -o ~{filebase}.unphased.bcf ~{bcf} \
			-x ID,QUAL,^INFO/ALLELE_A,^INFO/ALLELE_B,^FMT/GT,^FMT/BAF,^FMT/LRR && \
			bcftools index -f ~{filebase}.unphased.bcf
	>>>

	runtime {
		memory: memory + " GiB"
		disks: "local-disk " + disk_size + " HDD"
		cpu: threads
		docker: docker
		preemptible_attempts: preemptible_attempts
	}
}

task ExcludeVariants {
	input {
		File bcf
		File bcf_idx
		File sex_file

		File segdup_exclusion_annotations

		Int disk_size = 10
		String docker = "cwhelan/mocha:v1.0"
		Int threads = 1
		Int memory = 2
		Int preemptible_attempts = 3

	}

	String filebase = basename(bcf, ".bcf")
	String excluded_variants = "~{filebase}.xcl.bcf"

	output {
		File excluded_variants = "~{excluded_variants}"
		File excluded_variants_idx = "~{excluded_variants}.csi"
	}

	command <<<

		set -euo pipefail

		n=$(bcftools query -l ~{bcf} |wc -l); \
		ns=$((n*98/100)); \
		echo '##INFO=<ID=JK,Number=1,Type=Float,Description="Jukes Cantor">' | \
			bcftools annotate --no-version -Ou -a ~{segdup_exclusion_annotations} -c CHROM,FROM,TO,JK -h /dev/stdin ~{bcf} | \
			bcftools +fill-tags --no-version -Ou -t ^Y,MT,chrY,chrM -- -t NS,ExcHet | \
			bcftools +mochatools --no-version -Ou -- -x ~{sex_file} -G | \
			bcftools annotate --no-version -Ob -o ~{excluded_variants} \
			-i 'FILTER!="." && FILTER!="PASS" || JK<.02 || NS<'$ns' || ExcHet<1e-6 || AC_Sex_Test>6' \
			-x FILTER,^INFO/JK,^INFO/NS,^INFO/ExcHet,^INFO/AC_Sex_Test && \
			bcftools index -f ~{excluded_variants}
	>>>

	runtime {
		memory: memory + " GiB"
		disks: "local-disk " + disk_size + " HDD"
		cpu: threads
		docker: docker
		preemptible_attempts: preemptible_attempts
	}
}

task EagleChromosome {

	input {
		File unphased_bcf
		File unphased_bcf_idx

		File genetic_map_file

		File vcf_ref
		File vcf_ref_idx

		File excluded_variants
		File excluded_variants_idx

		String chromosome		

		Int disk_size = 100
		String docker = "cwhelan/mocha:v1.0"
		Int threads = 4
		Int memory = 8
		Int preemptible_attempts = 3

	}

	String filebase = basename(unphased_bcf, ".bcf")
	String outname = "~{filebase}.~{chromosome}.bcf"

	output {
		File phased_bcf = "~{outname}"
		File phased_bcf_idx = "~{outname}.csi"
	}


	command <<<
		set -euo pipefail

		eagle \
			--geneticMapFile ~{genetic_map_file} \
			--outPrefix ~{filebase}.~{chromosome} \
			--numThreads ~{threads} \
			--vcfRef ~{vcf_ref} \
			--vcfTarget ~{unphased_bcf} \
			--vcfOutFormat b \
			--noImpMissing \
			--outputUnphased \
			--vcfExclude ~{excluded_variants} \
			--chrom ~{chromosome} \
			--pbwtIters 3 && \
			bcftools index -f ~{outname}
	>>>

	runtime {
		memory: memory + " GiB"
		disks: "local-disk " + disk_size + " HDD"
		cpu: threads
		docker: docker
		preemptible_attempts: preemptible_attempts
	}

}

task ExtractUnphasedChromosomeData {
	input {
		File bcf
		File bcf_idx

		Array[String] phased_chromosomes

		Int disk_size = 10
		String docker = "cwhelan/mocha:v1.0"
		Int threads = 1
		Int memory = 2
		Int preemptible_attempts = 3

	}

	String filebase = basename(bcf, ".bcf")

	output {
		File out = "~{filebase}.other.bcf"
		File out_idx = "~{filebase}.other.bcf.csi"
	}

	command <<<
		set -euo pipefail

		bcftools view --no-version -Ob -o ~{out} ~{bcf} \
			-t ^~{sep="," phased_chromosomes} && \
			bcftools index -f ~{out}
	>>>

	runtime {
		memory: memory + " GiB"
		disks: "local-disk " + disk_size + " HDD"
		cpu: threads
		docker: docker
		preemptible_attempts: preemptible_attempts
	}
}


task GatherFinalVcf {
	input {
		Array[File] phased_chrom_files
		File other_chrom_file

		File ref_fasta

		Int disk_size = 10
		String docker = "cwhelan/mocha:v1.0"
		Int threads = 4
		Int memory = 2
		Int preemptible_attempts = 3

	}

	String filebase = basename(other_chrom_file, "other.bcf")

	output {
		File out = "~{filebase}.phased.bcf"
		File out_idx = "~{filebase}.phased.bcf.csi"
	}

	command <<<
		set -euo pipefail

		bcftools concat --no-version -Ou --threads ~{threads} ~{sep=" " phased_chrom_files} ~{other_chrom_file} | \
			bcftools +mochatools --no-version -Ob -o ~{out} -- -f ~{ref_fasta} && \
			bcftools index -f ~{out}
	>>>

	runtime {
		memory: memory + " GiB"
		disks: "local-disk " + disk_size + " HDD"
		cpu: threads
		docker: docker
		preemptible_attempts: preemptible_attempts
	}
}
