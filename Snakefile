import os
import certifi
os.environ['SSL_CERT_FILE'] = certifi.where()

configfile: "config.yaml"
all_samples= config["all_samples"]
reference = config["reference"]

rule targets:
	input:
		expand("analysis/{reference}/{sample}.dedup.bam",sample=all_samples, reference=reference),
		expand("analysis/{reference}/{sample}.flagstat.txt",sample=all_samples, reference=reference),
		expand("analysis/{reference}/coverage_plot.png", reference=reference),
		expand("analysis/{reference}/genotyped_cohort.vcf.gz", reference=reference),
		#expand("{sample}_haplotagged_15.bam",sample=all_samples),
		

rule prep_ref:
	output:
		"analysis/{reference}/ref.prepped"
	shell:
		"""
		bash prep_ref.sh {reference}
		"""

rule download_reads:
	output:
		"{sample}_1.fastq.gz",
		"{sample}_2.fastq.gz",
	threads: 6  # Number of CPU cores/threads
	resources:
		mem_mb=12000  # RAM in Megabytes (8GB)
	shell:
		"""
		module load SRA-Toolkit/3.0.0-rhel8
		fasterq-dump {wildcards.sample}
		bgzip {wildcards.sample}_1.fastq
		bgzip {wildcards.sample}_2.fastq
		"""
rule bwa_align:
	input:
		r1="{sample}_1.fastq.gz",
		r2="{sample}_2.fastq.gz", 
		ref="analysis/{reference}/ref.prepped"
	output:
		"analysis/{reference}/{sample}.dedup.bam"
	shell:
		"""
		bash alignment.sh {reference} {input.r1} {input.r2} {wildcards.sample}
		"""

rule alignment_qc:
	input:
		"analysis/{reference}/{sample}.dedup.bam"
	output:
		"analysis/{reference}/{sample}.flagstat.txt"	
	shell:
		"""
		bash qc_alignment.sh "{wildcards.sample}" {reference}
		"""

bams=" ".join(expand("analysis/{reference}/{sample}.dedup.bam",sample=all_samples, reference=reference))

rule plot_coverage:
	input:
		expand("analysis/{reference}/{sample}.dedup.bam",sample=all_samples, reference=reference)	
	output:
		"analysis/{reference}/coverage_plot.png"
	shell:
		"""
		plotCoverage -p 32 -b {bams} -o {output}
		"""

rule per_sample_variant_calling:
	input:
		"analysis/{reference}/{sample}.dedup.bam"
	threads: 6  # Number of CPU cores/threads
	resources:
		mem_mb=12000  # RAM in Megabytes (8GB)
	output:
		"analysis/{reference}/{sample}_spark.g.vcf.gz"
	shell:
		"""
		bash per_sample_variant_calling.sh {input} {output} {reference}
		"""

gvcfs=" ".join(expand("--variant analysis/{reference}/{sample}_spark.g.vcf.gz",sample=all_samples, reference=reference))
rule joint_variant_calling:
	input:
		expand("analysis/{reference}/{sample}_spark.g.vcf.gz",sample=all_samples, reference=reference)
	output:
		"analysis/{reference}/genotyped_cohort.vcf.gz"
	shell:
		"""
		singularity exec docker://broadinstitute/gatk:4.1.3.0 gatk CombineGVCFs \
		--reference {reference}.fasta {gvcfs} -O analysis/{reference}/cohort.g.vcf.gz
		singularity exec docker://broadinstitute/gatk:4.1.3.0 gatk GenotypeGVCFs \
		--reference {reference}.fasta -V analysis/{reference}/cohort.g.vcf.gz -O analysis/{reference}/genotyped_cohort.vcf.gz
		""" 

#rule calculate_IBD:
#	input:
#		"genotyped_cohort.vcf.gz"
#	output:
#		

rule phase_variants:
	input:
		combined="analysis/{reference}/genotyped_cohort.vcf.gz"
	output:
		"analysis/{reference}/{sample}_haplotagged.bam"	
	shell:
		"""
		bash phase_variants.sh {input.combined} {wildcards.sample} {reference}
		"""
rule extract_chr:
	input:
		"{sample}_haplotagged.bam"
	output:
		"{sample}_haplotagged_15.bam"
	shell:
		"""
		bash extract_15.sh {wildcards.sample}
		"""
