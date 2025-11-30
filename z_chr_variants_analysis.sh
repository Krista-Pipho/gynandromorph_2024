#!/bin/bash

sbatch partial_filtering.sh genotyped_cohort.vcf.gz
bgzip filtered_genotyped_cohort.vcf.gz.vcf_3
singularity exec -B $(pwd) docker://biocontainers/bcftools:v1.9-1-deb_cv1 bcftools index filtered_genotyped_cohort.vcf.gz.vcf_3.gz
singularity exec -B $(pwd) docker://biocontainers/bcftools:v1.9-1-deb_cv1 bcftools view filtered_genotyped_cohort.vcf.gz.vcf_3.gz --regions Hmel221001o_RagTag > z_chr.vcf
singularity exec -B $(pwd) docker://biocontainers/bcftools:v1.9-1-deb_cv1 bcftools stats z_chr.vcf > z_chr_stats.txt 
singularity exec -B $(pwd) docker://biocontainers/bcftools:v1.9-1-deb_cv1 bcftools view --min-alleles 3 z_chr.vcf > z_chr_multi.vcf

singularity exec -B $(pwd) docker://biocontainers/bcftools:v1.9-1-deb_cv1 bcftools query -f '%AF\n' z_chr_multi.vcf | sort | uniq -c
singularity exec -B $(pwd) docker://biocontainers/bcftools:v1.9-1-deb_cv1 bcftools query -f '%ALT\n' z_chr_multi.vcf | awk '{print split($1, a, ",")}' | sort | uniq -c

#singularity exec -B $(pwd) docker://biocontainers/bcftools:v1.9-1-deb_cv1 bcftools query -f '%AF\n' z_chr_multi_q800.vcf.recode.vcf | sort | uniq -c
#ingularity exec -B $(pwd) docker://biocontainers/bcftools:v1.9-1-deb_cv1 bcftools query -f '%ALT\n' z_chr_multi_q800.vcf.recode.vcf | awk '{print split($1, a, ",")}' | sort | uniq -c