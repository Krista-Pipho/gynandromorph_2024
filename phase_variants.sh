#!/bin/bash

module load bcftools/1.4

combined=$1
sample=$2
ref=$3

bcftools view -s $sample $combined > ${sample}.vcf
bgzip ${sample}.vcf
tabix ${sample}.vcf.gz

whatshap phase -o ${sample}_phased.vcf --reference $ref ${sample}.vcf.gz ${sample}.dedup.bam 
bgzip ${sample}_phased.vcf
tabix ${sample}_phased.vcf.gz

whatshap haplotag -o ${sample}_haplotagged.bam --reference $ref ${sample}_phased.vcf.gz ${sample}.dedup.bam

