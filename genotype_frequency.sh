

#!/bin/bash
#SBATCH --mem=50G
#SBATCH -c10
#SBATCH --partition=scavenger

module load VCFtools/0.1.17
vcftools --gzvcf filtered_genotyped_cohort.vcf.gz.vcf_4 --minQ 200 --recode --recode-INFO-all --out filtered_joint_quality_200.vcf
vcftools --gzvcf filtered_genotyped_cohort.vcf.gz.vcf_4 --minQ 800 --recode --recode-INFO-all --out filtered_joint_quality_800.vcf
pixi run plink --vcf filtered_joint_quality_50.vcf.recode.vcf --freqx --allow-extra-chr --out filtered_joint_quality_50
pixi run plink --vcf filtered_joint_quality_200.vcf.recode.vcf --freqx --allow-extra-chr --out filtered_joint_quality_200
pixi run plink --vcf filtered_joint_quality_800.vcf.recode.vcf --freqx --allow-extra-chr --out filtered_joint_quality_800
cat filtered_joint_quality_50.frqx | cut -f 5-7 | sort | uniq -c
cat filtered_joint_quality_200.frqx | cut -f 5-7 | sort | uniq -c
cat filtered_joint_quality_800.frqx | cut -f 5-7 | sort | uniq -c
