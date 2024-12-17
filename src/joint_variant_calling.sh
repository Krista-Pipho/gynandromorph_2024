#!/bin/bash
#SBATCH --mem=100G
#SBATCH --cpus-per-task 24
#SBATCH --partition scavenger

singularity exec -B /hpc/group/wraylab:/hpc/group/wraylab docker://broadinstitute/gatk:4.1.3.0 gatk CombineGVCFs --reference /hpc/group/wraylab/kp275/Butterfly/Heliconius_melpomene_melpomene_Hmel2.5.scaffolds.fa --variant ${1}_spark.g.vcf.gz --variant ${2}_spark.g.vcf.gz --variant ${3}_spark.g.vcf.gz --variant ${4}_spark.g.vcf.gz -O cohort.g.vcf.gz
singularity exec -B /hpc/group/wraylab:/hpc/group/wraylab docker://broadinstitute/gatk:4.1.3.0 gatk GenotypeGVCFs --reference /hpc/group/wraylab/kp275/Butterfly/Heliconius_melpomene_melpomene_Hmel2.5.scaffolds.fa -V cohort.g.vcf.gz -O genotyped_cohort.vcf.gz

