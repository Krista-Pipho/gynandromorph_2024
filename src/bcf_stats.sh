#!/bin/sh
#SBATCH --mem=16G

input_vcf=$1

module load bcftools/1.4

bcftools stats $input_vcf > vcf_stats.txt 
plot-vcfstats vcf_stats.txt -p /work/sm997/gynandromorph_2024/src/plots

