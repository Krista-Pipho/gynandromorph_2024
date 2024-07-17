#!/bin/sh
#SBATCH --mem=32G

input_vcf=$1
output_file="/work/sm997/gynandromorph_2024/testing/homozygous_vcfs/output.vcf"

module load bcftools/1.4

bcftools view --genotype hom $input_vcf -o $output_file
echo "bcftools command completed successfully. Output saved to '$output_file'."
