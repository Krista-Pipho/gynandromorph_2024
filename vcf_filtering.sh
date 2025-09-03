#!/bin/sh

input_vcf=$1
output_file="filtered_${input_vcf##*/}.vcf"
echo filtering $input_vcf ...

module load VCFtools/0.1.17
module load bcftools/1.4

# Keep only SNPs (remove indels)
bcftools filter -e 'TYPE!="snp"' $input_vcf > $output_file
bcftools stats $output_file > filtered_vcf_stats.txt
echo successfully removed indels and others, keeping only SNPs

#getting rid of missing positions
bcftools view --genotype ^miss $output_file -o ${output_file}_2
bcftools stats ${output_file}_2 > remove_missing_vcf_stats.txt
echo successfully removed missing position type

#remove fixed differences from reference
bcftools view --max-af 0.9:alt1 ${output_file}_2 -o ${output_file}_3
bcftools stats ${output_file}_3 > remove_fixeddifs_vcf_stats.txt
echo successfully removed fixed differences from reference

#remove multiallelic regions
bcftools view --max-alleles 2 ${output_file}_3 -o ${output_file}_4
bcftools stats ${output_file}_4 > remove_multiallelic_vcf_stats.txt
echo successfully removed multiallelic regions
