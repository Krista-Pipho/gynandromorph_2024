#!/bin/sh

input_vcf=$1
output_file="filtered_${input_vcf##*/}.vcf"
echo filtering $input_vcf ...

# Generate stats before filtering
singularity exec -B $(pwd) docker://biocontainers/bcftools:v1.9-1-deb_cv1 bcftools stats $input_vcf > before_filtering_stats.txt


# Keep only SNPs (remove indels)
singularity exec -B $(pwd) docker://biocontainers/bcftools:v1.9-1-deb_cv1 bcftools filter -e 'TYPE!="snp"' $input_vcf > $output_file
singularity exec -B $(pwd) docker://biocontainers/bcftools:v1.9-1-deb_cv1 bcftools stats $output_file > filtered_vcf_stats.txt
echo successfully removed indels and others, keeping only SNPs

#getting rid of missing positions
singularity exec -B $(pwd) docker://biocontainers/bcftools:v1.9-1-deb_cv1 bcftools view --genotype ^miss $output_file -o ${output_file}_2
singularity exec -B $(pwd) docker://biocontainers/bcftools:v1.9-1-deb_cv1 bcftools stats ${output_file}_2 > remove_missing_vcf_stats.txt
echo successfully removed missing position type

#remove fixed differences from reference
singularity exec -B $(pwd) docker://biocontainers/bcftools:v1.9-1-deb_cv1 bcftools view --max-af 0.9:alt1 ${output_file}_2 -o ${output_file}_3
singularity exec -B $(pwd) docker://biocontainers/bcftools:v1.9-1-deb_cv1 bcftools stats ${output_file}_3 > remove_fixeddifs_vcf_stats.txt
echo successfully removed fixed differences from reference

#remove multiallelic regions
singularity exec -B $(pwd) docker://biocontainers/bcftools:v1.9-1-deb_cv1 bcftools view --max-alleles 2 ${output_file}_3 -o ${output_file}_4
singularity exec -B $(pwd) docker://biocontainers/bcftools:v1.9-1-deb_cv1 bcftools stats ${output_file}_4 > remove_multiallelic_vcf_stats.txt
echo successfully removed multiallelic regions
