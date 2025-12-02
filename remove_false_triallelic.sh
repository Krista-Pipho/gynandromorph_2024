#singularity exec -B $(pwd) docker://biocontainers/bcftools:v1.9-1-deb_cv1 bcftools query -f '%AF\n' z_chr_snp_and_indel_q800.vcf | awk -F ',' '{sum=$1+$2; print sum}' | sort | uniq -c
singularity exec -B $(pwd) docker://biocontainers/bcftools:v1.9-1-deb_cv1 bcftools query -f '%AF\n' z_chr_snp_and_indel_q800.vcf | awk -F',' '{sum=$1+$2; print sum}' > alt_allele_freq_sum.txt
grep -v '^[[:space:]]*#' z_chr_snp_and_indel_q800.vcf > headerless_z_chr_snp_and_indel_q800.vcf
paste alt_allele_freq_sum.txt headerless_z_chr_snp_and_indel_q800.vcf > sum_headerless_z_chr_snp_and_indel_q800.vcf
grep -v '^[[:space:]]*1' sum_headerless_z_chr_snp_and_indel_q800.vcf > filtered_sum_headerless_z_chr_snp_and_indel_q800.vcf
cut -f2- -d$'\t' filtered_sum_headerless_z_chr_snp_and_indel_q800.vcf > filtered_headerless_z_chr_snp_and_indel_q800.vcf
grep '^[[:space:]]*#' z_chr_snp_and_indel_q800.vcf > header_z_chr_snp_and_indel_q800.vcf
cat header_z_chr_snp_and_indel_q800.vcf filtered_headerless_z_chr_snp_and_indel_q800.vcf > filtered_z_chr_snp_and_indel_q800.vcf
cat filtered_headerless_z_chr_snp_and_indel_q800.vcf | cut -f 1-2 | awk 'BEGIN{OFS="\t"}{sum=$2+1; print $1,$2,sum, "tri_allelic_site","black"}' > karyotype_trial