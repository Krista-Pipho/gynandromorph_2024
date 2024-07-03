#!/bin/sh
#SBATCH --mem=16G

module load bcftools/1.4

bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%TYPE\n' /work/sm997/gynandromorph_2024/raw_data/gynandromorph-combined.vcf.gz | grep 'S
NP' > snps.txt
