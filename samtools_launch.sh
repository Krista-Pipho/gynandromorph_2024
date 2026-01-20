#!/bin/bash
#SBATCH --mem=200G
#SBATCH -c64
#SBATCH --partition=scavenger

#for file in ../../*dedup.bam; do ./nQuire create -b $file -o ${file%dedup.bam}_nquire_output -x; 
#./nQuire denoise -o ${file%dedup.bam}_nquire_denoised ${file%dedup.bam}_nquire_output.bin;
#./nQuire view ${file%dedup.bam}_nquire_output.bin -a $file > ${file%dedup.bam}_nquire.txt;
#done
module load samtools
module load Java/23.0.1
#samtools merge --threads 64 -o ${1}/merged.bam ${1}/*.dedup.bam
#samtools view -H ${1}/merged.bam | grep '^@RG'
#java -jar /work/kp275/permuting_hifiasm/all_higher_values_with_telo/gynandromorph_2024/picard.jar AddOrReplaceReadGroups  -I ${1}/merged.bam -O ${1}/merged_rg.bam -RGID 1 -RGSM sample1 -RGPL illumina -RGLB lib1 -RGPU unit1
#samtools view -H ${1}/merged_rg.bam | grep '^@RG'
#samtools sort -@ 64 ${1}/merged_rg.bam -o ${1}/merged_sorted.bam
#samtools index ${1}/merged_sorted.bam
mkdir ${1}/hmr_contig_bams
for chr in Hmr{0,1}{0..9} Hmr20 Hmr21 HmrW; do samtools view -@ 64 -b ${1}/merged_sorted.bam $chr > ${1}/hmr_contig_bams/${chr}.bam; done