#!/bin/bash
#SBATCH --mem=200G
#SBATCH -c64
#SBATCH --partition=scavenger

#for file in ../../*dedup.bam; do ./nQuire create -b $file -o ${file%dedup.bam}_nquire_output -x; 
#./nQuire denoise -o ${file%dedup.bam}_nquire_denoised ${file%dedup.bam}_nquire_output.bin;
#./nQuire view ${file%dedup.bam}_nquire_output.bin -a $file > ${file%dedup.bam}_nquire.txt;
#done
module load samtools
#samtools sort all_samp.bam -o all_samp_sorted.bam
#samtools index all_samp_sorted.bam
for chr in Hmr{0,1}{0..9} Hmr20 Hmr21 HmrW; do samtools view -b all_samp_sorted.bam $chr > ${chr}.bam; done