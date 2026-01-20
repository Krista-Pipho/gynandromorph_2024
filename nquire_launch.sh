#!/bin/bash
#SBATCH --mem=200G
#SBATCH -c64
#SBATCH --partition=scavenger

for file in hmr_contig_bams/*.bam; do ./nQuire create -b $file -o ${file%.bam}_nquire; 
./nQuire denoise -o ${file%.bam}_nquire_denoised ${file%.bam}_nquire.bin;
#./nQuire view ${file%.bam}_nquire.bin -a $file > ${file%dedup.bam}_nquire.txt;
#done
./nQuire lrdmodel ${file%.bam}_nquire_denoised.bin;
done