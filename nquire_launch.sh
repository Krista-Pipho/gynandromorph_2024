#!/bin/bash
#SBATCH --mem=200G
#SBATCH -c64
#SBATCH --partition=scavenger

for file in hmr_contig_bams/*; do ./nQuire create -b $file -o ${file%.bam}_nquire; 
#./nQuire denoise -o ${file%dedup.bam}_nquire_denoised ${file%dedup.bam}_nquire_output.bin;
#./nQuire view ${file%.bam}_nquire.bin -a $file > ${file%dedup.bam}_nquire.txt;
#done
./nQuire lrdmodel ${file%.bam}_nquire.bin;
done