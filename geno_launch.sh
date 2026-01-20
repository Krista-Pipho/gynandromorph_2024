#!/bin/bash
#SBATCH --mem=200G
#SBATCH -c64
#SBATCH --partition=scavenger

# Make intermediate files for genomescope analysis using Jellyfish
# Find more information here. -m is kmer length, -s is memory, -t is threads
pixi run singularity exec -B $(pwd) docker://biodckrdev/jellyfish:2.2.3 jellyfish count -C -m 22 -s 1000000000 -t 60 ${1}/*.fastq -o ${1}/combined_gynandromorph.jf
pixi run singularity exec -B $(pwd) docker://biodckrdev/jellyfish:2.2.3 jellyfish histo -t 10 ${1}/combined_gynandromorph.jf > ${1}/combined_gynandromorph.histo

# Run genomescope
pixi run genomescope2 -i ${1}/combined_gynandromorph.histo -o ${1}/genomescope_combined -k 22 #kmer_size