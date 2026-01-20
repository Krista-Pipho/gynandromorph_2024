#!/bin/bash
#SBATCH --mem=200G
#SBATCH -c64
#SBATCH --partition=scavenger

#pixi run python ./Genotype_Likelihoods.py -p ploidy.list -d .1 -m .05 -r hmr_new_official_masked.fasta new_names.filelist
pixi run Rscript ./HMMploidy.R fileList=up_to_10.filelist wind=1000 maxPloidy=8