#!/bin/bash
#SBATCH --mem=200G
#SBATCH -c64
#SBATCH --partition=scavenger

snakemake --rulegraph | dot -Tpng > rulegraph.png
snakemake
