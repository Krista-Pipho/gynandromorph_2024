#!/bin/bash
#SBATCH --mem=200G
#SBATCH -c64
#SBATCH --partition=scavenger
#SBATCH --time=04-00:00:00

#snakemake --rulegraph | dot -Tpng > rulegraph.png
#snakemake
pixi run snakemake --cores 62 -k
