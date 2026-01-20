#!/bin/bash
#SBATCH --mem=200G
#SBATCH -c64
#SBATCH --partition=scavenger

#conda activate smudgeplot
FastK -v -t4 -k22 -M100 -T60 ${1}/*.fastq -N${1}/FastK_Table
smudgeplot hetmers -L 12 -t 60 -o ${1}/kmerpairs --verbose ${1}/FastK_Table
smudgeplot all -o ${1}/output ${1}/kmerpairs.smu