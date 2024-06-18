#!/bin/bash
#SBATCH --mem=200G
#SBATCH -c32
#SBATCH --mail-type=END
#SBATCH --mail-user=shriya.minocha@duke.edu

REF=$1
r1=$2
r2=$3
sample=$4

module load BWA/0.7.17
module load samtools/1.9

bwa index ${REF}
samtools faidx ${REF}

bwa mem -M -t 32 "${REF}" ${r1} ${r2} | samtools sort -@32 -o ${sample}.bam

samtools index ${sample}.bam -@32
