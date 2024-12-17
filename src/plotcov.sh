#!/bin/sh

module load Java/1.8.0_60
module load samtools/1.9

samtools index $1
samtools index $2
samtools index $3
samtools index $4
plotCoverage -b $1 $2 $3 $4 -o coverage_plot.png
