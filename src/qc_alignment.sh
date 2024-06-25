#!/bin/sh
#SBATCH --mem=16G

sample=$1
module load Java/1.8.0_60
module load samtools/1.9

#marking duplicates
echo "Processing sample: $sample"
java -jar -Xmx7g /hpc/group/wraylab/kp275/team_heliconius/picard.jar MarkDuplicates -INPUT $sample -OUTPUT ${sample}.dedup.bam -METRICS_FILE metrics.txt -VALIDATION_STRINGENCY LENIENT


file="${sample}.dedup.bam"


echo "Output file: $file"

#generating bam statistics
samtools flagstat $file
