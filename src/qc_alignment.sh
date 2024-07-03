#!/bin/sh
#SBATCH --mem=16G

leg=$1
module load Java/1.8.0_60
module load samtools/1.9

#marking duplicates
echo "Processing sample: $leg"
java -jar -Xmx7g /hpc/group/wraylab/kp275/team_heliconius/picard.jar MarkDuplicates -INPUT $leg -OUTPUT ${leg}.dedup.bam -METRICS_FILE metrics.txt -VALIDATION_STRINGENCY LENIENT


file="${leg}.dedup.bam"


echo "Output file: ${leg}.dedup.bam"

#generating bam statistics
samtools flagstat $leg
samtools index $leg
