#!/bin/sh
#SBATCH --mem=16G

leg=$1
module load Java/1.8.0_60
module load samtools/1.9

#marking duplicates
echo "Processing sample: $leg"
java -jar -Xmx7g picard.jar MarkDuplicates -INPUT ${leg}.bam -OUTPUT ${leg}.dedup.bam -METRICS_FILE metrics.txt -VALIDATION_STRINGENCY LENIENT

echo "Output file: ${leg}.dedup.bam"

#generating bam statistics
samtools flagstat ${leg}.dedup.bam
samtools index ${leg}.dedup.bam
