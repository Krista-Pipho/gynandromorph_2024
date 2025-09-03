#!/bin/bash

REF=$1
r1=$2
r2=$3
sample=$4

module load BWA/0.7.17
module load samtools/1.9
module load Java/1.8.0_60

# Set sample ID and Name. If the readgroup is left as defailt, downstream programs will complain
RG="@RG\tID:${sample}\tSM:${sample}\t" 

# Aligning reads
echo "Align sample: $sample"
bwa mem -M ${REF} ${r1} ${r2} -R ${RG} | samtools sort -o ${sample}.bam
samtools index ${sample}.bam

# Marking duplicates
echo "MarkDup sample: $sample"
java -jar -Xmx7g picard.jar MarkDuplicates -INPUT ${sample}.bam -OUTPUT ${sample}.dedup.bam -METRICS_FILE metrics.txt -VALIDATION_STRINGENCY LENIENT
samtools index ${sample}.dedup.bam
echo "Output file: ${sample}.dedup.bam"
