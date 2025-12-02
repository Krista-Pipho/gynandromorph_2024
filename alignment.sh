#!/bin/bash

REF=$1
r1=$2
r2=$3
sample=$4
cores=$5

# Set sample ID and Name. If the readgroup is left as defailt, downstream programs will complain
RG="@RG\tID:${sample}\tSM:${sample}" 

# Aligning reads
echo "Align sample: $sample"
singularity exec -B $(pwd) docker://biocontainers/bwa:v0.7.17_cv1 bwa mem -t ${cores} -M ${REF}.fasta ${r1} ${r2} -R ${RG} | samtools sort -o analysis/${REF}/${sample}.bam
samtools index analysis/${REF}/${sample}.bam


# Marking duplicates
echo "MarkDup sample: $sample"
java -jar -Xmx7g picard.jar MarkDuplicates -INPUT analysis/${REF}/${sample}.bam -OUTPUT analysis/${REF}/${sample}.dedup.bam -METRICS_FILE analysis/${REF}/metrics.txt -VALIDATION_STRINGENCY LENIENT
samtools index analysis/${REF}/${sample}.dedup.bam
echo "Output file: analysis/${REF}/${sample}.dedup.bam"
