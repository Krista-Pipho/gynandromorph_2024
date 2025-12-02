#!/bin/bash

module load BWA/0.7.17
module load samtools/1.9
module load Java/1.8.0_60

REF=${1}

bwa index ${REF}.fasta
samtools faidx ${REF}.fasta
java -jar picard.jar CreateSequenceDictionary -R ${REF}.fasta -O ${REF}.dict
touch analysis/${REF}/ref.prepped
