#!/bin/bash 

module load samtools/1.9
samtools index ${1}_haplotagged.bam
samtools view ${1}_haplotagged.bam Hmel215003o -b > ${1}_haplotagged_15.bam
samtools index ${1}_haplotagged_15.bam
