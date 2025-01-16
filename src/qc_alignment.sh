#!/bin/sh
#SBATCH --mem=16G

module load Java/1.8.0_60
module load samtools/1.9
# Get flagstats and create processed summary file
samtools flagstat ${1}.dedup.bam > ${1}.flagstat.txt
# Process summary by pulling out four numbers of interest and tagging them with sample of origin and identity of statistic
paste --delimiter="\t" <(echo $1) <(echo "total") <(cat ${1}.flagstat.txt | grep total | cut -d " " -f1) > ${1}.processed.flagstat.txt
paste --delimiter="\t" <(echo $1) <(echo "mapped") <(cat ${1}.flagstat.txt | grep "mapped (" | cut -d " " -f1) >> ${1}.processed.flagstat.txt
paste --delimiter="\t" <(echo $1) <(echo "paired") <(cat ${1}.flagstat.txt | grep proper | cut -d " " -f1) >> ${1}.processed.flagstat.txt
paste --delimiter="\t" <(echo $1) <(echo "duplicate") <(cat ${1}.flagstat.txt | grep duplicate | cut -d " " -f1) >> ${1}.processed.flagstat.txt

echo "Generated Statistics Summary ${1}.processed.flagstat.txt"

# Get insertsize metrics from GATK
#singularity exec -B /hpc/group/wraylab:/hpc/group/wraylab docker://broadinstitute/gatk:4.1.3.0 gatk CollectInsertSizeMetrics -R /hpc/group/wraylab/kp275/Butterfly/Heliconius_melpomene_melpomene_Hmel2.5.scaffolds.fa -I ${1}.dedup.bam -O ${1}.insert_size_metrics.txt -H ${1}_hist.txt

# Extract relevant data
#sed -n '/insert_size/,$ p' ${1}.insert_size_metrics.txt > ${1}.insert_size.txt
