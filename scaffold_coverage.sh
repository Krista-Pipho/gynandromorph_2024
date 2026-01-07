#!/bin/bash
#SBATCH --mem=50G
#SBATCH --cpus-per-task=24
#SBATCH --partition scavenger
#SBATCH --array=0-250

module load Java/1.8.0_60
module load samtools/1.9

assembly=$1

if [ ! -d "analysis/${assembly}/${assembly}_scaffold_coverage" ]; then
			mkdir analysis/${assembly}/${assembly}_scaffold_coverage
fi

readarray -t scaffoldArr < ${assembly}.fasta.fai

#cut is -d delimiter of space, -f 1 is field 1
scaffoldName=$(echo ${scaffoldArr[$SLURM_ARRAY_TASK_ID]} | cut -d ' ' -f 1)
scaffoldLen=$(echo ${scaffoldArr[$SLURM_ARRAY_TASK_ID]} | cut -d ' ' -f 2)
samplePoints=$((scaffoldLen/100))
ssamplePoints=${samplePoints%.*}


echo $scaffoldName
plotCoverage --skipZeros -b analysis/${assembly}/TN2L_S1.dedup.bam analysis/${assembly}/TN2R_S2.dedup.bam analysis/${assembly}/TN3L_S3.dedup.bam analysis/${assembly}/TN3R_S4.dedup.bam --plotFile analysis/${assembly}/${assembly}_scaffold_coverage/${scaffoldName}_coverage  --outRawCounts analysis/${assembly}/${assembly}_scaffold_coverage/${scaffoldName}_cov.tab --region $scaffoldName -n $samplePoints &> "analysis/${assembly}/${assembly}_scaffold_coverage/${scaffoldName}_coverageStats.txt"
#plotCoverage -b S1.bc2010.dedup.bam S2.dedup.bam S3.dedup.bam S4.dedup.bam --plotFile utg000001l_coverage --outRawCounts utg000001l_cov.tab --region utg000001l -n 7178

#means    
averages=""

while read sample mean std min p25 p50 p75 max; do
    averages="${averages} ${mean}"
done < analysis/${assembly}/${assembly}_scaffold_coverage/${scaffoldName}_coverageStats.txt

averages=$(echo "$averages" | awk '{print $2 " " $3 " " $4 " " $5}')
touch analysis/${assembly}/${assembly}_scaffold_coverage/scaffoldCov.txt | echo $scaffoldName $averages | tr ' ' '\t' >> analysis/${assembly}/${assembly}_scaffold_coverage/scaffoldCov.txt

