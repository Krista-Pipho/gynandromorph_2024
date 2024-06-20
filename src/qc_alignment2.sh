samples=("$1","$2","$3","$4")
module load Java/1.8.0_60
module load samtools/1.9

#marking duplicates
for sample in "${samples[@]}"; do 
echo "Processing sample: $sample" 
java -jar -Xmx7g /hpc/home/kp275/0_utils/picard.jar MarkDuplicates INPUT=${sample}.bam OUTPUT=${sample}.dedup.bam METRICS_FILE=metrics.txt  VALIDATION_STRINGENCY=LENIENT; done

#generating bam statistics
for file in */*.dedup.bam; do 
echo $file; samtools flagstat $file; done

#calculating read coverage
echo plotCoverage -b $1/$1.dedup.bam $2/$2.dedup.bam $3/$3.dedup.bam $4/$4.dedup.bam -o sample_coverage.png
