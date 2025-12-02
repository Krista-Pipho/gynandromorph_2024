#!/bin/bash
#SBATCH --mem=100G
#SBATCH --cpus-per-task 24
#SBATCH --partition scavenger

#replace read group information with the correct values

#singularity exec docker://broadinstitute/gatk:4.1.3.0 gatk AddOrReplaceReadGroups -I ${1}.dedup.bam -O ${1}.dedup.RG.bam -RGID $1 -RGLB lib1 -RGSM $1 -RGPL ILLUMINA -RGPU unit1

#samtools index ${1}.dedup.RG.bam

echo singularity exec docker://broadinstitute/gatk:4.1.3.0 gatk HaplotypeCallerSpark --input ./${1} --output ${2} --reference ${3}.fasta  -ERC GVCF
singularity exec docker://broadinstitute/gatk:4.1.3.0 gatk HaplotypeCallerSpark --input ./${1} --output ${2} --reference ${3}.fasta  -ERC GVCF
#echo changing sample name of $1
#sed -i s/Seq01/${1}/g ${1}_spark.g.vcf; done

echo indexing $2
tabix --force -p vcf ${2}

