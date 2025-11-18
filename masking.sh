#!/bin/bash
#SBATCH --mem=150G
#SBATCH -c62
#SBATCH --partition=scavenger

assembly="bcov_scaffold"

#mkdir ${assembly}
#pixi run singularity exec -B $(pwd) docker://dfam/tetools:latest BuildDatabase -name "${assembly}/${assembly}_db" ${assembly}.gfa
#pixi run singularity exec -B $(pwd) docker://dfam/tetools:latest RepeatModeler -database ${assembly}/${assembly}_db -engine ncbi -threads 60
pixi run singularity exec -B $(pwd) docker://dfam/tetools:latest RepeatMasker -pa 60 -gff -lib RM_3824175.FriNov142049232025/consensi.fa.classified -dir ${assembly}_MaskerOutput ${assembly}.gfa
