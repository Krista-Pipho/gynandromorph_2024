#!/bin/bash
#SBATCH --mem=100G
#SBATCH --cpus-per-task 24
#SBATCH --partition scavenger

#pixi run bamCoverage -b TN2L_S1.dedup.bam -o S1_coverage.bw
pixi run bigWigToBedGraph S1_coverage.bw S1_coverage.bedGraph

awk '{
    if ($4 >= 0 && $4 <= 1) {
        category = "gray"
    } else if ($4 >= 2 && $4 <= 15) {
        category = "white"
    } else if ($4 >= 16 && $4 <= 30) {
        category = "yellow"
    } else if ($4 >= 31 && $4 <= 45) {
        category = "orange"
    } else {
        category = "red"
    }
    print $1, $2, $3, $4, category
}' S1_coverage.bedGraph > S1_coverage_categorized.bedGraph