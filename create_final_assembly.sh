#!/bin/bash
#SBATCH --mem=100G
#SBATCH --cpus-per-task 24
#SBATCH --partition scavenger

grep -v -A 1 -f ids_to_remove.txt input.fasta > output.fasta

# Change chromosome naming scheme
sed s/"Hmel2"/"Hmr"/g $2 > $1
sed -i s/"o_RagTag"//g $1
# Add mitochondrial scaffold
#cat hmr.mito.ctg.fasta >> hmr.gfa
#sed -i s/">ctg000001l"/">HmrMt001"/g hmr.gfa
# Change scaffold numbering
sed -i s/"Hmr03003"/"Hmr03001"/g $1
sed -i s/"Hmr14004"/"Hmr14001"/g $1
sed -i s/"Hmr15003"/"Hmr15001"/g $1
sed -i s/"Hmr16002"/"Hmr16001"/g $1
sed -i s/"Hmr18002"/"Hmr18001"/g $1
sed -i s/"Hmr18003"/"Hmr18002"/g $1
sed -i s/"Hmr20003"/"Hmr20002"/g $1
sed -i s/"ptg000080l"/"HmrW001"/g $1
sed -i s/"ptg000083l"/"HmrW002"/g $1
