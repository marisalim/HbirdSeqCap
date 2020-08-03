#!/bin/bash 
#SBATCH -N 1 
#SBATCH -p LM 
#SBATCH --mem=128GB 
#SBATCH -t 2:00:00 
#SBATCH --mail-type=ALL 
#SBATCH --mail-user=marisa.lim@stonybrook.edu 

#echo commands to stdout
set -x 

module load python/2.7.11_gcc

python /pylon5/bi4iflp/mlim/SeqCapData/angsd_analysis/selection_tests/ccoru_MAFfiltercode.py

sed 's/\"//g' MAFfilteredgeno_ccoru97_out.csv | sed '/^\s*$/d' | cut -f2- > MAFfiltered_ccoru97_FINAL.txt
