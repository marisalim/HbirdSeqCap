#!/bin/bash 
#SBATCH -N 1 
#SBATCH -p LM 
#SBATCH --mem=128GB 
#SBATCH -t 24:00:00 
#SBATCH --mail-type=ALL 
#SBATCH --mail-user=marisa.lim@stonybrook.edu 

#echo commands to stdout
set -x 

module load fastqc/0.11.3

for file in /pylon1/bi4iflp/mlim/SeqCapData/ScrubReads_Results/runfastqc/*.fq; do	fastqc "$file" -o /pylon1/bi4iflp/mlim/SeqCapData/ScrubReads_Results/evaluation; done

rm /pylon1/bi4iflp/mlim/SeqCapData/ScrubReads_Results/evaluation/*fastqc.zip
