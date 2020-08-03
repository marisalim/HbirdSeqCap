#!/bin/bash 
#SBATCH -N 1 
#SBATCH -p LM 
#SBATCH --mem=128GB 
#SBATCH -t 48:00:00 
#SBATCH --mail-type=ALL 
#SBATCH --mail-user=marisa.lim@stonybrook.edu 

#echo commands to stdout
set -x 

module load python/2.7.11_gcc

python /pylon5/bi4iflp/mlim/SeqCapData/Novoalign_outs/cviol_removesites.py
