#!/bin/bash 
#SBATCH -N 1 
#SBATCH -p LM 
#SBATCH --mem=128GB 
#SBATCH -t 24:00:00 
#SBATCH --mail-type=ALL 
#SBATCH --mail-user=marisa.lim@stonybrook.edu 

#echo commands to stdout
set -x 

module load R/3.3.3-mkl

R --slave < /pylon5/bi4iflp/mlim/SeqCapData/angsd_analysis/Fst_dxy_theta/dxy/snpcallgenofiles_fordxy/make_nonvar_matrix.r
