#!/bin/bash 
#SBATCH -N 1 
#SBATCH -p LM 
#SBATCH --mem=128GB 
#SBATCH -t 48:00:00 
#SBATCH --mail-type=ALL 
#SBATCH --mail-user=marisa.lim@stonybrook.edu 

#echo commands to stdout
set -x 

module load perl/5.24.0-threads
module load spades/3.8.1

perl /pylon1/bi4iflp/mlim/SeqCapData/2-GenerateAssembliesPhylo spades -reads /pylon1/bi4iflp/mlim/SeqCapData/readsforspades -out /pylon1/bi4iflp/mlim/SeqCapData/readsforspades/spades_outs -np 1
