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
module load blat/v35
module load blast/2.2.31
source ~/.bashrc 

perl /pylon1/bi4iflp/mlim/SeqCapData/3-FindingTargetsV9 combineExon -t /pylon1/bi4iflp/mlim/SeqCapData/targeted_loci.fasta -a /pylon1/bi4iflp/mlim/SeqCapData/raw_assembly_dir -p 80 -b 1 -e 4 -f 500

