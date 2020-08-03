#!/bin/bash 
#SBATCH -N 1 
#SBATCH -p LM 
#SBATCH --mem=128GB 
#SBATCH -t 48:00:00 
#SBATCH --mail-type=ALL 
#SBATCH --mail-user=marisa.lim@stonybrook.edu 

#echo commands to stdout
set -x 

module load perl/5.18.4-threads

perl /pylon5/bi4iflp/mlim/SeqCapData/Novoalign_outs/SNPcleaner225.pl -u 3 -k 71 -d 214 -a 0 -H 0.0001 -h 0 -B /pylon5/bi4iflp/mlim/SeqCapData/Novoalign_outs/ccoru_sites_to_keep.bed -p /pylon5/bi4iflp/mlim/SeqCapData/Novoalign_outs/ccoru_sites_dumped.bed -v /pylon5/bi4iflp/mlim/SeqCapData/Novoalign_outs/ccoru.raw.vcf
