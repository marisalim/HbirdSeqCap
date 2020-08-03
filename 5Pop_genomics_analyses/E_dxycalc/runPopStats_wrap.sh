#!/bin/bash
#SBATCH -N 1
#SBATCH -p LM
#SBATCH --mem=128GB
#SBATCH -t 24:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=marisa.lim@stonybrook.edu

#echo commands to stdout
set -x

module load perl/5.18.4-threads

perl 6-PopStats Dxyp -f /pylon5/bi4iflp/mlim/SeqCapData/angsd_analysis/Fst_dxy_theta/dxy/snpcallgenofiles_fordxy/genofiles_forPopStats/cviol_dxy -u cviol_dxy

perl 6-PopStats Dxyp -f /pylon5/bi4iflp/mlim/SeqCapData/angsd_analysis/Fst_dxy_theta/dxy/snpcallgenofiles_fordxy/genofiles_forPopStats/ccoru_dxy -u ccoru_dxy
