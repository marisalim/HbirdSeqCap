#!/bin/bash 
#SBATCH -N 1 
#SBATCH -p LM 
#SBATCH --mem=128GB 
#SBATCH -t 48:00:00 
#SBATCH --mail-type=ALL 
#SBATCH --mail-user=marisa.lim@stonybrook.edu 

#echo commands to stdout
set -x 

source ~/.bashrc

# from Ke
# minInd = 70% of individuals (rounded) = 71
angsd -bam /pylon5/bi4iflp/mlim/SeqCapData/angsd_analysis/ccorubams.list -ref /pylon5/bi4iflp/mlim/SeqCapData/species_ref_files/Ccoruscans_reffasta/ref_ccoru_targetedRegionAndFlanking.fasta -anc /pylon5/bi4iflp/mlim/SeqCapData/species_ref_files/Ccoruscans_reffasta/ref_ccoru_targetedRegionAndFlanking.fasta -only_proper_pairs 0 -minMapQ 10 -minQ 20 -fold 1 -GL 1 -doGeno 32 -doPost 1 -doMaf 1 -doMajorMinor 1 -doSaf 1 -sites /pylon5/bi4iflp/mlim/SeqCapData/Novoalign_outs/ccoru.keep -rf /pylon5/bi4iflp/mlim/SeqCapData/Novoalign_outs/ccoru.rf -out /pylon5/bi4iflp/mlim/SeqCapData/angsd_analysis/ngsRelate/ccoru_relatetest_out -skipTriallelic 1 -minInd 71 -minIndDepth 3 -doCounts 1 -SNP_pval 1e-6 -minmaf 0.05 -doGlf 3
