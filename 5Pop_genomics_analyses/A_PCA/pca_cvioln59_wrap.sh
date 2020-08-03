#!/bin/bash 
#SBATCH -N 1 
#SBATCH -p LM 
#SBATCH --mem=128GB 
#SBATCH -t 12:00:00 
#SBATCH --mail-type=ALL 
#SBATCH --mail-user=marisa.lim@stonybrook.edu 

#echo commands to stdout
set -x 

source ~/.bashrc

# from Ke
# minInd = 70% of individuals (rounded) = 41
angsd -bam /pylon5/bi4iflp/mlim/SeqCapData/angsd_analysis/cviolbams_n59.list -ref /pylon5/bi4iflp/mlim/SeqCapData/species_ref_files/Cviolifer_reffasta/ref_cviol_targetedRegionAndFlanking.fasta -anc /pylon5/bi4iflp/mlim/SeqCapData/species_ref_files/Cviolifer_reffasta/ref_cviol_targetedRegionAndFlanking.fasta -only_proper_pairs 0 -minMapQ 10 -minQ 20 -fold 1 -GL 1 -doGeno 32 -doPost 1 -doMaf 1 -doMajorMinor 1 -doSaf 1 -sites /pylon5/bi4iflp/mlim/SeqCapData/Novoalign_outs/cviol.keep -rf /pylon5/bi4iflp/mlim/SeqCapData/Novoalign_outs/cviol.rf -out /pylon5/bi4iflp/mlim/SeqCapData/angsd_analysis/ngsToolsPCA/cviol_forPCAn59_out -skipTriallelic 1 -minInd 41 -minIndDepth 3 -doCounts 1
