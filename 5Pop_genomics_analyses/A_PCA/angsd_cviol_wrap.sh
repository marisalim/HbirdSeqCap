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

#testing
#angsd -bam /pylon5/bi4iflp/mlim/SeqCapData/angsd_analysis/cviolbams.list -ref /pylon5/bi4iflp/mlim/SeqCapData/species_ref_files/Cviolifer_reffasta/ref_cviol_targetedRegionAndFlanking.fasta -anc /pylon5/bi4iflp/mlim/SeqCapData/species_ref_files/Cviolifer_reffasta/ref_cviol_targetedRegionAndFlanking.fasta -only_proper_pairs 0 -minMapQ 0 -minQ 20 -fold 1 -GL 1 -doGeno 2 -doPost 1 -postCutoff 0.95 -SNP_pval 1e-6 -doMaf 2 -doMajorMinor 1 -doSaf 1 -sites /pylon5/bi4iflp/mlim/SeqCapData/Novoalign_outs/cviol.keep -rf /pylon5/bi4iflp/mlim/SeqCapData/Novoalign_outs/cviol.rf -out /pylon5/bi4iflp/mlim/SeqCapData/angsd_analysis/cviol_geno2_outindexed

# from Ke
# minInd = 80% of individuals (rounded) = 50
# setMinDepth = minInd * 5 = 250
angsd -bam /pylon5/bi4iflp/mlim/SeqCapData/angsd_analysis/cviolbams.list -ref /pylon5/bi4iflp/mlim/SeqCapData/species_ref_files/Cviolifer_reffasta/ref_cviol_targetedRegionAndFlanking.fasta -anc /pylon5/bi4iflp/mlim/SeqCapData/species_ref_files/Cviolifer_reffasta/ref_cviol_targetedRegionAndFlanking.fasta -only_proper_pairs 0 -minMapQ 10 -minQ 20 -fold 1 -GL 1 -doGeno 32 -doPost 1 -doMaf 1 -doMajorMinor 1 -doSaf 1 -sites /pylon5/bi4iflp/mlim/SeqCapData/Novoalign_outs/cviol.keep -rf /pylon5/bi4iflp/mlim/SeqCapData/Novoalign_outs/cviol.rf -out /pylon5/bi4iflp/mlim/SeqCapData/angsd_analysis/ngsToolsPCA/cviol_geno_ngstoolsPCA_out -skipTriallelic 1 -minInd 50 -setMinDepth 250 -doCounts 1 -SNP_pval 1e-3 
