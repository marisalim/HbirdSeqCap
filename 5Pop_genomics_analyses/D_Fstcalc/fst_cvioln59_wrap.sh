#!/bin/bash 
#SBATCH -N 1 
#SBATCH -p LM 
#SBATCH --mem=128GB 
#SBATCH -t 24:00:00 
#SBATCH --mail-type=ALL 
#SBATCH --mail-user=marisa.lim@stonybrook.edu 

#echo commands to stdout
set -x 

source ~/.bashrc

# minInd 70%
# pop1 n=6 * 0.7 = 4
angsd -bam /pylon5/bi4iflp/mlim/SeqCapData/angsd_analysis/Fst_dxy_theta/cviolbams_pop1.list -ref /pylon5/bi4iflp/mlim/SeqCapData/species_ref_files/Cviolifer_reffasta/ref_cviol_targetedRegionAndFlanking.fasta -anc /pylon5/bi4iflp/mlim/SeqCapData/species_ref_files/Cviolifer_reffasta/ref_cviol_targetedRegionAndFlanking.fasta -only_proper_pairs 0 -minMapQ 10 -minQ 20 -GL 1 -doSaf 1 -sites /pylon5/bi4iflp/mlim/SeqCapData/Novoalign_outs/cviol.keep -rf /pylon5/bi4iflp/mlim/SeqCapData/Novoalign_outs/cviol.rf -out /pylon5/bi4iflp/mlim/SeqCapData/angsd_analysis/Fst_dxy_theta/Fst/cviol_forFstpop1_out -minInd 4 -minIndDepth 3 -doCounts 1

# pop2 n=14 * 0.7 = 10
angsd -bam /pylon5/bi4iflp/mlim/SeqCapData/angsd_analysis/Fst_dxy_theta/cviolbams_pop2.list -ref /pylon5/bi4iflp/mlim/SeqCapData/species_ref_files/Cviolifer_reffasta/ref_cviol_targetedRegionAndFlanking.fasta -anc /pylon5/bi4iflp/mlim/SeqCapData/species_ref_files/Cviolifer_reffasta/ref_cviol_targetedRegionAndFlanking.fasta -only_proper_pairs 0 -minMapQ 10 -minQ 20 -GL 1 -doSaf 1 -sites /pylon5/bi4iflp/mlim/SeqCapData/Novoalign_outs/cviol.keep -rf /pylon5/bi4iflp/mlim/SeqCapData/Novoalign_outs/cviol.rf -out /pylon5/bi4iflp/mlim/SeqCapData/angsd_analysis/Fst_dxy_theta/Fst/cviol_forFstpop2_out -minInd 10 -minIndDepth 3 -doCounts 1

# pop3 n=4 * 0.7 = 3
angsd -bam /pylon5/bi4iflp/mlim/SeqCapData/angsd_analysis/Fst_dxy_theta/cviolbams_pop3.list -ref /pylon5/bi4iflp/mlim/SeqCapData/species_ref_files/Cviolifer_reffasta/ref_cviol_targetedRegionAndFlanking.fasta -anc /pylon5/bi4iflp/mlim/SeqCapData/species_ref_files/Cviolifer_reffasta/ref_cviol_targetedRegionAndFlanking.fasta -only_proper_pairs 0 -minMapQ 10 -minQ 20 -GL 1 -doSaf 1 -sites /pylon5/bi4iflp/mlim/SeqCapData/Novoalign_outs/cviol.keep -rf /pylon5/bi4iflp/mlim/SeqCapData/Novoalign_outs/cviol.rf -out /pylon5/bi4iflp/mlim/SeqCapData/angsd_analysis/Fst_dxy_theta/Fst/cviol_forFstpop3_out -minInd 3 -minIndDepth 3 -doCounts 1

# pop4 n=10 * 0.7 = 7
angsd -bam /pylon5/bi4iflp/mlim/SeqCapData/angsd_analysis/Fst_dxy_theta/cviolbams_pop4.list -ref /pylon5/bi4iflp/mlim/SeqCapData/species_ref_files/Cviolifer_reffasta/ref_cviol_targetedRegionAndFlanking.fasta -anc /pylon5/bi4iflp/mlim/SeqCapData/species_ref_files/Cviolifer_reffasta/ref_cviol_targetedRegionAndFlanking.fasta -only_proper_pairs 0 -minMapQ 10 -minQ 20 -GL 1 -doSaf 1 -sites /pylon5/bi4iflp/mlim/SeqCapData/Novoalign_outs/cviol.keep -rf /pylon5/bi4iflp/mlim/SeqCapData/Novoalign_outs/cviol.rf -out /pylon5/bi4iflp/mlim/SeqCapData/angsd_analysis/Fst_dxy_theta/Fst/cviol_forFstpop4_out -minInd 7 -minIndDepth 3 -doCounts 1

# pop5 n=15 * 0.7 = 11
angsd -bam /pylon5/bi4iflp/mlim/SeqCapData/angsd_analysis/Fst_dxy_theta/cviolbams_pop5.list -ref /pylon5/bi4iflp/mlim/SeqCapData/species_ref_files/Cviolifer_reffasta/ref_cviol_targetedRegionAndFlanking.fasta -anc /pylon5/bi4iflp/mlim/SeqCapData/species_ref_files/Cviolifer_reffasta/ref_cviol_targetedRegionAndFlanking.fasta -only_proper_pairs 0 -minMapQ 10 -minQ 20 -GL 1 -doSaf 1 -sites /pylon5/bi4iflp/mlim/SeqCapData/Novoalign_outs/cviol.keep -rf /pylon5/bi4iflp/mlim/SeqCapData/Novoalign_outs/cviol.rf -out /pylon5/bi4iflp/mlim/SeqCapData/angsd_analysis/Fst_dxy_theta/Fst/cviol_forFstpop5_out -minInd 11 -minIndDepth 3 -doCounts 1

# pop6 n=10 * 0.7 = 7
angsd -bam /pylon5/bi4iflp/mlim/SeqCapData/angsd_analysis/Fst_dxy_theta/cviolbams_pop6.list -ref /pylon5/bi4iflp/mlim/SeqCapData/species_ref_files/Cviolifer_reffasta/ref_cviol_targetedRegionAndFlanking.fasta -anc /pylon5/bi4iflp/mlim/SeqCapData/species_ref_files/Cviolifer_reffasta/ref_cviol_targetedRegionAndFlanking.fasta -only_proper_pairs 0 -minMapQ 10 -minQ 20 -GL 1 -doSaf 1 -sites /pylon5/bi4iflp/mlim/SeqCapData/Novoalign_outs/cviol.keep -rf /pylon5/bi4iflp/mlim/SeqCapData/Novoalign_outs/cviol.rf -out /pylon5/bi4iflp/mlim/SeqCapData/angsd_analysis/Fst_dxy_theta/Fst/cviol_forFstpop6_out -minInd 7 -minIndDepth 3 -doCounts 1
