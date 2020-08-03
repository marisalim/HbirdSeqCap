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
# pop1 n=9 * 0.7 = 6
angsd -bam /pylon5/bi4iflp/mlim/SeqCapData/angsd_analysis/Fst_dxy_theta/ccorubams_pop1.list -ref /pylon5/bi4iflp/mlim/SeqCapData/species_ref_files/Ccoruscans_reffasta/ref_ccoru_targetedRegionAndFlanking.fasta -anc /pylon5/bi4iflp/mlim/SeqCapData/species_ref_files/Ccoruscans_reffasta/ref_ccoru_targetedRegionAndFlanking.fasta -only_proper_pairs 0 -minMapQ 10 -minQ 20 -GL 1 -doSaf 1 -sites /pylon5/bi4iflp/mlim/SeqCapData/Novoalign_outs/ccoru.keep -rf /pylon5/bi4iflp/mlim/SeqCapData/Novoalign_outs/ccoru.rf -out /pylon5/bi4iflp/mlim/SeqCapData/angsd_analysis/Fst_dxy_theta/Fst/ccoru_forFstpop1_out -minInd 6 -minIndDepth 3 -doCounts 1

# pop2 n=11 * 0.7 = 8
angsd -bam /pylon5/bi4iflp/mlim/SeqCapData/angsd_analysis/Fst_dxy_theta/ccorubams_pop2.list -ref /pylon5/bi4iflp/mlim/SeqCapData/species_ref_files/Ccoruscans_reffasta/ref_ccoru_targetedRegionAndFlanking.fasta -anc /pylon5/bi4iflp/mlim/SeqCapData/species_ref_files/Ccoruscans_reffasta/ref_ccoru_targetedRegionAndFlanking.fasta -only_proper_pairs 0 -minMapQ 10 -minQ 20 -GL 1 -doSaf 1 -sites /pylon5/bi4iflp/mlim/SeqCapData/Novoalign_outs/ccoru.keep -rf /pylon5/bi4iflp/mlim/SeqCapData/Novoalign_outs/ccoru.rf -out /pylon5/bi4iflp/mlim/SeqCapData/angsd_analysis/Fst_dxy_theta/Fst/ccoru_forFstpop2_out -minInd 8 -minIndDepth 3 -doCounts 1

# pop3 n=4 * 0.7 = 3
angsd -bam /pylon5/bi4iflp/mlim/SeqCapData/angsd_analysis/Fst_dxy_theta/ccorubams_pop3.list -ref /pylon5/bi4iflp/mlim/SeqCapData/species_ref_files/Ccoruscans_reffasta/ref_ccoru_targetedRegionAndFlanking.fasta -anc /pylon5/bi4iflp/mlim/SeqCapData/species_ref_files/Ccoruscans_reffasta/ref_ccoru_targetedRegionAndFlanking.fasta -only_proper_pairs 0 -minMapQ 10 -minQ 20 -GL 1 -doSaf 1 -sites /pylon5/bi4iflp/mlim/SeqCapData/Novoalign_outs/ccoru.keep -rf /pylon5/bi4iflp/mlim/SeqCapData/Novoalign_outs/ccoru.rf -out /pylon5/bi4iflp/mlim/SeqCapData/angsd_analysis/Fst_dxy_theta/Fst/ccoru_forFstpop3_out -minInd 3 -minIndDepth 3 -doCounts 1

# pop4 n=39 * 0.7 = 27
angsd -bam /pylon5/bi4iflp/mlim/SeqCapData/angsd_analysis/Fst_dxy_theta/ccorubams_pop4.list -ref /pylon5/bi4iflp/mlim/SeqCapData/species_ref_files/Ccoruscans_reffasta/ref_ccoru_targetedRegionAndFlanking.fasta -anc /pylon5/bi4iflp/mlim/SeqCapData/species_ref_files/Ccoruscans_reffasta/ref_ccoru_targetedRegionAndFlanking.fasta -only_proper_pairs 0 -minMapQ 10 -minQ 20 -GL 1 -doSaf 1 -sites /pylon5/bi4iflp/mlim/SeqCapData/Novoalign_outs/ccoru.keep -rf /pylon5/bi4iflp/mlim/SeqCapData/Novoalign_outs/ccoru.rf -out /pylon5/bi4iflp/mlim/SeqCapData/angsd_analysis/Fst_dxy_theta/Fst/ccoru_forFstpop4_out -minInd 27 -minIndDepth 3 -doCounts 1

# pop5 n=8 * 0.7 = 6 
angsd -bam /pylon5/bi4iflp/mlim/SeqCapData/angsd_analysis/Fst_dxy_theta/ccorubams_pop5.list -ref /pylon5/bi4iflp/mlim/SeqCapData/species_ref_files/Ccoruscans_reffasta/ref_ccoru_targetedRegionAndFlanking.fasta -anc /pylon5/bi4iflp/mlim/SeqCapData/species_ref_files/Ccoruscans_reffasta/ref_ccoru_targetedRegionAndFlanking.fasta -only_proper_pairs 0 -minMapQ 10 -minQ 20 -GL 1 -doSaf 1 -sites /pylon5/bi4iflp/mlim/SeqCapData/Novoalign_outs/ccoru.keep -rf /pylon5/bi4iflp/mlim/SeqCapData/Novoalign_outs/ccoru.rf -out /pylon5/bi4iflp/mlim/SeqCapData/angsd_analysis/Fst_dxy_theta/Fst/ccoru_forFstpop5_out -minInd 6 -minIndDepth 3 -doCounts 1

# pop6 n=6 * 0.7 = 4
angsd -bam /pylon5/bi4iflp/mlim/SeqCapData/angsd_analysis/Fst_dxy_theta/ccorubams_pop6.list -ref /pylon5/bi4iflp/mlim/SeqCapData/species_ref_files/Ccoruscans_reffasta/ref_ccoru_targetedRegionAndFlanking.fasta -anc /pylon5/bi4iflp/mlim/SeqCapData/species_ref_files/Ccoruscans_reffasta/ref_ccoru_targetedRegionAndFlanking.fasta -only_proper_pairs 0 -minMapQ 10 -minQ 20 -GL 1 -doSaf 1 -sites /pylon5/bi4iflp/mlim/SeqCapData/Novoalign_outs/ccoru.keep -rf /pylon5/bi4iflp/mlim/SeqCapData/Novoalign_outs/ccoru.rf -out /pylon5/bi4iflp/mlim/SeqCapData/angsd_analysis/Fst_dxy_theta/Fst/ccoru_forFstpop6_out -minInd 4 -minIndDepth 3 -doCounts 1

# pop7 n=9 * 0.7 = 6
angsd -bam /pylon5/bi4iflp/mlim/SeqCapData/angsd_analysis/Fst_dxy_theta/ccorubams_pop7.list -ref /pylon5/bi4iflp/mlim/SeqCapData/species_ref_files/Ccoruscans_reffasta/ref_ccoru_targetedRegionAndFlanking.fasta -anc /pylon5/bi4iflp/mlim/SeqCapData/species_ref_files/Ccoruscans_reffasta/ref_ccoru_targetedRegionAndFlanking.fasta -only_proper_pairs 0 -minMapQ 10 -minQ 20 -GL 1 -doSaf 1 -sites /pylon5/bi4iflp/mlim/SeqCapData/Novoalign_outs/ccoru.keep -rf /pylon5/bi4iflp/mlim/SeqCapData/Novoalign_outs/ccoru.rf -out /pylon5/bi4iflp/mlim/SeqCapData/angsd_analysis/Fst_dxy_theta/Fst/ccoru_forFstpop7_out -minInd 6 -minIndDepth 3 -doCounts 1

# pop8 n=6 * 0.7 = 4
angsd -bam /pylon5/bi4iflp/mlim/SeqCapData/angsd_analysis/Fst_dxy_theta/ccorubams_pop8.list -ref /pylon5/bi4iflp/mlim/SeqCapData/species_ref_files/Ccoruscans_reffasta/ref_ccoru_targetedRegionAndFlanking.fasta -anc /pylon5/bi4iflp/mlim/SeqCapData/species_ref_files/Ccoruscans_reffasta/ref_ccoru_targetedRegionAndFlanking.fasta -only_proper_pairs 0 -minMapQ 10 -minQ 20 -GL 1 -doSaf 1 -sites /pylon5/bi4iflp/mlim/SeqCapData/Novoalign_outs/ccoru.keep -rf /pylon5/bi4iflp/mlim/SeqCapData/Novoalign_outs/ccoru.rf -out /pylon5/bi4iflp/mlim/SeqCapData/angsd_analysis/Fst_dxy_theta/Fst/ccoru_forFstpop8_out -minInd 4 -minIndDepth 3 -doCounts 1

# pop9 n=6 * 0.7 = 4
angsd -bam /pylon5/bi4iflp/mlim/SeqCapData/angsd_analysis/Fst_dxy_theta/ccorubams_pop9.list -ref /pylon5/bi4iflp/mlim/SeqCapData/species_ref_files/Ccoruscans_reffasta/ref_ccoru_targetedRegionAndFlanking.fasta -anc /pylon5/bi4iflp/mlim/SeqCapData/species_ref_files/Ccoruscans_reffasta/ref_ccoru_targetedRegionAndFlanking.fasta -only_proper_pairs 0 -minMapQ 10 -minQ 20 -GL 1 -doSaf 1 -sites /pylon5/bi4iflp/mlim/SeqCapData/Novoalign_outs/ccoru.keep -rf /pylon5/bi4iflp/mlim/SeqCapData/Novoalign_outs/ccoru.rf -out /pylon5/bi4iflp/mlim/SeqCapData/angsd_analysis/Fst_dxy_theta/Fst/ccoru_forFstpop9_out -minInd 4 -minIndDepth 3 -doCounts 1
