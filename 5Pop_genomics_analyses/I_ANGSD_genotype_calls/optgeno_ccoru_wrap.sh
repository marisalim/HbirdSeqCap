#!/bin/bash
#SBATCH -N 1
#SBATCH -p LM
#SBATCH --mem=128GB
#SBATCH -t 10:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=marisa.lim@stonybrook.edu

#echo commands to stdout
set -x

source ~/.bashrc

# 70% minInd = 71
#angsd -bam /pylon5/bi4iflp/mlim/SeqCapData/angsd_analysis/ccorubams.list -ref /pylon5/bi4iflp/mlim/SeqCapData/species_ref_files/Ccoruscans_reffasta/ref_ccoru_targetedRegionAndFlanking.fasta -anc /pylon5/bi4iflp/mlim/SeqCapData/species_ref_files/Ccoruscans_reffasta/ref_ccoru_targetedRegionAndFlanking.fasta -doGeno 2 -SNP_pval 1e-3 -postCutoff 0.75 -minInd 71 -geno_minDepth 5 -only_proper_pairs 0 -minMapQ 10 -minQ 20 -fold 1 -GL 1 -doPost 1 -doMaf 1 -doMajorMinor 1 -doSaf 1 -skipTriallelic 1 -doCounts 1 -sites /pylon5/bi4iflp/mlim/SeqCapData/Novoalign_outs/ccoru.keep -rf /pylon5/bi4iflp/mlim/SeqCapData/Novoalign_outs/ccoru.rf -out /pylon5/bi4iflp/mlim/SeqCapData/angsd_analysis/ANGSD_missingdat_check/ccoru_missingdatcheck_out

# minInd 70% = 68
#angsd -bam /pylon5/bi4iflp/mlim/SeqCapData/angsd_analysis/ccorubams_n97.list -ref /pylon5/bi4iflp/mlim/SeqCapData/species_ref_files/Ccoruscans_reffasta/ref_ccoru_targetedRegionAndFlanking.fasta -anc /pylon5/bi4iflp/mlim/SeqCapData/species_ref_files/Ccoruscans_reffasta/ref_ccoru_targetedRegionAndFlanking.fasta -doGeno 2 -SNP_pval 1e-6 -postCutoff 0.75 -hetBias_pval 0.001 -dosnpstat 1 -minInd 68 -geno_minDepth 5 -HWE_pval 1 -only_proper_pairs 0 -minMapQ 10 -minQ 20 -fold 1 -GL 1 -doPost 2 -doMaf 1 -doMajorMinor 1 -doSaf 1 -skipTriallelic 1 -doCounts 1 -sites /pylon5/bi4iflp/mlim/SeqCapData/Novoalign_outs/ccoru.keep -rf /pylon5/bi4iflp/mlim/SeqCapData/Novoalign_outs/ccoru.rf -out /pylon5/bi4iflp/mlim/SeqCapData/angsd_analysis/ANGSD_optgenocall_check/ccoru97_optgeno_out

# 10/12/17: testing different postCutoff values
# keeping SNP_pval 1e-3, removing hetBias_pval flag (too stringent)
# postCutoff: 0.75
#angsd -bam /pylon5/bi4iflp/mlim/SeqCapData/angsd_analysis/ccorubams_n97.list -ref /pylon5/bi4iflp/mlim/SeqCapData/species_ref_files/Ccoruscans_reffasta/ref_ccoru_targetedRegionAndFlanking.fasta -anc /pylon5/bi4iflp/mlim/SeqCapData/species_ref_files/Ccoruscans_reffasta/ref_ccoru_targetedRegionAndFlanking.fasta -doGeno 2 -SNP_pval 1e-3 -postCutoff 0.75 -dosnpstat 1 -minInd 68 -geno_minDepth 5 -HWE_pval 1 -only_proper_pairs 0 -minMapQ 10 -minQ 20 -fold 1 -GL 1 -doPost 2 -doMaf 1 -doMajorMinor 1 -doSaf 1 -skipTriallelic 1 -doCounts 1 -sites /pylon5/bi4iflp/mlim/SeqCapData/Novoalign_outs/ccoru.keep -rf /pylon5/bi4iflp/mlim/SeqCapData/Novoalign_outs/ccoru.rf -out /pylon5/bi4iflp/mlim/SeqCapData/angsd_analysis/ANGSD_optgenocall_check/ccoru97_optP75_out

# postCutoff: 0.85
#angsd -bam /pylon5/bi4iflp/mlim/SeqCapData/angsd_analysis/ccorubams_n97.list -ref /pylon5/bi4iflp/mlim/SeqCapData/species_ref_files/Ccoruscans_reffasta/ref_ccoru_targetedRegionAndFlanking.fasta -anc /pylon5/bi4iflp/mlim/SeqCapData/species_ref_files/Ccoruscans_reffasta/ref_ccoru_targetedRegionAndFlanking.fasta -doGeno 2 -SNP_pval 1e-3 -postCutoff 0.85 -dosnpstat 1 -minInd 68 -geno_minDepth 5 -HWE_pval 1 -only_proper_pairs 0 -minMapQ 10 -minQ 20 -fold 1 -GL 1 -doPost 2 -doMaf 1 -doMajorMinor 1 -doSaf 1 -skipTriallelic 1 -doCounts 1 -sites /pylon5/bi4iflp/mlim/SeqCapData/Novoalign_outs/ccoru.keep -rf /pylon5/bi4iflp/mlim/SeqCapData/Novoalign_outs/ccoru.rf -out /pylon5/bi4iflp/mlim/SeqCapData/angsd_analysis/ANGSD_optgenocall_check/ccoru97_optP85_out

# FINAL COMMAND
# postCutoff: 0.95
angsd -bam /pylon5/bi4iflp/mlim/SeqCapData/angsd_analysis/ccorubams_n97.list -ref /pylon5/bi4iflp/mlim/SeqCapData/species_ref_files/Ccoruscans_reffasta/ref_ccoru_targetedRegionAndFlanking.fasta -anc /pylon5/bi4iflp/mlim/SeqCapData/species_ref_files/Ccoruscans_reffasta/ref_ccoru_targetedRegionAndFlanking.fasta -doGeno 2 -SNP_pval 1e-3 -postCutoff 0.95 -dosnpstat 1 -minInd 68 -geno_minDepth 5 -only_proper_pairs 0 -minMapQ 10 -minQ 20 -fold 1 -GL 1 -doPost 2 -doMaf 1 -doMajorMinor 1 -doSaf 1 -skipTriallelic 1 -doCounts 1 -sites /pylon5/bi4iflp/mlim/SeqCapData/Novoalign_outs/ccoru.keep -rf /pylon5/bi4iflp/mlim/SeqCapData/Novoalign_outs/ccoru.rf -out /pylon5/bi4iflp/mlim/SeqCapData/angsd_analysis/ANGSD_optgenocall_check/ccoru97_optP95_out

# for dxy calculations - no snp calling, replace geno_minDepth with minIndDepth
#angsd -bam /pylon5/bi4iflp/mlim/SeqCapData/angsd_analysis/ccorubams_n97.list -ref /pylon5/bi4iflp/mlim/SeqCapData/species_ref_files/Ccoruscans_reffasta/ref_ccoru_targetedRegionAndFlanking.fasta -anc /pylon5/bi4iflp/mlim/SeqCapData/species_ref_files/Ccoruscans_reffasta/ref_ccoru_targetedRegionAndFlanking.fasta -doGeno 2 -minInd 68 -minIndDepth 5 -only_proper_pairs 0 -minMapQ 10 -minQ 20 -fold 1 -GL 1 -doPost 2 -doMaf 1 -doMajorMinor 1 -doSaf 1 -skipTriallelic 1 -doCounts 1 -sites /pylon5/bi4iflp/mlim/SeqCapData/Novoalign_outs/ccoru.keep -rf /pylon5/bi4iflp/mlim/SeqCapData/Novoalign_outs/ccoru.rf -out /pylon5/bi4iflp/mlim/SeqCapData/angsd_analysis/Fst_dxy_theta/dxy/nosnpcalls_fordxy/ccoru97_nosnpcalls_fordxy_outs
