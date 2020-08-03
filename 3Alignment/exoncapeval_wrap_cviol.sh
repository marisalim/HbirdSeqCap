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
module load samtools/1.3
module load bedtools/2.25.0

perl /pylon5/bi4iflp/mlim/SeqCapData/6-exonCaptureEvaluation Evaluation -genome /pylon5/bi4iflp/mlim/SeqCapData/species_ref_files/Cviolifer_reffasta/ref_cviol_targetedRegionAndFlanking.fasta -cleanDir /pylon5/bi4iflp/mlim/SeqCapData/ScrubReads_Results/files_to_align/Cviol_readstoalign/ -rawDir /pylon5/bi4iflp/mlim/SeqCapData/Lim/pre-clean/orginal/rawfq_cviol -bamDir /pylon5/bi4iflp/mlim/SeqCapData/Novoalign_outs/novo_cviol -resDir /pylon5/bi4iflp/mlim/SeqCapData/ExonCapEvaluation/Cviol_exoncapeval -bedFile /pylon5/bi4iflp/mlim/SeqCapData/species_ref_files/Cviol_bedfiles/combined_targeted_region.bed
