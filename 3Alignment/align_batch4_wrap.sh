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
module load gatk/3.6
module load picard/2.1.1
source ~/.bashrc

perl /pylon5/bi4iflp/mlim/SeqCapData/5-Alignment -f /pylon5/bi4iflp/mlim/SeqCapData/species_ref_files/Ccoruscans_reffasta/ref_ccoru_targetedRegionAndFlanking.fasta -r /pylon5/bi4iflp/mlim/SeqCapData/ScrubReads_Results/files_to_align/Ccoru_readstoalign/ccoru_bat3 -o /pylon5/bi4iflp/mlim/SeqCapData/Novoalign_outs/novo_ccoru -G $GATK_HOME/GenomeAnalysisTK.jar -P $PICARD_HOME/picard.jar -c 0 -k 0 -m 1 -i 250 -v 100 -l 151 -t 180 
