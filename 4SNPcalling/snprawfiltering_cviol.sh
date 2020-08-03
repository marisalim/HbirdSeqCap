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
module load samtools/1.3.1
module load bcftools/1.3.1

samtools mpileup -s -Q 20 -A -t DP,SP -d 100000 -E -ugf /pylon5/bi4iflp/mlim/SeqCapData/species_ref_files/Cviolifer_reffasta/ref_cviol_targetedRegionAndFlanking.fasta /pylon5/bi4iflp/mlim/SeqCapData/Novoalign_outs/novo_cviol/*.bam | bcftools call -c -p 0.1 - > /pylon5/bi4iflp/mlim/SeqCapData/Novoalign_outs/cviol.raw.vcf

