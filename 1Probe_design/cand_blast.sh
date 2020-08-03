#!/bin/bash 
#PBS -l nodes=1:ppn=15 
#PBS -l walltime=24:00:00 
#PBS -j oe 
#PBS -q batch 
#PBS -m e 
#PBS -M marisa.lim@stonybrook.edu 
#PBS -N blast_cand
#
set echo module load blast/2.2.31-all
#
blastp -query /crucible/bi4iflp/mlim/Calypte_genome_exoncheck/cand_blastthese.fasta -db /crucible/bi4iflp/mlim/Calypte_genome_exoncheck/Calypte_anna.pep -evalue 
1e-10 -outfmt 6 -out /crucible/bi4iflp/mlim/Calypte_genome_exoncheck/check_candblast.txt
