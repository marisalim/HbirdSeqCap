#!/bin/bash
#SBATCH -N 1
#SBATCH -p LM
#SBATCH --mem=128GB
#SBATCH -t 2:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=marisa.lim@stonybrook.edu

#echo commands to stdout
set -x

grep -w -v 'combined_Contig233' MAFfiltered_cviol59_FINAL.txt | grep -w -v 'combined_Contig11577' | grep -w -v 'combined_Contig11576' | grep -w -v 'combined_Contig11575' | gzip > cviol59_Filtered_noHBBHBE.geno.gz

grep -w -v 'combined_Contig233' MAFfiltered_ccoru97_FINAL.txt | grep -w -v 'combined_Contig11577' | grep -w -v 'combined_Contig11576' | grep -w -v 'combined_Contig11575' | gzip > ccoru97_Filtered_noHBBHBE.geno.gz
