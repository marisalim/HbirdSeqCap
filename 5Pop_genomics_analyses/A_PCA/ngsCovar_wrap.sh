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

# To do the PCA

# Ccoru
#ngsCovar -probfile ccoru_forPCAn99_out.geno -outfile ccoru_forPCAn99_out.covar -nind 99 -nsites 9426843 -call 0 -norm 0
ngsCovar -probfile ccoru_forPCAn97_out.geno -outfile ccoru_forPCAn97_out.covar -nind 97 -nsites 9420945 -call 0 -norm 0

# Cviol
ngsCovar -probfile cviol_forPCAn59_out.geno -outfile cviol_forPCAn59_out.covar -nind 59 -nsites 10244586 -call 0 -norm 0
