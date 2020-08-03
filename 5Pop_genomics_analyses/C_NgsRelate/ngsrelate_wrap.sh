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

ngsRelate -g cviol_relatetest_out.glf.gz -n 62 -f cviol_freq -z cviol_samp.list > cviol_outfile.res

ngsRelate -g ccoru_relatetest_out.glf.gz -n 102 -f ccoru_freq -z ccoru_samp.list > ccoru_outfile.res
