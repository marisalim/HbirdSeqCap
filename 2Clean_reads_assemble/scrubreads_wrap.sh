#!/bin/bash 
#SBATCH -N 1 
#SBATCH -p LM 
#SBATCH --mem=128GB 
#SBATCH -t 48:00:00 
#SBATCH --mail-type=ALL 
#SBATCH --mail-user=marisa.lim@stonybrook.edu 

#echo commands to stdout
set -x 

module load perl/5.24.0-threads
module load bowtie2/2.2.7
module load fastqc/0.11.3
module load flash/1.2.11
module load trimmomatic/0.36
source ~/.bashrc

perl /pylon1/bi4iflp/mlim/SeqCapData/1-ScrubReads_d cleanPE -f /pylon1/bi4iflp/mlim/SeqCapData/Lim/ -o /pylon1/bi4iflp/mlim/SeqCapData/ScrubReads_Results/ -t $TRIMMOMATIC_HOME/trimmomatic-0.36.jar -c /pylon1/bi4iflp/mlim/SeqCapData/e_coli_K12.fasta -d all -M 0 -l 150 -z 

