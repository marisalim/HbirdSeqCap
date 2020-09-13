#!/bin/bash
#SBATCH -p long-40core
#SBATCH --mem=10G
#SBATCH -t 20:00:00
#SBATCH -J lfmmK2

# K=2
for i in {1..5}; do ../LFMM_CL_v1.5/bin/LFMM -x ./inputs/cviol_Elev400m_allelefreqs_forLFMMREDO.txt -v ./inputs/Cviol_meanstandardized_elevationbins400mREDO.txt -K 2 -o ./lfmm_results/run_K2_i250000b25000_iter${i} -i 250000 -b 25000; done
