#!/bin/bash

#create geno probability file for calculating dissimilarity matrix in ngsDist as input for EEMs
# using genotype probabilities rather than calling genotypes
# same command as for calling genotypes, but use doGeno 8

# postCutoff: 0.95
angsd -bam ccorubams_n97.list -ref ref_ccoru_targetedRegionAndFlanking.fasta -anc ref_ccoru_targetedRegionAndFlanking.fasta -doGeno 8 -SNP_pval 1e-3 -postCutoff 0.95 -minInd 68 -geno_minDepth 5 -only_proper_pairs 0 -minMapQ 10 -minQ 20 -fold 1 -GL 1 -doPost 2 -doMaf 1 -doMajorMinor 1 -doSaf 1 -skipTriallelic 1 -doCounts 1 -sites ccoru.keep -rf ccoru.rf -out ccoru_forngsDist/ccoru97_genoforngsdist_out
