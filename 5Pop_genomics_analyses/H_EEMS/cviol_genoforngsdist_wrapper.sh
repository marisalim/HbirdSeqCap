#!/bin/bash

#create geno probability file for calculating dissimilarity matrix in ngsDist as input for EEMs
# using genotype probabilities rather than calling genotypes
# same command as for calling genotypes, but use doGeno 8

# postCutoff: 0.95
angsd -bam cviolbams_n59.list -ref ref_cviol_targetedRegionAndFlanking.fasta -anc ref_cviol_targetedRegionAndFlanking.fasta -doGeno 8 -SNP_pval 1e-3 -postCutoff 0.95 -minInd 41 -geno_minDepth 5 -only_proper_pairs 0 -minMapQ 10 -minQ 20 -fold 1 -GL 1 -doPost 2 -doMaf 1 -doMajorMinor 1 -doSaf 1 -skipTriallelic 1 -doCounts 1 -sites cviol.keep -rf cviol.rf -out cviol_forngsDist/cviol59_genoforngsdist_out
