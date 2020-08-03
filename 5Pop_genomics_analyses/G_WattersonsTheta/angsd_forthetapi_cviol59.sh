module load angsd/v0.921

# 22Oct18 version
bamfile_pop1=/groups/hologenomics/mlim/data/FstDxyThetaPi_CcoruCviolFix_21Oct18/cviolbams_pop1.list 
bamfile_pop2=/groups/hologenomics/mlim/data/FstDxyThetaPi_CcoruCviolFix_21Oct18/cviolbams_pop2.list 
bamfile_pop3=/groups/hologenomics/mlim/data/FstDxyThetaPi_CcoruCviolFix_21Oct18/cviolbams_pop3.list 
bamfile_pop4=/groups/hologenomics/mlim/data/FstDxyThetaPi_CcoruCviolFix_21Oct18/cviolbams_pop4.list 
bamfile_pop5=/groups/hologenomics/mlim/data/FstDxyThetaPi_CcoruCviolFix_21Oct18/cviolbams_pop5.list 
bamfile_pop6=/groups/hologenomics/mlim/data/FstDxyThetaPi_CcoruCviolFix_21Oct18/cviolbams_pop6.list 
reffile=/groups/hologenomics/data/sequencing/Marisa/Seqcap_bioinformatics/species_ref_files/Cviolifer_reffasta/ref_cviol_targetedRegionAndFlanking.fasta 
ancfile=/groups/hologenomics/data/sequencing/Marisa/Seqcap_bioinformatics/species_ref_files/Cviolifer_reffasta/ref_cviol_targetedRegionAndFlanking.fasta 
sitesfile=/groups/hologenomics/data/sequencing/Marisa/Seqcap_bioinformatics/Novoalign_outs/cviol.keep 
rf_file=/groups/hologenomics/data/sequencing/Marisa/Seqcap_bioinformatics/Novoalign_outs/cviol.rf 
outfile_pop1=/groups/hologenomics/mlim/data/FstDxyThetaPi_CcoruCviolFix_21Oct18/ThetaPi_CcoruCviolFix/FOLDED_forthetapi/cviol_fordxythetapop1_out 
outfile_pop2=/groups/hologenomics/mlim/data/FstDxyThetaPi_CcoruCviolFix_21Oct18/ThetaPi_CcoruCviolFix/FOLDED_forthetapi/cviol_fordxythetapop2_out 
outfile_pop3=/groups/hologenomics/mlim/data/FstDxyThetaPi_CcoruCviolFix_21Oct18/ThetaPi_CcoruCviolFix/FOLDED_forthetapi/cviol_fordxythetapop3_out 
outfile_pop4=/groups/hologenomics/mlim/data/FstDxyThetaPi_CcoruCviolFix_21Oct18/ThetaPi_CcoruCviolFix/FOLDED_forthetapi/cviol_fordxythetapop4_out 
outfile_pop5=/groups/hologenomics/mlim/data/FstDxyThetaPi_CcoruCviolFix_21Oct18/ThetaPi_CcoruCviolFix/FOLDED_forthetapi/cviol_fordxythetapop5_out 
outfile_pop6=/groups/hologenomics/mlim/data/FstDxyThetaPi_CcoruCviolFix_21Oct18/ThetaPi_CcoruCviolFix/FOLDED_forthetapi/cviol_fordxythetapop6_out

# minInd 70% pop1 n=6 * 0.7 = 4
angsd -bam $bamfile_pop1 -ref $reffile -anc $ancfile -only_proper_pairs 0 -minMapQ 10 -minQ 20 -GL 1 -doSaf 1 -doGeno 32 -doPost 1 -doMaf 1 -fold 1 -doMajorMinor 1 -skipTriallelic 1 -sites $sitesfile -rf $rf_file -minInd 4 -minIndDepth 3 -doCounts 1 -out $outfile_pop1
# pop2 n=14 * 0.7 = 10
angsd -bam $bamfile_pop2 -ref $reffile -anc $ancfile -only_proper_pairs 0 -minMapQ 10 -minQ 20 -GL 1 -doSaf 1 -doGeno 32 -doPost 1 -doMaf 1 -fold 1 -doMajorMinor 1 -skipTriallelic 1 -sites $sitesfile -rf $rf_file -minInd 10 -minIndDepth 3 -doCounts 1 -out $outfile_pop2
# pop3 n=4 * 0.7 = 3
angsd -bam $bamfile_pop3 -ref $reffile -anc $ancfile -only_proper_pairs 0 -minMapQ 10 -minQ 20 -GL 1 -doSaf 1 -doGeno 32 -doPost 1 -doMaf 1 -fold 1 -doMajorMinor 1 -skipTriallelic 1 -sites $sitesfile -rf $rf_file -minInd 3 -minIndDepth 3 -doCounts 1 -out $outfile_pop3
# pop4 n=15 * 0.7 = 11
angsd -bam $bamfile_pop4 -ref $reffile -anc $ancfile -only_proper_pairs 0 -minMapQ 10 -minQ 20 -GL 1 -doSaf 1 -doGeno 32 -doPost 1 -doMaf 1 -fold 1 -doMajorMinor 1 -skipTriallelic 1 -sites $sitesfile -rf $rf_file -minInd 11 -minIndDepth 3 -doCounts 1 -out $outfile_pop4
# pop5 n=10 * 0.7 = 7
angsd -bam $bamfile_pop5 -ref $reffile -anc $ancfile -only_proper_pairs 0 -minMapQ 10 -minQ 20 -GL 1 -doSaf 1 -doGeno 32 -doPost 1 -doMaf 1 -fold 1 -doMajorMinor 1 -skipTriallelic 1 -sites $sitesfile -rf $rf_file -minInd 7 -minIndDepth 3 -doCounts 1 -out $outfile_pop5
# pop6 n=10 * 0.7 = 7
angsd -bam $bamfile_pop6 -ref $reffile -anc $ancfile -only_proper_pairs 0 -minMapQ 10 -minQ 20 -GL 1 -doSaf 1 -doGeno 32 -doPost 1 -doMaf 1 -fold 1 -doMajorMinor 1 -skipTriallelic 1 -sites $sitesfile -rf $rf_file -minInd 7 -minIndDepth 3 -doCounts 1 -out $outfile_pop6
