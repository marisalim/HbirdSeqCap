module load angsd/v0.921

# 22Oct18 version
bamfile_pop1=/groups/hologenomics/mlim/data/FstDxyThetaPi_CcoruCviolFix_21Oct18/ccorubams_pop1.list 
bamfile_pop2=/groups/hologenomics/mlim/data/FstDxyThetaPi_CcoruCviolFix_21Oct18/ccorubams_pop2.list 
bamfile_pop3=/groups/hologenomics/mlim/data/FstDxyThetaPi_CcoruCviolFix_21Oct18/ccorubams_pop3.list 
bamfile_pop4=/groups/hologenomics/mlim/data/FstDxyThetaPi_CcoruCviolFix_21Oct18/ccorubams_pop4.list 
bamfile_pop5=/groups/hologenomics/mlim/data/FstDxyThetaPi_CcoruCviolFix_21Oct18/ccorubams_pop5.list 
bamfile_pop6=/groups/hologenomics/mlim/data/FstDxyThetaPi_CcoruCviolFix_21Oct18/ccorubams_pop6.list 
bamfile_pop7=/groups/hologenomics/mlim/data/FstDxyThetaPi_CcoruCviolFix_21Oct18/ccorubams_pop7.list 
reffile=/groups/hologenomics/data/sequencing/Marisa/Seqcap_bioinformatics/species_ref_files/Ccoruscans_reffasta/ref_ccoru_targetedRegionAndFlanking.fasta 
ancfile=/groups/hologenomics/data/sequencing/Marisa/Seqcap_bioinformatics/species_ref_files/Ccoruscans_reffasta/ref_ccoru_targetedRegionAndFlanking.fasta 
sitesfile=/groups/hologenomics/data/sequencing/Marisa/Seqcap_bioinformatics/Novoalign_outs/ccoru.keep 
rf_file=/groups/hologenomics/data/sequencing/Marisa/Seqcap_bioinformatics/Novoalign_outs/ccoru.rf 
pest_pop1=/groups/hologenomics/mlim/data/FstDxyThetaPi_CcoruCviolFix_21Oct18/ThetaPi_CcoruCviolFix/ccoru_fordxythetapop1_out.sfs 
pest_pop2=/groups/hologenomics/mlim/data/FstDxyThetaPi_CcoruCviolFix_21Oct18/ThetaPi_CcoruCviolFix/ccoru_fordxythetapop2_out.sfs 
pest_pop3=/groups/hologenomics/mlim/data/FstDxyThetaPi_CcoruCviolFix_21Oct18/ThetaPi_CcoruCviolFix/ccoru_fordxythetapop3_out.sfs 
pest_pop4=/groups/hologenomics/mlim/data/FstDxyThetaPi_CcoruCviolFix_21Oct18/ThetaPi_CcoruCviolFix/ccoru_fordxythetapop4_out.sfs 
pest_pop5=/groups/hologenomics/mlim/data/FstDxyThetaPi_CcoruCviolFix_21Oct18/ThetaPi_CcoruCviolFix/ccoru_fordxythetapop5_out.sfs 
pest_pop6=/groups/hologenomics/mlim/data/FstDxyThetaPi_CcoruCviolFix_21Oct18/ThetaPi_CcoruCviolFix/ccoru_fordxythetapop6_out.sfs 
pest_pop7=/groups/hologenomics/mlim/data/FstDxyThetaPi_CcoruCviolFix_21Oct18/ThetaPi_CcoruCviolFix/ccoru_fordxythetapop7_out.sfs 
outfile_pop1=/groups/hologenomics/mlim/data/FstDxyThetaPi_CcoruCviolFix_21Oct18/ThetaPi_CcoruCviolFix/ccoru_forthetapop1_out 
outfile_pop2=/groups/hologenomics/mlim/data/FstDxyThetaPi_CcoruCviolFix_21Oct18/ThetaPi_CcoruCviolFix/ccoru_forthetapop2_out 
outfile_pop3=/groups/hologenomics/mlim/data/FstDxyThetaPi_CcoruCviolFix_21Oct18/ThetaPi_CcoruCviolFix/ccoru_forthetapop3_out 
outfile_pop4=/groups/hologenomics/mlim/data/FstDxyThetaPi_CcoruCviolFix_21Oct18/ThetaPi_CcoruCviolFix/ccoru_forthetapop4_out 
outfile_pop5=/groups/hologenomics/mlim/data/FstDxyThetaPi_CcoruCviolFix_21Oct18/ThetaPi_CcoruCviolFix/ccoru_forthetapop5_out 
outfile_pop6=/groups/hologenomics/mlim/data/FstDxyThetaPi_CcoruCviolFix_21Oct18/ThetaPi_CcoruCviolFix/ccoru_forthetapop6_out 
outfile_pop7=/groups/hologenomics/mlim/data/FstDxyThetaPi_CcoruCviolFix_21Oct18/ThetaPi_CcoruCviolFix/ccoru_forthetapop7_out
# minInd 70% pop1 n=9 * 0.7 = 6
angsd -bam $bamfile_pop1 -ref $reffile -anc $ancfile -only_proper_pairs 0 -minMapQ 10 -minQ 20 -GL 1 -doSaf 1 -sites $sitesfile -rf $rf_file -minInd 6 -minIndDepth 3 -doCounts 1 -doThetas 1 -fold 1 -pest $pest_pop1 -out $outfile_pop1
# pop2 n=15 * 0.7 = 11
angsd -bam $bamfile_pop2 -ref $reffile -anc $ancfile -only_proper_pairs 0 -minMapQ 10 -minQ 20 -GL 1 -doSaf 1 -sites $sitesfile -rf $rf_file -minInd 11 -minIndDepth 3 -doCounts 1 -doThetas 1 -fold 1 -pest $pest_pop2 -out $outfile_pop2
# pop3 n=38 * 0.7 = 27
angsd -bam $bamfile_pop3 -ref $reffile -anc $ancfile -only_proper_pairs 0 -minMapQ 10 -minQ 20 -GL 1 -doSaf 1 -sites $sitesfile -rf $rf_file -minInd 27 -minIndDepth 3 -doCounts 1 -doThetas 1 -fold 1 -pest $pest_pop3 -out $outfile_pop3
# pop4 n=8 * 0.7 = 6
angsd -bam $bamfile_pop4 -ref $reffile -anc $ancfile -only_proper_pairs 0 -minMapQ 10 -minQ 20 -GL 1 -doSaf 1 -sites $sitesfile -rf $rf_file -minInd 6 -minIndDepth 3 -doCounts 1 -doThetas 1 -fold 1 -pest $pest_pop4 -out $outfile_pop4
# pop5 n=13 * 0.7 = 9
angsd -bam $bamfile_pop5 -ref $reffile -anc $ancfile -only_proper_pairs 0 -minMapQ 10 -minQ 20 -GL 1 -doSaf 1 -sites $sitesfile -rf $rf_file -minInd 9 -minIndDepth 3 -doCounts 1 -doThetas 1 -fold 1 -pest $pest_pop5 -out $outfile_pop5
# pop6 n=6 * 0.7 = 4
angsd -bam $bamfile_pop6 -ref $reffile -anc $ancfile -only_proper_pairs 0 -minMapQ 10 -minQ 20 -GL 1 -doSaf 1 -sites $sitesfile -rf $rf_file -minInd 4 -minIndDepth 3 -doCounts 1 -doThetas 1 -fold 1 -pest $pest_pop6 -out $outfile_pop6
# pop7 n=8 * 0.7 = 6
angsd -bam $bamfile_pop7 -ref $reffile -anc $ancfile -only_proper_pairs 0 -minMapQ 10 -minQ 20 -GL 1 -doSaf 1 -sites $sitesfile -rf $rf_file -minInd 6 -minIndDepth 3 -doCounts 1 -doThetas 1 -fold 1 -pest $pest_pop7 -out $outfile_pop7
