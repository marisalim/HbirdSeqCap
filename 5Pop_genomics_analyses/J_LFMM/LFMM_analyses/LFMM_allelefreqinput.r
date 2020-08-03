# Script for making input file for LFMM (allele frequency per elevation bins), calculating q-value for LFMM results, and analyzing outliers
# note that all genotypes used in this code have been filtered to remove low frequency SNPs (MAF < 5%)
# elevation bin allele frequencies calculated in plotAlleleFreqs.r

# Marisa Lim (c)2018

# Libraries
library(reshape2)
library(ggplot2)
library(cowplot)
library(qvalue)
library(ggmap)
library(ggpubr)
library(dplyr)
library(MASS)

# -- 1. format LFMM input file ----
# set up allele frequency and environmental files for LFMM
# Cviol 400m bins
# Ccoru 300m and 400m bins
# Allele frequencies calculated for SNPs in plotAlleleFreqs.r
# environmental file:
  # make vector of average elevation for each bin
  # use mean and std of all elevations per species to standardize the avg bin elevations

# rows are elevation bins, columns are SNP allele frequencies

# read in allele frequency by bin files
cviol_allelefreqs1 <- read.table('../../Allelefreqplots_byelevation/Cviol_Elevation_bins400m_snpbatch1_allelefreqsbybins.txt', header=TRUE)
cviol_allelefreqs2 <- read.table('../../Allelefreqplots_byelevation/Cviol_Elevation_bins400m_snpbatch2_allelefreqsbybins.txt', header=TRUE)
cviol_allelefreqs3 <- read.table('../../Allelefreqplots_byelevation/Cviol_Elevation_bins400m_snpbatch3_allelefreqsbybins.txt', header=TRUE)
cviol_allelefreqs4 <- read.table('../../Allelefreqplots_byelevation/Cviol_Elevation_bins400m_snpbatch4_allelefreqsbybins.txt', header=TRUE)
cviol_allelefreqs5 <- read.table('../../Allelefreqplots_byelevation/Cviol_Elevation_bins400m_snpbatch5_allelefreqsbybins.txt', header=TRUE)
cviol_allelefreqs6 <- read.table('../../Allelefreqplots_byelevation/Cviol_Elevation_bins400m_snpbatch6_allelefreqsbybins.txt', header=TRUE)
cviol_allelefreqs7 <- read.table('../../Allelefreqplots_byelevation/Cviol_Elevation_bins400m_snpbatch7_allelefreqsbybins.txt', header=TRUE)
cviol_allelefreqs8 <- read.table('../../Allelefreqplots_byelevation/Cviol_Elevation_bins400m_snpbatch8_allelefreqsbybins.txt', header=TRUE)
cviol_allelefreqs9 <- read.table('../../Allelefreqplots_byelevation/Cviol_Elevation_bins400m_snpbatch9_allelefreqsbybins.txt', header=TRUE)

cviol_allelefreqs_all <- rbind(cviol_allelefreqs1, cviol_allelefreqs2, cviol_allelefreqs3, cviol_allelefreqs4, cviol_allelefreqs5, cviol_allelefreqs6,
                               cviol_allelefreqs7, cviol_allelefreqs8, cviol_allelefreqs9)

#400m bins
ccoru_allelefreqs1 <- read.table('../../Allelefreqplots_byelevation/Ccoru_Elevation_bin400m_snpbatch1_allelefreqsbybins.txt', header=TRUE)
ccoru_allelefreqs2 <- read.table('../../Allelefreqplots_byelevation/Ccoru_Elevation_bin400m_snpbatch2_allelefreqsbybins.txt', header=TRUE)
ccoru_allelefreqs3 <- read.table('../../Allelefreqplots_byelevation/Ccoru_Elevation_bin400m_snpbatch3_allelefreqsbybins.txt', header=TRUE)
ccoru_allelefreqs4 <- read.table('../../Allelefreqplots_byelevation/Ccoru_Elevation_bin400m_snpbatch4_allelefreqsbybins.txt', header=TRUE)
ccoru_allelefreqs5 <- read.table('../../Allelefreqplots_byelevation/Ccoru_Elevation_bin400m_snpbatch5_allelefreqsbybins.txt', header=TRUE)
ccoru_allelefreqs6 <- read.table('../../Allelefreqplots_byelevation/Ccoru_Elevation_bin400m_snpbatch6_allelefreqsbybins.txt', header=TRUE)
ccoru_allelefreqs7 <- read.table('../../Allelefreqplots_byelevation/Ccoru_Elevation_bin400m_snpbatch7_allelefreqsbybins.txt', header=TRUE)
ccoru_allelefreqs8 <- read.table('../../Allelefreqplots_byelevation/Ccoru_Elevation_bin400m_snpbatch8_allelefreqsbybins.txt', header=TRUE)
ccoru_allelefreqs9 <- read.table('../../Allelefreqplots_byelevation/Ccoru_Elevation_bin400m_snpbatch9_allelefreqsbybins.txt', header=TRUE)
ccoru_allelefreqs10 <- read.table('../../Allelefreqplots_byelevation/Ccoru_Elevation_bin400m_snpbatch10_allelefreqsbybins.txt', header=TRUE)
ccoru_allelefreqs11 <- read.table('../../Allelefreqplots_byelevation/Ccoru_Elevation_bin400m_snpbatch11_allelefreqsbybins.txt', header=TRUE)

ccoru_allelefreqs_all <- rbind(ccoru_allelefreqs1, ccoru_allelefreqs2, ccoru_allelefreqs3, ccoru_allelefreqs4, ccoru_allelefreqs5,
                               ccoru_allelefreqs6, ccoru_allelefreqs7, ccoru_allelefreqs8, ccoru_allelefreqs9, ccoru_allelefreqs10,
                               ccoru_allelefreqs11)
# 300m bins
ccoru_allelefreqs1_300m <- read.table('../../Allelefreqplots_byelevation/Ccoru_Elevation_bin300m_snpbatch1_allelefreqsbybins.txt', header=TRUE)
ccoru_allelefreqs2_300m <- read.table('../../Allelefreqplots_byelevation/Ccoru_Elevation_bin300m_snpbatch2_allelefreqsbybins.txt', header=TRUE)
ccoru_allelefreqs3_300m <- read.table('../../Allelefreqplots_byelevation/Ccoru_Elevation_bin300m_snpbatch3_allelefreqsbybins.txt', header=TRUE)
ccoru_allelefreqs4_300m <- read.table('../../Allelefreqplots_byelevation/Ccoru_Elevation_bin300m_snpbatch4_allelefreqsbybins.txt', header=TRUE)
ccoru_allelefreqs5_300m <- read.table('../../Allelefreqplots_byelevation/Ccoru_Elevation_bin300m_snpbatch5_allelefreqsbybins.txt', header=TRUE)
ccoru_allelefreqs6_300m <- read.table('../../Allelefreqplots_byelevation/Ccoru_Elevation_bin300m_snpbatch6_allelefreqsbybins.txt', header=TRUE)
ccoru_allelefreqs7_300m <- read.table('../../Allelefreqplots_byelevation/Ccoru_Elevation_bin300m_snpbatch7_allelefreqsbybins.txt', header=TRUE)
ccoru_allelefreqs8_300m <- read.table('../../Allelefreqplots_byelevation/Ccoru_Elevation_bin300m_snpbatch8_allelefreqsbybins.txt', header=TRUE)
ccoru_allelefreqs9_300m <- read.table('../../Allelefreqplots_byelevation/Ccoru_Elevation_bin300m_snpbatch9_allelefreqsbybins.txt', header=TRUE)
ccoru_allelefreqs10_300m <- read.table('../../Allelefreqplots_byelevation/Ccoru_Elevation_bin300m_snpbatch10_allelefreqsbybins.txt', header=TRUE)
ccoru_allelefreqs11_300m <- read.table('../../Allelefreqplots_byelevation/Ccoru_Elevation_bin300m_snpbatch11_allelefreqsbybins.txt', header=TRUE)

ccoru_allelefreqs_all_300m <- rbind(ccoru_allelefreqs1_300m, ccoru_allelefreqs2_300m, ccoru_allelefreqs3_300m, ccoru_allelefreqs4_300m, ccoru_allelefreqs5_300m,
                               ccoru_allelefreqs6_300m, ccoru_allelefreqs7_300m, ccoru_allelefreqs8_300m, ccoru_allelefreqs9_300m, ccoru_allelefreqs10_300m,
                               ccoru_allelefreqs11_300m)


#format LFMM input
cviol_allelefreqs_for_lfmm <- acast(cviol_allelefreqs_all, thebin~theSNP, value='littlea_freq')
dim(cviol_allelefreqs_for_lfmm) # 5 SNPs were removed because bin(s) had all 9's, so now there are 91124-5 = 91119 SNPs left
write.table(cviol_allelefreqs_for_lfmm, 'cviol_Elev400m_allelefreqs_forLFMM.txt', row.names=FALSE, col.names=FALSE)

ccoru_allelefreq_for_lfmm <- acast(ccoru_allelefreqs_all, thebin~theSNP, value='littlea_freq')
dim(ccoru_allelefreq_for_lfmm) # no SNPs were removed, so there are still 215681 SNPs
write.table(ccoru_allelefreq_for_lfmm, 'ccoru_Elev400m_allelefreqs_forLFMM.txt', row.names=FALSE, col.names=FALSE)

ccoru_allelefreq_for_lfmm_300m <- acast(ccoru_allelefreqs_all_300m, thebin~theSNP, value='littlea_freq')
dim(ccoru_allelefreq_for_lfmm_300m) # 5 SNPs were removed, so there are 215676 SNPs
write.table(ccoru_allelefreq_for_lfmm_300m, 'ccoru_Elev300m_allelefreqs_forLFMM.txt', row.names=FALSE, col.names=FALSE)

# since some SNPs had to be removed because some bins only had missing genotype data, so no allele frequency
# could be calculated, need to make a new list of the order of SNPs and their names
# this should simply be extracted from the allele frequency $theSNP file
cviol_lfmm_snpnames <- colnames(cviol_allelefreqs_for_lfmm); length(cviol_lfmm_snpnames) #91119
ccoru_lfmm_snpnames <- colnames(ccoru_allelefreq_for_lfmm); length(ccoru_lfmm_snpnames) #215681
ccoru_lfmm_snpnames_300m <- colnames(ccoru_allelefreq_for_lfmm_300m); length(ccoru_lfmm_snpnames_300m) #215676

write.table(cviol_lfmm_snpnames, 'cviol_lfmmallelefreq_snpnames.txt', row.names=FALSE, col.names=FALSE)
write.table(ccoru_lfmm_snpnames, 'ccoru_lfmmallelefreq_snpnames.txt', row.names=FALSE, col.names=FALSE)
write.table(ccoru_lfmm_snpnames_300m, 'ccoru_lfmmallelefreq_snpnames_300m.txt', row.names=FALSE, col.names=FALSE)

# make elevation files
# calculate average and std of all elevations
# calculate average elevation per bin
# calculate standardized elevation per bin
ccoru_elev <- read.table('../../Pcadapt_outflank/ccoru_outflank_elevation.txt', header=TRUE)
cviol_elev <- read.table('../../Pcadapt_outflank/cviol_outflank_elevation.txt', header=TRUE)

par(mfrow=c(1,2))
hist(ccoru_elev$Elevation, breaks=5)
abline(v=mean(ccoru_elev$Elevation), col='tomato')
abline(v=median(ccoru_elev$Elevation), col='turquoise')
hist(cviol_elev$Elevation, breaks=5)
abline(v=mean(cviol_elev$Elevation), col='tomato')
abline(v=median(cviol_elev$Elevation), col='turquoise')
# median tracks the higher frequency bin better than mean
# distributions definitely not normal...

# ccoru 400 m bins, cviol 400 m bins
ccoru_elev_bins <- c(2000,2400,2800,3200,3600,4000)
cviol_elev_bins <- c(2200,2600,3000,3400,3800)

# ccoru 300 m bins
ccoru_elev_bins300 <- c(1950, 2250, 2550, 2850, 3150, 3450, 3750, 3950)

# this is for LFMM
ccoru_meanstandardized_elev <- sapply(ccoru_elev_bins, function(x){
  mean_elev <- mean(ccoru_elev$Elevation)
  std_elev <- sd(ccoru_elev$Elevation)
  standardized_elev <- (x - mean_elev)/std_elev
})
# ccoru_medianstandardized_elev <- sapply(ccoru_elev_bins, function(x){
#   median_elev <- median(ccoru_elev$Elevation)
#   std_elev <- sd(ccoru_elev$Elevation)
#   standardized_elev <- (x - median_elev)/std_elev
# })

ccoru_meanstandardized_elev_300m <- sapply(ccoru_elev_bins300, function(x){
  mean_elev <- mean(ccoru_elev$Elevation)
  std_elev <- sd(ccoru_elev$Elevation)
  standardized_elev <- (x - mean_elev)/std_elev
})
# ccoru_medianstandardized_elev_300m <- sapply(ccoru_elev_bins300, function(x){
#   median_elev <- median(ccoru_elev$Elevation)
#   std_elev <- sd(ccoru_elev$Elevation)
#   standardized_elev <- (x - median_elev)/std_elev
# })

cviol_meanstandardized_elev <- sapply(cviol_elev_bins, function(x){
  mean_elev <- mean(cviol_elev$Elevation)
  std_elev <- sd(cviol_elev$Elevation)
  standardized_elev <- (x - mean_elev)/std_elev
})
# cviol_medianstandardized_elev <- sapply(cviol_elev_bins, function(x){
#   median_elev <- median(cviol_elev$Elevation)
#   std_elev <- sd(cviol_elev$Elevation)
#   standardized_elev <- (x - median_elev)/std_elev
# })

jpeg('standardizedelev_400m_meanVSmedian.jpg', height=4, width=10, units='in', res=500)
par(mfrow=c(1,4))
hist(ccoru_meanstandardized_elev, col='turquoise', main='Ccoru standardized w/ mean')
hist(ccoru_medianstandardized_elev, col='tomato', main='Ccoru standardized w/ median')
hist(cviol_meanstandardized_elev, col='turquoise', main='Cviol standardized w/ mean')
hist(cviol_medianstandardized_elev, col='tomato', main='Cviol standardized w/ median')
dev.off()

par(mfrow=c(1,2))
hist(ccoru_meanstandardized_elev_300m, col='turquoise', main='Ccoru standardized w/ mean')
hist(ccoru_medianstandardized_elev_300m, col='tomato', main='Ccoru standardized w/ median')

write.table(ccoru_meanstandardized_elev, 'Ccoru_meanstandardized_elevationbins400m.txt', col.names=FALSE, row.names=FALSE)
write.table(ccoru_meanstandardized_elev_300m, 'Ccoru_meanstandardized_elevationbins300m.txt', col.names=FALSE, row.names=FALSE)
write.table(cviol_meanstandardized_elev, 'Cviol_meanstandardized_elevationbins400m.txt', col.names=FALSE, row.names=FALSE)

# -- 2. check p-value distribution, lambda, & calculate q-value ----
# this version calculates lambda
lfmm_qvalcalc <- function(mypath, myK, myiters, mysp){
  z.table = NULL

  #loop through each file and cbind the zscores
  n = 5
  for (i in 1:n){
    file.name = paste(mypath, myK, myiters, i, "_s1.", myK, ".zscore", sep="")
    z.table = cbind(z.table, read.table(file.name)[,1])
  }

  #calculate median for z.scores arosss different runs for each site
  z.score = apply(z.table, MARGIN = 1, median)

  #compute lambda (genome inflation factor), the closer to one the better
  # lambda < 1 is too conservative
  lambda = median(z.score^2) / 0.456

  #print lambda

  ap.values = pchisq(z.score^2 / lambda, df = 1, lower = F)
  # plot the p-values
  jpeg(paste('LFMM_adjustedpvalue_hist_', mysp, 'K', myK, '_', lambda, myiters, '.jpg', sep=''), height=6, width=6, units='in', res=600)
  hist(ap.values, col = "red", main=paste('K: ', myK, ', Lambda: ', round(lambda,2), sep=''))
  dev.off()
  # calculate q-values
  qobj <- qvalue(ap.values, fdr.level=0.05, pi0.method = "bootstrap")
  write.qvalue(qobj, file=paste(mysp, "_K",myK, '_', lambda, myiters,".txt", sep=''))
  summary(qobj)

  totalSNPs <- length(qobj$significant) #215,681 for ccoru
  sig_SNPs <- length(qobj[qobj$significant == TRUE]) #21,638
  nonsig_SNPs <- length(qobj[qobj$significant == FALSE]) #194,043

  mytable <- table(lambda, totalSNPs, sig_SNPs, nonsig_SNPs)
  return(mytable)

}
# this version you set the lambda value
lfmm_qvalcalc_chooselambda <- function(mypath, myK, myiters, mylambda, mysp){
  z.table = NULL

  #loop through each file and cbind the zscores
  n = 5 # number of LFMM iterations
  for (i in 1:n){
    file.name = paste(mypath, myK, myiters, i, "_s1.", myK, ".zscore", sep="")
    z.table = cbind(z.table, read.table(file.name)[,1])
  }

  #calculate median for z.scores arosss different runs for each site
  z.score = apply(z.table, MARGIN = 1, median)

  #compute lambda (genome inflation factor), the closer to one the better
  # lambda < 1 is too conservative
  lambda = mylambda

  #print lambda

  ap.values = pchisq(z.score^2 / lambda, df = 1, lower = F)
  # plot the p-values
  jpeg(paste('LFMM_adjustedpvalue_hist_', mysp, 'K', myK, '_', lambda, myiters, '.jpg', sep=''), height=6, width=6, units='in', res=600)
  hist(ap.values, col = "red", main=paste('K: ', myK, ', Lambda: ', round(lambda,2), sep=''))
  dev.off()
  # calculate q-values
  qobj <- qvalue(ap.values, fdr.level=0.05, pi0.method = "bootstrap")
  write.qvalue(qobj, file=paste(mysp, "_K",myK, '_', lambda, myiters,".txt", sep=''))
  summary(qobj)

  totalSNPs <- length(qobj$significant) #215,681
  sig_SNPs <- length(qobj[qobj$significant == TRUE]) #21,638
  nonsig_SNPs <- length(qobj[qobj$significant == FALSE]) #194,043

  mytable <- table(lambda, totalSNPs, sig_SNPs, nonsig_SNPs)
  return(mytable)

}

# Ccoru LFMM results 400m bin
lfmm_qvalcalc('./Ccoru_lfmm_allelefreqinputs/ccorubinned400m_K', '2', "_i250000b25000_run", 'Ccoru')
lfmm_qvalcalc_chooselambda('./Ccoru_lfmm_allelefreqinputs/ccorubinned400m_K', '2', '_i250000b25000_run', 1, 'Ccoru')
lfmm_qvalcalc_chooselambda('./Ccoru_lfmm_allelefreqinputs/ccorubinned400m_K', '2', '_i250000b25000_run', 0.8, 'Ccoru')

# ok, K=2 lambda result makes MUCH more sense than K3 lambda!
lfmm_qvalcalc('./Cviol_lfmm_allelefreqinputs/cviolbinned_K', '2', '_i250000b25000_run', 'Cviol')
lfmm_qvalcalc_chooselambda('./Cviol_lfmm_allelefreqinputs/cviolbinned_K', '2', '_i250000b25000_run', 1, 'Cviol')
lfmm_qvalcalc_chooselambda('./Cviol_lfmm_allelefreqinputs/cviolbinned_K', '2', '_i250000b25000_run', 0.8, 'Cviol')

# old thresholds
# # Ccoru LFMM results 400 m bin
# lfmm_qvalcalc('./Ccoru_lfmm_allelefreqinputs/ccorubinned_K', '1', "_i250000b25000_run", 'Ccoru')
# #lfmm_qvalcalc_chooselambda('./Ccoru_lfmm_allelefreqinputs/ccorubinned_K', '1', '_i250000b25000_run', 0.5, 'Ccoru')
# #lfmm_qvalcalc_chooselambda('./Ccoru_lfmm_allelefreqinputs/ccorubinned_K', '1', '_i250000b25000_run', 0.8, 'Ccoru')
# lfmm_qvalcalc_chooselambda('./Ccoru_lfmm_allelefreqinputs/ccorubinned_K', '1', '_i250000b25000_run', 1, 'Ccoru')

# # Ccoru LFMM results 300 m bin
# lfmm_qvalcalc('./Ccoru_lfmm_allelefreqinputs_300m/ccorubinned300m_K', '1', '_i250000b25000_run', 'Ccoru_300m')
# lfmm_qvalcalc_chooselambda('./Ccoru_lfmm_allelefreqinputs_300m/ccorubinned300m_K', '1', '_i250000b25000_run', 1, 'Ccoru_300m')
#
# lfmm_qvalcalc('./Ccoru_lfmm_allelefreqinputs_300m/ccorubinned300m_K', '2', '_i250000b25000_run', 'Ccoru_300m')
# lfmm_qvalcalc_chooselambda('./Ccoru_lfmm_allelefreqinputs_300m/ccorubinned300m_K', '2', '_i250000b25000_run', 1, 'Ccoru_300m')
# lfmm_qvalcalc_chooselambda('./Ccoru_lfmm_allelefreqinputs_300m/ccorubinned300m_K', '2', '_i250000b25000_run', 0.8, 'Ccoru_300m')

# # Cviol LFMM results 400 m bin
# lfmm_qvalcalc('./Cviol_lfmm_allelefreqinputs/cviolbinned_K', '3', '_i250000b25000_run', 'Cviol')
# lfmm_qvalcalc_chooselambda('./Cviol_lfmm_allelefreqinputs/cviolbinned_K', '3', '_i250000b25000_run', 1, 'Cviol')
# lfmm_qvalcalc_chooselambda('./Cviol_lfmm_allelefreqinputs/cviolbinned_K', '3', '_i250000b25000_run', 2000, 'Cviol')

# ------- Conclusions about version to use for enrichment analysis -----
# Ccoru: 400m
# K=1
# lambda=0.5 is too liberal now (even more outliers than I had with the genotype input K=1 lambda=0.22)
# lambda=1.17 is a bit too conservative now (Ke thinks 649 outliers is a bit low)
# will look at GO and KEGG for lambda=1 and 1.17
# K=2 - use this

# Ccoru: 300m
# K=2 less stringent than K1

# Cviol: 400m
# K=3 results are WEIRD
# K=2 results look much more normal (in terms of lambda values) - use this

# For Ccoru 400m bin K=2 and Cviol 400m bin K=2, the calculated lambda results are more conservative
# I will look at the lambda=1 results for outlier SNPs


# -- 3. Format LFMM qvalue df and annotate genes for enrichment analysis ----

### This part formats the annotation files ###
# Run this part and then use the annotate_snp() function I wrote below to annotate LFMM SNPs
# adding annotations to significant and non-significant SNPs because for GO enrichment analysis, I need a comparison to the significant SNPs
# load contig name text file
contignames_Annarefs <- read.table('../../ref.fasta_rename_compared_Annarefs.txt')
colnames(contignames_Annarefs) <- c('Contig', 'Attributes', 'Gene')
contignames_Annarefs$Contig <- gsub('>', '', contignames_Annarefs$Contig)

contignames_transcriptrefs <- read.table('../../ref.fast_rename_compared_transcriptrefs.txt')
colnames(contignames_transcriptrefs) <- c('Contig', 'Gene')
contignames_transcriptrefs$Attributes <- 'NA' #set this to NA, need this column so I can combine with the Annarefs
contignames_transcriptrefs$Contig <- gsub('>', '', contignames_transcriptrefs$Contig)
#dim(contignames_transcriptrefs)
# need to remove [Letter].ENSTGUP[ensembleID]. from gene refs that came from transcriptome
contignames_transcriptrefs_torename <- contignames_transcriptrefs[c(232,235:269),]
transcriptref_genename <- gsub(pattern='\\.', replacement=' ', x=contignames_transcriptrefs_torename$Gene)
transcriptref_genename2 <- strsplit(x=transcriptref_genename, split=' ')
transcriptref_genename3 <- sapply(transcriptref_genename2, function(x){
  paste(x[[3]], sep="")
})
contignames_transcriptrefs_torename$Gene <- transcriptref_genename3
#now, reglue this df together
contignames_transcriptrefs2 <- rbind(contignames_transcriptrefs[c(1:231),], contignames_transcriptrefs_torename[1,],
                                     contignames_transcriptrefs[c(233:234),], contignames_transcriptrefs_torename[c(2:36),],
                                     contignames_transcriptrefs[c(270:271),])
#some are just easier to change manually -_-
contignames_transcriptrefs2$Gene <- as.character(contignames_transcriptrefs2$Gene)
contignames_transcriptrefs2$Gene[contignames_transcriptrefs2$Gene == 'G.HBB'] <- 'HBB'
contignames_transcriptrefs2$Gene[contignames_transcriptrefs2$Gene == 'G.HBAD'] <- 'HBAD'
contignames_transcriptrefs2$Gene[contignames_transcriptrefs2$Gene == 'F.ENSTGUP00000018304.COX1.MITOCHONDRIAL'] <- 'COX1'
contignames_transcriptrefs2$Gene[contignames_transcriptrefs2$Gene == 'F.ENSTGUP00000018308.COX3.MITOCHONDRIAL'] <- 'COX3'
contignames_transcriptrefs2$Gene <- as.factor(contignames_transcriptrefs2$Gene)
#dim(contignames_transcriptrefs2) #should be 271

# Further edits to gene names
# e.g., some are [NAME]_CHICK which refer to the chicken genome, but do not match a HGNC human genome gene name automatically
# note: there are other genes that don't match HGNC naming system, but not sure how to systematically edit them..yet!
# will have to check as we go...
# e.g., some of the XP_[numbers] genes are searchable in genbank..but not all
# some of the CHICK gene names are a HGNC gene name with '_CHICK' added to the end
# for these, can just remove the '_CHICK' part
# while others are a Uniprot database name with '_CHICK' added to the end
# for these, remove '_CHICK' part and then have to search Uniprot database to get corresponding gene symbol
# the org.Gg.eg.db sort of works - but not for all entries, it does require removing the _CHICK part
# this is the annotation db for Gallus gallus
# Loop:
# 1. If '_CHICK' not in gene name, keep original gene symbol
# 2. If '_CHICK' in gene name:
# a. remove '_CHICK'
# b. If uniprot symbol not available, keep gene name without the '_CHICK'
# c. If uniprot symbol is available, use it for gene name
contignames_transcriptrefs_rename <- data.frame('Contig'='contig', 'Gene'='gene', 'Attributes'='attributes')
for(i in 1:nrow(contignames_transcriptrefs2)){
  tryCatch({
    therow <- contignames_transcriptrefs2[i,]
    print(i)
    if(length(grep('_CHICK', therow$Gene)) == 1){
      #print(contignames_transcriptrefs$Gene[i])
      # take out _CHICK
      noCHICK_gene <- sapply(therow$Gene, function(x) gsub('_CHICK', '', x))

      # search for corresponding gene name
      # use the try() function to catch errors since select() fails if no matching symbol for uniprot key
      # the next if statement will evaluate whether there was an error
      uniprot_results <- try(select(org.Gg.eg.db, keys=noCHICK_gene, columns=c('UNIPROT', 'SYMBOL'), keytype='UNIPROT'), TRUE)

      if(isTRUE(class(uniprot_results) == 'try-error')){
        # if there is no match, keep the name with _CHICK removed
        print('no match :(')
        therow$Gene <- noCHICK_gene
        contignames_transcriptrefs_rename <- rbind(contignames_transcriptrefs_rename, therow)

      }else{
        # else there is a match, replace uniprot name with gene symbol
        print('a match!')
        therow$Gene <- uniprot_results$SYMBOL
        contignames_transcriptrefs_rename <- rbind(contignames_transcriptrefs_rename, therow)
      }
    }else{
      # no CHICK in name, so just keep row unchanged
      print('no changes to report!')
      contignames_transcriptrefs_rename <- rbind(contignames_transcriptrefs_rename, therow)
    }
  }, error=function(e){cat('ERROR: ', conditionMessage(e), '\n')})
}
contignames_transcriptrefs_rename2 <- contignames_transcriptrefs_rename[-1,] #remove empty 1st column

# The transcriptome ref genes already have the gene annotation
# but the Anna's hbird genome ref genes have scaffold and attribute names
# now, have to go back to the .gff genome file to find the annotations for these
# Use file with attribute and annotation columns from findcandidates.py and Calypte_genome_exoncheck [was for probe design step]
# to match against the Attribute column in *_outliercontigs dfs
# read in attribute + annotation file
anna_genome_annots <- read.csv('../mrna_attribute_genename.csv', header=T)
#head(anna_genome_annots)
#dim(anna_genome_annots)
#length(grep('_CHICK', anna_genome_annots$Gene_name)) #1614 of 16,000
#length(grep('XP_', anna_genome_annots$Gene_name)) #570 of 16,000
# OK, now do the same _CHICK search as above for the anna refs:
anna_genome_annots_rename <- data.frame('Attributes'='attributes', 'Gene_name'='gene')
for(i in 1:nrow(anna_genome_annots)){
  tryCatch({
    therow <- anna_genome_annots[i,]
    print(i)
    if(length(grep('_CHICK', therow$Gene_name)) == 1){
      # take out _CHICK
      noCHICK_gene <- sapply(therow$Gene_name, function(x) gsub('_CHICK', '', x))

      # search for corresponding gene name
      # use the try() function to catch errors since select() fails if no matching symbol for uniprot key
      # the next if statement will evaluate whether there was an error
      uniprot_results <- try(select(org.Gg.eg.db, keys=noCHICK_gene, columns=c('UNIPROT', 'SYMBOL'), keytype='UNIPROT'), TRUE)

      if(isTRUE(class(uniprot_results) == 'try-error')){
        # if there is no match, keep the name with _CHICK removed
        print('no match :(')
        therow$Gene_name <- noCHICK_gene
        anna_genome_annots_rename <- rbind(anna_genome_annots_rename, therow)

      }else{
        # else there is a match, replace uniprot name with gene symbol
        print('a match!')
        therow$Gene_name <- uniprot_results$SYMBOL
        anna_genome_annots_rename <- rbind(anna_genome_annots_rename, therow)
      }
    }else{
      # no CHICK in name, so just keep row unchanged
      print('no changes to report!')
      anna_genome_annots_rename <- rbind(anna_genome_annots_rename, therow)
    }
  }, error=function(e){cat('ERROR: ', conditionMessage(e), '\n')})
}
anna_genome_annots_rename2 <- anna_genome_annots_rename[-1,] #remove empty 1st column
#  head(anna_genome_annots_rename2)
#  dim(anna_genome_annots_rename2)

# at the end of these gene name edits, there are still going to be a bunch that don't match the
# HGNC naming framework, but at least I can retain more than before..

### now, back to the task of adding the annotations to LFMM outlier dataset ###
# combine the Annarefs and transcript refs (transcript ref gene names have been edited above)
#head(contignames_Annarefs) #note that these have attribute names, no gene names. need to get gene names from mrna_attribute_genename.csv
#head(contignames_transcriptrefs_rename2) #note that these have gene names!
all_contignames <- rbind(contignames_transcriptrefs_rename2, contignames_Annarefs)
#dim(all_contignames) #12501 loci
#head(all_contignames)

annotate_snps <- function(mylfmmqvals, myloci, contig_annot, anna_annot,  myspecies){
  lfmmqvalfile <- read.table(mylfmmqvals, header=T)
  theloci <- read.table(myloci)
  #add loci to qval df
  lfmmqvalfile$loci <- theloci$V1

  #how many significant SNPs in this dataset?
  print(paste('There are ', nrow(lfmmqvalfile[lfmmqvalfile$significant == 1,]), ' significant outlier SNPs', sep=''))

  # LFMM outliers --
  # first make a column with just contig name in LFMM results dataframe from $loci column
  # then match those contigs with corresponding gene names from all_contignames df
  # then make a column with the attribute (R#) name
  # SNPs are linked to a Contig number: combined_Contig#_#
  # the gene annotation for transcriptome ref genes are in the all_contignames df
  # but the gene name for genes targeted with Anna's genome (rather than transcriptomes) have a scaffold#_attribute#_start_stop
  # rather than a gene name. To get the gene name, go to Anna's genome mRNA .gff file and extract annotation
  #testdat <- ccoru_K1_signif[1:100,]

  getLFMM_contigs <- function(LFMM_outlierfile){
    contig_list <- list()
    for(i in 1:nrow(LFMM_outlierfile)){
      thelocus <- strsplit(as.character(LFMM_outlierfile$loci[i]), '_')
      # Make column with contig #
      thecontig <- paste(thelocus[[1]][2])
      contig_list[[i]] <- thecontig
    }
    the_outliercontigs <- data.frame('Contig'=unlist(contig_list))
    lfmmdat1 <- cbind(LFMM_outlierfile, the_outliercontigs)

    #match Contig name gene annots
    lfmmdat2 <- merge(lfmmdat1, contig_annot, by='Contig')

    #since Anna refs don't have annotations
    #add in the annotations using attribute column
    lfmmdat3 <- merge(x=lfmmdat2, y=anna_annot, by='Attributes', all.x=TRUE)

    return(lfmmdat3)
  }

  print('Now, annotate LFMM contigs...')
  LFMM_outliercontigs <- getLFMM_contigs(lfmmqvalfile)

  # Since the Contig annotations and the attribute annotations are now in separate columns,
  # we need to put them all into a single column
  # genes with Contig annotations should have '<NA>' in the $Attribute column
  # separate the df by the $Gene_name column
  # then repaste together but with single gene name column
  contigannot <- LFMM_outliercontigs[LFMM_outliercontigs$Attributes == 'NA',]
  contigannot2 <- contigannot[,c(1:10)]
  attributeannot <- LFMM_outliercontigs[LFMM_outliercontigs$Attributes != 'NA',]
  attributeannot2 <- attributeannot[,c(1:9,11)]
  colnames(attributeannot2)[10] <- 'Gene'
  # stick them back together
  lfmm_annotated <- rbind(contigannot2, attributeannot2)
  print('Done annotating LFMM contigs!')

  # how many SNPs are not annotated?
  print(paste(length(which(lfmm_annotated$Gene == '')), ' SNPs are NOT annotated total', sep=''))
  print(paste(length(which(lfmm_annotated$significant == 1 & lfmm_annotated$Gene == '')),
              ' significant SNPs are not annotated', sep=''))

  # how many unique genes have significant outlier SNPs and are annotated?
  print(paste(length(unique(lfmm_annotated[lfmm_annotated$significant == 1 & lfmm_annotated$Gene != '',]$Gene)),
        ' unique genes with annotated significant SNPs',  sep=''))
  # let's save these annotated SNP files so i don't have to keep rerunning above code:
  write.table(lfmm_annotated, paste('LFMM_SNPs_annotated_', myspecies, '.txt', sep=''))

}

### Ccoru 400 m bins ###
# annotate_snps('Ccoru_K1_1.17272294534211_i250000b25000_run.txt', 'ccoru_lfmmallelefreq_snpnames.txt',
#               all_contignames, anna_genome_annots_rename2,
#               'CcoruK1lambda1.17')
# [1] "There are 649 significant outlier SNPs"
# [1] "11291 SNPs are NOT annotated total"
# [1] "31 significant SNPs are not annotated"
# [1] "472 unique genes with annotated significant SNPs"
# annotate_snps('Ccoru_K1_1_i250000b25000_run.txt', 'ccoru_lfmmallelefreq_snpnames.txt',
#               all_contignames, anna_genome_annots_rename2,
#               'CcoruK1lambda1')
# [1] "There are 2021 significant outlier SNPs"
# [1] "11291 SNPs are NOT annotated total"
# [1] "95 significant SNPs are not annotated"
# [1] "1263 unique genes with annotated significant SNPs"

annotate_snps('Ccoru_K2_1.31057743359868_i250000b25000_run.txt', 'ccoru_lfmmallelefreq_snpnames.txt',
              all_contignames, anna_genome_annots_rename2,
              'CcoruK2lambda1.31')
# [1] "There are 622 significant outlier SNPs"
# [1] "11291 SNPs are NOT annotated total"
# [1] "27 significant SNPs are not annotated"
# [1] "448 unique genes with annotated significant SNPs"
annotate_snps('Ccoru_K2_1_i250000b25000_run.txt', 'ccoru_lfmmallelefreq_snpnames.txt',
              all_contignames, anna_genome_annots_rename2,
              'CcoruK2lambda1')
# [1] "There are 3962 significant outlier SNPs"
# [1] "11291 SNPs are NOT annotated total"
# [1] "219 significant SNPs are not annotated"
# [1] "2136 unique genes with annotated significant SNPs"

### Ccoru 300 m bins ###
# annotate_snps('Ccoru_300m_K1_1.12827042913706_i250000b25000_run.txt', 'ccoru_lfmmallelefreq_snpnames_300m.txt',
#               all_contignames, anna_genome_annots_rename2, 'CcoruK1lambda1.13_300m')
# [1] "There are 559 significant outlier SNPs"
# [1] "11291 SNPs are NOT annotated total"
# [1] "23 significant SNPs are not annotated"
# [1] "411 unique genes with annotated significant SNPs"
# annotate_snps('Ccoru_300m_K1_1_i250000b25000_run.txt', 'ccoru_lfmmallelefreq_snpnames_300m.txt',
#               all_contignames, anna_genome_annots_rename2, 'CcoruK1lambda1_300m')
# [1] "There are 1373 significant outlier SNPs"
# [1] "11291 SNPs are NOT annotated total"
# [1] "71 significant SNPs are not annotated"
# [1] "933 unique genes with annotated significant SNPs"

### Cviol 400 m bins ###
# annotate_snps('Cviol_K3_3186.16325052632_i250000b25000_run.txt', 'cviol_lfmmallelefreq_snpnames.txt',
#               all_contignames, anna_genome_annots_rename2, 'CviolK3lambda3186')
# [1] "There are 408 significant outlier SNPs"
# [1] "5086 SNPs are NOT annotated total"
# [1] "24 significant SNPs are not annotated"
# [1] "312 unique genes with annotated significant SNPs"

# annotate_snps('Cviol_K3_2000_i250000b25000_run.txt', 'cviol_lfmmallelefreq_snpnames.txt',
#               all_contignames, anna_genome_annots_rename2, 'CviolK3lambda2000')
# [1] "There are 3807 significant outlier SNPs"
# [1] "5086 SNPs are NOT annotated total"
# [1] "186 significant SNPs are not annotated"
# [1] "2126 unique genes with annotated significant SNPs"

annotate_snps('Cviol_K2_1.74352771716009_i250000b25000_run.txt', 'cviol_lfmmallelefreq_snpnames.txt',
              all_contignames, anna_genome_annots_rename2, 'CviolK2lambda1.74')
# [1] "There are 526 significant outlier SNPs"
# [1] "5086 SNPs are NOT annotated total"
# [1] "26 significant SNPs are not annotated"
# [1] "380 unique genes with annotated significant SNPs"

annotate_snps('Cviol_K2_1_i250000b25000_run.txt', 'cviol_lfmmallelefreq_snpnames.txt',
              all_contignames, anna_genome_annots_rename2, 'CviolK2lambda1')
# [1] "There are 6433 significant outlier SNPs"
# [1] "Now, annotate LFMM contigs..."
# [1] "Done annotating LFMM contigs!"
# [1] "5086 SNPs are NOT annotated total"
# [1] "291 significant SNPs are not annotated"
# [1] "3099 unique genes with annotated significant SNPs"

# -- 4. GO and KEGG enrichment ----
# Note: in previous version, I used unique genes, but I think I'm supposed to be using
# top 500ish unique SNPs even if there are repeat genes (if I subset by top X number of SNPs)
# but since I have fewer SNPs with higher lambda, I'm going to look at all the significant SNPs

GO_KEGG_enrich <- function(myannotlfmm, myspecies){
  thedat <- read.table(myannotlfmm, header=TRUE)
  # remove rows with no $Gene annotation because won't be able to assign GO category to it
  theannotateddat <- thedat[thedat$Gene != '',]
  totalSNPs_NOT_annot <- nrow(thedat[thedat$Gene == '',])
  print(paste('Removed ', totalSNPs_NOT_annot, ' un-annotated SNPs from analysis', sep=''))
  sigSNPs_NOT_annot <- nrow(thedat[thedat$Gene == '' & thedat$significant == 1,])
  print(paste('Removed ', sigSNPs_NOT_annot, ' un-annotated significant SNPs from analysis', sep=''))

  # GO enrichment:
  # significant outiers $significant = 1
  # background set = all the non-sig outliers ($significant = 0)
  # genes that have both significant and non-significant SNPs retained in significant outlier category

  ######## CHECK THIS ##############
  # NOTE: nullp() function does NOT allow duplicate gene names (duplicate row names)
  theoutliers <- as.vector(droplevels(unique(theannotateddat[theannotateddat$significant == 1, ]$Gene)))
  bg0 <- as.vector(droplevels(unique(theannotateddat[theannotateddat$significant == 0, ]$Gene))) #all non-sig outliers

  # what to do about genes that are in both outlier and bg set? it means that there are multiple SNPs in a gene - some sig, some not
  # i guess remove the gene from the bg set because it does have sig SNPs
  # later, can go back to the specific SNPs and check out the sequence if interesting to follow up
  # this means adding a step to check if a gene is in both theoutliers and bg0.
  # if it is in both, remove it from the bg0 vector

  # find list of shared genes, remove these from bg0
  gene_overlaps <- intersect(theoutliers, bg0)
  bg1 <- bg0[! bg0 %in% gene_overlaps]
  #length(bg1)
  #length(intersect(theoutliers, bg1)) #this should now be 0 overlaps
  ######################################

  allgenes <- c(theoutliers, bg1)
  gene.vector <- as.integer(allgenes%in%theoutliers)
  names(gene.vector) <- allgenes
  thepwf=nullp(gene.vector, 'hg19', 'geneSymbol')

  # Wallenius approximation method
  GO.wall=goseq(thepwf, 'hg19', 'geneSymbol')

  # Edit if you prefer random sampling method...
  # Random sampling method
  #GO.samp=goseq(thepwf, 'hg19', 'geneSymbol', method='Sampling', repcnt=1000)
  #head(GO.samp)
  # compare Wallenius vs. sampling method
  # they are very similar (fall on the 0,1 line)
  #plot(log10(GO.wall[,2]), log10(GO.samp[match(GO.wall[,1], GO.samp[,1]),2]),
  #     xlab='log10(Wallenius p-values)', ylab='log10(Sampling p-values)',
  #     xlim=c(-3,0))
  #abline(0,1, col=3,lty=2)

  # making sense of the results after BH correction for multiple tests
  #The first column gives the name of the category, the second gives the p-value for the associated
  #category being over represented amongst DE genes. The third column gives the p-value for the associated
  #category being under represented amongst DE genes. ***The p-values have not been corrected for multiple
  #hypothesis testing.*** The fourth and fifth columns give the number of differentially expressed genes in the
  #category and total genes in the category respectively.
  #If any of the categories was a GO term, there will be two additional columns for the GO term and
  #its ontology.
  enriched.GO <- GO.wall$category[p.adjust(GO.wall$over_represented_pvalue, method='BH') <0.05]
  print(paste(length(enriched.GO), ' enriched GO terms', sep=''))
    # NONE pass the BH correction of p < 0.05 - what to do? (only 1 GO term pass with p < 0.1!)


  # KEGG analysis
  # the kegg categories refer to kegg pathway maps
  # create df of kegg pathways and map numbers
  # the mappings below for kegg are from KEGG.db which i believe is based (mostly/all??) on human pathway information
  kegg_path_maps <- as.data.frame(keggList('pathway'))
  kegg_path_maps$KEGG_map_cat <- rownames(kegg_path_maps)
  rownames(kegg_path_maps) <- NULL
  colnames(kegg_path_maps) <- c('KEGG_pathway', 'KEGG_map_cat')
  kegg_path_maps$KEGG_map_number <- sapply(kegg_path_maps$KEGG_map_cat, function(x) substr(x, 9, 13))

  KEGG <- goseq(thepwf, 'hg19', 'geneSymbol', test.cats='KEGG')

  # make df of the kegg pathways that have significant over_represented_pvalue
  enriched.kegg <- KEGG$category[p.adjust(KEGG$over_represented_pvalue, method='BH') <0.05]
  print(paste(length(enriched.kegg), ' enriched KEGG terms', sep=''))

}

GO_KEGG_enrich('LFMM_SNPs_annotated_CcoruK1lambda1.17.txt', 'CcoruK1lambda1.17')
GO_KEGG_enrich('LFMM_SNPs_annotated_CcoruK1lambda1.txt', 'CcoruK1lambda1')


# not finding any enriched categories...too spread over different categories?


# -- 5. Categorize top SNPs in previously known cand genes vs. new ----
# maybe one way to do this is to count the number of
# 1) previously identified candidate genes with significant SNPs (q-value < 0.05)
      # compare gene names with sig SNPs to list of candidate genes targeted in exon capture
# 2) new candidate genes with significant SNPs that seem to have some relation to elevation adaptation (q-value < 0.05)
      # the remaining genes that aren't part of 1). but will be more work to parse out their functions...
# Make plots for any that are interesting, but not all 4000...
    # minor allele frequency for each elevation bin
    # map of genotypes
    # genotype x elevation boxplots + wilcoxon rank sum test

# read in file with targeted candidate genes
cand_targs <- read.csv('C:/Users/mcwlim/Dropbox/Marisacompfiles/Seq cap files/mycandgenes.csv', header=T)

# ------- 5a) Ccoru 400m bins coeffs ----------
# read in the pearson corr coefficients for allele freqs per elev bin
# to merge with and compare to LFMM q-values
ccoru_coef1 <- read.table('/2nd pass ANGSD - remove outliers/Allelefreqplots_byelevation/Ccoru_Elevation_bin400m_snpbatch1_Pearsoncor_allelefreq-bin.txt', header=TRUE)
ccoru_coef2 <- read.table('/2nd pass ANGSD - remove outliers/Allelefreqplots_byelevation/Ccoru_Elevation_bin400m_snpbatch2_Pearsoncor_allelefreq-bin.txt', header=TRUE)
ccoru_coef3 <- read.table('/2nd pass ANGSD - remove outliers/Allelefreqplots_byelevation/Ccoru_Elevation_bin400m_snpbatch3_Pearsoncor_allelefreq-bin.txt', header=TRUE)
ccoru_coef4 <- read.table('/2nd pass ANGSD - remove outliers/Allelefreqplots_byelevation/Ccoru_Elevation_bin400m_snpbatch4_Pearsoncor_allelefreq-bin.txt', header=TRUE)
ccoru_coef5 <- read.table('/2nd pass ANGSD - remove outliers/Allelefreqplots_byelevation/Ccoru_Elevation_bin400m_snpbatch5_Pearsoncor_allelefreq-bin.txt', header=TRUE)
ccoru_coef6 <- read.table('/2nd pass ANGSD - remove outliers/Allelefreqplots_byelevation/Ccoru_Elevation_bin400m_snpbatch6_Pearsoncor_allelefreq-bin.txt', header=TRUE)
ccoru_coef7 <- read.table('/2nd pass ANGSD - remove outliers/Allelefreqplots_byelevation/Ccoru_Elevation_bin400m_snpbatch7_Pearsoncor_allelefreq-bin.txt', header=TRUE)
ccoru_coef8 <- read.table('/2nd pass ANGSD - remove outliers/Allelefreqplots_byelevation/Ccoru_Elevation_bin400m_snpbatch8_Pearsoncor_allelefreq-bin.txt', header=TRUE)
ccoru_coef9 <- read.table('/2nd pass ANGSD - remove outliers/Allelefreqplots_byelevation/Ccoru_Elevation_bin400m_snpbatch9_Pearsoncor_allelefreq-bin.txt', header=TRUE)
ccoru_coef10 <- read.table('/2nd pass ANGSD - remove outliers/Allelefreqplots_byelevation/Ccoru_Elevation_bin400m_snpbatch10_Pearsoncor_allelefreq-bin.txt', header=TRUE)
ccoru_coef11 <- read.table('/2nd pass ANGSD - remove outliers/Allelefreqplots_byelevation/Ccoru_Elevation_bin400m_snpbatch11_Pearsoncor_allelefreq-bin.txt', header=TRUE)

ccoru_coef_all <- rbind(ccoru_coef1, ccoru_coef2, ccoru_coef3, ccoru_coef4, ccoru_coef5, ccoru_coef6,
                        ccoru_coef7, ccoru_coef8, ccoru_coef9, ccoru_coef10, ccoru_coef11)
colnames(ccoru_coef_all)[2] <- 'loci'

# calculate q-value
qobj_ccoruPearson <- qvalue(ccoru_coef_all$p.value, fdr.level=0.05, pi0.method = "bootstrap")
summary(qobj_ccoruPearson)
ccoru_coef_all <- cbind(ccoru_coef_all, data.frame('Pearsoncor_qvalue'=qobj_ccoruPearson$qvalues))
head(ccoru_coef_all)
colnames(ccoru_coef_all)[4] <- 'Pearsoncor_pvalue'
# note: only 1 SNP survives after FDR correction for the Pearson's correlation results. The correlation is low-power, so let's focus on LFMM results

# ------- 5b) K=2 results, lambda = 1 ----
ccoru_annot_snps_lambda1 <- read.table('LFMM_SNPs_annotated_CcoruK2lambda1.txt', header=TRUE)
ccoru_coef_annot_lambda1 <- merge(ccoru_coef_all, ccoru_annot_snps_lambda1[,c(4,7,9,10)], by='loci')
ccoru_coef_sorted_lambda1 <- ccoru_coef_annot_lambda1[order(ccoru_coef_annot_lambda1$qvalue),]
head(ccoru_coef_sorted_lambda1)
# check a few genes...
# head(ccoru_coef_sorted_lambda1[ccoru_coef_sorted_lambda1$significant == 1,], 10) # look at top 10 significant lfmm snps that were also correlated with elev bin by allele freq
# length(ccoru_coef_sorted_lambda1[ccoru_coef_sorted_lambda1$significant == 1,]$Gene)
# head(ccoru_coef_sorted_lambda1[ccoru_coef_sorted_lambda1$Gene == 'EPAS1' & ccoru_coef_sorted$significant == 1,], 10) #none of the top corr'd are also signif...
# SENP1 <- ccoru_coef_sorted_lambda1[ccoru_coef_sorted_lambda1$Gene == 'SENP1' & ccoru_coef_sorted_lambda1$significant == 1,] #this has 2 hits [this gene was picked up in OutFLANK analysis]
# head(ccoru_coef_sorted_lambda1[ccoru_coef_sorted_lambda1$Gene == 'PPARA',], 10) #nope
# head(ccoru_coef_sorted_lambda1[ccoru_coef_sorted_lambda1$Gene == 'HIF1AN',], 10) #no
# RYR2_ccoru <- ccoru_coef_sorted_lambda1[ccoru_coef_sorted_lambda1$Gene == 'RYR2' & ccoru_coef_sorted_lambda1$significant == 1,]

# note that one of these SNPs is in a gene called Sep-11 (this is not a formatting error, that's how it's named in the C. anna gff file too)

# looking at which SNP survives FDR for Pearson's correlation results, continue with just LFMM
#ccoru_400mK2lamb1_05qval_Pearson <- ccoru_coef_sorted_lambda1[ccoru_coef_sorted_lambda1$Pearsoncor_qvalue < 0.05,]

ccoru_400mK2lamb1_05qval <- ccoru_coef_sorted_lambda1[ccoru_coef_sorted_lambda1$qvalue < 0.05,]
dim(ccoru_400mK2lamb1_05qval) #3962 significant SNPs left after FDR (meaning that all were kept...)
length(unique(ccoru_400mK2lamb1_05qval$Gene)) # SNPs on 2137 unique genes
head(ccoru_400mK2lamb1_05qval)

# there are close to 4000 SNPs with q-values < 0.05
# so the top 100 are extremely significant, but can still look at the other 3900...

# -- compare genes from LFMM with sig and non-sig SNPs to targeted cand genes ------
# make column for contig number and SNP site >> use this to check if SNP results in dN or dS amino acid change
ccoru_loci_split <- strsplit(as.character(ccoru_400mK2lamb1_05qval$loci), '_')
ccoru_400mK2lamb1_05qval$contig_num <- sapply(ccoru_loci_split, function(x){
  paste(x[[1]], x[[2]], sep='_')
})
ccoru_400mK2lamb1_05qval$SNP_site <- sapply(ccoru_loci_split, function(x){
  x[[3]]
})

colnames(cand_targs) <- 'Gene'
#ccoru_previouslyknowngenes <- merge(x=ccoru_400mK2lamb1_05qval, y=cand_targs, by.y='Gene')
  # this misses the genes that have odd names e.g., ACE_CHICK, NRG1_CHICK, but they are candidates!
  # this is why there are some candidate gene SNPs in the control genes file
  # however, since I know that combined_Contig1-271 are my candidate genes and 272+ are control genes,
  # it is better to subset by contig number instead of trying to match the gene name (and the variations)
# create new col to order contig numbers
ccoru_400mK2lamb1_05qval$contig_index <- gsub('combined_Contig','',ccoru_400mK2lamb1_05qval$contig_num)
# put in order by contig number
ccoru_400mK2lamb1_05qval_ordered0 <- ccoru_400mK2lamb1_05qval[order(as.numeric(ccoru_400mK2lamb1_05qval$contig_index)),]
# subset to get all contigs 1-271
ccoru_previouslyknowngenes <- ccoru_400mK2lamb1_05qval_ordered0[as.numeric(ccoru_400mK2lamb1_05qval_ordered0$contig_index) <= 271,]

print(paste('There are ', nrow(ccoru_previouslyknowngenes), ' significant SNPs in ', length(unique(ccoru_previouslyknowngenes$Gene)), ' previously identified candidate genes.', sep=''))
print(paste('There are ', nrow(ccoru_400mK2lamb1_05qval) - nrow(ccoru_previouslyknowngenes), ' significant SNPs in ',
            length(unique(ccoru_400mK2lamb1_05qval$contig_num)) - length(unique(ccoru_previouslyknowngenes$contig_num)),
            ' contigs that are either new candidates for elevation adaptation, or are unrelated -- must check gene function for relevance',
            sep='')) #note that these are contigs not genes - because a number of the contigs have no gene ID
# [1] "There are 757 significant SNPs in 174 previously identified candidate genes."
# [1] "There are 3205 significant SNPs in 2176 contigs that are either new candidates for elevation adaptation, or are unrelated -- must check gene function for relevance"

# check same but for the non-sig SNPs
ccoru_nonsig_SNPs <- ccoru_annot_snps_lambda1[ccoru_annot_snps_lambda1$significant == 0,]
dim(ccoru_nonsig_SNPs)
# create new col to order contig numbers
ccoru_nonsig_SNPs$contig_index <- gsub('Contig','',ccoru_nonsig_SNPs$Contig)
# put in order by contig number
ccoru_nonsig_SNPs_ordered0 <- ccoru_nonsig_SNPs[order(as.numeric(ccoru_nonsig_SNPs$contig_index)),]
# subset to get all contigs 1-271
ccoru_cand_nonsig <- ccoru_nonsig_SNPs_ordered0[as.numeric(ccoru_nonsig_SNPs_ordered0$contig_index) <= 271,]
dim(ccoru_cand_nonsig)
print(paste('There are ', nrow(ccoru_cand_nonsig), ' non-significant SNPs in ', length(unique(ccoru_cand_nonsig$Gene)), ' previously identified candidate genes.', sep=''))
print(paste('There are ', nrow(ccoru_nonsig_SNPs) - nrow(ccoru_cand_nonsig), ' non-significant SNPs in ',
            length(unique(ccoru_nonsig_SNPs$Contig)) - length(unique(ccoru_cand_nonsig$Contig)),
            ' contigs that have no previous association with elevation adaptation',
            sep=''))
# [1] "There are 40048 non-significant SNPs in 258 previously identified candidate genes."
# [1] "There are 171671 non-significant SNPs in 10997 contigs that have no previous association with elevation adaptation"

# -- make SNP files --------
# Format .sites and .rf files to run these sites through ANGSD for PCA and Fst calcs
# these do not differentiate between SNPs on candidate vs. control genes
# .sites file needs the combined_Contig# and site #
# .rf file needs only unique combined_Contig#:
ccoru_400mK2lamb1_05qval_ordered <- ccoru_400mK2lamb1_05qval[order(ccoru_400mK2lamb1_05qval$contig_num),]
ccoru_sites <- ccoru_400mK2lamb1_05qval_ordered[,c(9,10)]
ccoru_rf <- data.frame('contig'=paste(c(unique(ccoru_sites$contig_num)), ':', sep=''))

ccoru_sites_file <- file('ccoru_400mK2LFMMSNPs.sites', 'wb')
write.table(ccoru_sites, file=ccoru_sites_file, col.names=FALSE, row.names=FALSE, quote=FALSE, sep='\t')
close(ccoru_sites_file)

ccoru_rf_file <- file('ccoru_400mK2LFMMSNPs.rf', 'wb')
write.table(ccoru_rf, file=ccoru_rf_file, col.names=FALSE, row.names=FALSE, quote=FALSE)
close(ccoru_rf_file)

# use the following to calculate Fst between candidate vs. control genes
# sig SNPs on candidate genes
ccoru_previouslyknowngenes_ordered <- ccoru_previouslyknowngenes[order(ccoru_previouslyknowngenes$contig_num),]
ccoru_sites_candSNPs <- ccoru_previouslyknowngenes_ordered[,c(9,10)]
ccoru_rf_candSNPs <- data.frame('contig'=paste(c(unique(ccoru_sites_candSNPs$contig_num)), ':', sep=''))

ccoru_candsites_file <- file('ccoru_400mK2LFMM_candSNPs.sites', 'wb')
write.table(ccoru_sites_candSNPs, file=ccoru_candsites_file, col.names=FALSE, row.names=FALSE, quote=FALSE, sep='\t')
close(ccoru_candsites_file)

ccoru_candrf_file <- file('ccoru_400mK2LFMM_candSNPs.rf', 'wb')
write.table(ccoru_rf_candSNPs, file=ccoru_candrf_file, col.names=FALSE, row.names=FALSE, quote=FALSE)
close(ccoru_candrf_file)

# sig SNPs on control genes
control_sigSNPs_ccoru <- data.frame('loci'=setdiff(ccoru_400mK2lamb1_05qval$loci, ccoru_previouslyknowngenes$loci))
dim(control_sigSNPs_ccoru)
ccoru_loci_split2 <- strsplit(as.character(control_sigSNPs_ccoru$loci), '_')
control_sigSNPs_ccoru$contig_num <- sapply(ccoru_loci_split2, function(x){
  paste(x[[1]], x[[2]], sep='_')
})
control_sigSNPs_ccoru$SNP_site <- sapply(ccoru_loci_split2, function(x){
  x[[3]]
})

ccoru_controlSNPs_ordered <- control_sigSNPs_ccoru[order(control_sigSNPs_ccoru$contig_num),]
ccoru_sites_controlSNPs <- ccoru_controlSNPs_ordered[,c(2,3)]
ccoru_rf_controlSNPs <- data.frame('contig'=paste(c(unique(ccoru_sites_controlSNPs$contig_num)), ':', sep=''))

ccoru_controlsites_file <- file('ccoru_400mK2LFMM_controlSNPs.sites', 'wb')
write.table(ccoru_sites_controlSNPs, file=ccoru_controlsites_file, col.names=FALSE, row.names=FALSE, quote=FALSE, sep='\t')
close(ccoru_controlsites_file)

ccoru_controlrf_file <- file('ccoru_400mK2LFMM_controlSNPs.rf', 'wb')
write.table(ccoru_rf_controlSNPs, file=ccoru_controlrf_file, col.names=FALSE, row.names=FALSE, quote=FALSE)
close(ccoru_controlrf_file)

# -- Zchr ccoru --------
# use this txt file in a grepf search of Z chr contigs from Z chr folder
# want to make sure none of these are on the Z chr
write.table(ccoru_400mK2lamb1_05qval[,c(8:10)], 'C:/Users/mcwlim/Desktop/ccoru_Zchrcheck.txt', quote=FALSE, row.names=FALSE, sep='\t')

# Ok, now, how many LFMM SNPs are in previously ID'd candidate genes that are also on the Z chr?
Zchr_ccorulfmmSNPs <- read.table('ccoru_Zcontigs_removefromLFMM.txt', na.strings=c('', 'NA'), sep='\t')
colnames(Zchr_ccorulfmmSNPs) <- c('Gene', 'Contig', 'SNP_site')
length(intersect(Zchr_ccorulfmmSNPs$Gene, cand_targs$Gene))
# 0 Z chr contigs with LFMM SNPs are in previously ID'd candidate genes

# ------- 5c) Cviol 400m bins coeffs ----------
# read in the pearson corr coefficients for allele freqs per elev bin
# to merge with and compare to LFMM q-values
cviol_coef1 <- read.table('/2nd pass ANGSD - remove outliers/Allelefreqplots_byelevation/Cviol_Elevation_bins400m_snpbatch1_Pearsoncor_allelefreq-bin.txt', header=TRUE)
cviol_coef2 <- read.table('/2nd pass ANGSD - remove outliers/Allelefreqplots_byelevation/cviol_Elevation_bins400m_snpbatch2_Pearsoncor_allelefreq-bin.txt', header=TRUE)
cviol_coef3 <- read.table('/2nd pass ANGSD - remove outliers/Allelefreqplots_byelevation/cviol_Elevation_bins400m_snpbatch3_Pearsoncor_allelefreq-bin.txt', header=TRUE)
cviol_coef4 <- read.table('/2nd pass ANGSD - remove outliers/Allelefreqplots_byelevation/cviol_Elevation_bins400m_snpbatch4_Pearsoncor_allelefreq-bin.txt', header=TRUE)
cviol_coef5 <- read.table('/2nd pass ANGSD - remove outliers/Allelefreqplots_byelevation/cviol_Elevation_bins400m_snpbatch5_Pearsoncor_allelefreq-bin.txt', header=TRUE)
cviol_coef6 <- read.table('/2nd pass ANGSD - remove outliers/Allelefreqplots_byelevation/cviol_Elevation_bins400m_snpbatch6_Pearsoncor_allelefreq-bin.txt', header=TRUE)
cviol_coef7 <- read.table('/2nd pass ANGSD - remove outliers/Allelefreqplots_byelevation/cviol_Elevation_bins400m_snpbatch7_Pearsoncor_allelefreq-bin.txt', header=TRUE)
cviol_coef8 <- read.table('/2nd pass ANGSD - remove outliers/Allelefreqplots_byelevation/cviol_Elevation_bins400m_snpbatch8_Pearsoncor_allelefreq-bin.txt', header=TRUE)
cviol_coef9 <- read.table('/2nd pass ANGSD - remove outliers/Allelefreqplots_byelevation/cviol_Elevation_bins400m_snpbatch9_Pearsoncor_allelefreq-bin.txt', header=TRUE)

cviol_coef_all <- rbind(cviol_coef1, cviol_coef2, cviol_coef3, cviol_coef4, cviol_coef5, cviol_coef6,
                        cviol_coef7, cviol_coef8, cviol_coef9)
colnames(cviol_coef_all)[2] <- 'loci'

# calculate q-value
qobj_cviolPearson <- qvalue(cviol_coef_all$p.value, fdr.level=0.05, pi0.method = "bootstrap")
summary(qobj_cviolPearson)
cviol_coef_all <- cbind(cviol_coef_all, data.frame('Pearsoncor_qvalue'=qobj_cviolPearson$qvalues))
head(cviol_coef_all)
colnames(cviol_coef_all)[4] <- 'Pearsoncor_pvalue'
# note: no SNP survives after FDR correction for the Pearson's correlation results. The correlation is low-power, so let's focus on LFMM results

# ------- 5d) K=2 results, lambda = 1 ----
cviol_annot_snps_lambda1 <- read.table('LFMM_SNPs_annotated_CviolK2lambda1.txt', header=TRUE)
cviol_coef_annot_lambda1 <- merge(cviol_coef_all, cviol_annot_snps_lambda1[,c(4,7,9,10)], by='loci')
cviol_coef_sorted_lambda1 <- cviol_coef_annot_lambda1[order(cviol_coef_annot_lambda1$qvalue),]
head(cviol_coef_sorted_lambda1)
# check a few genes...
# cviol_coef_sorted_lambda1[cviol_coef_sorted_lambda1$significant == 1,] # look at top 10 significant lfmm snps that were also correlated with elev bin by allele freq
# length(cviol_coef_sorted_lambda1[cviol_coef_sorted_lambda1$significant == 1,]$Gene)
# head(cviol_coef_sorted_lambda1[cviol_coef_sorted_lambda1$Gene == 'EPAS1' & cviol_coef_sorted$significant == 1,], 10) #none of the top corr'd are also signif...
# cviol_coef_sorted_lambda1[cviol_coef_sorted_lambda1$Gene == 'SENP1' & cviol_coef_sorted_lambda1$significant ==1,] #no
# cviol_coef_sorted_lambda1[cviol_coef_sorted_lambda1$Gene == 'PPARA' & cviol_coef_sorted_lambda1$significant==1,] # no
# head(cviol_coef_sorted_lambda1[cviol_coef_sorted_lambda1$Gene == 'HIF1AN',], 10) #no
# RYR2 <- cviol_coef_sorted_lambda1[cviol_coef_sorted_lambda1$Gene == 'RYR2' & cviol_coef_sorted_lambda1$significant==1,]

cviol_400mK2lamb1_05qval <- cviol_coef_sorted_lambda1[cviol_coef_sorted_lambda1$qvalue < 0.05,]
dim(cviol_400mK2lamb1_05qval) # 6433 significant SNPs left after FDR

# -- compare genes from LFMM with sig and non-sig SNPs to targeted cand genes -----
# make column for contig number and SNP site >> use this to check if SNP results in dN or dS amino acid change
cviol_loci_split <- strsplit(as.character(cviol_400mK2lamb1_05qval$loci), '_')
cviol_400mK2lamb1_05qval$contig_num <- sapply(cviol_loci_split, function(x){
  paste(x[[1]], x[[2]], sep='_')
})
cviol_400mK2lamb1_05qval$SNP_site <- sapply(cviol_loci_split, function(x){
  x[[3]]
})

colnames(cand_targs) <- 'Gene'
#cviol_previouslyknowngenes <- merge(x=cviol_400mK2lamb1_05qval, y=cand_targs, by.y='Gene')
# this misses the genes that have odd names e.g., ACE_CHICK, NRG1_CHICK, but they are candidates!
# this is why there are some candidate gene SNPs in the control genes file
# however, since I know that combined_Contig1-271 are my candidate genes and 272+ are control genes,
# it is better to subset by contig number instead of trying to match the gene name (and the variations)
# create new col to order contig numbers
cviol_400mK2lamb1_05qval$contig_index <- gsub('combined_Contig','',cviol_400mK2lamb1_05qval$contig_num)
# put in order by contig number
cviol_400mK2lamb1_05qval_ordered0 <- cviol_400mK2lamb1_05qval[order(as.numeric(cviol_400mK2lamb1_05qval$contig_index)),]
# subset to get all contigs 1-271
cviol_previouslyknowngenes <- cviol_400mK2lamb1_05qval_ordered0[as.numeric(cviol_400mK2lamb1_05qval_ordered0$contig_index) <= 271,]

print(paste('There are ', nrow(cviol_previouslyknowngenes), ' significant SNPs in ', length(unique(cviol_previouslyknowngenes$Gene)), ' previously identified candidate genes.', sep=''))
print(paste('There are ', nrow(cviol_400mK2lamb1_05qval) - nrow(cviol_previouslyknowngenes), ' significant SNPs in ',
            length(unique(cviol_400mK2lamb1_05qval$contig_num)) - length(unique(cviol_previouslyknowngenes$contig_num)),
            ' contigs that are either new candidates for elevation adaptation, or are unrelated -- must check gene function for relevance',
            sep='')) #note that these are contigs not genes - because a number of the contigs have no gene ID
# [1] "There are 1263 significant SNPs in 188 previously identified candidate genes."
# [1] "There are 5170 significant SNPs in 3198 contigs that are either new candidates for elevation adaptation, or are unrelated -- must check gene function for relevance"

# check same but for the non-sig SNPs
cviol_nonsig_SNPs <- cviol_annot_snps_lambda1[cviol_annot_snps_lambda1$significant == 0,]
dim(cviol_nonsig_SNPs)
# create new col to order contig numbers
cviol_nonsig_SNPs$contig_index <- gsub('Contig','',cviol_nonsig_SNPs$Contig)
# put in order by contig number
cviol_nonsig_SNPs_ordered0 <- cviol_nonsig_SNPs[order(as.numeric(cviol_nonsig_SNPs$contig_index)),]
# subset to get all contigs 1-271
cviol_cand_nonsig <- cviol_nonsig_SNPs_ordered0[as.numeric(cviol_nonsig_SNPs_ordered0$contig_index) <= 271,]
dim(cviol_cand_nonsig)
print(paste('There are ', nrow(cviol_cand_nonsig), ' non-significant SNPs in ', length(unique(cviol_cand_nonsig$Gene)), ' previously identified candidate genes.', sep=''))
print(paste('There are ', nrow(cviol_nonsig_SNPs) - nrow(cviol_cand_nonsig), ' non-significant SNPs in ',
            length(unique(cviol_nonsig_SNPs$Contig)) - length(unique(cviol_cand_nonsig$Contig)),
            ' contigs that have no previous association with elevation adaptation',
            sep=''))
# [1] "There are 70008 non-significant SNPs in 10861 contigs that have no previous association with elevation adaptation"
# [1] "There are 14678 non-significant SNPs in 255 previously identified candidate genes."

# -- make SNP files -----
# Format .sites and .rf files to run these sites through ANGSD for PCA and Fst calcs
# .sites file needs the combined_Contig# and site #
# .rf file needs only unique combined_Contig#:
cviol_400mK2lamb1_05qval_ordered <- cviol_400mK2lamb1_05qval[order(cviol_400mK2lamb1_05qval$contig_num),]
cviol_sites <- cviol_400mK2lamb1_05qval_ordered[,c(9,10)]
cviol_rf <- data.frame('contig'=paste(c(unique(cviol_sites$contig_num)), ':', sep=''))

cviol_sites_file <- file('cviol_400mK2LFMMSNPs.sites', 'wb')
write.table(cviol_sites, file=cviol_sites_file, col.names=FALSE, row.names=FALSE, quote=FALSE, sep='\t')
close(cviol_sites_file)

cviol_rf_file <- file('cviol_400mK2LFMMSNPs.rf', 'wb')
write.table(cviol_rf, file=cviol_rf_file, col.names=FALSE, row.names=FALSE, quote=FALSE)
close(cviol_rf_file)

# use the following to calculate Fst between candidate vs. control genes
# sig SNPs on candidate genes (1148 SNPs, 165 genes)
cviol_previouslyknowngenes_ordered <- cviol_previouslyknowngenes[order(cviol_previouslyknowngenes$contig_num),]
cviol_sites_candSNPs <- cviol_previouslyknowngenes_ordered[,c(9,10)]
cviol_rf_candSNPs <- data.frame('contig'=paste(c(unique(cviol_sites_candSNPs$contig_num)), ':', sep=''))

cviol_candsites_file <- file('cviol_400mK2LFMM_candSNPs.sites', 'wb')
write.table(cviol_sites_candSNPs, file=cviol_candsites_file, col.names=FALSE, row.names=FALSE, quote=FALSE, sep='\t')
close(cviol_candsites_file)

cviol_candrf_file <- file('cviol_400mK2LFMM_candSNPs.rf', 'wb')
write.table(cviol_rf_candSNPs, file=cviol_candrf_file, col.names=FALSE, row.names=FALSE, quote=FALSE)
close(cviol_candrf_file)

# sig SNPs on control genes (5285 SNPs, 3221 contigs)
control_sigSNPs_cviol <- data.frame('loci'=setdiff(cviol_400mK2lamb1_05qval$loci, cviol_previouslyknowngenes$loci))
dim(control_sigSNPs_cviol)
cviol_loci_split2 <- strsplit(as.character(control_sigSNPs_cviol$loci), '_')
control_sigSNPs_cviol$contig_num <- sapply(cviol_loci_split2, function(x){
  paste(x[[1]], x[[2]], sep='_')
})
control_sigSNPs_cviol$SNP_site <- sapply(cviol_loci_split2, function(x){
  x[[3]]
})

cviol_controlSNPs_ordered <- control_sigSNPs_cviol[order(control_sigSNPs_cviol$contig_num),]
cviol_sites_controlSNPs <- cviol_controlSNPs_ordered[,c(2,3)]
cviol_rf_controlSNPs <- data.frame('contig'=paste(c(unique(cviol_sites_controlSNPs$contig_num)), ':', sep=''))

cviol_controlsites_file <- file('cviol_400mK2LFMM_controlSNPs.sites', 'wb')
write.table(cviol_sites_controlSNPs, file=cviol_controlsites_file, col.names=FALSE, row.names=FALSE, quote=FALSE, sep='\t')
close(cviol_controlsites_file)

cviol_controlrf_file <- file('cviol_400mK2LFMM_controlSNPs.rf', 'wb')
write.table(cviol_rf_controlSNPs, file=cviol_controlrf_file, col.names=FALSE, row.names=FALSE, quote=FALSE)
close(cviol_controlrf_file)

# -- Zchr cviol ----
# use this txt file in a grepf search of Z chr contigs from Z chr folder
# want to make sure none of these are on the Z chr
write.table(cviol_400mK2lamb1_05qval[,c(8:10)], 'C:/Users/mcwlim/Desktop/cviol_Zchrcheck.txt', quote=FALSE, row.names=FALSE, sep='\t')

# Ok, now, how many LFMM SNPs are in previously ID'd candidate genes that are also on the Z chr?
Zchr_cviollfmmSNPs <- read.table('cviol_Zcontigs_removefromLFMM.txt')
colnames(Zchr_cviollfmmSNPs) <- c('Gene', 'Contig', 'SNP_site')
length(intersect(Zchr_cviollfmmSNPs$Gene, cand_targs$Gene))
# 0 Z chr contigs with LFMM SNPs are in previously ID'd candidate genes

# ------- 5e) Any overlapping genes for Ccoru vs. Cviol? ---------
ccoruVScviol_previous <- intersect(ccoru_previouslyknowngenes$Gene, cviol_previouslyknowngenes$Gene)
length(ccoruVScviol_previous)
# 124

ccoruVScviol_all <- intersect(ccoru_400mK2lamb1_05qval$Gene, cviol_400mK2lamb1_05qval$Gene)
length(ccoruVScviol_all)
# 828

# on novel candidates
ccoruVScviol_novel <- setdiff(ccoruVScviol_all, ccoruVScviol_previous)
length(ccoruVScviol_novel)
write.csv(ccoru_previouslyknowngenes, 'Ccoru_SNPsinKnownCandGenes.csv', row.names=FALSE, quote=FALSE)
write.csv(cviol_previouslyknowngenes, 'Cviol_SNPsinKnownCandGenes.csv', row.names=FALSE, quote=FALSE)

write.table(ccoruVScviol_previous, 'Ccoru-Cviol_previousIDcandgenes.txt', quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(ccoruVScviol_novel, 'Ccoru-Cviol_novelcandgenes.txt', quote=FALSE, row.names=FALSE, col.names=FALSE)

# ------- 5f) plot contingency table as bar plot for cline vs. no cline/cand vs. control/violifer & coruscans ---------
# (following code from Liliana from Figure 1: https://www.nature.com/articles/s41559-018-0727-8)

countsToCases <- function(x, countcol = "SNP_freq") {
  # Get the row indices to pull from x
  idx <- rep.int(seq_len(nrow(x)), x[[countcol]])
  # Drop count column#
  x[[countcol]] <- NULL
  # Get the rows from x
  x[idx, ]
}
# LFMM - SNP is clinal or not clinal
# species - for 2 highland species
# SNP_type - SNPs are on candidate or control genes
# SNP_freq - number of SNPs that are clinal-candidate, clinal-control, not clinal-candidate, not clinal-control
# gene_freq - number of genes with SNPs that are clinal-candidate, clinal-control, not clinal-candidate, not clinal-control
dcounts <- data.frame(
  LFMM=c('cline', 'no cline', 'cline', 'no cline', 'cline', 'no cline', 'cline', 'no cline'),
  species=c(rep('Coe. violifer',4), rep('Col. coruscans',4)),
  SNP_type=c(rep('candidate',2), rep('control',2), rep('candidate',2), rep('control',2)),
  SNP_freq=c(1263, 14678, 5170, 70008, 757, 40048, 3205, 171671),
  gene_freq=c(188, 225, 3198, 10861, 174, 258, 2176, 10997))

dcounts_v2 <- data.frame(
  LFMM=c('yes_cline', 'no_cline', 'yes_cline', 'no_cline', 'yes_cline', 'no_cline', 'yes_cline', 'no_cline'),
  species=c(rep('Coe. violifer',4), rep('Col. coruscans',4)),
  SNP_type=c(rep('candidate',2), rep('control',2), rep('candidate',2), rep('control',2)),
  SNP_freq=c(1263, 14678, 5170, 70008, 757, 40048, 3205, 171671),
  gene_freq=c(188, 225, 3198, 10861, 174, 258, 2176, 10997))

dcases<-countsToCases(dcounts)

# -- by SNP ------
# x-axis = cline vs. no cline and colors = candidate vs. control
dcounts$factorC <-with(dcounts, interaction(LFMM, species))
dc3<-as.data.frame(dcounts %>% group_by(factorC) %>% mutate(fraction=SNP_freq/sum(SNP_freq)))
dc3$SNP_type<- factor(dc3$SNP_type, levels = c("candidate", "control"))
ggplot() + geom_bar(aes(y = fraction, x = factorC, fill = SNP_type),
                    data = dc3, stat="identity")+
  labs(x="Cline vs. no cline", y="proportion") +
  guides(fill=guide_legend(title = "SNP type"))+
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size=15))

dc3$species<- factor(dc3$species, levels = c('Coe. violifer', 'Col. coruscans'))
ggplot() + geom_bar(aes(y = fraction, x = LFMM, fill = SNP_type),
                    data = dc3, stat="identity")+
  labs(x="Cline vs. no cline", y="proportion") +
  guides(fill=guide_legend(title = "SNP type"))+
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size=15))+
  facet_grid(~ species, margins=F)
ggsave('C:/Users/mcwlim/Dropbox/Marisacompfiles/Seq cap files/Pop_GEA_chapter/Pop_GEA_figs/cand-control_cline-not/SNPfacetCviolCcoru_cline-no_cand-control.jpg',
       h=5, w=6, dpi=600)

# x-axis = candidate vs. control and colors = cline vs. no cline
dcounts$factorC2 <-with(dcounts, interaction(SNP_type, species))
dc4<-as.data.frame(dcounts %>% group_by(factorC2) %>% mutate(fraction=SNP_freq/sum(SNP_freq)))
dc4$LFMM<- factor(dc4$LFMM, levels = c("cline", "no cline"))
ggplot() + geom_bar(aes(y = fraction, x = factorC2, fill = LFMM),
                    data = dc4, stat="identity")+
  labs(x="Candidate vs. control", y="proportion") +
  guides(fill=guide_legend(title = "Cline type"))+
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size=15))

dc4$species<- factor(dc4$species, levels = c('Coe. violifer', 'Col. coruscans'))
ggplot() + geom_bar(aes(y = fraction, x = SNP_type, fill = LFMM),
                    data = dc4, stat="identity")+
  labs(x="Candidate vs. control", y="proportion") +
  guides(fill=guide_legend(title = "Cline type"))+
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size=15))+
  facet_grid(~ species, margins=F)
ggsave('C:/Users/mcwlim/Dropbox/Marisacompfiles/Seq cap files/Pop_GEA_chapter/Pop_GEA_figs/cand-control_cline-not/SNPfacetCviolCcoru_cand-control_cline-no.jpg',
       h=5, w=6, dpi=600)

# -- by gene ----
#(by gene instead of by SNP shows the diffs between cand and control better)
# x-axis = cline vs. no cline and colors = candidate vs. control
dcounts$factorC <-with(dcounts, interaction(LFMM, species))
dc5<-as.data.frame(dcounts %>% group_by(factorC) %>% mutate(fraction=gene_freq/sum(gene_freq)))
dc5$species<- factor(dc5$species, levels = c('Coe. violifer', 'Col. coruscans'))
ggplot() + geom_bar(aes(y = fraction, x = LFMM, fill = SNP_type),
                    data = dc5, stat="identity")+
  labs(x="Cline vs. no cline", y="proportion") +
  guides(fill=guide_legend(title = "gene type"))+
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size=15))+
  facet_grid(~ species, margins=F)
ggsave('C:/Users/mcwlim/Dropbox/Marisacompfiles/Seq cap files/Pop_GEA_chapter/Pop_GEA_figs/cand-control_cline-not/genefacetCviolCcoru_cline-no_cand-control.jpg',
       h=5, w=6, dpi=600)

# USE THIS ONE!!
# x-axis = candidate vs. control and colors = cline vs. no cline
dcounts$factorC2 <-with(dcounts, interaction(SNP_type, species))
dc6<-as.data.frame(dcounts %>% group_by(factorC2) %>% mutate(fraction=gene_freq/sum(gene_freq)))
dc6$species<- factor(dc6$species, levels = c('Coe. violifer', 'Col. coruscans'))
ggplot() + geom_bar(aes(y = fraction, x = SNP_type, fill = LFMM),
                    data = dc6, stat="identity")+
  labs(x="Candidate vs. control", y="proportion") +
  guides(fill=guide_legend(title = "Cline type"))+
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size=15))+
  facet_grid(~ species, margins=F)
ggsave('C:/Users/mcwlim/Dropbox/Marisacompfiles/Seq cap files/Pop_GEA_chapter/Pop_GEA_figs/cand-control_cline-not/genefacetCviolCcoru_cand-control_cline-no.jpg',
       h=5, w=6, dpi=600)

# -- models by gene ------

# poisson distributed errors
mod1_poi <- glm(gene_freq ~ species*LFMM*SNP_type, data=dcounts, family=poisson())
mod2_poi <- update(mod1_poi, ~.-species:LFMM:SNP_type) #gene_freq ~ species + LFMM + SNP_type + species:LFMM + species:SNP_type + LFMM:SNP_type
mod3_poi <- update(mod2_poi, ~.-LFMM:SNP_type) #gene_freq ~ species + LFMM + SNP_type + species:LFMM + species:SNP_type
mod4_poi <- update(mod3_poi, ~. -species:LFMM -species:SNP_type -LFMM:SNP_type) #gene_freq ~ species + LFMM + SNP_type

anova(mod1_poi, mod2_poi, mod3_poi, mod4_poi, test='Chisq')

# slightly overdispersed poisson
# same coefficient estimates as poisson but adjusted standard errors
mod1_qpoi <- glm(gene_freq ~ species*LFMM*SNP_type, data=dcounts, family=quasipoisson())
mod2_qpoi <- update(mod1_qpoi, ~.-species:LFMM:SNP_type) #gene_freq ~ species + LFMM + SNP_type + species:LFMM + species:SNP_type + LFMM:SNP_type
mod3_qpoi <- update(mod2_qpoi, ~.-LFMM:SNP_type) #gene_freq ~ species + LFMM + SNP_type + species:LFMM + species:SNP_type
mod4_qpoi <- update(mod3_qpoi, ~. -species:LFMM -species:SNP_type -LFMM:SNP_type) #gene_freq ~ species + LFMM + SNP_type

anova(mod1_qpoi, mod2_qpoi, mod3_qpoi, mod4_qpoi)

## quasi-poisson - change order so in reference to yes_cline
mod1_qpoi2 <- glm(gene_freq ~ species*LFMM*SNP_type, data=dcounts_v2, family=quasipoisson())
mod2_qpoi2 <- update(mod1_qpoi2, ~.-species:LFMM:SNP_type) #gene_freq ~ species + LFMM + SNP_type + species:LFMM + species:SNP_type + LFMM:SNP_type
mod3_qpoi2 <- update(mod2_qpoi2, ~.-LFMM:SNP_type) #gene_freq ~ species + LFMM + SNP_type + species:LFMM + species:SNP_type
mod4_qpoi2 <- update(mod3_qpoi2, ~. -species:LFMM -species:SNP_type -LFMM:SNP_type) #gene_freq ~ species + LFMM + SNP_type

anova(mod1_qpoi2, mod2_qpoi2, mod3_qpoi2, mod4_qpoi2)
summary(mod2_qpoi2)



# negative binomial - i think this is not running because theta gets too high and becomes NaN...but in previous tests, the thetas were too high so Liliana didn't think neg binomial was good fit to the data
mod1_nb <- glm.nb(gene_freq ~ species*LFMM*SNP_type, data=dcounts)
mod2_nb <- update(mod1_nb, ~.-species:LFMM:SNP_type) #gene_freq ~ species + LFMM + SNP_type + species:LFMM + species:SNP_type + LFMM:SNP_type
mod3_nb <- update(mod2_nb, ~.-LFMM:SNP_type) #gene_freq ~ species + LFMM + SNP_type + species:LFMM + species:SNP_type
mod4_nb <- update(mod3_nb, ~. -species:LFMM -species:SNP_type -LFMM:SNP_type) #gene_freq ~ species + LFMM + SNP_type

sink('C:/Users/mcwlim/Dropbox/Marisacompfiles/Seq cap files/Pop_GEA_chapter/Pop_GEA_figs/cand-control_cline-not/mod_poissondistr_summaryresults.txt')
print(summary(mod1_poi))
print(summary(mod2_poi))
print(summary(mod3_poi))
print(summary(mod4_poi))
sink()

sink('C:/Users/mcwlim/Dropbox/Marisacompfiles/Seq cap files/Pop_GEA_chapter/Pop_GEA_figs/cand-control_cline-not/mod_quasipoissondistr_summaryresults.txt')
print(summary(mod1_qpoi))
print(summary(mod2_qpoi))
print(summary(mod3_qpoi))
print(summary(mod4_qpoi))
sink()

# sink('C:/Users/mcwlim/Dropbox/Marisacompfiles/Seq cap files/Pop_GEA_chapter/Pop_GEA_figs/cand-control_cline-not/mod_negbinom_summaryresults.txt')
# print(summary(mod1_nb))
# print(summary(mod2_nb))
# print(summary(mod3_nb))
# print(summary(mod4_nb))
# sink()

# prelim interpretation
# if species is not sig/only sorta sig: whether SNP is clinal does not really depend on which species you're looking at
# if SNP_type is sig: whether SNP is clinal DOES depend on whether the SNP is on a candidate vs. control gene
# if LFMM is sig: well, this ought to always be sig because whether a SNP is clinal SHOULD depend on whether there is cline or not...




# -- 6. Plot genotypes (0, 1, 2) for specific SNPs for each individual -----

# Read in geno files
cviol_geno <- read.table('../../Allelefreqplots_byelevation/cvioln59_genotypes.txt')
ccoru_geno <- read.table('../../Allelefreqplots_byelevation/ccorun97_genotypes.txt')
# Read in loci
cviol_loci <- read.table('../../Pcadapt_outflank/cviol59_outflank.loci')
ccoru_loci <- read.table('../../Pcadapt_outflank/ccoru97_outflank.loci')
# read in elev bin, lat/longs (also has southern winter/summer for ccoru and geo subpop info)
annot_cvioln59 <- read.table('../../ngsTools PCA - v2/cviol_n59.clst', sep="\t", header=T) ## need to add elev bin to this!!
colnames(annot_cvioln59)[9] <- 'elevbin'
annot_ccorun97 <- read.table('../../ngsTools PCA - v2/ccoru_n97.clst', sep='\t', header=T)
colnames(annot_ccorun97)[11] <- 'elevbin'

# ------- Map genotypes ------

genotype_map <- function(thegeno, theloci, theannotfile, theSNP, thespecies_gene){

  # 1. Format data
  # make loci rowname, then transpose
  rownames(thegeno) <- theloci$V1
  thegeno2 <- t(thegeno)

  # combine annot file columns with loci and geno
  thedat <- cbind(theannotfile, thegeno2)

  #find the column number to subset data by and subset by SNP
  annot_cols <- ncol(theannotfile) # different for cviol vs. ccoru, because ccoru also has the season info
  SNPdat <- thedat[, c(1:annot_cols, which(colnames(thedat)==theSNP))]
  # remove rows with 9s (missing genotype data)
  SNPdat <- SNPdat[SNPdat[annot_cols + 1] != 9, ]

  colnames(SNPdat)[annot_cols + 1] <- 'Genotype'

  # 2. Plot maps
  # purpose is to look at spatial distribution of genotypes
  soam2 <- get_map(location =c(lon= -75, lat=-8), zoom=5, maptype="hybrid", color="bw")
  # subset by elev bin
  ggmap(soam2, extent="normal", maprange=T) +
    geom_jitter(data=SNPdat,
                aes(x=Long, y=Lat, col=as.character(Genotype)),
                size=4, alpha=0.7, height=1, width=0.5) +
    facet_wrap(~elevbin, ncol=3) +
    xlab('Longitude (degrees)') + ylab('Latitude (degrees)') +
    ggtitle(paste(thespecies_gene, ' ', theSNP, sep='')) +
    scale_colour_discrete(name='Genotype')
  ggsave(paste('../SNPoutlierGenotype_byelevbin_map/', thespecies_gene, '_', theSNP, 'elevbingenoplot.jpg', sep=''), height=10, width=10, units='in', dpi=600)

}

genotype_map(ccoru_geno, ccoru_loci, annot_ccorun97[ ,c(1,4,6:11)], 'combined_Contig174_6567', 'Ccoru_SENP1')
genotype_map(ccoru_geno, ccoru_loci, annot_ccorun97[ ,c(1,4,6:11)], 'combined_Contig174_3984', 'Ccoru_SENP1')

genotype_map(cviol_geno, cviol_loci, annot_cvioln59[ ,c(1,4,5:8)], '....', 'Cviol') # decide on SNPs to test


# ------- plot genotypes by elevation boxplots & Wilcoxon rank sum test ----------
# loop through all SNPs and get wilcoxan test stats - then decide which to plot
# First, loop through SNPs in [species]_previouslyknowngenes
# Second, do wilcoxan test
genotype_elev_wilcox_ALLSNPS <- function(thegeno, theloci, looploci, theannotfile, thespecies){

  # 1. Format data
  # make loci rowname, then transpose
  rownames(thegeno) <- theloci$V1
  thegeno2 <- t(thegeno)

  # combine annot file columns with loci and geno
  thedat <- cbind(theannotfile, thegeno2)

  df <- data.frame('Gene'='mygene', 'SNP'='SNP', 'Elev-01_p.adj'=1, 'Elev-02_p.adj'=1,
                   'Elev-12_p.adj'=1, 'Elev-01_p.signif'='ns',
                   'Elev-02_p.signif'='ns', 'Elev-12_p.signif'='ns', 'Pearson_coef_littlea'=1)
  # loop through each SNP
  #find the column number to subset data by and subset by SNP
  annot_cols <- ncol(theannotfile) # different for cviol vs. ccoru, because ccoru also has the season info
  for(i in 1:nrow(looploci)){
    cat(i)
    SNPcol <- which(colnames(thedat) == looploci$loci[i])
    SNPdat <- thedat[, c(1:annot_cols, SNPcol)]

    theSNP <- looploci$loci[i]
    mygene <- looploci$Gene[i]
    myPearsoncor <- looploci$Pearson_coef_littlea[i]

    # remove rows with 9s (missing genotype data)
    SNPdat <- SNPdat[SNPdat[annot_cols + 1] != 9, ]
    colnames(SNPdat)[annot_cols + 1] <- 'Genotype'

    # a non-parametric comparison of means (per genotype category)
    # use Wilcoxon rank sum test (Mann Whitney U )
    wilcox_results <- compare_means(Elevation ~ Genotype, data=SNPdat,
                                    method = 'wilcox.test') # code below plots this type of test
    wilcox_results2 <- data.frame('Gene'=mygene, 'SNP'=theSNP, 'Elev-01_p.adj'=wilcox_results$p.adj[1], 'Elev-02_p.adj'=wilcox_results$p.adj[2],
                                  'Elev-12_p.adj'=wilcox_results$p.adj[3], 'Elev-01_p.signif'=wilcox_results$p.signif[1],
                                  'Elev-02_p.signif'=wilcox_results$p.signif[2], 'Elev-12_p.signif'=wilcox_results$p.signif[3],
                                  'Pearson_coef_littlea'=myPearsoncor)
    df <- rbind(df, wilcox_results2)
  }
  df1 <- df[-1,] # if there are NA's, it's likely because the SNP doesn't have all 3 genotypes, e.g. only 0 and 1, so can only do 0-1 comparison
  write.csv(df1, paste('Wilcoxon_test_', thespecies, '.csv', sep=''), row.names=FALSE, quote=FALSE)
}
genotype_elev_wilcox_ALLSNPS(ccoru_geno, ccoru_loci, ccoru_previouslyknowngenes, annot_ccorun97, 'Ccoru')
genotype_elev_wilcox_ALLSNPS(cviol_geno, cviol_loci, cviol_previouslyknowngenes, annot_cvioln59, 'Cviol')
# Third, now read in these files and subset by the rows with * or **, indicating signif difference between genotypes and elevation
ccoru_wilcoxresults <- read.csv('Wilcoxon_test_Ccoru.csv', header=TRUE)
cviol_wilcoxresults <- read.csv('Wilcoxon_test_Cviol.csv', header=TRUE)

dim(ccoru_wilcoxresults) #757
dim(cviol_wilcoxresults) #1263

# let's just check the ones with * or ** signif
ccoru_wilcoxresults_subset <- ccoru_wilcoxresults[unique(c(as.vector(grep('\\*', ccoru_wilcoxresults[,5])),
                                                               as.vector(grep('\\*', ccoru_wilcoxresults[,6])),
                                                               as.vector(grep('\\*', ccoru_wilcoxresults[,7])))),]

cviol_wilcoxresults_subset <- cviol_wilcoxresults[unique(c(as.vector(grep('\\*', cviol_wilcoxresults[,5])),
                                                    as.vector(grep('\\*', cviol_wilcoxresults[,6])),
                                                    as.vector(grep('\\*', cviol_wilcoxresults[,7])))),]
dim(ccoru_wilcoxresults_subset) #369 SNPs
dim(cviol_wilcoxresults_subset) #794 SNPs
length(unique(ccoru_wilcoxresults_subset$Gene)) #120
length(unique(cviol_wilcoxresults_subset$Gene)) #149

write.csv(ccoru_wilcoxresults_subset, '../../../../../Pop_GEA_chapter/TableS10_Ccoru_wilcoxTOPresults.csv', row.names=FALSE, quote=FALSE)
write.csv(cviol_wilcoxresults_subset, '../../../../../Pop_GEA_chapter/TableS9_Cviol_wilcoxTOPresults.csv', row.names=FALSE, quote=FALSE)

# Fourth, now loop through these SNPs that have * or ** and plot the boxplots
# genotype_elev_wilcoxboxplot <- function(thegeno, theloci, looploci, theannotfile, thespecies){
#
#   # 1. Format data
#   # make loci rowname, then transpose
#   rownames(thegeno) <- theloci$V1
#   thegeno2 <- t(thegeno)
#
#   # combine annot file columns with loci and geno
#   thedat <- cbind(theannotfile, thegeno2)
#
#   # loop through each SNP
#   #find the column number to subset data by and subset by SNP
#   annot_cols <- ncol(theannotfile) # different for cviol vs. ccoru, because ccoru also has the season info
#   for(i in 1:nrow(looploci)){
#     cat(i)
#     SNPcol <- which(colnames(thedat) == looploci$SNP[i])
#     SNPdat <- thedat[, c(1:annot_cols, SNPcol)]
#
#     theSNP <- looploci$SNP[i]
#     mygene <- looploci$Gene[i]
#
#     # remove rows with 9s (missing genotype data)
#     SNPdat <- SNPdat[SNPdat[annot_cols + 1] != 9, ]
#     colnames(SNPdat)[annot_cols + 1] <- 'Genotype'
#
#     # ggpubr package currently (as of May 2018) doesn't have a way to plot the adjusted
#     # p-values [https://github.com/kassambara/ggpubr/issues/65]
#     # so, instead, i'll plot the p.signif on boxplots - check that these reflect the
#     # adjusted p-values, rather than the raw p-values
#
#     mycomparisons <- list(c('0', '1'), c('1', '2'), c('0', '2'))
#     ggboxplot(data=SNPdat, x='Genotype', y='Elevation',
#               color='Genotype', palette='jco', show.legend=FALSE) +
#       stat_compare_means(comparisons=mycomparisons, label='p.signif') +
#       ggtitle(paste(mygene, ' ', theSNP, ' geno x elev', sep=''))
#     ggsave(paste('../SNPoutlierGenotype_byelevbin_map/', thespecies, '_', mygene, '_', theSNP, '_genoXelev_wilcox.jpg', sep=''),
#            height=5, width=6, units='in', dpi=600)
#   }
# }
# genotype_elev_wilcoxboxplot(ccoru_geno, ccoru_loci, ccoru_wilcoxresults_subset, annot_ccorun97, 'Ccoru')
# genotype_elev_wilcoxboxplot(cviol_geno, cviol_loci, cviol_wilcoxresults_subset, annot_cvioln59, 'Cviol')

# this is only for specifc SNPs, also makes the boxplots
genotype_elev_wilcox <- function(thegeno, theloci, theannotfile, theSNP, thespecies_gene){

  # 1. Format data
  # make loci rowname, then transpose
  rownames(thegeno) <- theloci$V1
  thegeno2 <- t(thegeno)

  # combine annot file columns with loci and geno
  thedat <- cbind(theannotfile, thegeno2)

  #find the column number to subset data by and subset by SNP
  annot_cols <- ncol(theannotfile) # different for cviol vs. ccoru, because ccoru also has the season info
  SNPdat <- thedat[, c(1:annot_cols, which(colnames(thedat)==theSNP))]
  # remove rows with 9s (missing genotype data)
  SNPdat <- SNPdat[SNPdat[annot_cols + 1] != 9, ]

  colnames(SNPdat)[annot_cols + 1] <- 'Genotype'

  # 2. Plots

  # statistically test for differences between the boxplots
  # yeah, these don't look parametric,
  # so we'll go with a non-parametric comparison of means (per genotype category)
  # ggplot(data=SENP1_combined_Contig174_6567) +
  #   geom_histogram(aes(x=Elevation), stat='bin') + facet_wrap(~Elevation_bin_400m, ncol=3)

  # use Wilcoxon rank sum test (Mann Whitney U )
  # wilcox_results <- compare_means(Elevation ~ Genotype, data=SNPdat,
                                                      # method = 'wilcox.test') # code below plots this type of test
  #write.table(wilcox_results, paste(thespecies_gene, '_', theSNP, 'wilcoxtest.txt', sep=''))

  # ggpubr package currently (as of May 2018) doesn't have a way to plot the adjusted
  # p-values [https://github.com/kassambara/ggpubr/issues/65]
  # so, instead, i'll plot the p.signif on boxplots - check that these reflect the
  # adjusted p-values, rather than the raw p-values

  mycomparisons <- list(c('0', '1'), c('1', '2'), c('0', '2'))
  ggboxplot(data=SNPdat, x='Genotype', y='Elevation',
            color='Genotype', palette='jco', show.legend=FALSE) +
    stat_compare_means(comparisons=mycomparisons, label='p.signif') +
    ggtitle(paste(thespecies_gene, ' ', theSNP, ' geno x elev', sep=''))
  ggsave(paste('../SNPoutlierGenotype_byelevbin_map/', thespecies_gene, '_', theSNP, '_genoXelev_wilcox.jpg', sep=''),
         height=5, width=6, units='in', dpi=600)
}
# genotype_elev_wilcox(ccoru_geno, ccoru_loci, annot_ccorun97[ ,c(1,4,6:11)], 'combined_Contig174_6567', 'Ccoru_SENP1')
# genotype_elev_wilcox(ccoru_geno, ccoru_loci, annot_ccorun97[ ,c(1,4,6:11)], 'combined_Contig174_3984', 'Ccoru_SENP1')
# genotype_elev_wilcox(ccoru_geno, ccoru_loci, annot_ccorun97[,c(1,4,6:11)], 'combined_Contig198_4637', 'Ccoru_ALOX5')
# genotype_elev_wilcox(ccoru_geno, ccoru_loci, annot_ccorun97[,c(1,4,6:11)], 'combined_Contig198_8465', 'Ccoru_ALOX5')
# genotype_elev_wilcox(cviol_geno, cviol_loci, annot_cvioln59[,c(1,4:9)], 'combined_Contig198_364', 'Cviol_ALOX5')


# -- 7. Plot allele frequencies by elevation bin ------
# set up data: plot same SNPs as for Wilcoxon box plots
# SNPs are in ccoru_wilcoxresults_subset or cviol_wilcoxresults_subset
# allele frequencies
# need to load the allele freqs
ccoru_af1 <- read.table('/2nd pass ANGSD - remove outliers/Allelefreqplots_byelevation/Ccoru_Elevation_bin400m_snpbatch1_allelefreqsbybins.txt', header=TRUE)
ccoru_af2 <- read.table('/2nd pass ANGSD - remove outliers/Allelefreqplots_byelevation/Ccoru_Elevation_bin400m_snpbatch2_allelefreqsbybins.txt', header=TRUE)
ccoru_af3 <- read.table('/2nd pass ANGSD - remove outliers/Allelefreqplots_byelevation/Ccoru_Elevation_bin400m_snpbatch3_allelefreqsbybins.txt', header=TRUE)
ccoru_af4 <- read.table('/2nd pass ANGSD - remove outliers/Allelefreqplots_byelevation/Ccoru_Elevation_bin400m_snpbatch4_allelefreqsbybins.txt', header=TRUE)
ccoru_af5 <- read.table('/2nd pass ANGSD - remove outliers/Allelefreqplots_byelevation/Ccoru_Elevation_bin400m_snpbatch5_allelefreqsbybins.txt', header=TRUE)
ccoru_af6 <- read.table('/2nd pass ANGSD - remove outliers/Allelefreqplots_byelevation/Ccoru_Elevation_bin400m_snpbatch6_allelefreqsbybins.txt', header=TRUE)
ccoru_af7 <- read.table('/2nd pass ANGSD - remove outliers/Allelefreqplots_byelevation/Ccoru_Elevation_bin400m_snpbatch7_allelefreqsbybins.txt', header=TRUE)
ccoru_af8 <- read.table('/2nd pass ANGSD - remove outliers/Allelefreqplots_byelevation/Ccoru_Elevation_bin400m_snpbatch8_allelefreqsbybins.txt', header=TRUE)
ccoru_af9 <- read.table('/2nd pass ANGSD - remove outliers/Allelefreqplots_byelevation/Ccoru_Elevation_bin400m_snpbatch9_allelefreqsbybins.txt', header=TRUE)
ccoru_af10 <- read.table('/2nd pass ANGSD - remove outliers/Allelefreqplots_byelevation/Ccoru_Elevation_bin400m_snpbatch10_allelefreqsbybins.txt', header=TRUE)
ccoru_af11 <- read.table('/2nd pass ANGSD - remove outliers/Allelefreqplots_byelevation/Ccoru_Elevation_bin400m_snpbatch11_allelefreqsbybins.txt', header=TRUE)

ccoru_af_all <- rbind(ccoru_af1, ccoru_af2, ccoru_af3, ccoru_af4, ccoru_af5, ccoru_af6, ccoru_af7, ccoru_af8,
                      ccoru_af9, ccoru_af10, ccoru_af11)
colnames(ccoru_af_all)[3] <- 'SNP'
ccoru_af_dat <- merge(ccoru_wilcoxresults_subset[,c(1,2,9)], ccoru_af_all, by='SNP')

# need to load the allele freqs
cviol_af1 <- read.table('/2nd pass ANGSD - remove outliers/Allelefreqplots_byelevation/Cviol_Elevation_bins400m_snpbatch1_allelefreqsbybins.txt', header=TRUE)
cviol_af2 <- read.table('/2nd pass ANGSD - remove outliers/Allelefreqplots_byelevation/Cviol_Elevation_bins400m_snpbatch2_allelefreqsbybins.txt', header=TRUE)
cviol_af3 <- read.table('/2nd pass ANGSD - remove outliers/Allelefreqplots_byelevation/Cviol_Elevation_bins400m_snpbatch3_allelefreqsbybins.txt', header=TRUE)
cviol_af4 <- read.table('/2nd pass ANGSD - remove outliers/Allelefreqplots_byelevation/Cviol_Elevation_bins400m_snpbatch4_allelefreqsbybins.txt', header=TRUE)
cviol_af5 <- read.table('/2nd pass ANGSD - remove outliers/Allelefreqplots_byelevation/Cviol_Elevation_bins400m_snpbatch5_allelefreqsbybins.txt', header=TRUE)
cviol_af6 <- read.table('/2nd pass ANGSD - remove outliers/Allelefreqplots_byelevation/Cviol_Elevation_bins400m_snpbatch6_allelefreqsbybins.txt', header=TRUE)
cviol_af7 <- read.table('/2nd pass ANGSD - remove outliers/Allelefreqplots_byelevation/Cviol_Elevation_bins400m_snpbatch7_allelefreqsbybins.txt', header=TRUE)
cviol_af8 <- read.table('/2nd pass ANGSD - remove outliers/Allelefreqplots_byelevation/Cviol_Elevation_bins400m_snpbatch8_allelefreqsbybins.txt', header=TRUE)
cviol_af9 <- read.table('/2nd pass ANGSD - remove outliers/Allelefreqplots_byelevation/Cviol_Elevation_bins400m_snpbatch9_allelefreqsbybins.txt', header=TRUE)

cviol_af_all <- rbind(cviol_af1, cviol_af2, cviol_af3, cviol_af4, cviol_af5, cviol_af6, cviol_af7, cviol_af8,
                      cviol_af9)
colnames(cviol_af_all)[3] <- 'SNP'
cviol_af_dat <- merge(cviol_wilcoxresults_subset[,c(1,2,9)], cviol_af_all, by='SNP')

# Version 0: try to subset by positive, flatish, and negative trends
myfacetIDs <- c('r > 0.80', '0.50 < r < 0.80', '-0.50 > r > -0.80', 'r < -0.80')
myplotfacets_ccoru <- rep(myfacetIDs[1], length(ccoru_af_dat[,1]))
myplotfacets_ccoru[ccoru_af_dat$Pearson_coef_littlea >= 0.80] <- myfacetIDs[1]
myplotfacets_ccoru[ccoru_af_dat$Pearson_coef_littlea < 0.80 & ccoru_af_dat$Pearson_coef_littlea > 0.50] <- myfacetIDs[2]
myplotfacets_ccoru[ccoru_af_dat$Pearson_coef_littlea < -0.50 & ccoru_af_dat$Pearson_coef_littlea > -0.80] <- myfacetIDs[3]
myplotfacets_ccoru[ccoru_af_dat$Pearson_coef_littlea <= -0.80] <- myfacetIDs[4]
ccoru_af_dat$myplotfacets <- myplotfacets_ccoru

myplotfacets_cviol <- rep(myfacetIDs[1], length(cviol_af_dat[,1]))
myplotfacets_cviol[cviol_af_dat$Pearson_coef_littlea >= 0.80] <- myfacetIDs[1]
myplotfacets_cviol[cviol_af_dat$Pearson_coef_littlea < 0.80 & cviol_af_dat$Pearson_coef_littlea > 0.50] <- myfacetIDs[2]
myplotfacets_cviol[cviol_af_dat$Pearson_coef_littlea < -0.50 & cviol_af_dat$Pearson_coef_littlea > -0.80] <- myfacetIDs[3]
myplotfacets_cviol[cviol_af_dat$Pearson_coef_littlea <= -0.80] <- myfacetIDs[4]
cviol_af_dat$myplotfacets <- myplotfacets_cviol

ggplot(data=ccoru_af_dat, aes(x=thebin, y=littlea_freq, group=SNP), col='black') +
  facet_wrap(~myplotfacets) +
  geom_line(alpha=0.5) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab('Minor allele frequency') + xlab('Elevation bin (meters)') +
  ggtitle('Ccoru SNP MAF x elevation bin')
ggsave('C:/Users/mcwlim/Dropbox/Marisacompfiles/Seq cap files/Pop_GEA_chapter/Pop_GEA_figs/Ccoru_SNPMAFxElevbin_Pearsoncorr.jpg',
       height=6, width=6, units='in', dpi=600)
ggplot(data=cviol_af_dat, aes(x=thebin, y=littlea_freq, group=SNP), col='black') +
  facet_wrap(~myplotfacets) +
  geom_line(alpha=0.5) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab('Minor allele frequency') + xlab('Elevation bin (meters)') +
  ggtitle('Cviol SNP MAF x elevation bin')
ggsave('C:/Users/mcwlim/Dropbox/Marisacompfiles/Seq cap files/Pop_GEA_chapter/Pop_GEA_figs/Cviol_SNPMAFxElevbin_Pearsoncorr.jpg',
       height=6, width=6, units='in', dpi=600)


# plot function for minor allele with elev bins for specific SNPs
# plot_allelefrequency <- function(num_bins, SNPlist, myallelefreqfile, myspecies){
  # Version 1: plots and saves all SNPs to 1 plot
  # only plot of minor allele frequencies (littlea_freq)
  # version with all genes faceted (but now there are too many to see properly)
  # ggplot(data=myallelefreqfile, aes(x=thebin, y=littlea_freq, group=SNP)) +
  #   facet_wrap(~Gene) +
  #   geom_line() +
  #   #geom_point(size=4, alpha=0.6) +
  #   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  #   ylab('Minor allele frequency') + xlab('Elevation bin (meters)') +
  # ggtitle(paste(myspecies, ' SNP MAF x elevation bin', sep=''))
  # ggsave(paste('./LFMM_outlier_allelefreqplots/', myspecies, '_SNPMAF_elevbin_plot.jpg', sep=''),
  #        height=10, width=10, dpi=600)

  #for(i in 1:nrow(SNPlist)){
      # Version 2: plots and saves each SNP separately
      #dattoplot <- myallelefreqfile[myallelefreqfile$SNP == SNPlist$SNP[i], ]
      # only plot of minor allele frequencies (littlea_freq)
      # ggplot(data=dattoplot, aes(x=thebin, y=littlea_freq, group=SNP)) +
      #   geom_line() +
      #   geom_point(size=4, alpha=0.6) +
      #   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      #   ggtitle(paste(myspecies, '_', dattoplot$Gene[1], ' allele frequencies', '\n',
      #                 dattoplot$SNP[1], sep='')) +
      #   ylab('Minor allele frequency') + xlab('Elevation bin (meters)')
      # ggsave(paste('./LFMM_outlier_allelefreqplots/', myspecies, '_', dattoplot$Gene[1],dattoplot$SNP[1], '_plot.jpg', sep=''),
      #        height=6, width=10, dpi=600)
  #}
# }
# plot_allelefrequency(6, ccoru_wilcoxresults_subset, ccoru_af_dat, 'Ccoru')
# plot_allelefrequency(5, cviol_wilcoxresults_subset, cviol_af_dat, 'Cviol')

# ------- 7a) Ccoru 400m K2 lambda 1 allele frequency example plots specific genes -----
#SENP1
SENP1_af <- rbind(ccoru_af_all[ccoru_af_all$theSNP == 'combined_Contig174_6567', ],
                  ccoru_af_all[ccoru_af_all$theSNP == 'combined_Contig174_3984', ])
#RYR2
RYR2_af_ccoru <- rbind(ccoru_af_all[ccoru_af_all$theSNP == 'combined_Contig201_80046', ],
                       ccoru_af_all[ccoru_af_all$theSNP == 'combined_Contig201_80055', ],
                       ccoru_af_all[ccoru_af_all$theSNP == 'combined_Contig201_80053', ],
                       ccoru_af_all[ccoru_af_all$theSNP == 'combined_Contig201_86414', ],
                       ccoru_af_all[ccoru_af_all$theSNP == 'combined_Contig201_26562', ],
                       ccoru_af_all[ccoru_af_all$theSNP == 'combined_Contig201_80027', ],
                       ccoru_af_all[ccoru_af_all$theSNP == 'combined_Contig201_48638', ],
                       ccoru_af_all[ccoru_af_all$theSNP == 'combined_Contig201_45212', ],
                       ccoru_af_all[ccoru_af_all$theSNP == 'combined_Contig201_32896', ],
                       ccoru_af_all[ccoru_af_all$theSNP == 'combined_Contig201_77710', ],
                       ccoru_af_all[ccoru_af_all$theSNP == 'combined_Contig201_63323', ],
                       ccoru_af_all[ccoru_af_all$theSNP == 'combined_Contig201_30768', ],
                       ccoru_af_all[ccoru_af_all$theSNP == 'combined_Contig201_38067', ],
                       ccoru_af_all[ccoru_af_all$theSNP == 'combined_Contig201_56922', ],
                       ccoru_af_all[ccoru_af_all$theSNP == 'combined_Contig201_76638', ],
                       ccoru_af_all[ccoru_af_all$theSNP == 'combined_Contig201_68597', ],
                       ccoru_af_all[ccoru_af_all$theSNP == 'combined_Contig201_69754', ])
#ALOX5
ALOX5_af_ccoru <- rbind(ccoru_af_all[ccoru_af_all$theSNP == 'combined_Contig198_4637',],
                        ccoru_af_all[ccoru_af_all$theSNP == 'combined_Contig198_8465',])


plot_allelefrequency(SENP1_af, 'Elevation (400 m bins)', 'Ccoru400mbin', 'SENP1', 'bottom')
plot_allelefrequency(RYR2_af_ccoru[c(1:36),], 'Elevation (400 m bins)', 'Ccoru400mbin_set1', 'RYR2', 'bottom')
plot_allelefrequency(RYR2_af_ccoru[c(37:72),], 'Elevation (400 m bins)', 'Ccoru400mbin_set2', 'RYR2', 'bottom')
plot_allelefrequency(RYR2_af_ccoru[c(73:102),], 'Elevation (400 m bins)', 'Ccoru400mbin_set3', 'RYR2', 'bottom')
plot_allelefrequency(ALOX5_af_ccoru, 'Elevation (400 m bins)', 'Ccoru400mbin', 'ALOX5', 'bottom')

# ------- 7b) Cviol 400m K2 lambda 1 allele frequency example plots specific genes -----

#RYR2
RYR2_af_cviol <- rbind(cviol_af_all[cviol_af_all$theSNP == 'combined_Contig201_47345', ],
                       cviol_af_all[cviol_af_all$theSNP == 'combined_Contig201_6678', ],
                       cviol_af_all[cviol_af_all$theSNP == 'combined_Contig201_37216', ],
                       cviol_af_all[cviol_af_all$theSNP == 'combined_Contig201_76336', ])
#ALOX5
ALOX5_af_cviol <- cviol_af_all[cviol_af_all$theSNP == 'combined_Contig198_364',]


plot_allelefrequency(RYR2_af_cviol, 'Elevation (400 m bins)', 'Cviol400mbin', 'RYR2', 'bottom')
plot_allelefrequency(ALOX5_af_cviol, 'Elevation (400 m bins)', 'Cviol400mbin', 'ALOX5', 'bottom')






# -- 8. For given SNP, plot allele freq x elev bin & genotype x elevation boxplot -----
SNPplots <- function(thegeno, theloci, theannotfile, theSNP, thespecies_gene, thespecies, alleledat){

  alleledat2 <- alleledat[alleledat$SNP == theSNP,]
  allele_elevbin <- ggplot(data=alleledat2, aes(x=thebin, y=littlea_freq, group=SNP), col='black') +
    geom_line(linetype=2, alpha=0.5) +
    geom_point() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ylab('Minor allele frequency') + xlab('Elevation bin (meters)')

  # 1. Format data
  # make loci rowname, then transpose
  rownames(thegeno) <- theloci$V1
  thegeno2 <- t(thegeno)

  # combine annot file columns with loci and geno
  thedat <- cbind(theannotfile, thegeno2)

  #find the column number to subset data by and subset by SNP
  annot_cols <- ncol(theannotfile) # different for cviol vs. ccoru, because ccoru also has the season info
  SNPdat <- thedat[, c(1:annot_cols, which(colnames(thedat)==theSNP))]
  # remove rows with 9s (missing genotype data)
  SNPdat <- SNPdat[SNPdat[annot_cols + 1] != 9, ]

  colnames(SNPdat)[annot_cols + 1] <- 'Genotype'

  # 2. Plots
  # use Wilcoxon rank sum test (Mann Whitney U )
  # ggpubr package currently (as of May 2018) doesn't have a way to plot the adjusted
  # p-values [https://github.com/kassambara/ggpubr/issues/65]
  # so, instead, i'll plot the p.signif on boxplots - check that these reflect the
  # adjusted p-values, rather than the raw p-values

  mycomparisons <- list(c('0', '1'), c('1', '2'), c('0', '2'))
  geno_elev <- ggboxplot(data=SNPdat, x='Genotype', y='Elevation',
            color='Genotype', palette='jco', show.legend=FALSE) +
    stat_compare_means(comparisons=mycomparisons, label='p.signif')

  title <- ggdraw() +
    draw_label(paste(thespecies_gene, '\n', theSNP, sep=''),
               fontface = 'bold')
  theplot <- plot_grid(allele_elevbin, geno_elev, ncol=2, align = 'h', labels=c('a)', 'b)'))
  theplot2 <- plot_grid(title, theplot, ncol = 1, rel_heights = c(0.1, 1))
  save_plot(paste('C:/Users/mcwlim/Dropbox/Marisacompfiles/Seq cap files/Pop_GEA_chapter/Pop_GEA_figs/',thespecies, thespecies_gene, theSNP, '.jpg', sep=''),
            theplot2, base_height=5, base_width=8)
}
SNPplots(ccoru_geno, ccoru_loci, annot_ccorun97[ ,c(1,4,6:11)], 'combined_Contig174_6567', 'SENP1', 'Ccoru', ccoru_af_dat)
SNPplots(ccoru_geno, ccoru_loci, annot_ccorun97[ ,c(1,4,6:11)], 'combined_Contig174_3984', 'SENP1', 'Ccoru', ccoru_af_dat)
SNPplots(ccoru_geno, ccoru_loci, annot_ccorun97[ ,c(1,4,6:11)], 'combined_Contig270_8121', 'COX1', 'Ccoru', ccoru_af_dat)
SNPplots(cviol_geno, cviol_loci, annot_cvioln59[ ,c(1,4:9)], 'combined_Contig63_1750', 'WNT7B', 'Cviol', cviol_af_dat)
SNPplots(cviol_geno, cviol_loci, annot_cvioln59[ ,c(1,4:9)], 'combined_Contig261_2356', 'ASH2L', 'Cviol', cviol_af_dat)

# --------- plot overlapping genes in wilcox subset ----------
### to do:
  # make these all on one plot somehow
  # let's start by just making the allele freq plots (more complicated to figure out good way to combine boxplots)
# SMURF2, RYR2
smurf2_cviol <- cviol_af_dat[cviol_af_dat$Gene == 'SMURF2',]
smurf2_ccoru <- ccoru_af_dat[ccoru_af_dat$Gene == 'SMURF2',]

ryr2_cviol <- cviol_af_dat[cviol_af_dat$Gene == 'RYR2',]
ryr2_ccoru <- ccoru_af_dat[ccoru_af_dat$Gene == 'RYR2',]

smurf2cviol_allele_elevbin <- ggplot(data=smurf2_cviol, aes(x=thebin, y=littlea_freq, group=SNP), col='black') +
  geom_line(linetype=2, alpha=0.5) +
  geom_point(aes(col=SNP), size=3) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab('Minor allele frequency') + xlab('Elevation bin (meters)') +
  theme(legend.position='none')
smurf2ccoru_allele_elevbin <- ggplot(data=smurf2_ccoru, aes(x=thebin, y=littlea_freq, group=SNP), col='black') +
  geom_line(linetype=2, alpha=0.5) +
  geom_point(aes(col=SNP), size=3) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab('Minor allele frequency') + xlab('Elevation bin (meters)') +
  theme(legend.position = 'none')

ryr2cviol_allele_elevbin <- ggplot(data=ryr2_cviol, aes(x=thebin, y=littlea_freq, group=SNP), col='black') +
  geom_line(linetype=2, alpha=0.5) +
  geom_point(aes(col=SNP), size=3) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab('Minor allele frequency') + xlab('Elevation bin (meters)') +
  theme(legend.position = 'none')
ryr2ccoru_allele_elevbin <- ggplot(data=ryr2_ccoru, aes(x=thebin, y=littlea_freq, group=SNP), col='black') +
  geom_line(linetype=2, alpha=0.5) +
  geom_point(aes(col=SNP), size=3) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab('Minor allele frequency') + xlab('Elevation bin (meters)') +
  theme(legend.position = 'none')

title1 <- ggdraw() + draw_label('SMURF2', fontface = 'italic')
title2 <- ggdraw() + draw_label('RYR2', fontface = 'italic')
theplot <- plot_grid(smurf2cviol_allele_elevbin, smurf2ccoru_allele_elevbin,
                     ncol=2, align = 'h', labels=c('a)', 'b)'))
theplot2 <- plot_grid(title1, theplot, ncol = 1, rel_heights = c(0.1, 1))
theplot3 <- plot_grid(ryr2cviol_allele_elevbin, ryr2ccoru_allele_elevbin,
                      ncol=2, align = 'h', labels=c('c)', 'd)'))
theplot4 <- plot_grid(title2, theplot3, ncol=1, rel_heights=c(0.1, 1))
theplot5 <- plot_grid(theplot2, theplot4, ncol=1)
save_plot(paste('C:/Users/mcwlim/Dropbox/Marisacompfiles/Seq cap files/Pop_GEA_chapter/Pop_GEA_figs/cviolccoru_smurf2ryr2_allelefreqplots.jpg', sep=''),
          theplot5, base_height=8, base_width=8)


# -- 9. Get GO terms for the Wilcoxon subset of genes & novel candidates --------
# use panther database
# wilcox subset
write.table(paste(shQuote(sort(unique(cviol_wilcoxresults_subset$Gene)), type='cmd'), collapse=', '),
            'C:/Users/mcwlim/Dropbox/Marisacompfiles/Seq cap files/Pop_GEA_chapter/Pop_GEA_figs/cviol_wilcoxGenes_forpanther.txt',
          row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(paste(shQuote(sort(unique(ccoru_wilcoxresults_subset$Gene)), type='cmd'), collapse=', '),
            'C:/Users/mcwlim/Dropbox/Marisacompfiles/Seq cap files/Pop_GEA_chapter/Pop_GEA_figs/ccoru_wilcoxGenes_forpanther.txt',
            row.names=FALSE, col.names=FALSE, quote=FALSE)
# the novel candidates - get some kind of pie or barplot of GO terms for very basic functional information
write.table(paste(shQuote(sort(setdiff(unique(cviol_400mK2lamb1_05qval$Gene), unique(cviol_previouslyknowngenes$Gene))), type='cmd'), collapse=', '),
            'C:/Users/mcwlim/Dropbox/Marisacompfiles/Seq cap files/Pop_GEA_chapter/Pop_GEA_figs/cviol_Novelcands_forpanther.txt',
            row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(paste(shQuote(sort(setdiff(unique(ccoru_400mK2lamb1_05qval$Gene), unique(ccoru_previouslyknowngenes$Gene))), type='cmd'), collapse=', '),
            'C:/Users/mcwlim/Dropbox/Marisacompfiles/Seq cap files/Pop_GEA_chapter/Pop_GEA_figs/ccoru_Novelcands_forpanther.txt',
            row.names=FALSE, col.names=FALSE, quote=FALSE)

# ---------------- make barplots: Wilcox previously ID'd candidates --------------------
cviol_Wilcoxcand_pathway <- read.table('C:/Users/mcwlim/Dropbox/Marisacompfiles/Seq cap files/Pop_GEA_chapter/Pop_GEA_figs/pantherChart_cviolWilcoxgenes_pathway.txt',
                                      sep='\t')
cviol_Wilcoxcand_BP <- read.table('C:/Users/mcwlim/Dropbox/Marisacompfiles/Seq cap files/Pop_GEA_chapter/Pop_GEA_figs/pantherChart_cviolWilcoxgenes_bioprocess.txt',
                                 sep='\t')
ccoru_Wilcoxcand_pathway <- read.table('C:/Users/mcwlim/Dropbox/Marisacompfiles/Seq cap files/Pop_GEA_chapter/Pop_GEA_figs/pantherChart_ccoruWilcoxgenes_pathway.txt',
                                      sep='\t')
ccoru_Wilcoxcand_BP<- read.table('C:/Users/mcwlim/Dropbox/Marisacompfiles/Seq cap files/Pop_GEA_chapter/Pop_GEA_figs/pantherChart_ccoruWilcoxgenes_bioprocess.txt',
                                sep='\t')

levels(cviol_Wilcoxcand_BP$V2) <- gsub(' \\(', '\n(', levels(cviol_Wilcoxcand_BP$V2))
levels(ccoru_Wilcoxcand_BP$V2) <- gsub(' \\(', '\n(', levels(ccoru_Wilcoxcand_BP$V2))

cviol_bp_plot2 <- ggplot(cviol_Wilcoxcand_BP) +
  geom_bar(stat='identity', aes(y=V3, x=V2)) +
  geom_text(aes(x=V2, y=V3, label=V3), hjust=0, col='black', size=4.5) +
  coord_flip() + xlab('') + ylab('Number of genes')
cviol_pathway_plot2 <- ggplot(cviol_Wilcoxcand_pathway) +
  geom_bar(stat='identity', aes(y=V3, x=V2)) +
  geom_text(aes(x=V2, y=V3, label=V3), hjust=0, col='black', size=5) +
  coord_flip() + xlab('') + ylab('Number of genes')
ccoru_bp_plot2 <- ggplot(ccoru_Wilcoxcand_BP) +
  geom_bar(stat='identity', aes(y=V3, x=V2)) +
  geom_text(aes(x=V2, y=V3, label=V3), hjust=0, col='black', size=4.5) +
  coord_flip() + xlab('') + ylab('Number of genes')
ccoru_pathway_plot2 <- ggplot(ccoru_Wilcoxcand_pathway) +
  geom_bar(stat='identity', aes(y=V3, x=V2)) +
  geom_text(aes(x=V2, y=V3, label=V3), hjust=0, col='black', size=5) +
  coord_flip() + xlab('') + ylab('Number of genes')

bpplots2 <- plot_grid(cviol_bp_plot2, ccoru_bp_plot2, ncol=1, align='h', labels=c('a)', 'b)'))
pathwayplots2 <- plot_grid(cviol_pathway_plot2, ccoru_pathway_plot2, ncol=2, align='h', labels=c('a)', 'b)'))

save_plot('C:/Users/mcwlim/Dropbox/Marisacompfiles/Seq cap files/Pop_GEA_chapter/Pop_GEA_figs/CviolCcoru_Wilcoxcand_bpbar.jpg',
          bpplots2, base_height=10, base_width=10)
save_plot('C:/Users/mcwlim/Dropbox/Marisacompfiles/Seq cap files/Pop_GEA_chapter/Pop_GEA_figs/CviolCcoru_Wilcoxcand_pathwaysbar.jpg',
          pathwayplots2, base_height=15, base_width=25)
# ---------------- make barplots: novel candidates --------------------
# note: have to consolidate paths for cviol and ccoru, too many for plot viz

cviol_novelcand_pathway <- read.table('C:/Users/mcwlim/Dropbox/Marisacompfiles/Seq cap files/Pop_GEA_chapter/Pop_GEA_figs/pantherChart_cviolnovelcands_pathways.txt',
                                      sep='\t', header=T)
cviol_novelcand_BP <- read.table('C:/Users/mcwlim/Dropbox/Marisacompfiles/Seq cap files/Pop_GEA_chapter/Pop_GEA_figs/pantherChart_cviolnovelcands_bioprocess.txt',
                                      sep='\t')
ccoru_novelcand_pathway <- read.table('C:/Users/mcwlim/Dropbox/Marisacompfiles/Seq cap files/Pop_GEA_chapter/Pop_GEA_figs/pantherChart_ccorunovelcands_pathway.txt',
                                      sep='\t', header=T)
ccoru_novelcand_BP<- read.table('C:/Users/mcwlim/Dropbox/Marisacompfiles/Seq cap files/Pop_GEA_chapter/Pop_GEA_figs/pantherChart_ccorunovelcands_bioprocess.txt',
                                      sep='\t')

levels(cviol_novelcand_BP$V2) <- gsub(' \\(', '\n(', levels(cviol_novelcand_BP$V2))
#levels(cviol_novelcand_pathway$pathway) <- gsub(' \\(', '\n(', levels(cviol_novelcand_pathway$pathway))
levels(ccoru_novelcand_BP$V2) <- gsub(' \\(', '\n(', levels(ccoru_novelcand_BP$V2))
#levels(ccoru_novelcand_pathway$pathway) <- gsub(' \\(', '\n(', levels(ccoru_novelcand_pathway$pathway))

cviol_bp_plot <- ggplot(cviol_novelcand_BP) +
  geom_bar(stat='identity', aes(y=V3, x=V2)) +
  geom_text(aes(x=V2, y=V3, label=V3), hjust=0, col='black', size=4.5) +
  coord_flip() + xlab('') + ylab('Number of genes')
cviol_pathway_plot <- ggplot(cviol_novelcand_pathway) +
  geom_bar(stat='identity', aes(y=count, x=pathway)) +
  geom_text(aes(x=pathway, y=count, label=count), hjust=0, col='black', size=5) +
  coord_flip() + xlab('') + ylab('Number of genes')
ccoru_bp_plot <- ggplot(ccoru_novelcand_BP) +
  geom_bar(stat='identity', aes(y=V3, x=V2)) +
  geom_text(aes(x=V2, y=V3, label=V3), hjust=0, col='black', size=4.5) +
  coord_flip() + xlab('') + ylab('Number of genes')
ccoru_pathway_plot <- ggplot(ccoru_novelcand_pathway) +
  geom_bar(stat='identity', aes(y=count, x=pathway)) +
  geom_text(aes(x=pathway, y=count, label=count), hjust=0, col='black', size=5) +
  coord_flip() + xlab('') + ylab('Number of genes')

# this is by species, but due to diff in heights, doesn't look as nice
#plot_grid(cviol_bp_plot, cviol_pathway_plot, ncol=2)
#plot_grid(ccoru_bp_plot, ccoru_pathway_plot, ncol=2)

bpplots <- plot_grid(cviol_bp_plot, ccoru_bp_plot, ncol=1, align='h', labels=c('a)', 'b)'))
pathwayplots <- plot_grid(cviol_pathway_plot, ccoru_pathway_plot, ncol=2, align='h', labels=c('a)', 'b)'))

save_plot('C:/Users/mcwlim/Dropbox/Marisacompfiles/Seq cap files/Pop_GEA_chapter/Pop_GEA_figs/CviolCcoru_novelcand_bpbar.jpg',
          bpplots, base_height=10, base_width=10)
save_plot('C:/Users/mcwlim/Dropbox/Marisacompfiles/Seq cap files/Pop_GEA_chapter/Pop_GEA_figs/CviolCcoru_novelcand_pathwaysbar.jpg',
          pathwayplots, base_height=15, base_width=25)

# -- 10. Are Z chr genes in novel or known candidate gene set? --------------
Zchr_cviol <- read.table('cviol_Zcontigs_removefromLFMM.txt'); colnames(Zchr_cviol) <- c('Gene', 'Contig', 'SNP')
Zchr_ccoru <- read.table('ccoru_Zcontigs_removefromLFMM.txt'); colnames(Zchr_ccoru) <- c('Gene', 'Contig', 'SNP')
length(unique(Zchr_cviol$Gene)) #82 unique genes on Zchr (excluding the n/a for the contigs that don't have known gene symbol)
length(unique(Zchr_ccoru$Gene)) #66 unique genes on Zchr (excluding the n/a for the contigs that don't have known gene symbol)

# compare Zchr with _previouslyknowngenes
intersect(unique(Zchr_cviol$Gene), unique(cviol_previouslyknowngenes$Gene))
intersect(unique(Zchr_ccoru$Gene), unique(ccoru_previouslyknowngenes$Gene))
## 0 of the Z chr genes have significant SNPs on previously identified candidate genes (based on our search for candidate genes)

# compare Zchr with the setdiff novel cands
Zchr_novelcands_cviol <- intersect(unique(Zchr_cviol$Gene), setdiff(unique(cviol_400mK2lamb1_05qval$Gene), unique(cviol_previouslyknowngenes$Gene)))
length(Zchr_novelcands_cviol) # All 82 Z chr genes are from the novel cand set for Cviol
Zchr_novelcands_ccoru <- intersect(unique(Zchr_ccoru$Gene), setdiff(unique(ccoru_400mK2lamb1_05qval$Gene), unique(ccoru_previouslyknowngenes$Gene)))
length(Zchr_novelcands_ccoru) #All 66 Z chr genes are from the novel cand set for Ccoru

write.table(paste(shQuote(sort(Zchr_novelcands_cviol), type='cmd'), collapse=', '),
            'C:/Users/mcwlim/Dropbox/Marisacompfiles/Seq cap files/Pop_GEA_chapter/Pop_GEA_figs/cviol_Zchrgenelist.txt',
            row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(paste(shQuote(sort(Zchr_novelcands_ccoru), type='cmd'), collapse=', '),
            'C:/Users/mcwlim/Dropbox/Marisacompfiles/Seq cap files/Pop_GEA_chapter/Pop_GEA_figs/ccoru_Zchrgenelist.txt',
            row.names=FALSE, col.names=FALSE, quote=FALSE)
