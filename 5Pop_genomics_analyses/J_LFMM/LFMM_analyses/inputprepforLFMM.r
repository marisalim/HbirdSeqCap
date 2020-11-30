# 1. calculate allele frequencies
# 2. create new input files for LFMM
# 3. analyze lambda and p-value distributions for redo LFMM results

# Load libraries
library(reshape2)
# qvalue isn't available from CRAN for R v4, so you have to install from bioconductor
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("qvalue")
library(qvalue)

# ------ Prep inputs for LFMM analysis ----
# ---- original files ----
# example code for one species
cviol_geno <- read.table('cvioln59_genotypes.txt')
cviol_loci <- read.table('cviol59_outflank.loci')
cviol_bins <- read.table('cvioln59_lat_elevbins.txt', header=TRUE)

# ---- Calculate allele frequencies ----
# this can't run on my computer with all the snps - memory crashes
# so, run in batches
calc_allelefrequency <- function(genofile, binfile, mybin, myspecies, snpbatch){
  
  #add to bin file to geno file
  #binfile first 2 cols = index and specimen, rest are bins
  targetbin <- data.frame(binfile[, mybin])
  colnames(targetbin) <- mybin
  geno_bin_files <- cbind(t(genofile), targetbin)
  
  # number snps and bin types to loop through
  n_snps <- nrow(genofile)
  bin_cols <- ncol(geno_bin_files)
  
  bin_df <- data.frame('bin_type'='bin_type', 'thebin'='thebin', 'theSNP'='theSNP', 'bigA_freq'=1, 'littlea_freq'=1)
  # cor_df <- data.frame('bin_type'='bin_type', 'theSNP'='theSNP', 'Pearson_coef_littlea'=1, 'p-value'=1)
  #loop thru each SNP column
  for(k in 1:n_snps){
    if(k %% 500 == 0){
      cat(paste0('iterations: ', k, '\n')) # print out iterations for every 500
    }
    
    # get snp
    snpcol <- geno_bin_files[ , c(k, bin_cols)]
    
    uniq_bins <- as.character(unique(snpcol[order(snpcol[2]),2]))
    # then, subset bins
    for(j in 1:length(uniq_bins)){
      binsnps <- snpcol[snpcol[2] == uniq_bins[j],]
      # subset for each bin and remove 9's (ignore missing data in calculations)
      if(sum(binsnps[1] == 9) == nrow(binsnps)){
        # if there are no 0s, 1s, or 2s, then can't calculate allele frequency
        # remove SNP from dataset because 1) no information about how LFMM handles missing allele frequency and 
        # 2) can't calculate correlation coefficient for snp with missing bin frequencies
        littleA <- NA
        bigA <- NA
        
      }else{
        # otherwise, remove 9s and add up the rest
        binforfreqcalc <- binsnps[binsnps[1] != 9,]
        alt_allele_total <- sum(binforfreqcalc[1])
        total_indivs <- nrow(binforfreqcalc) * 2
        
        #calculate allele frequency (p^2 + 2pq + q^2 =1 for AA, Aa, aa)
        littleA <- alt_allele_total/total_indivs
        bigA <- 1-littleA
      }
      allele_freqs <- data.frame('bin_type'=mybin, 'thebin'= uniq_bins[j], 'theSNP'=colnames(snpcol)[1], 'bigA_freq'=bigA, 'littlea_freq'=littleA)
      
      # now, save results
      bin_df <- rbind(bin_df, allele_freqs)
    }
    
    # remove any rows with NAs (SNPs that have bins with missing allele frequencies)
    snptocor <- bin_df[bin_df$theSNP == colnames(snpcol)[1],]
    if(sum(is.na(snptocor$bigA_freq)) != 0){
      # there is at least 1 row of NA frequencies - need this section to remove NAs!
      # remove from bin_df and don't add to cor_df
      bin_df <- bin_df[bin_df$theSNP != colnames(snpcol)[1],]
      print(paste(colnames(snpcol)[1], ' SNP has been removed', sep=''))
    }
  }
  
  bin_df2 <- bin_df[-1,]
  write.table(bin_df2, paste(myspecies, '_', mybin, '_', snpbatch, '_allelefreqsbybins.txt', sep=''), 
              row.names=FALSE, quote=FALSE)
}

# run function, note snps that are removed because can't calculate allele frequency
calc_allelefrequency(cviol_geno[1:10000,], cviol_bins_new, 'Elevation_bins400m', 'Cviol', 'snpbatch1')
calc_allelefrequency(cviol_geno[10001:20000,], cviol_bins, 'Elevation_bins400m', 'Cviol', 'snpbatch2')
calc_allelefrequency(cviol_geno[20001:30000,], cviol_bins, 'Elevation_bins400m', 'Cviol', 'snpbatch3')
calc_allelefrequency(cviol_geno[30001:40000,], cviol_bins, 'Elevation_bins400m', 'Cviol', 'snpbatch4')
calc_allelefrequency(cviol_geno[40001:50000,], cviol_bins, 'Elevation_bins400m', 'Cviol', 'snpbatch5')
calc_allelefrequency(cviol_geno[50001:60000,], cviol_bins, 'Elevation_bins400m', 'Cviol', 'snpbatch6')
calc_allelefrequency(cviol_geno[60001:70000,], cviol_bins, 'Elevation_bins400m', 'Cviol', 'snpbatch7')
calc_allelefrequency(cviol_geno[70001:80000,], cviol_bins, 'Elevation_bins400m', 'Cviol', 'snpbatch8')
calc_allelefrequency(cviol_geno[80001:91124,], cviol_bins, 'Elevation_bins400m', 'Cviol', 'snpbatch9')
  
# ---- format input files for LFMM ---- 
# set up allele frequency and environmental files for LFMM
# Cviol 400m bins
# environmental file:
# make vector of average elevation for each bin
# use mean and std of all elevations per species to standardize the avg bin elevations

# rows are elevation bins, columns are SNP allele frequencies

# read in allele frequency by bin files
cviol_allelefreqs1 <- read.table('./Cviol_Elevation_bins400m_snpbatch1_allelefreqsbybins.txt', header=TRUE)
cviol_allelefreqs2 <- read.table('./Cviol_Elevation_bins400m_snpbatch2_allelefreqsbybins.txt', header=TRUE)
cviol_allelefreqs3 <- read.table('./Cviol_Elevation_bins400m_snpbatch3_allelefreqsbybins.txt', header=TRUE)
cviol_allelefreqs4 <- read.table('./Cviol_Elevation_bins400m_snpbatch4_allelefreqsbybins.txt', header=TRUE)
cviol_allelefreqs5 <- read.table('./Cviol_Elevation_bins400m_snpbatch5_allelefreqsbybins.txt', header=TRUE)
cviol_allelefreqs6 <- read.table('./Cviol_Elevation_bins400m_snpbatch6_allelefreqsbybins.txt', header=TRUE)
cviol_allelefreqs7 <- read.table('./Cviol_Elevation_bins400m_snpbatch7_allelefreqsbybins.txt', header=TRUE)
cviol_allelefreqs8 <- read.table('./Cviol_Elevation_bins400m_snpbatch8_allelefreqsbybins.txt', header=TRUE)
cviol_allelefreqs9 <- read.table('./Cviol_Elevation_bins400m_snpbatch9_allelefreqsbybins.txt', header=TRUE)

cviol_allelefreqs_all <- rbind(cviol_allelefreqs1, cviol_allelefreqs2, cviol_allelefreqs3, cviol_allelefreqs4, cviol_allelefreqs5, cviol_allelefreqs6,
                               cviol_allelefreqs7, cviol_allelefreqs8, cviol_allelefreqs9)

#format LFMM input
cviol_allelefreqs_for_lfmm <- acast(cviol_allelefreqs_all, thebin~theSNP, value='littlea_freq')
dim(cviol_allelefreqs_for_lfmm) # 0 SNPs were removed because bin(s) had all 9's, so there are still 91124 snps
write.table(cviol_allelefreqs_for_lfmm, 'cviol_Elev400m_allelefreqs_forLFMM.txt', row.names=FALSE, col.names=FALSE)

cviol_lfmm_snpnames <- colnames(cviol_allelefreqs_for_lfmm); length(cviol_lfmm_snpnames) #91124
write.table(cviol_lfmm_snpnames, 'cviol_lfmmallelefreq_snpnames.txt', row.names=FALSE, col.names=FALSE)

# ---- make elevation files ----
# calculate average and std of all elevations
# calculate average elevation per bin
# calculate standardized elevation per bin
cviol_elev <- read.table('./cviolLFMM_elevs.txt', sep='\t', header=TRUE)

hist(cviol_elev$Elevation, breaks=30)
abline(v=mean(cviol_elev$Elevation), col='tomato')
abline(v=median(cviol_elev$Elevation), col='turquoise')
# median tracks the higher frequency bin better than mean

cviol_elev_bins <- c(2200,2600,3000,3400,3800)
# this is for LFMM
cviol_meanstandardized_elev <- sapply(cviol_elev_bins, function(x){
  mean_elev <- mean(cviol_elev$Elevation)
  std_elev <- sd(cviol_elev$Elevation)
  standardized_elev <- (x - mean_elev)/std_elev
})

hist(cviol_meanstandardized_elev, col='turquoise', main='Cviol standardized w/ mean')
write.table(cviol_meanstandardized_elev, 'Cviol_meanstandardized_elevationbins400m.txt', col.names=FALSE, row.names=FALSE)

# now you have input files for LFMM!
# run command line LFMM and then move on to post-LFMM analysis below

# ------ Post-LFMM analysis ----
# ---- lambda and p-value distribution ----
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
  jpeg(paste('./LFMM_adjustedpvalue_hist_', mysp, 'K', myK, '_', round(lambda, 3), '_', myiters, '.jpg', sep=''), height=6, width=6, units='in', res=600)
  hist(ap.values, col = "red", main=paste('K: ', myK, ', Lambda: ', round(lambda,2), sep=''))
  dev.off()
  # calculate q-values
  qobj <- qvalue(ap.values, fdr.level=0.05, pi0.method = "bootstrap")
  write.qvalue(qobj, file=paste('./qval_', mysp, "_K",myK, '_', round(lambda,3), myiters,".txt", sep=''))
  capture.output(summary(qobj), file=paste('./qvalstats_', mysp, '_K', myK, '_', round(lambda, 3), myiters, '.txt', sep=''))
  
  totalSNPs <- length(qobj$significant) 
  sig_SNPs <- length(qobj[qobj$significant == TRUE]) 
  nonsig_SNPs <- length(qobj[qobj$significant == FALSE])
  
  mytable <- table(lambda, totalSNPs, sig_SNPs, nonsig_SNPs)
  write.table(mytable, paste('./lfmmstats_', mysp, '_K', myK, '_', round(lambda, 3), myiters, '.txt', sep=''))
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
  jpeg(paste('./LFMM_adjustedpvalue_hist_', mysp, 'K', myK, '_', round(lambda, 3), '_', myiters, '.jpg', sep=''), height=6, width=6, units='in', res=600)
  hist(ap.values, col = "red", main=paste('K: ', myK, ', Lambda: ', round(lambda,2), sep=''))
  dev.off()
  # calculate q-values
  qobj <- qvalue(ap.values, fdr.level=0.05, pi0.method = "bootstrap")
  write.qvalue(qobj, file=paste('./qval_', mysp, "_K",myK, '_', lambda, myiters,".txt", sep=''))
  capture.output(summary(qobj), file=paste('./qvalstats_', mysp, '_K', myK, '_', round(lambda, 3), myiters, '.txt', sep=''))
  
  totalSNPs <- length(qobj$significant) 
  sig_SNPs <- length(qobj[qobj$significant == TRUE]) 
  nonsig_SNPs <- length(qobj[qobj$significant == FALSE])
  
  mytable <- table(lambda, totalSNPs, sig_SNPs, nonsig_SNPs)
  write.table(mytable, paste('./lfmmstats_', mysp, '_K', myK, '_', round(lambda, 3), myiters, '.txt', sep=''))
}

# K=2
lfmm_qvalcalc('./lfmm_results/run_K', '2', '_i250000b25000_iter', 'Cviol')
lfmm_qvalcalc_chooselambda('./lfmm_results/run_K', '2', '_i250000b25000_iter', 1, 'Cviol')
lfmm_qvalcalc_chooselambda('./lfmm_results/run_K', '2', '_i250000b25000_iter', 2, 'Cviol') # too conservative
lfmm_qvalcalc_chooselambda('./lfmm_results/run_K', '2', '_i250000b25000_iter', 0.5, 'Cviol') # distribution flatted out, but maybe too liberal now

# using K=2, lambda=1
