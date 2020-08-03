# allele frequency vs. elevation plot
# use genotype calls (0,1,2 format) from ANGSD -doGeno 2 to calculate allele frequencies
# for LFMM outlier SNPs of interest
# Bin elevation by 300 or 400 meter intervals

# note that all genotypes used in this code have been filtered to remove low frequency SNPs (MAF < 5%)

# libraries for plots
library(ggplot2)
library(reshape2)
library(cowplot)

# ----------- BELOW: code to calculate allele frequencies for specific SNPs (and plot) and ALL SNPs -------
# ---- Read in data ----
# Read in geno files
cviol_geno <- read.table('cvioln59_genotypes.txt')
ccoru_geno <- read.table('ccorun97_genotypes.txt')

# Read in loci
cviol_loci <- read.table('../Pcadapt_outflank/cviol59_outflank.loci')
ccoru_loci <- read.table('../Pcadapt_outflank/ccoru97_outflank.loci')

# Read in elevation and latitude bin files
cviol_bins <- read.table('cvioln59_lat_elevbins.txt', header=TRUE)
ccoru_bins <- read.table('ccorun97_lat_elevbins.txt', header=TRUE)

# ---- Format data for SPECIFIC SNP allele freq calculation ----
# add loci back to geno files
cviol_geno2 <- cbind(cviol_loci, cviol_geno); colnames(cviol_geno2)[1] <- 'SNP'
ccoru_geno2 <- cbind(ccoru_loci, ccoru_geno); colnames(ccoru_geno2)[1] <- 'SNP'

# Read in LFMM outlier SNP files (created in lfmm_qvalue_calc.r script)
cviol_outliers_epas1 <- read.table('cviol_epas1_snps.txt', header=TRUE)
ccoru_outliers_epas1 <- read.table('ccoru_epas1_snps.txt', header=TRUE)
ccoru_outliers_calhm1 <- read.table('ccoru_calhm1_snps.txt', header=TRUE)
cviol_outliers_sbds <- read.table('cviol_sbds_snps.txt', header=TRUE)

# ---- SPECIFIC SNPs: Subset indivs by bin, calculate allele frequency, plot and calc Pearson's cor coef ----

calcplot_allelefrequency <- function(mySNPfile, genofile, binfile, myspecies, thegene, legendposition){
  # loop through outlier SNPs
  bin_df <- data.frame('bin_type'='bin_type', 'thebin'='thebin', 'theSNP'='theSNP', 'bigA_freq'=1, 'littlea_freq'=1)
  for(k in 1:nrow(mySNPfile)){
    theSNP <- grep(mySNPfile$loci[k], genofile$SNP)

    mygenos <- data.frame(t(genofile[theSNP, c(2:ncol(genofile))]))
    colnames(mygenos) <- 'genotypes'

    #add to bin file
    geno_bin_files <- cbind(mygenos, binfile)

    # number types to loop through
    # there are always 3 columns: genotypes, Sample_order_index, and Specimen
    # the rest are bin categories (last columns)
    bin_names <- c(colnames(geno_bin_files[4:ncol(geno_bin_files)]))
    bin_cols <- c(4:ncol(geno_bin_files))
    bin_type_num <- ncol(geno_bin_files) - 3

    # first, subset bin groups
    for(i in 1:bin_type_num){
      bincol <- geno_bin_files[, c(1,bin_cols[i])]
      uniq_bins <- as.character(unique(bincol[order(bincol[2]),2]))

      for(j in 1:length(uniq_bins)){
        # subset for each bin and remove 9's (ignore missing data in calculations)
        binforfreqcalc <- bincol[bincol[2] == uniq_bins[j] & bincol$genotypes != 9,]
        alt_allele_total <- sum(binforfreqcalc$genotypes)
        total_indivs <- nrow(binforfreqcalc) * 2

        #calculate allele frequency (p^2 + 2pq + q^2 =1 for AA, Aa, aa)
        littleA <- alt_allele_total/total_indivs
        bigA <- 1-littleA
        allele_freqs <- data.frame('bin_type'=bin_names[i], 'thebin'= uniq_bins[j], 'theSNP'=mySNPfile$loci[k], 'bigA_freq'=bigA, 'littlea_freq'=littleA)

        # now, save results
        bin_df <- rbind(bin_df, allele_freqs)
      }
    }
  }
  bin_df2 <- bin_df[-1,]
  write.table(bin_df2, paste(myspecies, '_', thegene, '_allelefreqsbybins.txt', sep=''),
              row.names=FALSE, quote=FALSE)

  # make plots and calculate Pearson's correlation coeff
  # need to count number of bin types
  num_bintypes_plot <- length(unique(bin_df2$bin_type))
  bintypes <- unique(bin_df2$bin_type)

  cor_df <- data.frame('bin_type'='bin_type', 'theSNP'='theSNP', 'Pearson_coef_littlea'=1, 'p-value'=1)
  for(m in 1:num_bintypes_plot){
    bintoplot <- bin_df2[bin_df2$bin_type == bintypes[m],]
    bin_df3 <- melt(bintoplot, id.vars=c('bin_type', 'thebin', 'theSNP'))
    ggplot(bin_df3, aes(x=thebin, y=value, col=theSNP, group=theSNP)) +
      geom_line() +
      geom_point(size=4, alpha=0.6) +
      facet_grid(.~variable) + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position=legendposition) +
      ggtitle(paste(myspecies, '_', thegene, ' allele frequencies', sep='')) +
      ylab('Allele frequency') + xlab(bintypes[m])
    ggsave(paste('./allelefreq_plots/', myspecies, '_', thegene, '_', bintypes[m], '_plot.jpg', sep=''),
           height=6, width=10, dpi=600)

    uniq_snps <- unique(bintoplot$theSNP)
    for(n in 1:length(uniq_snps)){
      # Pearson's correlation. Change bin to numeric.
      snptocor <- bintoplot[bintoplot$theSNP == uniq_snps[n],]
      cor_coef_littlea <- cor.test(as.numeric(snptocor$thebin), snptocor$littlea_freq)
      cor_df <- rbind(cor_df, data.frame('bin_type'=snptocor$bin_type[1], 'theSNP'=snptocor$theSNP[1],
                                         'Pearson_coef_littlea'=cor_coef_littlea$estimate, 'p-value'=cor_coef_littlea$p.value))

    }
  }
  cor_df2 <- cor_df[-1,]
  write.table(cor_df2, paste(myspecies, '_', thegene, '_Pearsoncor_allelefreq-bin.txt', sep=''), row.names=FALSE)
}

calcplot_allelefrequency(cviol_outliers_epas1, cviol_geno2, cviol_bins, 'Cviol', 'EPAS1', 'right')
# break up ccoru epas1 snps because there are too many on 1 plot
calcplot_allelefrequency(ccoru_outliers_epas1[1:6,], ccoru_geno2, ccoru_bins, 'Ccoru', 'EPAS1_set1', 'right')
calcplot_allelefrequency(ccoru_outliers_epas1[7:13,], ccoru_geno2, ccoru_bins, 'Ccoru', 'EPAS1_set2', 'right')
calcplot_allelefrequency(ccoru_outliers_epas1[14:20,], ccoru_geno2, ccoru_bins, 'Ccoru', 'EPAS1_set3', 'right')
calcplot_allelefrequency(ccoru_outliers_epas1[21:27,], ccoru_geno2, ccoru_bins, 'Ccoru', 'EPAS1_set4', 'right')
calcplot_allelefrequency(ccoru_outliers_epas1[28:33,], ccoru_geno2, ccoru_bins, 'Ccoru', 'EPAS1_set5', 'right')

# these were the top genes with outlier snps (smallest q-values; though not necessarily biologically meaningful)
calcplot_allelefrequency(ccoru_outliers_calhm1, ccoru_geno2, ccoru_bins, 'Ccoru', 'CALHM1', 'right')
calcplot_allelefrequency(cviol_outliers_sbds, cviol_geno2, cviol_bins, 'Cviol', 'SBDS', 'right')

# ---- Format data for ALL SNP allele freq calculation ----
# note: for comparison with LFMM results, not really needed to do this for sgeof or ccoel for ALL SNPs. Only needed for specific SNPs section

# make geno files with SNPs as rownames
rownames(cviol_geno) <- cviol_loci$V1
rownames(ccoru_geno) <- ccoru_loci$V1

dim(cviol_geno)
dim(ccoru_geno)

# ---- ALL SNPs: Subset indivs by bin, calculate allele frequency, calc Pearson's cor coef ----

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
  cor_df <- data.frame('bin_type'='bin_type', 'theSNP'='theSNP', 'Pearson_coef_littlea'=1, 'p-value'=1)
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
    # Pearson's correlation. Change bin to numeric.
    # remove any rows with NAs (SNPs that have bins with missing allele frequencies)
    # have to make sure there's a way to still match up the SNP loci for later analysis
    snptocor <- bin_df[bin_df$theSNP == colnames(snpcol)[1],]
    if(sum(is.na(snptocor$bigA_freq)) == 0){ # if 0 is TRUE, then there are no frequencies with 0 and corr can be calculated
      cor_coef_littlea <- cor.test(as.numeric(snptocor$thebin), snptocor$littlea_freq)
      cor_df <- rbind(cor_df, data.frame('bin_type'=snptocor$bin_type[1], 'theSNP'=snptocor$theSNP[1],
                                         'Pearson_coef_littlea'=cor_coef_littlea$estimate, 'p-value'=cor_coef_littlea$p.value))
    }else{
      # there is at least 1 row of NA frequencies
      # remove from bin_df and don't add to cor_df
      bin_df <- bin_df[bin_df$theSNP != colnames(snpcol)[1],]
      print(paste(colnames(snpcol)[1], ' SNP has been removed', sep=''))
    }
  }

  bin_df2 <- bin_df[-1,]
  write.table(bin_df2, paste(myspecies, '_', mybin, '_', snpbatch, '_allelefreqsbybins.txt', sep=''),
              row.names=FALSE, quote=FALSE)
  cor_df2 <- cor_df[-1,]
  write.table(cor_df2, paste(myspecies, '_', mybin, '_', snpbatch, '_Pearsoncor_allelefreq-bin.txt', sep=''),
              row.names=FALSE, quote=FALSE)
}

calc_allelefrequency(cviol_geno[1:10000,], cviol_bins, 'Elevation_bins400m', 'Cviol', 'snpbatch1')
    # "combined_Contig20_17118 SNP has been removed"
calc_allelefrequency(cviol_geno[10001:20000,], cviol_bins, 'Elevation_bins400m', 'Cviol', 'snpbatch2')
    # "combined_Contig264_1522 SNP has been removed"
calc_allelefrequency(cviol_geno[20001:30000,], cviol_bins, 'Elevation_bins400m', 'Cviol', 'snpbatch3')
calc_allelefrequency(cviol_geno[30001:40000,], cviol_bins, 'Elevation_bins400m', 'Cviol', 'snpbatch4')
calc_allelefrequency(cviol_geno[40001:50000,], cviol_bins, 'Elevation_bins400m', 'Cviol', 'snpbatch5')
calc_allelefrequency(cviol_geno[50001:60000,], cviol_bins, 'Elevation_bins400m', 'Cviol', 'snpbatch6')
calc_allelefrequency(cviol_geno[60001:70000,], cviol_bins, 'Elevation_bins400m', 'Cviol', 'snpbatch7')
    # "combined_Contig8287_422 SNP has been removed"
calc_allelefrequency(cviol_geno[70001:80000,], cviol_bins, 'Elevation_bins400m', 'Cviol', 'snpbatch8')
calc_allelefrequency(cviol_geno[80001:91124,], cviol_bins, 'Elevation_bins400m', 'Cviol', 'snpbatch9')
    # "combined_Contig10793_215 SNP has been removed"
    # "combined_Contig12216_852 SNP has been removed"

# 400 m bins
calc_allelefrequency(ccoru_geno[1:20000,], ccoru_bins, 'Elevation_bin400m', 'Ccoru', 'snpbatch1')
calc_allelefrequency(ccoru_geno[20001:40000,], ccoru_bins, 'Elevation_bin400m', 'Ccoru', 'snpbatch2')
calc_allelefrequency(ccoru_geno[40001:60000,], ccoru_bins, 'Elevation_bin400m', 'Ccoru', 'snpbatch3')
calc_allelefrequency(ccoru_geno[60001:80000,], ccoru_bins, 'Elevation_bin400m', 'Ccoru', 'snpbatch4')
calc_allelefrequency(ccoru_geno[80001:100000,], ccoru_bins, 'Elevation_bin400m', 'Ccoru', 'snpbatch5')
calc_allelefrequency(ccoru_geno[100001:120000,], ccoru_bins, 'Elevation_bin400m', 'Ccoru', 'snpbatch6')
calc_allelefrequency(ccoru_geno[120001:140000,], ccoru_bins, 'Elevation_bin400m', 'Ccoru', 'snpbatch7')
calc_allelefrequency(ccoru_geno[140001:160000,], ccoru_bins, 'Elevation_bin400m', 'Ccoru', 'snpbatch8')
calc_allelefrequency(ccoru_geno[160001:180000,], ccoru_bins, 'Elevation_bin400m', 'Ccoru', 'snpbatch9')
calc_allelefrequency(ccoru_geno[180001:200000,], ccoru_bins, 'Elevation_bin400m', 'Ccoru', 'snpbatch10')
calc_allelefrequency(ccoru_geno[200001:215681,], ccoru_bins, 'Elevation_bin400m', 'Ccoru', 'snpbatch11')
    # no SNPs removed for Ccoru at Elevation_bin400m bin type

# 300 m bins
calc_allelefrequency(ccoru_geno[1:20000,], ccoru_bins, 'Elevation_bin300m', 'Ccoru', 'snpbatch1')
  # "combined_Contig14_9293 SNP has been removed"
calc_allelefrequency(ccoru_geno[20001:40000,], ccoru_bins, 'Elevation_bin300m', 'Ccoru', 'snpbatch2')
  # "combined_Contig127_2650 SNP has been removed"
calc_allelefrequency(ccoru_geno[40001:60000,], ccoru_bins, 'Elevation_bin300m', 'Ccoru', 'snpbatch3')
calc_allelefrequency(ccoru_geno[60001:80000,], ccoru_bins, 'Elevation_bin300m', 'Ccoru', 'snpbatch4')
calc_allelefrequency(ccoru_geno[80001:100000,], ccoru_bins, 'Elevation_bin300m', 'Ccoru', 'snpbatch5')
calc_allelefrequency(ccoru_geno[100001:120000,], ccoru_bins, 'Elevation_bin300m', 'Ccoru', 'snpbatch6')
calc_allelefrequency(ccoru_geno[120001:140000,], ccoru_bins, 'Elevation_bin300m', 'Ccoru', 'snpbatch7')
calc_allelefrequency(ccoru_geno[140001:160000,], ccoru_bins, 'Elevation_bin300m', 'Ccoru', 'snpbatch8')
calc_allelefrequency(ccoru_geno[160001:180000,], ccoru_bins, 'Elevation_bin300m', 'Ccoru', 'snpbatch9')
calc_allelefrequency(ccoru_geno[180001:200000,], ccoru_bins, 'Elevation_bin300m', 'Ccoru', 'snpbatch10')
calc_allelefrequency(ccoru_geno[200001:215681,], ccoru_bins, 'Elevation_bin300m', 'Ccoru', 'snpbatch11')
  # "combined_Contig11474_955 SNP has been removed"
  # "combined_Contig11474_957 SNP has been removed"
  # "combined_Contig11474_964 SNP has been removed"
