# script to calculate delta K (by Evanno method) from likelihoods given by NgsAdmix

# ---- evanno calc function
K1dat_cviol <- read.table('K1_likelihood_cviol.txt')
K2dat_cviol <- read.table('K2_likelihood_cviol.txt')
K3dat_cviol <- read.table('K3_likelihood_cviol.txt')
K4dat_cviol <- read.table('K4_likelihood_cviol.txt')
K5dat_cviol <- read.table('K5_likelihood_cviol.txt')

K1dat_ccoru <- read.table('K1_likelihood_ccoru.txt')
K2dat_ccoru <- read.table('K2_likelihood_ccoru.txt')
K3dat_ccoru <- read.table('K3_likelihood_ccoru.txt')
K4dat_ccoru <- read.table('K4_likelihood_ccoru.txt')
K5dat_ccoru <- read.table('K5_likelihood_ccoru.txt')

cviolKdat <- cbind(K1dat_cviol, K2dat_cviol, K3dat_cviol, K4dat_cviol, K5dat_cviol)
ccoruKdat <- cbind(K1dat_ccoru, K2dat_ccoru, K3dat_ccoru, K4dat_ccoru, K5dat_ccoru)

evanno <- function(spKdat, mysp, mysptitle){
  colnames(spKdat) <- c('K1', 'K2', 'K3', 'K4', 'K5')

  # create empty matrix and then fill in values
  # this is for a matrix where [1:5, 1:10] = L values from NgsAdmix
  # [2:4, 11:20] = L''(k) calculated below
  # [2:4, 21] = row averages for L''(k) calculated below
  # [2:4, 22] = row st dev for L calculated below
  # [2:4, 23] = delta K calculated below
  num_K = 5
  num_iter = 10

  sp_matrix <- matrix(nrow=num_K,ncol=((num_iter*2)+3))
  # transpose K data
  t_spKdat <- t(spKdat)
  sp_matrix[1:num_K,1:num_iter] <- t_spKdat # add K data
  # calculate values for the L''(k) columns
  # iterate over columns (num_iter i) and rows (num_K j)

  Lkval_matrix <- matrix(nrow=num_K-2, ncol=num_iter)
  for (i in 1:num_iter){
    for(j in 1:(num_K-2)){
      Lkval <- abs(t_spKdat[j+2,i] - 2*t_spKdat[j+1,i] + t_spKdat[j,i])
      Lkval_matrix[j,i] <- Lkval
    }
  }
  sp_matrix[2:4,11:20] <- Lkval_matrix # add L''(k) data

  # calculate average L''(k)
  sp_matrix[2:4, 21] <- apply(Lkval_matrix, 1, mean)

  # calculate and add standard deviation L
  sp_matrix[2:4, 22] <- apply(t_spKdat[2:4,], 1, sd)

  # calculate delta K
  sp_matrix[1:5, 23] <- apply(sp_matrix[2:4,21:22], 1,
                                 function(x) sp_matrix[,21]/sp_matrix[,22])[,1]

  # plot delta K
  jpeg(paste(mysp, '_deltaK_plot.jpg', sep=''), height=6, width=6, units='in', res=600)
  plot(sp_matrix[2:4,23], pch=20, xlab='K', ylab='Delta K (Evanno method)')
  lines(sp_matrix[2:4,23])
  title(mysptitle)
  dev.off()
}

evanno(cviolKdat, 'cviol', 'C. violifer')
evanno(ccoruKdat, 'ccoru', 'C. coruscans')
