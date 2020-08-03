# Create variant + non-variant sites geno files
# This is for calculating dxy with Ke's 6-PopStats script

setwd('/pylon5/bi4iflp/mlim/SeqCapData/angsd_analysis/Fst_dxy_theta/dxy/snpcallgenofiles_fordxy')

# unique non-variant sites dimensions
#C. coeligena 9,237,986 sites x 35 individuals
#C. violifer 9,384,035 sites x 59 individuals
#C. coruscans 8,606,434 sites x 97 individuals
#S. geoffroyi 9,592,855 sites x 28 individuals

#1 Create the non-variant matrices of zeros
ccoel_nonvar <- matrix(0L, ncol=35, nrow=9237986)
cviol_nonvar <- matrix(0L, ncol=59, nrow=9384035)
ccoru_nonvar <- matrix(0L, ncol=97, nrow=8606434)
sgeof_nonvar <- matrix(0L, ncol=28, nrow=9592855)

#2 Add contig # and site #
# First 2 columns are contig # and site #
ccoel_names <- matrix(c(rep('nonvar_contig',9237986), seq(1:9237986)), ncol=2,nrow=9237986)
cviol_names <- matrix(c(rep('nonvar_contig',9384035), seq(1:9384035)), ncol=2,nrow=9384035)
ccoru_names <- matrix(c(rep('nonvar_contig',8606434), seq(1:8606434)), ncol=2,nrow=8606434)
sgeof_names <- matrix(c(rep('nonvar_contig',9592855), seq(1:9592855)), ncol=2,nrow=9592855)

#3 cbind the matrices together
ccoel_nonvar2 <- cbind(ccoel_names, ccoel_nonvar)
cviol_nonvar2 <- cbind(cviol_names, cviol_nonvar)
ccoru_nonvar2 <- cbind(ccoru_names, ccoru_nonvar)
sgeof_nonvar2 <- cbind(sgeof_names, sgeof_nonvar)

#4 Save matrix as .gz file
write.table(ccoel_nonvar2, file='ccoel_nonvar0mat.txt',sep='\t',col.names=FALSE,row.names=FALSE)
system('gzip ccoel_nonvar0mat.txt')
write.table(cviol_nonvar2, file='cviol_nonvar0mat.txt',sep='\t',col.names=FALSE,row.names=FALSE)
system('gzip cviol_nonvar0mat.txt')
write.table(ccoru_nonvar2, file='ccoru_nonvar0mat.txt',sep='\t',col.names=FALSE,row.names=FALSE)
system('gzip ccoru_nonvar0mat.txt')
write.table(sgeof_nonvar2, file='sgeof_nonvar0mat.txt',sep='\t',col.names=FALSE,row.names=FALSE)
system('gzip sgeof_nonvar0mat.txt')



# test code on small dataset
#ccoel_nonvar <- matrix(0L, ncol=10, nrow=10)
#ccoel_names <- matrix(c(rep('nonvar_contig',10), rep('1',10)), ncol=2,nrow=10)
#ccoel_nonvar2 <- cbind(ccoel_names, ccoel_nonvar)
#class(ccoel_nonvar2)
#write.csv(ccoel_nonvar2, 'ccoelnonvar.csv')
#system('gzip ccoelnonvar.csv')
