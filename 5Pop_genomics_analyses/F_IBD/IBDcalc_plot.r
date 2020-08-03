# calculate pairwise geodesic distances to plot with pairwise genetic distance
# to test for isolation by distance

library(Imap)
library(ggplot2)
library(cowplot)

# read in the lat/long filesc
annotpath <- './2nd pass ANGSD - remove outliers/ngsTools PCA - v2/'
annot_cvioln59 <- read.table(paste(annotpath,'cviol_n59.clst', sep=''), sep="\t", header=T)
annot_ccorun97 <- read.table(paste(annotpath,'ccoru_n97.clst', sep=''), sep='\t', header=T)

# calculate pairwise distances
# using functions that rely on gdist() from Imap package to calculate distance matrix
# modified from here: https://eurekastatistics.com/calculating-a-distance-matrix-for-geographic-points-using-r/
ReplaceLowerOrUpperTriangle <- function(m, triangle.to.replace){
  # If triangle.to.replace="lower", replaces the lower triangle of a square matrix with its upper triangle.
  # If triangle.to.replace="upper", replaces the upper triangle of a square matrix with its lower triangle.

  if (nrow(m) != ncol(m)) stop("Supplied matrix must be square.")
  if      (tolower(triangle.to.replace) == "lower") tri <- lower.tri(m)
  else if (tolower(triangle.to.replace) == "upper") tri <- upper.tri(m)
  else stop("triangle.to.replace must be set to 'lower' or 'upper'.")
  m[tri] <- t(m)[tri]
  return(m)
}
GeoDistanceInMetresMatrix <- function(df.geopoints){
  # Returns a matrix (M) of distances between geographic points.
  # M[i,j] = M[j,i] = Distance between (df.geopoints$lat[i], df.geopoints$lon[i]) and
  # (df.geopoints$lat[j], df.geopoints$lon[j]).
  # The row and column names are given by df.geopoints$name.

  GeoDistanceInMetres <- function(g1, g2){
    # Returns a vector of distances. (But if g1$index > g2$index, returns zero.)
    # The 1st value in the returned vector is the distance between g1[[1]] and g2[[1]].
    # The 2nd value in the returned vector is the distance between g1[[2]] and g2[[2]]. Etc.
    # Each g1[[x]] or g2[[x]] must be a list with named elements "index", "lat" and "lon".
    # E.g. g1 <- list(list("index"=1, "lat"=12.1, "lon"=10.1), list("index"=3, "lat"=12.1, "lon"=13.2))
    DistM <- function(g1, g2){
      require("Imap")
      return(ifelse(g1$index > g2$index, 0, gdist(lat.1=g1$DEC_LAT, lon.1=g1$DEC_LONG, lat.2=g2$DEC_LAT, lon.2=g2$DEC_LONG, units="m")))
    }
    return(mapply(DistM, g1, g2))
  }

  n.geopoints <- nrow(df.geopoints)

  # The index column is used to ensure we only do calculations for the upper triangle of points
  df.geopoints$index <- 1:n.geopoints

  # Create a list of lists
  list.geopoints <- by(df.geopoints[,c("index", "DEC_LAT", "DEC_LONG")], 1:n.geopoints, function(x){return(list(x))})

  # Get a matrix of distances (in metres)
  mat.distances <- ReplaceLowerOrUpperTriangle(outer(list.geopoints, list.geopoints, GeoDistanceInMetres), "lower")

  # Set the row and column names
  rownames(mat.distances) <- df.geopoints$name
  colnames(mat.distances) <- df.geopoints$name

  return(mat.distances)
}

# units should be in km, since dividing by 1000
geodist_cviol <- round(GeoDistanceInMetresMatrix(annot_cvioln59[,c(4,6,7)]) / 1000, 2)
geodist_ccoru <- round(GeoDistanceInMetresMatrix(annot_ccorun97[,c(4,7,8)]) / 1000, 2)

# read in genetic distances from NGSdist
# using the NGSdist matrixes calculated for EEMS, deleted blank 1st line and 2nd line that had sample size
# the first column is Ind#
IBDpath <- './2nd pass ANGSD - remove outliers/Isolation_by_Dist/'
gendist_cvioln59 <- read.table(paste(IBDpath,'cviol59_forIBD.txt', sep=''), sep='\t')[-1]
gendist_ccorun97 <- read.table(paste(IBDpath,'ccoru97_forIBD.txt', sep=''), sep='\t')[-1]

# plot
jitterfactor <- 100
jpeg(paste(IBDpath, 'IBD_4sp_jittered.jpg', sep=''), height=8, width=8, units='in', res=600)
par(mfrow=c(2,2))
plot(gendist_ccorun97[lower.tri(gendist_ccorun97)]~
       jitter(geodist_ccoru[lower.tri(geodist_ccoru)], factor=jitterfactor), pch=20,
     ylab='Genetic distance', xlab='Geographic distance (km)', main='C. coruscans n=97')
plot(gendist_cvioln59[lower.tri(gendist_cvioln59)]~ 
     jitter(geodist_cviol[lower.tri(geodist_cviol)], factor=jitterfactor), pch=20,
     ylab='Genetic distance', xlab='Geographic distance (km)', main='C. violifer n=59')
dev.off()

# check the Cviol individuals with pairwise geographic disance ~700km
cviolcheck_geo <- as.vector(geodist_cviol[lower.tri(geodist_cviol)])
cviolcheck_gen <- as.vector(gendist_cvioln59[lower.tri(gendist_cvioln59)])

which((cviolcheck_geo > 680) & (cviolcheck_geo < 790))
cviolcheck_geo[18:23]; cviolcheck_gen[18:23]
  # top cluster
cviolcheck_geo[75:80]; cviolcheck_gen[75:80]
  # lower cluster - bottom
cviolcheck_geo[131:136]; cviolcheck_gen[131:136]
  # lower cluster - top
cviolcheck_geo[851:856]; cviolcheck_gen[851:856]
  # middle cluster
# how to figure out which pairs of individuals these are??
cviol_geodist_withsamps <- cbind(as.character(annot_cvioln59$Specimen), as.character(annot_cvioln59$Geographic_subpopulations), geodist_cviol)
cviol_gendist_withsamps <- cbind(as.character(annot_cvioln59$Specimen), as.character(annot_cvioln59$Geographic_subpopulations), gendist_cvioln59)
View(cviol_geodist_withsamps) #search for 706.85
View(cviol_gendist_withsamps)

# pairwise dist between 19:21,23:25 and 1:3, 18 have the 706.85
# the indexes are 1:3, 18:24, because 22 doesn't exist in FID, but the row count obviously is continuous
annot_cvioln59[c(1:3,18:24),]
annot_cvioln59[annot_cvioln59$Geographic_subpopulations == 'pop1',]
annot_cvioln59[annot_cvioln59$Geographic_subpopulations == 'pop3',]
