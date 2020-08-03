# cost distance analysis for pop gen samples

library(raster)
library(gdistance)
library(spatialEco)
library(Imap)
library(ade4)
library(ggplot2)
library(cowplot)

# the transition() function assumes that the highest raster cell value is a barrier
# the low raster values are higher conductance
# so, for the raster layer values for high/low suitability,
# need to invert the raster layer so that high suitability = low raster values
# and low suitability = high raster values
# you can test this by looking at the shortestPath() in inverted vs. not inverted layers

# from https://dyerlab.github.io/applied_population_genetics/ecological-distance.html
# One way to estimate distances in R, is through the use of the gdistance library. In this approach, we define a transition object based upon:
# 1. The cost distance raster. By default, the transition() function works on conductance, which is the inverse of resistance. In these examples, we have used a single raster, though use of RasterBrick objects is just as appropriate.
# 2. The function by which we estimate the pixel-to-pixel distances.
# 3. The neighborhood around each pixel that we look at during each step. Options for this include:
# - A von Neumann neighborhood (directions=4) consists of the pixels immediately above, below, and on both sides of the target pixel.
# - A Moore's neighborhood (directions=8) consisting of all pixels surrounding the target pixel (e.g., even the diagonal ones)
# - A Knight & One-Cell Queen move (directions=16), which adds the next layer of cells outside of Moore's neighborhood.
#
# Once estimated, the transition object must be corrected for if you are either using a
# large extent based upon Latitude and Longitude datum (e.g., they are not unit-wise consistent in area),
# or you have used a direction option other than the von Neuman Neighborhood.

# Calculate cost distance
# these look at pixels with Moore's neighborhood definition (directions=8)
# ----- example -----
## Create cost surface where "land" exists in the middle
cost <- raster(nrow=100, ncol=100,
               xmn=0, xmx=100, ymn=0, ymx=100, crs="+proj=utm")
cost[] <- 1
cost[cellFromRowColCombine(cost, 50:55,20:80)] <- 10000

## Produce transition matrices, and correct because 8 directions
trCost <- transition(1/cost, mean, directions=8)
trCost <- geoCorrection(trCost, type="c")

## Create three points (representing three points in time series)
pts <- cbind(x=c(20, 60, 40), y=c(80, 60, 20))

## inverted raster version
r.cost <- raster.invert(cost)
trCost_invert <- transition(1/r.cost, mean, directions=8)
trCost_invert <- geoCorrection(trCost_invert, type="c")

## Display results
par(mfrow=c(1,2))
plot(cost)
plot(SpatialPoints(pts), add=TRUE, pch=20, col="red")
text(pts[,1]+2, pts[,2]+2, 1:nrow(pts))
plot(shortestPath(trCost, pts[1,], pts[2,], output="SpatialLines"), add=TRUE)
plot(shortestPath(trCost, pts[2,], pts[3,], output="SpatialLines"), add=TRUE)
plot(shortestPath(trCost, pts[1,], pts[3,], output='SpatialLines'), add=TRUE)

plot(r.cost)
plot(SpatialPoints(pts), add=TRUE, pch=20, col="red")
text(pts[,1]+2, pts[,2]+2, 1:nrow(pts))
plot(shortestPath(trCost_invert, pts[1,], pts[2,], output="SpatialLines"), add=TRUE)
plot(shortestPath(trCost_invert, pts[2,], pts[3,], output="SpatialLines"), add=TRUE)
plot(shortestPath(trCost_invert, pts[1,], pts[3,], output='SpatialLines'), add=TRUE)

# units are meters, divide by 1000 to get km
costDistance(trCost_invert, pts)/1000

# ----- Cviol -----
cviol_sdm <- raster('Coeligena_violifer.tif')
cviol_sdm_invert <- raster.invert(cviol_sdm)
cviol_pts <- read.table('../SeqCapBioinformatics/ANGSD analyses/2nd pass ANGSD - remove outliers/ngsTools PCA - v2/cviol_n59.clst', header=TRUE)[,c('Specimen', 'DEC_LAT', 'DEC_LONG')]
cviol_pts_matrix <- as.matrix(cviol_pts[,c(3,2)])

trCost_invert_cviol <- transition(1/cviol_sdm_invert, mean, directions=8)
trCost_invert_cviol <- geoCorrection(trCost_invert_cviol, type="c")

cviol_costdist <- costDistance(trCost_invert_cviol, cviol_pts_matrix)/1000
write.table(as.matrix(cviol_costdist), 'cviol_costdist_km.txt', quote=FALSE)

# ----- Ccoru -----
ccoru_sdm <- raster('Colibri_coruscans.tif')
ccoru_sdm_invert <- raster.invert(ccoru_sdm)
ccoru_pts <- read.table('../SeqCapBioinformatics/ANGSD analyses/2nd pass ANGSD - remove outliers/ngsTools PCA - v2/ccoru_n97.clst', header=TRUE)[,c('Specimen', 'DEC_LAT', 'DEC_LONG')]
ccoru_pts_matrix <- as.matrix(ccoru_pts[,c(3,2)])

trCost_invert_ccoru <- transition(1/ccoru_sdm_invert, mean, directions=8)
trCost_invert_ccoru <- geoCorrection(trCost_invert_ccoru, type="c")

ccoru_costdist <- costDistance(trCost_invert_ccoru, ccoru_pts_matrix)/1000
write.table(as.matrix(ccoru_costdist), 'ccoru_costdist_km.txt', quote=FALSE)


# ----- plot genetic distance vs. cost distance -------

# read in cost distance matrices
costdist_ccoru <- as.matrix(read.table('ccoru_costdist_km.txt'))
costdist_cviol <- as.matrix(read.table('cviol_costdist_km.txt'))

# read in genetic distances from NGSdist
# using the NGSdist matrixes calculated for EEMS, deleted blank 1st line and 2nd line that had sample size
# the first column is Ind#
IBDpath <- '../SeqCapBioinformatics/ANGSD analyses/2nd pass ANGSD - remove outliers/Isolation_by_Dist/'
gendist_cvioln59 <- read.table(paste(IBDpath,'cviol59_forIBD.txt', sep=''), sep='\t')[-1]
gendist_ccorun97 <- read.table(paste(IBDpath,'ccoru97_forIBD.txt', sep=''), sep='\t')[-1]

# plot
jitterfactor <- 100
jpeg('gendistXcostdist_2sp_jittered.jpg', height=8, width=8, units='in', res=600)
par(mfrow=c(2,2))
plot(gendist_ccorun97[lower.tri(gendist_ccorun97)]~
       jitter(costdist_ccoru[lower.tri(costdist_ccoru)], factor=jitterfactor), pch=20,
     ylab='Genetic distance', xlab='Cost distance (km)', main='C. coruscans n=97')
plot(gendist_cvioln59[lower.tri(gendist_cvioln59)]~
       jitter(costdist_cviol[lower.tri(costdist_cviol)], factor=jitterfactor), pch=20,
     ylab='Genetic distance', xlab='Cost distance (km)', main='C. violifer n=59')
dev.off()

# ----- plot genetic distance vs. cost distance and eucldean dist on same plot ----
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
geodist_cviol <- round(GeoDistanceInMetresMatrix(cviol_pts) / 1000, 2)
geodist_ccoru <- round(GeoDistanceInMetresMatrix(ccoru_pts) / 1000, 2)

#overlapping
jitterfactor <- 100
jpeg('gendistXcostdistXeuclid_2sp_jittered.jpg', height=8, width=8, units='in', res=600)
par(mfrow=c(2,2))
plot(gendist_cvioln59[lower.tri(gendist_cvioln59)]~
       jitter(costdist_cviol[lower.tri(costdist_cviol)], factor=jitterfactor), pch=20,
     ylab='Genetic distance', xlab='Cost distance (km)', main='C. violifer n=59',
     xlim=c(0,1300))
points(gendist_cvioln59[lower.tri(gendist_cvioln59)]~
       jitter(geodist_cviol[lower.tri(geodist_cviol)], factor=jitterfactor), pch=23, col='red')
plot(gendist_ccorun97[lower.tri(gendist_ccorun97)]~
       jitter(costdist_ccoru[lower.tri(costdist_ccoru)], factor=jitterfactor), pch=20,
     ylab='Genetic distance', xlab='Cost distance (km)', main='C. coruscans n=97',
     xlim=c(0,1600))
points(gendist_ccorun97[lower.tri(gendist_ccorun97)]~
         jitter(geodist_ccoru[lower.tri(geodist_ccoru)], factor=jitterfactor), pch=23, col='red')

dev.off()

# ----- Mantel test ------
#https://stats.idre.ucla.edu/r/faq/how-can-i-perform-a-mantel-test-in-r/
# The test consists of calculating the correlation of the entries in the matrices,
# then permuting the matrices and calculating the same test statistic under each permutation
# and comparing the original test statistic to the distribution of test statistics from the
# permutations to generate a p-value. The number of permutations defines the precision with
# which the p-value can be calculated.  The function to perform the Mantel test is mantel.rtest
# and the required arguments are the two distance matrices.
# The number of permutations can also be specified by the user, but is by default 99.

## Genetic distance vs. Geodesic distance
mantel.rtest(dist(gendist_ccorun97), dist(geodist_ccoru), nrepet=10000)
# Observation: 0.04907599
#
# Based on 10000 replicates
# Simulated p-value: 0.1506849
# Alternative hypothesis: greater
#
# Std.Obs  Expectation     Variance
# 1.0264600635 0.0007746536 0.0022142888

mantel.rtest(dist(gendist_cvioln59), dist(geodist_cviol), nrepet=10000)
# Observation: 0.7462659
#
# Based on 10000 replicates
# Simulated p-value: 9.999e-05
# Alternative hypothesis: greater
#
# Std.Obs   Expectation      Variance
# 18.9888314843 -0.0003861461  0.0015461086

## Genetic distance vs. Cost distance
mantel.rtest(dist(gendist_ccorun97), dist(costdist_ccoru), nrepet=10000)
# Observation: 0.03390982
#
# Based on 10000 replicates
# Simulated p-value: 0.230077
# Alternative hypothesis: greater
#
# Std.Obs   Expectation      Variance
# 0.6396318993 -0.0009232507  0.0029656737

mantel.rtest(dist(gendist_cvioln59), dist(costdist_cviol), nrepet=10000)
# Observation: 0.6952883
#
# Based on 10000 replicates
# Simulated p-value: 9.999e-05
# Alternative hypothesis: greater
#
# Std.Obs  Expectation     Variance
# 16.233564179  0.000612716  0.001831202
