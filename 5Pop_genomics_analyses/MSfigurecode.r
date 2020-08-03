# title: "Hbird pop gen & cline plots
# author: "Marisa Lim"
# date: "August 2, 2020"
# output: html_document

# Load libraries:
# for sample distribution figures
library(raster)
library(sf)
library(maptools)
library(RColorBrewer)
library(scales) #for alpha = transparent color points

# for PCA figures
library(ggplot2)
library(ggrepel)
library(gridExtra)
library(cowplot)

# for admixture figures
# library(raster) # already loaded for sample distribution maps
library(mapplots)
require(marmap)

# for include_graphics()
library(knitr)
# for reading images and plotting
library(grid)
library(magick)

# for IBD figures
# library(raster) # already loaded
library(gdistance)
library(Imap)
library(ade4)
# library(ggplot2) # already loaded
# library(cowplot) # already loaded

# ----- Main figures -----

# ------ Sample distribution maps ------
# read in data
# cviol_sdm <- raster('Coeligena_violifer.tif')
# ccoru_sdm <- raster('Colibri_coruscans.tif')

peru_waterways <- st_read('Hidrografia/Hidrografia/hidro1.shp')
peru_elevation <- raster('Peru_elevation_fromFanny/HDR.ADF')

annotpath <- '2nd pass ANGSD - remove outliers/NgsTools PCA - v2/'
annot_cvioln59 <- read.table(paste(annotpath,'cviol_n59.clst', sep=''), sep="\t", header=T)
annot_ccorun97 <- read.table(paste(annotpath,'ccoru_n97.clst', sep=''), sep='\t', header=T)

# species range maps
cviol_dichroura_range <- readShapePoly('Cviolifer_subsp_distributions/species_22726758_Cvdichroura/species_22726758/species_22726758.shp')
cviol_albicaudata_range <- readShapePoly('Cviolifer_subsp_distributions/species_22726764_Cvalbicaudata/species_22726764/species_22726764.shp')
cviol_osculans_range <- readShapePoly('Cviolifer_subsp_distributions/species_22726770_Cvosculans/species_22726770/species_22726770.shp')
ccoru_range <- readShapePoly('species_22687114_Ccoruscans/species_22687114/species_22687114.shp')

# ------ make masked rasters to Peru ------
# get peru country outline
peru <- getData('GADM', country='Peru', level=0) #country outline only (no depts, counties, etc.)
peru_adm <- getData('GADM', country='Peru', level=1)

# Run this to mask raster and save new raster - run once!
#ccoru_masked <- mask(ccoru_sdm, peru)
#cviol_masked <- mask(cviol_sdm, peru)

#writeRaster(ccoru_masked, 'ccoru_Perumasked.tif')
#writeRaster(cviol_masked, 'cviol_perumasked.tif')

crs(cviol_dichroura_range) <- crs(peru)
crs(cviol_albicaudata_range) <- crs(peru)
crs(cviol_osculans_range) <- crs(peru)
crs(ccoru_range) <- crs(peru)

cviol_d_range_cropped <- crop(cviol_dichroura_range, peru)
cviol_a_range_cropped <- crop(cviol_albicaudata_range, peru)
cviol_o_range_cropped <- crop(cviol_osculans_range, peru)
ccoru_range_cropped <- crop(ccoru_range, peru)

# ------ plot samples on masked raster maps -------
ccoru_masked_raster <- raster('ccoru_Perumasked.tif')
cviol_masked_raster <- raster('cviol_perumasked.tif')

#set up sample colors by geo subpop
cviolcols <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#c51b7d")
ccorucols <- c("#a6cee3", "#fdbf6f", "#ff7f00", "#cab2d6", "#33a02c", "#c51b7d", "#6a3d9a")

myplotcolors_cviol <- rep(cviolcols[1], length(annot_cvioln59[,1]))
myplotcolors_cviol[annot_cvioln59$Geographic_subpopulations == 'pop2'] <- cviolcols[2]
myplotcolors_cviol[annot_cvioln59$Geographic_subpopulations == 'pop3'] <- cviolcols[3]
myplotcolors_cviol[annot_cvioln59$Geographic_subpopulations == 'pop4'] <- cviolcols[4]
myplotcolors_cviol[annot_cvioln59$Geographic_subpopulations == 'pop5'] <- cviolcols[5]
myplotcolors_cviol[annot_cvioln59$Geographic_subpopulations == 'pop6'] <- cviolcols[6]
annot_cvioln59$myplotcolors <- myplotcolors_cviol

myplotcolors_ccoru <- rep(ccorucols[1], length(annot_ccorun97[,1]))
myplotcolors_ccoru[annot_ccorun97$Geographic_subpopulations == 'pop2'] <- ccorucols[2]
myplotcolors_ccoru[annot_ccorun97$Geographic_subpopulations == 'pop3'] <- ccorucols[3]
myplotcolors_ccoru[annot_ccorun97$Geographic_subpopulations == 'pop4'] <- ccorucols[4]
myplotcolors_ccoru[annot_ccorun97$Geographic_subpopulations == 'pop5'] <- ccorucols[5]
myplotcolors_ccoru[annot_ccorun97$Geographic_subpopulations == 'pop6'] <- ccorucols[6]
myplotcolors_ccoru[annot_ccorun97$Geographic_subpopulations == 'pop7'] <- ccorucols[7]
annot_ccorun97$myplotcolors <- myplotcolors_ccoru

# ------ Plots with elevation & elevation bins ------
# plot elevation by sorted individuals
annot_cvioln59$FID <- factor(annot_cvioln59$FID,
                             levels=annot_cvioln59$FID[order(annot_cvioln59$Elevation)])
annot_ccorun97$FID <- factor(annot_ccorun97$FID,
                             levels=annot_ccorun97$FID[order(annot_ccorun97$Elevation)])

elevcolors <- c('#8c510a', '#d8b365', '#f6e8c3', '#c7eae5', '#5ab4ac', '#01665e')

# use this if want to plot as shapes too: pch=as.integer(annot_ccorun97$myplotshapes)
par(mfrow=c(1,2), oma=c(3,2,2,3))
jpeg('Cviol_sample-elev_greysmap_21june.jpg', height=6, width=6, units='in', res=600)
plot(peru_elevation, xlab='Longitude (degrees)', ylab='Latitude (degrees)',
     col=brewer.pal(9, 'Greys'), main='Cviol')
plot(peru, add=T)
plot(peru_waterways$geometry, add=T, col='turquoise')
points(annot_cvioln59[,c(7,6)], pch=21,
       bg=alpha(annot_cvioln59$myplotcolors, 0.7), cex=2)
dev.off()

jpeg('Ccoru_sample-elev_greysmap_21june.jpg', height=6, width=6, units='in', res=600)
plot(peru_elevation, xlab='Longitude (degrees)', ylab='Latitude (degrees)',
     col=brewer.pal(9,'Greys'), main='Ccoru')
plot(peru, add=T)
plot(peru_waterways$geometry, add=T, col='turquoise')
points(annot_ccorun97[,c(8,7)], pch=21,
       bg=alpha(annot_ccorun97$myplotcolors, 0.7), cex=2)
dev.off()
par(mfrow=c(1,1))

ggplot(annot_cvioln59, aes(x=FID, y=Elevation)) +
  geom_bar(stat='identity', aes(fill=Geographic_subpopulations)) +
  scale_fill_manual(values=cviolcols) +
  xlab('Sample') +
  ylab('Elevation (m)') +
  coord_cartesian(ylim=c(min(annot_cvioln59$Elevation), max(annot_cvioln59$Elevation))) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  labs(fill='Geopop')
ggsave('Fig1a_cviolelevbins_geopopcolor.jpg', height=5, width=7, units='in', dpi=500)

ggplot(annot_ccorun97, aes(x=FID, y=Elevation)) +
  geom_bar(stat='identity', aes(fill=Geographic_subpopulations)) +
  scale_fill_manual(values=ccorucols) +
  xlab('Sample') +
  ylab('Elevation (m)') +
  coord_cartesian(ylim=c(min(annot_ccorun97$Elevation), max(annot_ccorun97$Elevation))) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  labs(fill='Geopop')
ggsave('Fig1b_ccoruelevbins_geopopcolor.jpg', height=5, width=7, units='in', dpi=500)

# ----- EEMS maps -----

# Coeligena violifer EEMS map
cvioleems <- image_read(paste(eems_plots, '/cviol_eems_OUT_rejratetesting3-mrates01.png', sep=''))

# Colibri coruscans EEMS map
ccorueems <- image_read(paste(eems_plots, '/ccoru_eems_OUT_rejratetesting3-mrates01.png', sep=''))


# ----- Admixture plots -----
# input file paths
FULLpath <- '/2nd pass ANGSD - remove outliers/NgsAdmix - v2/NGSadmix_FULLDAT/'
#for attribute file inputs
pass2path_pca <- '2nd pass ANGSD - remove outliers/ngsTools PCA - v2/'

# output file path
popgenoutpath = './'

admix_ccoruK1_97 <- t(as.matrix(read.table(paste(FULLpath,'Ccoru_admix_n97FULL/ccoru97FULLDAT_K1-1.qopt',sep=''))))

admix_ccoruK2_97 <- t(as.matrix(read.table(paste(FULLpath,'Ccoru_admix_n97FULL/ccoru97FULLDAT_K2-3.qopt',sep=''))))
admix_cviolK2 <- t(as.matrix(read.table(paste(FULLpath,'Cviol_admix_n59FULL/cviol59FULLDAT_K2-7.qopt',sep=''))))

admix_ccoruK3_97 <- t(as.matrix(read.table(paste(FULLpath,'Ccoru_admix_n97FULL/ccoru97FULLDAT_K3-7.qopt',sep=''))))
admix_cviolK3 <- t(as.matrix(read.table(paste(FULLpath,'Cviol_admix_n59FULL/cviol59FULLDAT_K3-4.qopt',sep=''))))

pop_cviol <-read.table(paste(pass2path_pca, "cviol_n59.clst",sep=''),header=T)
pop_ccoru_97 <- read.table(paste(pass2path_pca, 'ccoru_n97.clst', sep=''), header=T)

plotadmix_geoorder_K123 <- function(plottype, annot, mysp, myadmix, outpath){
  # order by geographic subpops
  myadmix<-myadmix[,order(annot$Geographic_subpopulations)]

  annot2<-annot[order(annot$Geographic_subpopulations),]

  if(plottype=='K1'){
    # jpeg(paste(outpath, mysp, '_K1barplot.jpg',sep=''), height=3, width=12, units='in', res=600)
    barplot(myadmix,col=c('grey13'),space=0.05, width=0.5,
            ylab=paste(mysp, " admixture K1", sep=''), names.arg=annot2$Geographic_subpopulations, las=2)
    # dev.off()
  }
  else if(plottype=='K2'){
    # jpeg(paste(outpath, mysp, '_K2barplot.jpg',sep=''), height=3, width=12, units='in', res=600)
    barplot(myadmix,col=c('grey13', 'white'),space=0.05, width=0.5,
            ylab=paste(mysp, " admixture K2", sep=''), names.arg=annot2$Geographic_subpopulations, las=2)
    # dev.off()
  }
  else if(plottype=='K3'){
    # jpeg(paste(outpath, mysp, '_K3barplot.jpg',sep=''), height=3, width=12, units='in', res=600)
    barplot(myadmix,col=c('grey13', 'grey', 'white'),space=0.05, width=0.5,
            ylab=paste(mysp, " admixture K3", sep=''), names.arg=annot2$Geographic_subpopulations, las=2)
    # dev.off()
  }
}

# ------ Coeligena violifer admixture plots K=2 & 3 ------
plotadmix_geoorder_K123('K2', pop_cviol, 'Cviol', admix_cviolK2, popgenoutpath)
plotadmix_geoorder_K123('K3', pop_cviol, 'Cviol', admix_cviolK3, popgenoutpath)

# ------ Colibri coruscans admixture plots K=1, 2, & 3 ------
plotadmix_geoorder_K123('K1', pop_ccoru_97, 'Ccoru', admix_ccoruK1_97, popgenoutpath)
plotadmix_geoorder_K123('K2', pop_ccoru_97, 'Ccoru', admix_ccoruK2_97, popgenoutpath)
plotadmix_geoorder_K123('K3', pop_ccoru_97, 'Ccoru', admix_ccoruK3_97, popgenoutpath)

# ----- PCA plots -----
# Annotation file is in plink cluster format
# input file path
pass2path <- '/2nd pass ANGSD - remove outliers/NgsTools PCA - v2/'

# output file path
outpathpass2 <- './'

covar_cvioln59 <- read.table(paste(pass2path, 'cviol_forPCAn59_out.covar', sep=''), stringsAsFact=F)
covar_ccorun97 <- read.table(paste(pass2path, 'ccoru_forPCAn97_out.covar', sep=''), stringsAsFact=F)

annot_cvioln59 <- read.table(paste(pass2path,'cviol_n59.clst', sep=''), sep="\t", header=T)
annot_ccorun97 <- read.table(paste(pass2path,'ccoru_n97.clst', sep=''), sep='\t', header=T)

# plot all species together:
# Eigenvalues
eig_cviol <- eigen(covar_cvioln59, symm=TRUE);
eig_cviol$val <- eig_cviol$val/sum(eig_cviol$val);
cat(signif(eig_cviol$val, digits=3)*100,"\n");

eig_ccoru <- eigen(covar_ccorun97, symm=TRUE)
eig_ccoru$val <- eig_ccoru$val/sum(eig_ccoru$val)
cat(signif(eig_ccoru$val, digits=3)*100,'\n')

# distribution of PCs
barplot((eig_cviol$val)*100)
barplot((eig_ccoru$val)*100)

# Plot
PC_cviol <- as.data.frame(eig_cviol$vectors)
colnames(PC_cviol) <- gsub("V", "PC", colnames(PC_cviol))
# PC_cviol$Pop <- factor(annot_cvioln59$CLUSTER)
PC_cviol$GeoPop <- factor(annot_cvioln59$Geographic_subpopulations)
PC_cviol$Specimen <- factor(annot_cvioln59$Specimen)
PC_cviol$Lat <- factor(annot_cvioln59$DEC_LAT)
PC_cviol$Long <- factor(annot_cvioln59$DEC_LONG)
PC_cviol$elevbin <- factor(annot_cvioln59$Elevation_bins_400m)

PC_ccoru <- as.data.frame(eig_ccoru$vectors)
colnames(PC_ccoru) <- gsub("V", "PC", colnames(PC_ccoru))
# PC_ccoru$Pop <- factor(annot_ccorun97$CLUSTER)
PC_ccoru$GeoPop <- factor(annot_ccorun97$Geographic_subpopulations)
PC_ccoru$Specimen <- factor(annot_ccorun97$Specimen)
PC_ccoru$Lat <- factor(annot_ccorun97$DEC_LAT)
PC_ccoru$Long <- factor(annot_ccorun97$DEC_LONG)
PC_ccoru$season <- factor(annot_ccorun97$season)
PC_ccoru$sex <- factor(annot_ccorun97$Gender)
PC_ccoru$elevbin <- factor(annot_ccorun97$Elevation_bin_400m)

write.csv(PC_cviol, 'PCA_plotdata_cviol.csv', row.names=F, quote=F)
write.csv(PC_ccoru, 'PCA_plotdata_ccoru.csv', row.names=F, quote=F)

write.csv(eig_cviol$val, 'PCA_eigenvals_cviol.csv', row.names=F, quote=F)
write.csv(eig_ccoru$val, 'PCA_eigenvals_ccoru.csv', row.names=F, quote=F)

# Parse components to analyze
# 1-2 means PC1 and PC2
comp <- as.numeric(strsplit('1-2', "-", fixed=TRUE)[[1]])
x_axis = paste("PC",comp[1],sep="")
y_axis = paste("PC",comp[2],sep="")

cviolPCAplot <- ggplot() + geom_point(data=PC_cviol, aes_string(x=x_axis, y=y_axis, fill='GeoPop'), size=5, alpha=0.7, pch=21) +
  ggtitle(expression(paste(italic('Coeligena violifer'), ' (n=59)', sep=''))) +
  xlab(paste("PC",comp[1]," (",signif(eig_cviol$val[comp[1]], digits=3)*100,"%)",sep='')) +
  ylab(paste("PC",comp[2]," (",signif(eig_cviol$val[comp[2]], digits=3)*100,"%)",sep='')) +
  scale_fill_manual(values=cviolcols) +
  theme_cowplot()
ggsave('cviolPCA_21june.jpeg', height=5, width=5, units='in', dpi=600)

ccoruPCAplot <- ggplot() + geom_point(data=PC_ccoru, aes_string(x=x_axis, y=y_axis, fill='GeoPop'), size=5, alpha=0.7, pch=21) +
  ggtitle(expression(paste(italic('Colibri coruscans'), ' (n=97)', sep=''))) +
  xlab(paste("PC",comp[1]," (",signif(eig_ccoru$val[comp[1]], digits=3)*100,"%)",sep='')) +
  ylab(paste("PC",comp[2]," (",signif(eig_ccoru$val[comp[2]], digits=3)*100,"%)",sep='')) +
  scale_fill_manual(values=ccorucols)+
  theme_cowplot()
ggsave('ccoruPCA_21june.jpeg', height=5, width=5, units='in', dpi=600)

# jpeg('cviol_PCA.jpeg', height=5, width=5, units='in', res=600)
plot_grid(cviolPCAplot, ccoruPCAplot)
# dev.off()

# ----- SNP figures -----

# ------ Figure 3: ------
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
dcases<-countsToCases(dcounts)
# x-axis = candidate vs. control and colors = cline vs. no cline
dcounts$factorC2 <-with(dcounts, interaction(SNP_type, species))
dc6<-as.data.frame(dcounts %>% group_by(factorC2) %>% mutate(fraction=gene_freq/sum(gene_freq)))
dc6$species<- factor(dc6$species, levels = c('Coe. violifer', 'Col. coruscans'))
ggplot() + geom_bar(aes(y = fraction, x = SNP_type, fill = LFMM), col='black', data = dc6, stat="identity")+
  labs(x="Candidate vs. control", y="Proportion") +
  guides(fill=guide_legend(title = "Cline type"))+
  scale_fill_manual(values=c('grey', 'white')) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size=15))+
  facet_grid(~ species, margins=F)
ggsave('Fig3_genefacetCviolCcoru_cand-control_cline-no.jpg', h=5, w=6, dpi=600)

# ----- Supplementary figures -----

# ------ Figure S1: Plots with SDM & range map ------
```{r, out.width='90%',, cache=TRUE, echo=FALSE, message=FALSE, results='hide'}
par(mfrow=c(1,2), oma=c(3,3,2,2))
# 'Greys' doesn't have the best contrast, but I don't want extra colors...
# jpeg('Cviol_sample-SDM_map.jpg', height=6, width=6, units='in', res=600)
plot(cviol_masked_raster, xlab='Longitude (degrees)', ylab='Latitude (degrees)',
     col=brewer.pal(9,'Greys'), main='Cviol')
plot(peru, add=T)
plot(cviol_d_range_cropped, add=T, col=alpha('gold', 0.4), lwd=0.01)
plot(cviol_a_range_cropped, add=T, col=alpha('#02818a', 0.4), lwd=0.01)
plot(cviol_o_range_cropped, add=T, col=alpha('black', 0.4), lwd=0.01)
points(annot_cvioln59[,c(7,6)], pch=21,
       bg=alpha(annot_cvioln59$myplotcolors, 0.7), cex=1.3)
# dev.off()

# jpeg('Ccoru_sample-SDM_map.jpg', height=6, width=6, units='in', res=600)
plot(ccoru_masked_raster, xlab='Longitude (degrees)', ylab='Latitude (degrees)',
     col=brewer.pal(9,'Greys'), main='Ccoru')
plot(peru, add=T)
plot(ccoru_range_cropped, add=T, col=alpha('blue', 0.4), lwd=0.01)
points(annot_ccorun97[,c(8,7)], pch=21,
       bg=alpha(annot_ccorun97$myplotcolors, 0.7), cex=1.3)
# dev.off()
par(mfrow=c(1,1))

# ----- Isolation by distance ------
# read in cost distance matrices
costdist_ccoru <- as.matrix(read.table('ccoru_costdist_km.txt'))
costdist_cviol <- as.matrix(read.table('cviol_costdist_km.txt'))

# read in genetic distances from NGSdist
# using the NGSdist matrixes calculated for EEMS, deleted blank 1st line and 2nd line that had sample size
# the first column is Ind#
IBDpath <- '2nd pass ANGSD - remove outliers/Isolation_by_Dist/'
gendist_cvioln59 <- read.table(paste(IBDpath,'cviol59_forIBD.txt', sep=''), sep='\t')[-1]
gendist_ccorun97 <- read.table(paste(IBDpath,'ccoru97_forIBD.txt', sep=''), sep='\t')[-1]

cviol_pts <- read.table('/2nd pass ANGSD - remove outliers/ngsTools PCA - v2/cviol_n59.clst', header=TRUE)[,c('Specimen', 'DEC_LAT', 'DEC_LONG')]
ccoru_pts <- read.table('/2nd pass ANGSD - remove outliers/ngsTools PCA - v2/ccoru_n97.clst', header=TRUE)[,c('Specimen', 'DEC_LAT', 'DEC_LONG')]

# ------ Figure S2: Plot genetic distance vs. cost distance & genetic distance vs. eucldean distance ------
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


cviol_geo <- ggplot() + geom_point(aes(x=geodist_cviol[lower.tri(geodist_cviol)],
                                       y=gendist_cvioln59[lower.tri(gendist_cvioln59)]), size=1.5, alpha=0.5) +
  xlab('Geodesic distance (km)') + ylab('Genetic distance')+
  annotate(geom='text', x=1000, y=0.09, label='r = 0.75 \n    p = 0.0001', size=4)
cviol_cost <- ggplot() + geom_point(aes(x=costdist_cviol[lower.tri(costdist_cviol)],
                                        y=gendist_cvioln59[lower.tri(gendist_cvioln59)]), size=1.5, alpha=0.5) +
  xlab('Cost distance (km)') + ylab('Genetic distance')+
  annotate(geom='text', x=350, y=0.09, label='r = 0.70 \n    p = 0.0001', size=4)

ccoru_geo <- ggplot() + geom_point(aes(x=geodist_ccoru[lower.tri(geodist_ccoru)],
                                       y=gendist_ccorun97[lower.tri(gendist_ccorun97)]), size=1.5, alpha=0.5) +
  xlab('Geodesic distance (km)') + ylab('Genetic distance') +
  annotate(geom='text', x=1400, y=0.087, label='r = 0.05 \np = 0.15', size=4)
ccoru_cost <- ggplot() + geom_point(aes(x=costdist_ccoru[lower.tri(costdist_ccoru)],
                                        y=gendist_ccorun97[lower.tri(gendist_ccorun97)]), size=1.5, alpha=0.5) +
  xlab('Cost distance (km)') + ylab('Genetic distance')+
  annotate(geom='text', x=500, y=0.087, label='r = 0.03 \np = 0.23', size=4)

jpeg('FigS2_gendistXcostdistXeuclid_4plots_2sp.jpg', height=8, width=10, units='in', res=600)
plot_grid(cviol_geo, cviol_cost,
          ccoru_geo, ccoru_cost,
          labels=c('a', '' , 'b', ' '),
          nrow=2, ncol=2)
dev.off()

# ----- Figure S5-8: GO biological process and Panther pathway bar plots -----
# these are after fixing gene lists Nov 2019
# Previously ID'd candidate genes
cviol_Wilcoxcand_pathway <- read.table('cviol_wilcox_pathway.txt', sep='\t', header=F)
cviol_Wilcoxcand_BP <- read.table('cviol_wilcox_BP.txt', sep='\t', header=F)
ccoru_Wilcoxcand_pathway <- read.table('ccoru_wilcox_pathway.txt', sep='\t', header=F)
ccoru_Wilcoxcand_BP<- read.table('ccoru_wilcox_BP.txt', sep='\t', header=F)

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

save_plot('17nov19_CviolCcoru_Wilcoxcand_bpbar.jpg', bpplots2, base_height=10, base_width=10)
save_plot('17nov19_CviolCcoru_Wilcoxcand_pathwaysbar.jpg', pathwayplots2, base_height=15, base_width=25)

# Novel candidate genes
cviol_novelcand_pathway <- read.table('cviol_novel_pathway.txt', sep='\t', header=F)
cviol_novelcand_BP <- read.table('cviol_novel_BP.txt', sep='\t', header=F)
ccoru_novelcand_pathway <- read.table('ccoru_novel_pathway.txt', sep='\t', header=F)
ccoru_novelcand_BP<- read.table('ccoru_novel_BP.txt', sep='\t', header=F)

levels(cviol_novelcand_BP$V2) <- gsub(' \\(', '\n(', levels(cviol_novelcand_BP$V2))
levels(ccoru_novelcand_BP$V2) <- gsub(' \\(', '\n(', levels(ccoru_novelcand_BP$V2))

cviol_bp_plot <- ggplot(cviol_novelcand_BP) +
  geom_bar(stat='identity', aes(y=V3, x=V2)) +
  geom_text(aes(x=V2, y=V3, label=V3), hjust=0, col='black', size=4.5) +
  coord_flip() + xlab('') + ylab('Number of genes')
cviol_pathway_plot <- ggplot(cviol_novelcand_pathway) +
  geom_bar(stat='identity', aes(y=V3, x=V2)) +
  geom_text(aes(x=V2, y=V3, label=V3), hjust=0, col='black', size=5) +
  coord_flip() + xlab('') + ylab('Number of genes')
ccoru_bp_plot <- ggplot(ccoru_novelcand_BP) +
  geom_bar(stat='identity', aes(y=V3, x=V2)) +
  geom_text(aes(x=V2, y=V3, label=V3), hjust=0, col='black', size=4.5) +
  coord_flip() + xlab('') + ylab('Number of genes')
ccoru_pathway_plot <- ggplot(ccoru_novelcand_pathway) +
  geom_bar(stat='identity', aes(y=V3, x=V2)) +
  geom_text(aes(x=V2, y=V3, label=V3), hjust=0, col='black', size=5) +
  coord_flip() + xlab('') + ylab('Number of genes')

bpplots <- plot_grid(cviol_bp_plot, ccoru_bp_plot, ncol=1, align='h', labels=c('a)', 'b)'))
pathwayplots <- plot_grid(cviol_pathway_plot, ccoru_pathway_plot, ncol=2, align='h', labels=c('a)', 'b)'))

save_plot('17nov19_CviolCcoru_novelcand_bpbar.jpg', bpplots, base_height=10, base_width=10)
save_plot('17nov19_CviolCcoru_novelcand_pathwaysbar.jpg', pathwayplots, base_height=15, base_width=25)