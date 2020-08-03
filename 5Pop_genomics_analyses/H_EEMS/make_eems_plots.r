# plot EEMs results

#library(devtools)
#install_github('dipetkov/eems/plotting/rEEMSplots')
library(rEEMSplots)
library(RColorBrewer)
library(rgdal)
library(rworldmap)
library(rworldxtra)

#set Ghostscript executable
Sys.setenv(R_GSCMD="C:/Program Files (x86)/gs/gs9.23/bin/gswin32c.exe")

# --- Check the input genetic distances with heatmap ----
cviol_eems <- as.matrix(read.table('cviol_eems_in/cviol59.diffs'))
ccoru_eems <- as.matrix(read.table('ccoru_eems_in/ccoru97.diffs'))

jpeg('Cviol_pairwiseGenDist.jpg', height=5, width=5, units='in', res=500)
heatmap(cviol_eems, symm=T, main='C. violifer pairwise genetic distance')
dev.off()
jpeg('Ccoru_pairwiseGenDist.jpg', height=5, width=5, units='in', res=500)
heatmap(ccoru_eems, symm=T, main='C. coruscans pairwise genetic distance')
dev.off()
# red = 0 differences, that's the 1:1 line
# yellowish = more differences
# whitish = most different

# --- Make eems plots ----
eems_outs_path <- './26july18_acceptance_rate_testing/' #use the OUT_rejratetesting[#] mcmc and plot paths - these are based on proposal acceptance rates (see output stats in slurm files)
map_world <- getMap()
map_peru <- map_world[which(map_world@data$SOVEREIGNT=='Peru'), ]

# mcmcpath_cviol = paste(eems_outs_path, 'cviol_eems_OUT4e6deme100', sep='')
# plotpath_cviol = paste(eems_outs_path, 'cviol_eems_OUT4e6deme100', sep='')
mcmcpath_cviol = paste(eems_outs_path, 'cviol_eems_OUT_rejratetesting3', sep='')
plotpath_cviol = paste(eems_outs_path, 'cviol_eems_OUT_rejratetesting3', sep='')
coord_cviol <- read.table('./cviol_eems_in/cviol59.coord')

# mcmcpath_ccoru = paste(eems_outs_path, 'ccoru_eems_OUT4e6deme100', sep='')
# plotpath_ccoru = paste(eems_outs_path, 'ccoru_eems_OUT4e6deme100', sep='')
mcmcpath_ccoru = paste(eems_outs_path, 'ccoru_eems_OUT_rejratetesting3', sep='')
plotpath_ccoru = paste(eems_outs_path, 'ccoru_eems_OUT_rejratetesting3', sep='')
coord_ccoru <- read.table('./ccoru_eems_in/ccoru97.coord')

# --- with grids and demes -----
eems.plots(mcmcpath_cviol, plotpath_cviol, longlat=TRUE,
           m.plot.xy={points(coord_cviol, col='black', pch=20, cex=1.5)
             plot(map_peru, col=NA, add=T)},
           q.plot.xy={points(coord_cviol, col='black', pch=20, cex=1.5)
             plot(map_peru, col=NA, add=T)},
           eems.colors=brewer.pal(11,'RdYlBu'),
           add.grid=T, add.outline=T,
           col.grid = "gray90",
           lwd.grid = 2,
           col.outline = "gray90",
           lwd.outline = 2,
           add.demes = TRUE,
           col.demes = "red",
           pch.demes = 5,
           min.cex.demes = 0.5,
           max.cex.demes = 1.5)

eems.plots(mcmcpath_ccoru, plotpath_ccoru, longlat=TRUE,
           m.plot.xy={points(coord_ccoru, col='black', pch=20, cex=1.5)
             plot(map_peru, col=NA, add=T)},
           q.plot.xy={points(coord_ccoru, col='black', pch=20, cex=1.5)
             plot(map_peru, col=NA, add=T)},
           eems.colors=brewer.pal(11,'RdYlBu'),
           add.grid=T, add.outline=T,
           col.grid = "gray90",
           lwd.grid = 2,
           col.outline = "gray90",
           lwd.outline = 2,
           add.demes = TRUE,
           col.demes = "red",
           pch.demes = 5,
           min.cex.demes = 0.5,
           max.cex.demes = 1.5)

# the color scale works for m and q rate plots, but doesn't change from default for
# the posterior probability plots for some reason...

# --- no grids or demes ----
eems.plots(mcmcpath_cviol, plotpath_cviol, longlat=TRUE,
           m.plot.xy={points(coord_cviol, col='black', pch=20, cex=1.5)
             plot(map_peru, col=NA, add=T)},
           q.plot.xy={points(coord_cviol, col='black', pch=20, cex=1.5)
             plot(map_peru, col=NA, add=T)},
           eems.colors=brewer.pal(11,'RdYlBu'),
           add.outline = TRUE,
           col.outline = "gray90",
           lwd.outline = 2)

eems.plots(mcmcpath_ccoru, plotpath_ccoru, longlat=TRUE,
           m.plot.xy={points(coord_ccoru, col='black', pch=20, cex=1.5)
             plot(map_peru, col=NA, add=T)},
           q.plot.xy={points(coord_ccoru, col='black', pch=20, cex=1.5)
             plot(map_peru, col=NA, add=T)},
           eems.colors=brewer.pal(11,'RdYlBu'),
           add.outline = TRUE,
           col.outline = "gray90",
           lwd.outline = 2)
