# Isolation by Distance script description

Required software:
- ANGSD NGSdist
- R

## Calculate pairwise differences between each individual

Use NGSdist matrix that was the input for [EEMs analysis](../H_EEMS/readme.md).

## Geographic (geodesic) distance

Geodesic distance = shortest possible line between two points on a sphere or other curved surface) between individuals
	- note: geodesic is not the same as Euclidean distance because Euclidean dist is straight lines, not curved

Requires:
- genetic covariance matrix
- lat/long coordinates

Calculations done in R: `IBDcalc_plot.r`

## Cost distance

Use `gdistance` R package to calculate least cost distances using species distribution models from Javier Fajardo's [paper](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0114367).

Requires:
- genetic covariance matrix
- lat/long coordinates
- SDM files (e.g., tif file)

Calculations done in R: `calc_cost_dist.r`

Paper figure code in [`MSfigurecode.rmd`](../MSfigurecode.rmd)
