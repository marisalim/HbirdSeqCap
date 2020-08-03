# EEMS script descriptions

Required software
- [EEMS](https://github.com/dipetkov/eems)
  - [documentation](https://github.com/dipetkov/eems/blob/master/Documentation/EEMS-doc.pdf)
- angsd [NGSdist](https://github.com/fgvieira/ngsDist)

The input data:
- `datapath.diffs` = matrix of average pairwise differences (indiv x indiv matrix of genetic dissimilarity)
- `datapath.coord` = list of sample coordinates (2 columns: long, lat)
- `datapath.outer` = list of habitat boundary points (draw a polygon around geographical area that encompasses all samples and choose several points (long, lat) to put into this list)
-

## `datapath.diffs`: Generate genetic dissimilarity matrix with NGSdist

1. generate genotype probabilities

- use same ANGSD as was used to call genotypes for LFMM, but use genotype probabilities `-doGeno 8` instead of 2
- does not need the MAF filter since this is an analysis across all markers and it not a test for natural selection
- may need to add NGSdist to bashrc path to run

Wrapper scripts:
```
ccoru_genoforngsdist_wrapper.sh
cviol_genoforngsdist_wrapper.sh
```

2. Remove the first column so just the dissimilarity matrix is left.
```
cut -f1 --complement ccoru97_foreems_out | tail -n+3 > ccoru97.diffs
cut -f1 --complement cviol59_foreems_out | tail -n+3 > cviol59.diffs
```

## `datapath.coord`

- labeled as <species><sample number>.coord e.g., cviol59.coord
- make sure long, lat coordinates are in decimal degrees

## `datapath.outer`
- did this in QGIS
  - draw polygon, extract feature data to get lat/longs
  - save in csv and plot on map again to make sure they are correct
  - EEMS wants the coordinates going counterclockwise. QGIS defaults to list the first to the last vertex (a closed ring), so you need to reverse the order of points. The first vertex should also be the last vertex so the outline is a **closed ring**

## Run EEMS

Set up EEMS parameters and paths:
```
params-ccoru_chain1.ini
params-cviol_chain1.ini
```

EEMS command:
```
# these were run on HPC, basic command is
runeems_snps --params params-ccoru_chain1.ini --seed 120
runeems_snps --params params-cviol_chain1.ini --seed 120
```

Note: settings were optimized following recommendation by EEMS authors to aim for acceptance proportion rate of 20-30% though 10-40% is ok in practice (find the rates in EEMS output files).

## Plot results

- plot maps with `make_eems_plots.r`
