# ngsAdmix script descriptions

Required software:
- [NgsAdmix](http://www.popgen.dk/software/index.php/NgsAdmix)

## Create inputs for NgsAdmix

Run ANGSD command used for [PCA](../A_PCA/readme.md), but add `-doGlf 2` to save genotype likelihoods in beagle format.
  - use (70% \* samples) for `minInd` and `minIndDepth = 3`

Wrapper scripts:
```
sbatch admix_ccoru_wrap.sh
sbatch admix_cviol_wrap.sh
```

## Run NgsAdmix

Command:
```
NGSadmix -likes {input} -K {myK} -P 4 -minMaf 0.05 -o {output}
```
- `-likes` = from ANGSD output above
- `-K` = number of *a priori* populations you want to test (tested K=[1,2,3,4,5])
- `-minMaf` = minimum minor allele frequency filter

Ran each iteration (K) 10 times.

The NGSadmix command is run from a python script. Running the python scripts on the HPC with wrapper scripts.

Python scripts:
```
runadmix_ccoru.py
runadmix_cviol.py
```

Wrapper scripts:
```
sbatch runadmix_ccoru_wrap.sh
sbatch runadmix_cviol_wrap.sh
```

## Get likelihood values from each iteration per K to calculate delta K
The delta K value is how we assess which K value is the best estimate of number of populations - this is called the Evanno method.

The Evanno method was used for the first pass analysis of the data only, as it was not very informative. For the 2nd pass (final), I defined subpopulations by locality and by comparing to the other pop gen results (PCA, Fst, Dxy, EEMs).

1. Extract likihoods from log files:

Example: (done for each K per species)
```
grep best *K1*.log | awk '{print $2}' | tr -d 'like=' > K1_likelihood_sgeof.txt
grep best *K2*.log | awk '{print $2}' | tr -d 'like=' > K2_likelihood_sgeof.txt
# etc. to K5
```

The likelihood txt files are input for this R script that calculates delta K and plots the results: `evanno_method.r`

## Plot admixture results
- plotted as bar plots
