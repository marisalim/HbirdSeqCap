# Fst calculation script descriptions

This is a pair-wise calculation between each population, here defined by locality.

Required software:
- angsd
  - [Fst calc](http://popgen.dk/angsd/index.php/Fst)

Follow steps under [`Two Populations real data`](http://popgen.dk/angsd/index.php/Fst#Two_Populations_real_data) header:

1. ANGSD
- Fst global estimate for 2 populations, do not do the sliding window part
- run ANGSD command for [PCA](../A_PCA/readme.md) with `minInd` = 70% of samples and `minIndDepth=3`
- do not need these flags: `doGeno`, `doPost`, `doMaf`, `doMajorMinor`, `skipTriallelic`, `-fold` (don't need genotype calls for Fst)

Wrapper scripts:
```
sbatch fst_ccorun97_wrap.sh
sbatch fst_cvioln59_wrap.sh
```

2. realsFS
- in the output, you want 2nd number = `Fst.Weight` value for Fst

Wrapper scripts:
```
calcfst_ccorun97_wrap.sh
calcfst_cvioln59_wrap.sh
```

Copy/paste the `Fst.Weight` value for each pair-wise combination of subpopulations to create Fst matrix.
