# NgsRelate script descriptions

Required software:
- [NgsRelate](https://github.com/ANGSD/NgsRelate)

## Create input files for NgsRelate
1. Run ANGSD command used for [PCA](../A_PCA/readme.md), but add these flags
  - `snp_pval 1e-6` = 'Remove sites with a pvalue larger than'
  - `minmaf 0.05` = 'Remove sites with MAF below' (MAF = minor allele frequency)
  - `doGlf 3` = 'Output the log genotype likelihoods to a file, beagle binary format'

Wrapper scripts:
```
sbatch relate_ccoru_wrap.sh
sbatch relate_cviol_wrap.sh
```

2. Create file with allele frequencies
```
zcat cviol_relatetest_out.mafs.gz | cut -f7 | sed 1d > cviol_freq
zcat ccoru_relatetest_out.mafs.gz | cut -f7 | sed 1d > ccoru_freq
```

## Run NgsRelate
NgsRelate command:
```
ngsRelate -g <output .glf.gz file> -n 38 -f <allele freq file> -z <sample ID list> > <output file>
```

Notes on flags
- `-g` = genotype file (made above, step 1)
- `-n` = number of samples in .glf.gz file
- `-f` = allele frequencies file (made above, step 2
- `-z` = name of file with sample IDs

Wrapper script:
```
ngsrelate_wrap.sh
```

## Calculate [relatedness](http://www.popgen.dk/angsd/index.php/Relatedness)
relatedness coefficient (r) = (`k1`/2) + `k2`
relatedness (theta) = r/2

We are interested in theta. Theta and r are pair-wise calculations between each individual in the dataset.

Based on data from the Nachman lab on house mice, theta = 0.25 are full sibs, so theta < 0.25 is ok to keep. The following commands extract sample IDs ($3=individual a and $4=individual b) and calculate r and theta using the output from NgsRelate.

Repeat for all species:
```
cat cviol_outfile.res | tr -d $'\r' > cviol_outfile2.res
awk '{print $3 "\t" $4 "\t" ($7/2 + $8) "\t" ($7/2 + $8)/2}' cviol_outfile2.res > cviol_rtheta.txt
```

Filter the theta columns to check if any values are >0.25. This file has 4 columns: 1) individual a, 2) individual b, 3) r, and 4) theta.
```
awk '{if ($4 >= 0.25) print $1 "\t" $2 "\t" $3 "\t" $4}' cviol_rtheta.txt
```
