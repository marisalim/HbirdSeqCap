# ANGSD genotype calling script descriptions

Unlike genotype probabilities, genotype calling is even more sensitive. Consequently, I tested the sensitivity of several paramters used by the ANGSD program to call genotypes.

Steps:
1. run ANGSD command to call genotypes (see scripts below)
2. Check MAFs output for number of retained SNPs
```
zcat <output file name.mafs.gz> | tail -n+2 | wc -l
```
3. Check which SNPs are filtered out
	- can check the candidate genes I know I don’t want filtered out for # SNPs retained between versions (mafs file wc -l after using grep to search for specific contigs)
	- look at geno file to see what kinds of genotypes (het vs. homozygote) are filtered out

## Baseline test
These scripts contain commented out previous commands, the baseline to the final command used in analysis:

```
sbatch optgeno_ccoru_wrap.sh
sbatch optgeno_cviol_wrap.sh
```

## Final genotype call commands per species after sensitivity analyses

Final parameter choice: `SNP_pval` = 1e-3 and	`postCutOff` = 0.95; `hetBias` and `-C` not used

Wrapper scripts for final analysis in paper:
```
sbatch optgeno_cviol_wrap.sh
sbatch optgeno_ccoru_wrap.sh
```
