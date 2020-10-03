# ANGSD genotype calling script descriptions

Unlike genotype probabilities, genotype calling is even more sensitive. Consequently, I tested the sensitivity of several parameters used by the ANGSD program to call genotypes.

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

Full command:
```
angsd -bam <bam list file> -ref <reference fasta file> -anc <reference fasta file> -doGeno 2 -SNP_pval 1e-3 -postCutoff 0.95 -dosnpstat 1 -minInd <70% of individuals> -geno_minDepth 5 -only_proper_pairs 0 -minMapQ 10 -minQ 20 -fold 1 -GL 1 -doPost 2 -doMaf 1 -doMajorMinor 1 -doSaf 1 -skipTriallelic 1 -doCounts 1 -sites <.keep file> -rf <.rf file> -out <output file name>
```

Final parameter choice: `SNP_pval` = 1e-3 and	`postCutOff` = 0.95; `hetBias` and `-C` not used

Wrapper scripts for final analysis in paper:
```
sbatch optgeno_cviol_wrap.sh
sbatch optgeno_ccoru_wrap.sh
```
