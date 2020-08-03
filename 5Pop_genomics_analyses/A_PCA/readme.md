# Principal Component Analysis (PCA) script descriptions

Required software:
- [ngsTools](https://github.com/mfumagalli/ngsTools)
  - includes tools like angsd, ngsDist, and ngsPopGen >> add these to ~/.bashrc path
  - `PATH="$PATH:~/.local/bin/ngsTools/angsd`, etc.

## ngsTools PCA - round 1

*This was the first pass using all samples/species. Look at the results and remove any super outliers.*

1. Make inputs for PCA

ANGSD command:
```
angsd -bam <path to bam .list file> -ref <path to species reference fasta> -anc <path to species reference fasta> -only_proper_pairs 0 -minMapQ 10 -minQ 20 -fold 1 -GL 1 -doGeno 32 -doPost 1 -doMaf 1 -doMajorMinor 1 -doSaf 1 -sites <path to .keep file> -rf <path to .rf file> -out <output path> -skipTriallelic 1 -minInd <80% of individuals> -setMinDepth <minInd * 5> -doCounts 1 -SNP_pval 1e-3
```

Wrapper scripts:
```
sbatch angsd_ccoru_wrap.sh
sbatch angsd_cviol_wrap.sh
```

2. To do the PCA, we need the number of sites a.k.a., the number of contigs. You can find the number site, `$NSITES`, by doing: `zcat ALL.mafs.gz | tail -n+2 | wc -l`. These will be a subset of the total number of sites in the .keep file: `wc -l *.keep`.

3. PCA with ngsTools ngsCovar (from ngsPopGen)

PCA command: (this was for ccoel)
```
ngsCovar -probfile <path to output file from angsd command above>  -outfile <output file> -nind <number of samples> -nsites <number of sites from step 2 above>  -call 0 -norm 0
```

Wrapper script: `sbatch ngsCovar_wrap.sh`

## ngsTools PCA - round 2

*After initial PCA above, plus additional first pass analyses (SFS, admixture, ngsrelate), we identified outlier samples or samples that were too closely-related to others. These were excluded from downstream analysis and I had to redo the PCA. We adjusted parameter settings for the angsd step.*

1. Make inputs for PCA
ANGSD command:
```
angsd -bam <path to bam .list file> -ref <path to species reference fasta> -anc <path to species reference fasta> -only_proper_pairs 0 -minMapQ 10 -minQ 20 -fold 1 -GL 1 -doGeno 32 -doPost 1 -doMaf 1 -doMajorMinor 1 -doSaf 1 -sites <path to .keep file> -rf <path to .rf file> -out <output path> -skipTriallelic 1 -minInd <70% of individuals> -minIndDepth 3 -doCounts 1
```

Notes on flags:
- using `-minIndDepth` = 3, instead of `-setMinDepth`.

Wrapper scripts:
```
sbatch pca_ccorun97_wrap.sh
sbatch pca_cvioln59_wrap.sh
```

2. Get number of sites (same directions as above)
- the number of sites was much higher than previous approach and about the same as total number of sites in .keep file. It seems that removing outliers makes a difference.

3. PCA (same as above)

Wrapper script: `sbatch ngsCovar_wrap.sh`

## Visualize PCA plots

Use .covar files from ngsTools PCA output. Code in [`MSfigurecode.rmd`](../MSfigurecode.rmd)
