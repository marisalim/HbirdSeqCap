# LFMM selection test script descriptions

Required software:
- [LFMM](http://membres-timc.imag.fr/Olivier.Francois/lfmm/software.htm) - Linux command line version
- R

## LFMM with allele frequency input:

Note: I also tried LFMM with genotype input, but the results look too 'noisy'. Ke suggested using allele frequencies calculated by elevational band groups as input for LFMM instead of genotypes. This should remove some of the noise. Based on Megan Rixey’s work there should be very good correlation between the LFMM top outlier SNPs and allele frequency by elevation band for that SNP. For the analysis with genotypes, the correlations for the top SNPs with elevation were not that good.

What we need:
- use elevational bins (don’t necessarily have to be equal units, just shouldn’t be drastically different)
- calculate allele frequencies (with respect to alternative allele (q allele)) for each bin. The input should be rows = populations and columns = SNPs
- try K=1 and K=2 (note: allele frequency method doesn’t work well for large K values. In fact, it’s recommended by authors to choose K based on lambda value, rather than from structure methods)

### 1. Define elevation bins
- Using 400 m bins for Cviol and Ccoru (try to divide samples amongst bins as evenly as possible, tried a few different bin sizes before deciding to use 400 m)

### 2. Calculate allele frequencies for all SNPs per bin and format LFMM input file
- Use `plotAlleleFreqs.r` script and the ALL SNPs section to calculate fequencies
- Format LFMM allele frequency and elevation files in `LFMM_allelefreqinput.r` script
- Check if mean or median changes elevation enviro file:
	- median pushes distr to right (higher elevations) compared to mean, we used mean

### 3. Run LFMM

Run LFMM for each species, for K1-3 with 5 iterations each - I actually ran these for the different bin sizes as well. These batch scripts have example code for running LFMM; it does require a lot of parameter sensitivity testing - refer to LFMM documentation for help:

Example batch scripts:
```
lfmm_iter1K1.sh
lfmm_iter1K2.sh
lfmm_iter1K3.sh

lfmm_iter2K1.sh
lfmm_iter2K2.sh
lfmm_iter2K3.sh

lfmm_iter3K1.sh
lfmm_iter3K2.sh
lfmm_iter3K3.sh

lfmm_iter4K1.sh
lfmm_iter4K2.sh
lfmm_iter4K3.sh

lfmm_iter5K1.sh
lfmm_iter5K2.sh
lfmm_iter5K3.sh
```

**Used `LFMM_allelefreqinput.r` script to calculate p-value, q-value, lambda, do enrichment analysis, compare allele frequency x elevation, an other post-LFMM analyses:**

#### 4. Check p-value distribution and lambda; calculate q-values
- check a few different lambdas surrounding the calculated lambda to check out how conservative or liberal the lambda is

#### 5. GO and KEGG enrichment
- tested, but wasn't able to pick up enriched categories because the genes were spread over too many GO and KEGG categories

#### 6. Compare top SNPs with allele frequency ~ elevation bin correlation coefficients (pearson’s) and plots

#### 7. count number of genes with significant SNPs that are previously id’d candidate genes vs. genes with significant SNPs that are either new or possibly unrelated to elevation but still show a cline..
- and then make the allele frequency and genotype plots for a few examples of SNPs in interesting genes that we have a better idea of what they do and how they’re related to elevation adaptation

#### 8. to check the geographic location of genotypes, plot them on a map; subset by elevation bin

#### 9. do Wilcoxon rank sum test for comparing means of individuals with 0, 1, or 2 genotype
- groups are by genotype, but the y-axis of the box plots is elevation, so if certain genotypes are more common at certain elevations, there may be significant difference in genotype group means
