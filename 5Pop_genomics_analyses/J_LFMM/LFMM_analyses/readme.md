# LFMM selection test script descriptions

Required software:
- [LFMM](http://membres-timc.imag.fr/Olivier.Francois/lfmm/software.htm) - Linux command line version
- R

## LFMM with allele frequency input:

Note: I also tried LFMM with genotype input, but the results look too 'noisy'. Ke suggested using allele frequencies calculated by elevational band groups as input for LFMM instead of genotypes. This should remove some of the noise. Based on Megan Rixey’s work (former postdoc in Nachman lab), there should be very good correlation between the LFMM top outlier SNPs and allele frequency by elevation band for that SNP. For the analysis with genotypes, the correlations for the top SNPs with elevation were not that good. 

What we need:
- use elevational bins (don’t necessarily have to be equal units, just shouldn’t be drastically different)
- calculate allele frequencies (with respect to alternative allele (q allele)) for each bin. The input should be rows = populations and columns = SNPs
- try K=1, K=2, or K=3 (note: allele frequency method doesn’t work well for large K values. In fact, it’s recommended by authors to choose K based on lambda value, rather than from structure methods). Did K=1 and 2 for Ccoru, K=2 and 3 for Cviol.

### 1. Define elevation bins
- Using 400 m bins for Cviol and Ccoru (try to divide samples amongst bins as evenly as possible, tried a few different bin sizes before deciding to use 400 m). Drop samples from bins if they are too genetically distinct (based on pop gen structure analyses), while ensuring you have enough samples per bin to calculate allele frequency. 
  
### 2. Calculate allele frequencies for all SNPs per bin and format LFMM input file

- Use `inputprepforLFMM.r` script to calculate allele frequencies and prep inputs for LFMM

### 3. Run LFMM

Run LFMM for each species, for K1-3 with 5 iterations each - I actually ran these for the different bin sizes as well. I used batch scripts for running LFMM [example for allele frequency input](./LFMM_fromallelefreq); it does require a lot of parameter sensitivity testing - refer to LFMM documentation for help. The command is the same, the input files change. Note that the latest version of LFMM (v1.5) does not require the `-D` flag anymore (it's auto-detected).

### 4. Check lambda and p-value distributions for LFMM outputs

- Use `inputprepforLFMM.r` to calculate lambda and plot p-value distributions
