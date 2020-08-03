# Dxy script descriptions

Required software:
- ANGSD ngsStat for dxy (in ngsTools in [ngsPopGen](https://github.com/mfumagalli/ngsPopGen)) and realSFS for theta

# Dxy

## a. Call SNPs for all individuals together (already did this in [SNP filter optimization tests](../I_ANGSD_genotype_calls/readme.md)
- use the files that have SNP_pval 1e-3 and postCutOff 0.95)
- These were the wrapper [scripts](../I_ANGSD_genotype_calls/):
```
sbatch optgeno_cviol_wrap.sh
sbatch optgeno_ccoru_wrap.sh
```

## b. Get nonvariable sites

Use above ANGSD script, but take out these filters below to get nonvariable sites, also add minIndDespth=5 instead of geno_minDepth:
- SNP_pval
- postCutOff
- HWE_pval
- dosnpstat
- geno_minDepth

Use last line of wrapper [scripts](../I_ANGSD_genotype_calls/):
```
# for dxy calculations - no snp calling, replace geno_minDepth with minIndDepth
sbatch optgeno_cviol_wrap.sh
sbatch optgeno_ccoru_wrap.sh
```

## c. Now, for all the sites that are in step b geno file, but not step a geno file, add to the step a geno file
- The chr and position names don't matter here
- Just make sure for each individual, the geno value is 0
- You can even do all this with all the same chr and position, just need to get the right number of sites for this category

First, I need to compare the geno files to figure out which sites are in step b geno file (SNPs NOT called) (non-variant sites), but NOT in step a geno file (SNPs called - variant sites).
- Let’s create text files that just have contig and site numbers concatenated so that I can compare these between geno files:
```
#non-variant sites
zcat cviol59_nosnpcalls_fordxy_outs.geno.gz | cut -f 1,2 | awk '{print $1 "_" $NF}' > non-variantsites_cviol59.txt
zcat ccoru97_nosnpcalls_fordxy_outs.geno.gz | cut -f 1,2 | awk '{print $1 "_" $NF}' > non-variantsites_ccoru97.txt

#variant sites
zcat cviol59_optP95_out.geno.gz | cut -f 1,2 | awk '{print $1 "_" $NF}' > variantsites_cviol59.txt
zcat ccoru97_optP95_out.geno.gz | cut -f 1,2 | awk '{print $1 "_" $NF}' > variantsites_ccoru97.txt
```

Now, I need to compare the variantsites_* and non-variantsites_* text files (in dxy directory):
```
comm -3 ./nosnpcalls_fordxy/non-variantsites_cviol59.txt ./snpcallgenofiles_fordxy/variantsites_cviol59.txt > sitediff_cviol59.txt
#same thing for rest of species
```
- lines that are NOT indented in sitediff_* file are only in non-variantsites_*
	- I want to keep just these!
- lines that ARE indented in sitediff_* file are only in variantsites_*
```
grep 'combined_Contig1_246' ./nosnpcalls_fordxy/non-variantsites_cviol59.txt
grep 'combined_Contig1_246' ./snpcallgenofiles_fordxy/variantsites_cviol59.txt
```

Need to filter sitediff_* file to get just the non-indented lines. The easiest way to do this (I failed to use bash to remove lines with tabs) is subtract number lines with tabs from total number of lines in the file.
- Total number of lines: `wc -l [filename]`
- Number of lines with tabs: `< [filename] tr -dc \\t | wc -c`

Then, because the contig and site numbers themselves don’t matter, just get the dimensions of how many sites are non-variant. Create new geno file with variant sites + block of non-variant sites from above (call these all one contig name and one site name). Give all the non-variant sites 0’s for the genotype.
- In snpcallgenofiles_fordxy directory, run `make_nonvar_matrix.r`
- Run:
```
sbatch makenonvarmatrix_wrap.sh
```
- This script should 1) make matrix of 0’s for the non-variant sites, 2) make matrix of contig and site (sequential, can’t be the same number) labels for non-variant sites, 3) cbind these matrices, and 4) save them as .gz files. Then, you need to 5) get rid of the “” in the nonvar0mat files and 6) use bash to cat the variant sites geno files with the non-variant sites nonvar0mat files. Use sed command:
```
# the s/ searches for the pattern to replace, the /g tells the command to do this globally over the entire file rather than stopping at the first match
zcat cviol_nonvar0mat.txt.gz | sed ‘s/\”//g’ | gzip > cviol_nonvar0mat_noquotes.txt.gz
zcat cviol_nonvar0mat_noquotes.txt.gz | less -S
cat cviol59_optP95_out.geno.gz cviol_nonvar0mat_noquotes.txt.gz > cviol59_fordxy.geno.gz
```

Delete the nonvar0mat.txt.gz files when done (the files that have the “”).

## d. Need to subset the geno files by population (subset keeping first 2 columns, and then the correct columns afterwards for the individuals that belong to each population.

Subset by population designations used for Fst, theta, original dxy calc. Have to figure out order of individuals for each population (currently in order of all indiv/sp bam files, but check with the sp pop bam files).

In snpcallgenofiles_fordxy directory, make the geno files subsetted for each population:
```
zcat cviol59_fordxy.geno.gz | cut -f 1,2,21-26 | gzip > cviol_pop1_fordxy.SNP.gz
zcat cviol59_fordxy.geno.gz | cut -f 1,2,6-19 | gzip > cviol_pop2_fordxy.SNP.gz
zcat cviol59_fordxy.geno.gz | cut -f 1,2,3-5,20 | gzip > cviol_pop3_fordxy.SNP.gz
zcat cviol59_fordxy.geno.gz | cut -f 1,2,27-32,53-56 | gzip > cviol_pop4_fordxy.SNP.gz
zcat cviol59_fordxy.geno.gz | cut -f 1,2,43-52,57-61 | gzip > cviol_pop5_fordxy.SNP.gz
zcat cviol59_fordxy.geno.gz | cut -f 1,2,33-42 | gzip > cviol_pop6_fordxy.SNP.gz

zcat ccoru97_fordxy.geno.gz | cut -f 1,2,47,72-79 | gzip > ccoru_pop1_fordxy.SNP.gz
zcat ccoru97_fordxy.geno.gz | cut -f 1,2,63-67,93-99 | gzip > ccoru_pop2_fordxy.SNP.gz
zcat ccoru97_fordxy.geno.gz | cut -f 1,2,68-71 | gzip > ccoru_pop3_fordxy.SNP.gz
zcat ccoru97_fordxy.geno.gz | cut -f 1,2,3-8,22-46,54-56,89-92 | gzip > ccoru_pop4_fordxy.SNP.gz
zcat ccoru97_fordxy.geno.gz | cut -f 1,2,9-14,48 | gzip > ccoru_pop5_fordxy.SNP.gz
zcat ccoru97_fordxy.geno.gz | cut -f 1,2,20-21,49,86-88 | gzip > ccoru_pop6_fordxy.SNP.gz
zcat ccoru97_fordxy.geno.gz | cut -f 1,2,15-19,50-53 | gzip > ccoru_pop7_fordxy.SNP.gz
zcat ccoru97_fordxy.geno.gz | cut -f 1,2,57-62 | gzip > ccoru_pop8_fordxy.SNP.gz
zcat ccoru97_fordxy.geno.gz | cut -f 1,2,80-85 | gzip > ccoru_pop9_fordxy.SNP.gz
```

Move all pop*_fordxy.SNP.gz files to genofiles_forPopStats directory
- change extension to .SNP.gz instead of .geno.gz for PopStat script

Divide by species in to cviol_dxy and ccoru_dxy directories
	- check that they have same number of sites:
	`for i in cviol_pop*; do zcat $i | wc -l; done`

`gunzip` files (otherwise, script does not automatically gunzip and you’ll get tons of error messages)

## e. Then use [`6-PopStats`](../../CGRLScripts/6-PopStats) script to do the dxy calculations

In dxy directory, add the 6-PopStats script. From here: (will run through for all species)
```
sbatch runPopStats_wrap.sh
```

The output files will be in dxy directory - *.global.txt and also results for Dxy for each site between populations.
