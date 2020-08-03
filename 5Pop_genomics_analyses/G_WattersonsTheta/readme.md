# Watterson's theta script descriptions

0. Generate folded saf files as input for SFS global estimate, by subpopulation:
```
sbatch angsd_forthetapi_cviol59.sh
sbatch angsd_forthetapi_ccoru97.sh
# repeat for other species
```

1. find a ‘global estimate’ of the SFS

Obtain the maximum likelihood estimate of the SFS using the realSFS program command. Wrapper script:
```
sbatch theta_realSFS_wrap.sh # from Watterson_theta directory
```

2. Calculate the thetas for each site. Generate the .thetas.idx files and then make the .thetas.idx.pestPG files.
```
sbatch cvioln59_doThetas_wrap.sh
sbatch ccorun97_doThetas_wrap.sh
```

3. The output from the above command are two files out.thetas.gz and out.thetas.idx. A formal description of these files can be found in the doc/formats.pdf in the angsd package. It is possible to extract the logscale persite thetas using the ./thetaStat print program with this command:
```
thetaStat print out.thetas.idx 2>/dev/null |head   # but this is log-scaled value, so you have to take exponent first for avg theta calculation
```
Use thetaStat print output, take exponent and then calculate average with Ke’s script.

Wrapper script:
```
sbatch thetaStat_wrap.sh # outputs go to thetafolder directory
```

4. Above gives per site theta, so take average of column with Watterson theta to get average Watterson theta across all sites. Use Ke’s [`average_theta.pl`](../../CGRLScripts/average_theta.pl) script:
```
#first, have to remove header from *.thetas files
tail -n+2 [.thetas file name] > ./test/[.thetas file name 2]
#in thetafolder dir:
perl average_theta.pl test/
#output results are printed to screen
```
