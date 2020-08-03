# Before running selection tests, the SNPs with low frequency 
# (minor allele frequency (MAF) of 5%) have to be filtered out
# Remove sites with allele frequencies < 0.05

# Steps:
# 1. Read in mafs.gz files 
# 2. Filter by 7th column of maf file (this is the minor allele since I used -doMajorMinor 1)
# 3. Keep all sites with MAF > 0.05
# 4. Match these sites to 2nd column of geno file and keep just the sites that match

import pandas as pd
import csv

# Remove header from angsd mafs file first - loop can't skip that header line
# There are also duplicate site #s between contigs, so I can't just use the site number to filter otherwise
# in the genofile loop, I retain more sites than I should
# So, I concatenated the contig and site number IDs so that each one will be unique
# I'll remove the concatenated column from geno file after filtering (since I keep the separate contig and site columns in geno file)
##1. MAF file changes:
#zcat cviol59_optP95_out.mafs.gz | tail -n +2 | awk '{print $1"_"$2"\t"$7}' > cviol59_optP95_out_FILTER.mafs
##2. GENOfile changes:
#zcat cviol59_optP95_out.geno.gz | awk '{print $1"_"$2}' > cviol59_contigsite.txt
#gunzip cviol59_optP95_out.geno.gz
#paste cviol59_contigsite.txt cviol59_optP95_out.geno > cviol59_optP95_out_FILTER.geno

# FILE NAMES
# File1 (input mafs file): 'cviol59_optP95_out_FILTER.mafs'
# File2 (input geno file): 'cviol59_optP95_out_FILTER.geno'
# File3 (output filtered geno file): 'MAFfiltered_cviol59_out.csv'

maf_threshold = 0.05
MAF_list = []
ContigSite_list = []

maf_counter1 = 0
maf_counter2 = 0
with open('cviol59_optP95_out_FILTER.mafs', 'r') as maf:
	for line in maf:
		Contig_site = line.split('\t')[0]
		MAFcol = line.split('\t')[1]
		
		if float(MAFcol) >= float(maf_threshold):
			#print 'high freq SNP','\t', MAFcol, '\t', Site
			MAF_list.append(MAFcol)
			ContigSite_list.append(Contig_site)
			maf_counter1 = maf_counter1 + 1
			
		else:
			#print 'low freq SNP','\t', MAFcol, '\t', Site
			maf_counter2 = maf_counter2 + 1
			continue

ContigSite_list.sort() #sort the sites numerically

# check how many sites are kept, how many should be removed
print 'Length of Site_list: ', len(ContigSite_list)
print 'MAF FILE - Number of sites kept after MAF filter: ', maf_counter1
print 'MAF FILE - Number of sites removed: ', maf_counter2

# if genofile site is in the Site_list, then keep that line, else move to next line
newgenofile = []
counter = 0
counter2 = 0
with open('cviol59_optP95_out_FILTER.geno', 'r') as genofile:
	for line in genofile:
		thesite = line.split('\t')[0]
		if thesite in ContigSite_list:
			#print 'yes', thesite
			#print line
			counter = counter + 1
			newgenofile.append(line)
		else:
			#print 'no', thesite
			counter2 = counter2 + 1
			continue

# check again how many sites are kept, how many should be removed
print 'GENO FILE - Number of sites kept after MAF filter: ', counter
print 'GENO FILE - Number of sites removed: ', counter2

mygeno = {'geno': newgenofile}
df_geno = pd.DataFrame(mygeno)
df_geno.to_csv('MAFfilteredgeno_cviol59_out.csv', index=False, header=False)






