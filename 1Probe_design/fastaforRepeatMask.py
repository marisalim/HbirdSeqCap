#!/usr/bin/env python

# Create fasta file for Repeat Masker
# These sequences have been through the following filters: 1) >150bp, 2) GC% between 35-70%
# 1. Read in GC content file and subset to 35-70%
# 2. Save fasta file with SeqID and sequences 

import pandas as pd

myfile = '/pylon1/bi4iflp/mlim/Calypte_genome_exoncheck/GCcontentabove150bp.csv'
data = pd.read_csv(myfile)
# if GC% 35-70% keep row
gcdat = data[(data.GCcontent >= 35) & (data.GCcontent <= 70)]
gcdat.to_csv('above150bp_35-70GC.csv', index=False)

myfile2 = '/pylon1/bi4iflp/mlim/Calypte_genome_exoncheck/above150bp_35-70GC.csv'
data2 = pd.read_csv(myfile2)

list_seq = data2['Seq']
list_name = data2['SeqID']

outfile = open('forRMasker3.fa', 'w')
for i in range(len(list_seq)):
	outfile.write('>' + list_name[i] + '\n' + list_seq[i] + '\n')
outfile.close()
