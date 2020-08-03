#!/usr/bin/env python

# After running RepeatMasker, you will have .out files that show the query sequence id for each identified repeat sequence. 
# We want to remove these for downstream analysis.

# Step 1: Read in files to remove repeats from above150bp_35-70GC.csv
# Step 2: Read in RepeatMasker .out file (concatenated file)
# Step 3: If SeqID == query sequence from RepeatMasker, remove 
# Step 4: Save new .csv

import pandas as pd
import csv
from Bio import SeqIO 
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict

# Step 1
mycsvfile = '/pylon1/bi4iflp/mlim/Calypte_genome_exoncheck/above150bp_35-70GC.csv'
data = pd.read_csv(mycsvfile)
#myfafile = SeqIO.to_dict(SeqIO.parse(open('forRMasker3.fa'), 'fasta')) # this is a dictionary

# Step 2:
RMout = '/pylon1/bi4iflp/mlim/Calypte_genome_exoncheck/RM3.out'
repeats = []
for line in open(RMout, 'r'):
	repeat_seqid = line.split()[4]
#	print repeat_seqid
	repeats.append(repeat_seqid)	
#print repeats
print 'Length of repeat seqs: ', len(repeats)

# Step 3/4: 
# Run this ONCE and delete if run again because it appends:
#csvtobed = '/pylon1/bi4iflp/mlim/Calypte_genome_exoncheck/above150bp_35-70GC.csv'
#with open(csvtobed, 'r') as infile, open('order_above150bp_35-70GC.csv', 'a') as outfile:
#	fieldnames = ['SeqID', 'CDS_length', 'GCcontent', 'Seq']
#	writer = csv.DictWriter(outfile, fieldnames=fieldnames, delimiter='\t')
#	writer.writeheader()
#	for row in csv.DictReader(infile):
#		writer.writerow(row)

# BEFORE this, modify reordered csv file. First save as .bed, then cut out last column and remove header, resave as new file
# cat order_above150bp_35-70GC.csv > order_above150bp_35-70GC.bed
# cut -f1,2,3,4 order_above150bp_35-70GC.bed | tail -n +2 > order_above150bp_35-70GC.bed1

# Need to make csv a dictionary
mydat = defaultdict(list)	
with open('order_above150bp_35-70GC.bed1') as f:
	for line in f:
		seqid, cdslen, gccont, seq  = line.split()
		mydat[seqid].append((int(cdslen), float(gccont), seq))

		
# Find matches and remove		
print len(mydat) #47512
for key in repeats:
	if key in mydat:
		del mydat[key]
print len(mydat) #44319

# Save dictionary as new csv file
SeqID = mydat.keys()
CDS_length = []
GCcontent = []
Seqs = []

list_values = [v for v in mydat.values()]
for listovals in list_values:
	for vals in listovals:
		cdslen = vals[0]
		gccont = vals[1]
		myseqs = vals[2]
	
		CDS_length.append(cdslen)
		GCcontent.append(gccont)
		Seqs.append(myseqs)
		
mycsvdict = {'SeqID': SeqID, 'CDS_length': CDS_length, 'GCcontent': GCcontent, 'Seqs': Seqs}
	
df = pd.DataFrame(mycsvdict)
df.to_csv('norepeats.csv', index=False)	
	

# if you want a fasta file of this, look at fastaforRepeatMask.py for scripts		



 


