#!/usr/bin/env python

# 1. Remove all CDS <150bp
# 2. Grab sequence
#		- Need sequence ID, start val, and end val to search for CDS in Annas hbird genome fasta file.
# 3. Count G's and C's. Calculate GC content = G's + C's/sequence length

import pandas as pd
import csv
from Bio import SeqIO 
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
 
# Step 1:
myfile = '/crucible/bi4iflp/mlim/Calypte_genome_exoncheck/outfile.csv'
data = pd.read_csv(myfile)
# if CDS_length >= 150bp, keep row
abovedat = data[data['CDS_length'] >= 150]
#abovedat.to_csv('above150bp.csv', index=False)

# reorder columns
csvtobed = '/crucible/bi4iflp/mlim/Calypte_genome_exoncheck/above150bp.csv'
#with open(csvtobed, 'r') as infile, open('order_above150bp.csv', 'a') as outfile:
#	fieldnames = ['SeqID', 'Start_val', 'End_val', 'Attributes', 'CDS_length']
#	writer = csv.DictWriter(outfile, fieldnames=fieldnames, delimiter='\t')
#	writer.writeheader()
#	for row in csv.DictReader(infile):
#		writer.writerow(row)

# Step 2: 
# BEFORE this, modify reordered csv file. First save as .bed, then cut out last column and remove header, resave as new file
# cat order_above150bp.csv > order_above150bp.bed
# cut -f1,2,3,4 order_above150bp.bed | tail -n +2 > order_above150bp.bed1
# test file: 'scaff290.bed'
positions = defaultdict(list)	
with open('order_above150bp.bed1') as f:
	for line in f:
		name, start, stop, name2 = line.split()
		positions[name].append((int(start), int(stop), name2))

# test file: 'scaffold290.fa'		
records = SeqIO.to_dict(SeqIO.parse(open('Calypte_anna.fa'), 'fasta'))

short_seq_records = []
for name in positions:
	for (start, stop, name2) in positions[name]:
		long_seq_record = records[name] #Note: name must be in the fasta file 
		long_seq = long_seq_record.seq
		alphabet = long_seq.alphabet
		short_seq = str(long_seq)[start-1:stop]		
		short_seq_record = SeqRecord(Seq(short_seq, alphabet), id=name + '_start' + str(start) + '_' + name2, description=name)
		short_seq_records.append(short_seq_record)
#		print short_seq_records
		
with open('output.fasta', 'w') as f:
	SeqIO.write(short_seq_records, f, 'fasta')

# Step 3	
# for that CDS, count G's and C's/sequence length * 100 for GC content
SeqID = []
GCcontent = []
CDS_length = []
theseq = [] # save this as fasta output for RepeatMasker

trimmed_seqs = '/crucible/bi4iflp/mlim/Calypte_genome_exoncheck/output.fasta'
for myseq in SeqIO.parse(open(trimmed_seqs), 'fasta'):
	name = myseq.id
	sequence = myseq.seq
	cds_length = len(sequence)
	CDS_length.append(cds_length)
	theseq.append(sequence)

#	print name
#	print cds_length
#	print sequence
	gCount = 0
	cCount = 0
	for base in sequence:
		#print base
		if base == 'G':
			gCount = gCount + 1
		elif base == 'C':
			cCount = cCount + 1 
#	print 'Gs: ', gCount
#	print 'Cs: ', cCount
	gccont = round((float(gCount) + float(cCount))/cds_length * 100, 2)
#	print 'GC content ', gccont

	# save GC% + seqid (want values to be appended)
	SeqID.append(name)
	GCcontent.append(gccont)	
	
#print SeqID
#print GCcontent

mydict = {'SeqID':SeqID, 'GCcontent':GCcontent, 'CDS_length': CDS_length, 'Seq': theseq}
#print mydict

df = pd.DataFrame(mydict)
df.to_csv('GCcontentabove150bp.csv', index=False) # this cvs file has SeqID (scaffold, attribute ID, start value (because otherwise scaffold and attr ID are the same for a lot of CDSs)), GC%, and CDS length




