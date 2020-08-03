# After adding the buffer to random markers, now need to:
#1. re-check the GC content, filter to 35-70%
#2. remove sequences with Ns

import pandas as pd
import csv
from collections import defaultdict
import re
from Bio import SeqIO 
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# BEFORE this, modify for_gc_recalc.csv file.
# First save as .bed cat for_gc_recalc.csv > for_gc_recalc.bed
# tail -n +2 for_gc_recalc.bed > for_gc_recalc.bed1 # to get rid of header
mydat2 = defaultdict(list) 
with open('for_gc_recalc.bed1') as f:
	for line in f:
		end, geneid, seqid, start = line.split()
		mydat2[seqid].append((int(start), int(end), geneid)) 
		
records = SeqIO.to_dict(SeqIO.parse(open('Calypte_anna.fa'), 'fasta')) 

short_seq_records = [] 
for name in mydat2:
	for (start, stop, name2) in mydat2[name]:
		long_seq_record = records[name] #Note: name must be in the fasta file
		long_seq = long_seq_record.seq
		alphabet = long_seq.alphabet
		short_seq = str(long_seq)[start-1:stop]
		short_seq_record = SeqRecord(Seq(short_seq, alphabet), id=name + '_start' + str(start) + '_end' + str(stop) + '_' + name2, description=name)
		short_seq_records.append(short_seq_record)
#		print short_seq_records

with open('recalc_GC.fasta', 'w') as f:
	SeqIO.write(short_seq_records, f, 'fasta')

# recalculate GC and Ns
SeqID = [] 
GeneID = [] 
GCcontent = [] 
CDS_length = [] 
Startval = [] 
Endval = [] 
Nscount = []

recalc_gc = '/pylon1/bi4iflp/mlim/Calypte_genome_exoncheck/recalc_GC.fasta' 

for myseq in SeqIO.parse(open(recalc_gc), 'fasta'):
	name = myseq.id
		
	sequence = myseq.seq
	cds_length = len(sequence) # it looks like some sequences don't exist in the fasta file for the start/stops I made
	#if cds_length == 0:
	#	print (name, geneid) # not too many, let's set these aside
	
	if cds_length != 0:
		#theseq.append(sequence)

		gCount = 0
		cCount = 0
		nCount = 0
		for base in sequence:
			#print base
			if base == 'G':
				gCount = gCount + 1
			elif base == 'C':
				cCount = cCount + 1
			elif base == 'N':
				nCount = nCount + 1
#		print 'Gs: ', gCount 
#		print 'Cs: ', cCount 
		gccont = round((float(gCount) + float(cCount))/cds_length * 100, 2) 
#		print 'GC content ', gccont
		ncont = round((float(nCount))/cds_length * 100, 2)	
#		print 'Ns: ', nCount
	
		# Filter by GC content and Ns
		# Want GC% between 35-70%
		# Want to remove any sequences with Ns (so Ns count should be 0)
		if gccont >= 35 and gccont <= 70 and ncont == 0:
			name2 = name.split('_')[0]
			SeqID.append(name2) 
			CDS_length.append(cds_length)
			
			# parse out start and end vals
			thestart = name.split('_')[1]
			thestart2 = re.sub('start', '', thestart)
			Startval.append(thestart2)
	
			theend = name.split('_')[2]
			theend2 = re.sub('end','',theend)
			Endval.append(theend2)
	
			# parse out gene ID
			geneid = name.split('_')[3]
			GeneID.append(geneid)
			
			GCcontent.append(gccont) 
			Nscount.append(ncont)
				
	#else:
		#print 'This sequence has length 0. Skip.'
	
finaldict = {'SeqID': SeqID, 'GeneID': GeneID, 'Start': Startval, 'End': Endval, 'GC content': GCcontent, 'CDS_length': CDS_length, 'Ns': Nscount}
print 'Total CDS length after removing len=0 sequences:', sum(CDS_length) # 5181504
	
df = pd.DataFrame(finaldict) 
df.to_csv('randmarks_withbuffer_newfilt.csv', '\t', index=False)


	




		





