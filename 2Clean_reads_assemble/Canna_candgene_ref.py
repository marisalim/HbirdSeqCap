# this script is for making reference fasta for candidate genes designed from Anna's hbird genome
# 1. get gene name
# 2. extract sequences
# 3. add 39 N's to end of CDSs (not last CDS of gene though? or does it matter?)
# 4. concatenate CDSs for a gene
# 5. write as fasta file

import pandas as pd
import csv
from collections import defaultdict
import re
from Bio import SeqIO 
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Do the following in bash before running this script
#cat annacands_withbuffer.csv > annacands_withbuffer.bed
#tail -n +2 annacands_withbuffer.bed > annacands_withbuffer.bed1

annacands = defaultdict(list)
with open('annacands_withbuffer.bed1') as f:
	for line in f:
		CDS_ID, CDS_len, End, Gene_name, Seq_ID, Start = line.split()
		annacands[Seq_ID].append((Start, End, Gene_name, Seq_ID))

records = SeqIO.to_dict(SeqIO.parse(open('Calypte_anna.fa'), 'fasta'))

annacands_short_seq_records = []
for name in annacands:
	for (start, stop, genename, seqid) in annacands[name]:
		long_seq_record = records[name]
		long_seq = long_seq_record.seq 		
		alphabet = long_seq.alphabet
		short_seq = str(long_seq)[int(start) - 1:int(stop)]	
		short_seq_record = SeqRecord(Seq(short_seq, alphabet), id=genename, description=seqid + '_' + str(genename) + '_start' + str(start) + '_stop' + str(stop))
		annacands_short_seq_records.append(short_seq_record)

with open('annacands.fasta', 'w') as f:
	SeqIO.write(annacands_short_seq_records, f, 'fasta')

# NOW, you have a fasta file of all the candidate gene CDS sequences
# next, add 39 N's to the end of each CDS
sequences = '/pylon5/bi4iflp/mlim/SeqCapData/mkref_Canna_targs/annacands.fasta'
#sequences = '/pylon5/bi4iflp/mlim/SeqCapData/mkref_Canna_targs/testanna.fasta'
sequences2 = list(SeqIO.parse(sequences, 'fasta'))

for theseqs in sequences2:
	theseqs.seq = theseqs.seq + ('N' * 39)
SeqIO.write(sequences2, '/pylon5/bi4iflp/mlim/SeqCapData/mkref_Canna_targs/annacands_withNs.fasta', 'fasta')
	
# NOW, you have a fasta file with all the candidate gene CDSs with Ns added
# next, you have to concatenate all the CDSs for a given gene
# don't know how to concatenate gene CDSs within the same fasta file, so just going to manually concatenate for now





	
	
