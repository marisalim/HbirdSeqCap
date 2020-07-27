# this script is for making reference fasta for random markers designed from Anna's hbird genome
# 1. get gene name
# 2. extract sequences
# 3. write as fasta file

import pandas as pd
import csv
from collections import defaultdict
import re
from Bio import SeqIO 
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

#cat randmarks_wo_negcontrols.csv > randmarks_wo_negcontrols.bed
#tail -n +2 randmarks_wo_negcontrols.bed > randmarks_wo_negcontrols.bed1	

randmarks = defaultdict(list)
with open('randmarks_wo_negcontrols.bed1') as f:
	for line in f:
		CDS_length, End, GCcontent, GeneID, Ns, SeqID, Start = line.split()
		randmarks[SeqID].append((Start, End, GeneID, SeqID))

records = SeqIO.to_dict(SeqIO.parse(open('Calypte_anna.fa'), 'fasta'))

randmarks_short_seq_records = []
for name in randmarks:
	for (start, stop, genename, seqid) in randmarks[name]:
		long_seq_record = records[name]
		long_seq = long_seq_record.seq
		alphabet = long_seq.alphabet
		short_seq = str(long_seq)[int(start)-1:int(stop)]
		short_seq_record = SeqRecord(Seq(short_seq, alphabet), id=genename, description=seqid + '_' + str(genename) + '_start' + str(start) + '_stop' + str(stop))
		randmarks_short_seq_records.append(short_seq_record)

with open('randommarkers.fasta', 'w') as f:
	SeqIO.write(randmarks_short_seq_records, f, 'fasta')



	
	
	
	
	
