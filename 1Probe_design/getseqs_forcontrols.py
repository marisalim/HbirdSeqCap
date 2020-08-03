#This script is for designing qPCR primers for checking capture experiment. 
#It pulls out the sequences I want to send to company for primer design.
# get sequences for the pos and neg control sequences

import pandas as pd
import csv
from collections import defaultdict
import re
from Bio import SeqIO 
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

poscontrols = defaultdict(list) 
with open('poscontrols_captarg_adjusted_startstop.bed') as f:
	for line in f:
		seqid, oldstart, oldend, CDSlen, gccont, newstart, newend, newCDSlen = line.split()
		poscontrols[seqid].append((int(newstart), int(newend))) 
		
negcontrols = defaultdict(list) 
with open('negcontrols.bed') as f:
	for line in f:
		CDSlen, end, gccont, geneid, ns, seqid, start = line.split()
		negcontrols[seqid].append((int(start), int(end), geneid)) 

records = SeqIO.to_dict(SeqIO.parse(open('Calypte_anna.fa'), 'fasta')) 
		
pos_short_seq_records = [] 
for name in poscontrols:
	for (start, stop) in poscontrols[name]:
		long_seq_record = records[name] #Note: name must be in the fasta file
		long_seq = long_seq_record.seq
		alphabet = long_seq.alphabet
		short_seq = str(long_seq)[start-1:stop]
		short_seq_record = SeqRecord(Seq(short_seq, alphabet), id=name + '_start' + str(start) + '_end' + str(stop), description=name)
		pos_short_seq_records.append(short_seq_record)

with open('poscontrols.fasta', 'w') as f:
	SeqIO.write(pos_short_seq_records, f, 'fasta')

neg_short_seq_records = [] 
for name in negcontrols:
	for (start, stop, geneid) in negcontrols[name]:
		long_seq_record = records[name] #Note: name must be in the fasta file
		long_seq = long_seq_record.seq
		alphabet = long_seq.alphabet
		short_seq = str(long_seq)[start-1:stop]
		short_seq_record = SeqRecord(Seq(short_seq, alphabet), id=name + '_start' + str(start) + '_end' + str(stop) + '_' + geneid, description=name)
		neg_short_seq_records.append(short_seq_record)

with open('negcontrols.fasta', 'w') as f:
	SeqIO.write(neg_short_seq_records, f, 'fasta')

	




		






