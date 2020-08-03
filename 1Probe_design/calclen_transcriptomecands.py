# calculate sequence length of candidate genes pulled from my transcriptome data

from Bio import SeqIO 
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

CDS_length = []
#theseq = []

for myseq in SeqIO.parse(open('cands_in_transcriptomes.fa'), 'fasta'):
	name = myseq.id
	sequence = myseq.seq
	cds_length = len(sequence)
	print cds_length
	CDS_length.append(cds_length)
#	theseq.append(sequence)
	
print 'Target size of cands_in_transcriptomes: ', sum(CDS_length) #27114


