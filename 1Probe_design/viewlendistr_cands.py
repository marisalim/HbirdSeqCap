# plot length distribution of all candidate genes - from transcriptomes and from Anna's hbird genome
import pandas as pd

import matplotlib # must be here or get errors with histogram
matplotlib.use('Agg') # must be here or get errors with histogram
import matplotlib.pyplot as plt

from Bio import SeqIO 
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Candidate genes from transcriptomes:
CDS_length = []
genename = []

for myseq in SeqIO.parse(open('cands_in_transcriptomes.fa'), 'fasta'):
	name = myseq.id
	sequence = myseq.seq
	cds_length = len(sequence)
#	print cds_length
	CDS_length.append(cds_length)
	genename.append(name)

plt.figure(0)
plt.hist(CDS_length, bins=100, color='b', alpha=0.5)
plt.title("Candidate genes from transcriptomes length histogram")
plt.xlabel("Gene Length")
plt.ylabel("Frequency")
plt.savefig('Candgenes_transcriptomelen_hist.png')
print 'Number of candidate genes from transcriptomes:', len(genename)
print 'Range for cands in transcriptomes:', (min(CDS_length), max(CDS_length))

# Candidate genes from Anna's hbird genome:
fortarget = '/pylon1/bi4iflp/mlim/Calypte_genome_exoncheck/onlycandgenes.csv'
targlen_cand = pd.read_csv(fortarget, header=None)
#print targlen_cand
targlen_cand.columns=['Attributes', 'Gene_name', 'CDS_length', 'End', 'Seqid', 'Start']
#print targlen_cand

# Calculate the sum of CDS_length for all the candidate genes in the file
thecdslens = targlen_cand['CDS_length']

plt.figure(1)
plt.hist(thecdslens, bins=100, color='g', alpha=0.5)
plt.title("Candidate genes from Anna's hbird genome length histogram")
plt.xlabel("Gene Length")
plt.ylabel("Frequency")
plt.savefig('Candgenes_Annalen_hist.png')
print "Number of CDSs for candidate genes from Anna's hbird genome:", len(targlen_cand['Attributes'])
print "Range for cands in Anna's hbird genome, the fragments:", (min(thecdslens), max(thecdslens))

targlen_cand2 = targlen_cand.groupby(by=['Gene_name'])['CDS_length'].sum()
#print targlen_cand2

print "Number of candidate genes from Anna's hbird genome:", len(list(set(targlen_cand['Gene_name'])))

plt.figure(2)
plt.hist(targlen_cand2, bins=100, color='r', alpha=0.5)
plt.title("Candidate genes from Anna's hbird genome length concat histogram")
plt.xlabel("Gene Length")
plt.ylabel("Frequency")
plt.savefig('Candgenes_Annalenconcat_hist.png')

print "Range for cands in Anna's hbird genome, the total length:", (min(targlen_cand2), max(targlen_cand2))
