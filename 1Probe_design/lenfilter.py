# Select random markers based on length filter
# Input for this script are sequences that have already gone through these filters:
# 1. <150bp removed
# 2. 35-70%GC content kept

# Step 1: length filter
# Step 2: calculate total length

import pandas as pd
import csv
from collections import defaultdict

# Step 1:
# Input will be norepeats.csv 
mycsvfile = '/pylon1/bi4iflp/mlim/Calypte_genome_exoncheck/norepeats.csv'
data = pd.read_csv(mycsvfile)

# Make new column with just gene ID - Separate SeqID by '_' to get the gene id (order is scaffold_start_geneID)
geneID = [i.split('_',2)[2] for i in data['SeqID']]
data['GeneID'] = geneID

#print list(data.columns.values)
#data.to_csv('forlenfilter.csv', index=False)

# make list of unique geneIDs
uniqgenes = list(set(geneID))
#print uniqgenes
print 'Unique genes:', len(uniqgenes)
#print len(geneID)

# Make csv data into dictionary
# Run this ONCE and delete if run again because it appends:
#csvtobed = '/pylon1/bi4iflp/mlim/Calypte_genome_exoncheck/forlenfilter.csv'
#with open(csvtobed, 'r') as infile, open('order_forlenfilter.csv', 'a') as outfile:
#	fieldnames = ['SeqID', 'GeneID', 'CDS_length', 'GCcontent', 'Seqs']
#	writer = csv.DictWriter(outfile, fieldnames=fieldnames, delimiter='\t')
#	writer.writeheader()
#	for row in csv.DictReader(infile):
#		writer.writerow(row)

# BEFORE this, modify reordered csv file. First save as .bed, then cut out last column and remove header, resave as new file
# cat order_forlenfilter.csv > order_forlenfilter.bed
# tail -n +2 order_forlenfilter.bed > order_forlenfilter.bed1 # to get rid of header

# Need to make csv a dictionary
mydat = defaultdict(list)	
with open('order_forlenfilter.bed1') as f:
	for line in f:
		seqid, geneid, cdslen, gccont, seq  = line.split()
		mydat[seqid].append((seqid, geneid, int(cdslen), float(gccont), seq))

#print len(mydat)		
		
# Go through mydat, and if cdslen > 600, trim to 600bp
# save new dataframe
# then filter by GeneID to get 1 cds/gene (the longest one)
filtercdslen = []
filterseq = []
theseqid = []
thegeneid = []
thegccont = []
oldcdslen = []
Nscount = []

for vals in mydat.values():
	for thevals in vals:
		cdslen = thevals[2]
		myseqs = thevals[4]
		seqid = thevals[0]
		geneid = thevals[1]
		gccont = thevals[3]
		
		# For each of these, if any are >600bp, trim first 1-601
		if cdslen > int(600):
			#print cdslen
			newseq = myseqs[0:600]
			#print len(newseq)
			newcdslen = len(newseq)
					
		else: 
			newseq = myseqs
			#print len(newseq)
			newcdslen = cdslen
			
		nCount = 0
		for base in newseq:
			#print base
			if base == 'N':
				nCount = nCount + 1

#		print 'Ns: ', nCount
		ncont = round((float(nCount))/newcdslen * 100, 2)	
		
		Nscount.append(ncont)
		filtercdslen.append(newcdslen)
		filterseq.append(newseq)
		oldcdslen.append(cdslen) # to make sure that the right cds is getting trimmed
		theseqid.append(seqid)
		thegeneid.append(geneid)
		thegccont.append(gccont)
		
#print len(filtercdslen)
#print len(filterseq)

print 'Number of sequences:', len(filterseq)

myfilterdict = {'SeqID': theseqid, 'GeneID': thegeneid, 'GCcontent': thegccont, 'Filt_cdslen': filtercdslen, 'Filt_seq': filterseq, 'Old_cdslen': oldcdslen, 'Ns': Nscount}

df = pd.DataFrame(myfilterdict)
df.to_csv('len_filtered.csv', index=False)
						
# then randomly choose the longest CDS for the mRNA to keep 
# find duplicates in GeneID column and keep row that has largest CDS_length
myfile2 = '/pylon1/bi4iflp/mlim/Calypte_genome_exoncheck/len_filtered.csv'
subdat = pd.read_csv(myfile2)
# remove sequences with Ns
subdat2 = subdat[(subdat.Ns == 0)]
#subdat2.to_csv('check.csv', index=False) #just checking that it got rid of the sequences with Ns

nodupdat = subdat2.groupby('GeneID', group_keys=False).apply(lambda x: x.ix[x.Filt_cdslen.idxmax()]) 
nodupdat.to_csv('myrandommarkers.csv', index=False)
	
# Step 2:	
# Calculate total length 
# Something like this, where targlen_random is a csv file (taken from findcandidates.py script) -
targlen = nodupdat['Filt_cdslen']
numgenes = len(nodupdat)
thetotal = sum(targlen)
print 'Number of genes:', numgenes #12826
print 'Random marker target length:', thetotal #4,157,242

print 'Target length if add 50bp upstream of each gene:', (50*numgenes + thetotal) #4798542
print 'Target length if add 100bp upstream of each gene:', (100*numgenes + thetotal) #5439842

