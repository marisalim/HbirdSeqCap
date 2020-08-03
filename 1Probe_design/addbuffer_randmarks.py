# script to add 100bp (or less depending on start position) to CDS for random markers pulled from Anna's hbird genome

import pandas as pd
import csv
from collections import defaultdict
import re

# Step 1:
# Make csv data into dictionary
csvtobed = '/pylon1/bi4iflp/mlim/Calypte_genome_exoncheck/myrandommarkers.csv'
randmark = pd.read_csv(csvtobed, header=None)
randmark.to_csv('randmarkstobed.csv', '\t', index=False)

# BEFORE this, modify csv file. First save as .bed
# cat randmarkstobed.csv > randmarkstobed.bed
# tail -n +3 randmarkstobed.bed > randmarkstobed.bed1 # to get rid of header, have to use +3 because there's an extra line of numbers at the top

# Need to make csv a dictionary
mydat = defaultdict(list)	
with open('randmarkstobed.bed1') as f:
	for line in f:
		filtcdslen, filtseq, gccont, geneid, ns, old_cdslen, seqid  = line.split()
		mydat[seqid].append((seqid, filtcdslen, filtseq, gccont, geneid))

# Add buffer bp to start positions
mynewstarts = []
geneID = []
CDS_len = []
seqids = []
endvals = [] # need to calculate by adding filtcdslen + start values

for vals in mydat.values():
	for thevals in vals:
		seqid = thevals[0]
		filtcdslen = thevals[1]
		geneid = thevals[4]
		
		geneID.append(geneid)
		
		# 1. Get start position
		seqid2 = seqid.split('_')[0]
		seqids.append(seqid2)
		
		thestart = seqid.split('_')[1]
		# get rid of 'start'
		thestart2 = re.sub('start', '', thestart)
		
		# 2. get end vals		
		endval = int(thestart2) + int(filtcdslen) - int(1)
		endvals.append(endval)
		
		# 3. new start position
		newstartval = int(thestart2) - int(100)
		
		#print (thestart2, newstartval)
		
		# Start position must be at 1, not 0! If 0, then the sequence gets completely lost
		if newstartval < 0:
			newstartval = int(thestart2) - int(thestart2) + int(1)
			
		elif newstartval == 0:
			newstartval = newstartval + int(1)
			
		else:
			newstartval = newstartval
		
		# calculate new CDS lengths
		new_cds_length = int(endval) - int(newstartval)
		
		mynewstarts.append(newstartval)
		CDS_len.append(new_cds_length)
		
print 'Total CDS len:', sum(CDS_len) # 5,421,833	

#print geneID

mynewdict = {'SeqID': seqids, 'GeneID': geneID, 'CDS_len': CDS_len, 'Start': mynewstarts, 'End': endvals}

# for recalculate GC content
mynewdict2 = {'SeqID': seqids, 'Start': mynewstarts, 'End': endvals, 'GeneID': geneID}
df = pd.DataFrame(mynewdict2)
df.to_csv('for_gc_recalc.csv', '\t', index=False)

