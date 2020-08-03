# script to add 100bp (or less depending on start position) to CDS for candidate genes pulled from Anna's hbird genome

import pandas as pd
import csv
from collections import defaultdict

# Step 1:
# Make csv data into dictionary
csvtobed = '/pylon1/bi4iflp/mlim/Calypte_genome_exoncheck/onlycandgenes.csv'
targlen_cand = pd.read_csv(csvtobed, header=None)
targlen_cand.columns=['Attributes', 'Gene_name', 'CDS_length', 'End', 'Seqid', 'Start']

targlen_cand.to_csv('labeledcandgenes_annas.csv', '\t', index=False)

# BEFORE this, modify csv file. First save as .bed
# cat labeledcandgenes_annas.csv > labeledcandgenes_annas.bed
# tail -n +2 labeledcandgenes_annas.bed > labeledcandgenes_annas.bed1 # to get rid of header

# Need to make csv a dictionary
mydat = defaultdict(list)	
with open('labeledcandgenes_annas.bed1') as f:
	for line in f:
		cdsid, genename, cdslen, end, seqid, start  = line.split()
		myline = (genename, cdsid, cdslen, end, seqid, start)
		mydat[genename].append(myline) 
		
# Add buffer bp to start positions
mynewstarts = []
gene_name = []
CDS_id = []
myendval = []
CDS_len = []
seqid = []

for vals in mydat.values():	
	for thevals in vals:	
		genename = thevals[0]
		cdsid = thevals[1]
		myseqid = thevals[4]
		endval = thevals[3]
		startval = thevals[5]
		newstartval = int(startval) - int(100)
		
		#print (startval, newstartval)

		if newstartval < 0:
			newstartval2 = int(startval) - int(startval) + 1
			
		else:
			newstartval2 = newstartval
		
		# calculate new CDS lengths
		new_cds_length = int(endval) - int(newstartval2)
		
		mynewstarts.append(newstartval2)
		gene_name.append(genename)
		seqid.append(myseqid)
		CDS_id.append(cdsid)
		myendval.append(endval)
		CDS_len.append(new_cds_length)
		
print 'Total CDS len:', sum(CDS_len) #864774

mynewdict = {'Seq_ID':seqid,  'Start': mynewstarts, 'End': myendval, 'CDS_ID': CDS_id, 'Gene_name': gene_name, 'CDS_len': CDS_len}

df = pd.DataFrame(mynewdict)
df.to_csv('annacands_withbuffer.csv', '\t', index=False)


# fixed this, but this was how to check for and remove duplicates
#print 'CHECK FOR DUPLICATE LINES'
#df2 = df.drop_duplicates()

#print sum(df2['CDS_len'])

#df2.to_csv('annacands_withbuffer.csv', '\t', index=False)





