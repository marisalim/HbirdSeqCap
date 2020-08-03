#!/usr/bin/env python

# Calculate the length of CDS from Anna's hummingbird genome
# Use subsetted .gff file (made from: awk '$3 == "CDS"' Calypte_anna.gff)

import pandas

import matplotlib # must be here or get errors with histogram
matplotlib.use('Agg') # must be here or get errors with histogram
import matplotlib.pyplot as plt

# --- Open file ---
# open CDS_Calypte_anna.gff
myfile = '/crucible/bi4iflp/mlim/Calypte_genome_exoncheck/CDS_Calypte_anna.gff'
#myfile = '/crucible/bi4iflp/mlim/Calypte_genome_exoncheck/subset_test.gff'
f = open(myfile, 'r')

# --- Calculate CDS length ---
# want to subtract column 5 - column 4
# save as new df with seqid, length, attribute columns 
SeqID = []
CDS_length = []
Starts = []
Ends = []
Attributes = []
mylines = f.readlines()
for aline in mylines:
	columns = aline.split()
	endvals = columns[4]
	startvals = columns[3]
	seqids = columns[0]
	attr = columns[8]
	parentattr = attr.split('_')[1].replace(';','')
	
	cds_length = int(endvals) - int(startvals)
#	print '%s - %s = %s' % (endvals, startvals, cds_length)
	
	# create lists for values you want to keep
	SeqID.append(seqids)
	CDS_length.append(cds_length)
	Starts.append(startvals)
	Ends.append(endvals)
	Attributes.append(parentattr)
	
#print SeqID
#print CDS_length
#print Attributes

#print 'Min: ', min(CDS_length) # 0
#print 'Max: ', max(CDS_length) # 8779

# --- Create dictionary of lists ---
mydict = {'SeqID': SeqID, 'CDS_length': CDS_length, 'Start_val': Starts, 'End_val': Ends, 'Attributes': Attributes}
#print mydict

# --- Save dictionary as dataframe ---
## commented out because done, just uncomment to run again
df = pandas.DataFrame(mydict)
df.to_csv('outfile.csv', index=False)

# --- Plot histogram of CDS lengths ---
#plt.figure(0)
#plt.hist(CDS_length, bins=100, color='b', alpha=0.5)
#plt.title("CDS length histogram")
#plt.xlabel("Length (bp)")
#plt.ylabel("Frequency")
#plt.savefig('CDS_hist.png')

#below1000_cds = [i for i in CDS_length if i < 1000]
#plt.figure(1)
#plt.hist(below1000_cds, bins=100, color='r', alpha=0.5)
#plt.title("CDS length histogram")
#plt.xlabel("Length (bp)")
#plt.ylabel("Frequency")
#plt.savefig('CDSbelow1000_hist.png')

#bt400_600_cds = [i for i in CDS_length if i >= 400 and i <= 600]
#plt.figure(2)
#plt.hist(bt400_600_cds, bins=100, color='g', alpha=0.5)
#plt.title("CDS length histogram")
#plt.xlabel("Length (bp)")
#plt.ylabel("Frequency")
#plt.savefig('CDS400_600_hist.png')

#above1000_cds = [i for i in CDS_length if i > 1000]
#plt.figure(3)
#plt.hist(above1000_cds, bins=100, color='m', alpha=0.5)
#plt.title("CDS length histogram")
#plt.xlabel("Length (bp)")
#plt.ylabel("Frequency")
#plt.savefig('CDSabove1000_hist.png')
