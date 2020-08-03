#!/usr/bin/env python

# Use outfile.csv and mRNA_Calypte_anna.gff
# this code is going to do the same thing as get_annots.py, except it uses outfile.csv rather than final_nodupas_above500bp.csv
# Because for the candidate genes, I want all the exon sequence, not just 1/mRNA

# Add gene annotation column to outfile.csv by matching cds parent ID to mRNA ID column
# and keeping the Function attribute from mRNA

import pandas
import pandas as pd

cdsfile = '/pylon1/bi4iflp/mlim/Calypte_genome_exoncheck/outfile.csv'
mrnafile = '/pylon1/bi4iflp/mlim/Calypte_genome_exoncheck/mRNA_Calypte_anna.gff'
#mrnafile = '/pylon1/bi4iflp/mlim/Calypte_genome_exoncheck/test_mrna.gff'

cds_dat = pd.read_csv(cdsfile)
mrna_dat = open(mrnafile, 'r')

# get mRNA ID and gene annotation
Mrna_id = []
Function_attr = []
m_start = []
m_end = []
mylines = mrna_dat.readlines()
for aline in mylines:
	columns = aline.split()
	m_start = columns[3]
	m_end = columns[4]
	attr_id = columns[8].split(';')[0].split('_')[1]
	attr_func = columns[8].split(';')[2].split('"')[1]
	
	Mrna_id.append(attr_id)
	Function_attr.append(attr_func)
	
#print Mrna_id
#print Function_attr

mydict = {'Attributes': Mrna_id, 'Gene_name': Function_attr}

df = pandas.DataFrame(mydict)
#df.to_csv('mrna_id_genename.csv', index=False)

# Now, compare and merge into one file
mrnafile2 = '/pylon1/bi4iflp/mlim/Calypte_genome_exoncheck/mrna_id_genename.csv'
mrna_dat2 = pd.read_csv(mrnafile2)

merged = mrna_dat2.merge(cds_dat, on='Attributes', how='outer').fillna('')
#merged.to_csv('merged_forcandidates.csv', index=False)

# Now, match to candidate gene list and total the CDS_length column
matchtothis = '/pylon1/bi4iflp/mlim/Calypte_genome_exoncheck/merged_forcandidates.csv'
#matchtothis = '/pylon1/bi4iflp/mlim/Calypte_genome_exoncheck/testcand.csv'
cand_genes = '/pylon1/bi4iflp/mlim/Calypte_genome_exoncheck/mycandgenes.csv'
# Note: mycandgenes.csv should have a header called 'Gene_name' - be sure to add

match_dat = pd.read_csv(matchtothis)
cands = pd.read_csv(cand_genes)
mygenes = cands['Gene_name']
#print mygenes

out_csv = 'onlycandgenes.csv'
for gene in mygenes:
	cand_dat = match_dat[match_dat['Gene_name'] == gene]
#	print(cand_dat)
	cand_dat.to_csv(out_csv, index=False, header=False, mode='a')
	
# Note: since onlycandgenes.csv is made by appending rows, every time you run the last part of this script
# you'll want to rm onlycandgenes.csv so you don't keep appending every time the script is run

# Note: this doesn't get all of the genes, some of the names might not match perfectly so will need to check by hand
fortarget = '/pylon1/bi4iflp/mlim/Calypte_genome_exoncheck/onlycandgenes.csv'
targlen_cand = pd.read_csv(fortarget, header=None)
#print targlen_cand
targlen_cand.columns=['Attributes', 'Gene_name', 'CDS_length', 'End', 'Seqid', 'Start']
#print targlen_cand

# Calculate the sum of CDS_length for all the candidate genes in the file
thecdslens = targlen_cand['CDS_length']
thetotal = sum(thecdslens)
print thetotal
