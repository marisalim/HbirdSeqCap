# script to plot out cds and gccont distribution of random markers (following lenfilter.py script)

import pandas as pd

import matplotlib # must be here or get errors with histogram
matplotlib.use('Agg') # must be here or get errors with histogram
import matplotlib.pyplot as plt

randommarkfile = '/pylon1/bi4iflp/mlim/Calypte_genome_exoncheck/myrandommarkers.csv'
data = pd.read_csv(randommarkfile)

gctoplot = data['GCcontent']

plt.figure(0)
plt.hist(gctoplot, bins=100, color='b', alpha=0.5)
plt.title("GC% histogram")
plt.xlabel("GC%")
plt.ylabel("Frequency")
plt.savefig('RandmarksGCcont_hist.png')

cdstoplot = data['Filt_cdslen']

plt.figure(1)
plt.hist(cdstoplot, bins=100, color='r', alpha=0.5)
plt.title("CDS length histogram")
plt.xlabel("CDS length")
plt.ylabel("Frequency")
plt.savefig('RandmarksCDSlen_hist.png')

Nsfile = '/pylon1/bi4iflp/mlim/Calypte_genome_exoncheck/len_filtered.csv'
data1 = pd.read_csv(Nsfile)
nstoplot = data1['Ns']
above0Ns = [i for i in nstoplot if i > 0]
print len(above0Ns) #244

plt.figure(2)
plt.hist(above0Ns, bins=100, color='g', alpha=0.5)
plt.title("%N's histogram above 0")
plt.xlabel("%N's")
plt.ylabel("Frequency")
plt.savefig('RandmarksNumNs_hist.png')

Nsfile1 = data['Ns']
plt.figure(3)
plt.hist(Nsfile1, bins=100, color='g', alpha=0.5)
plt.title("%N's histogram final random marker set")
plt.xlabel("%N's")
plt.ylabel("Frequency")
plt.savefig('RandmarksNumNsFINAL_hist.png')






