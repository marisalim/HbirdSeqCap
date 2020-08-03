#!/usr/bin/env python

# View GC% distribution of CDSs > 150bp

import pandas as pd

import matplotlib # must be here or get errors with histogram
matplotlib.use('Agg') # must be here or get errors with histogram
import matplotlib.pyplot as plt

gcfile = '/crucible/bi4iflp/mlim/Calypte_genome_exoncheck/GCcontentabove150bp.csv'
data = pd.read_csv(gcfile)

gctoplot = data['GCcontent']

plt.figure(0)
plt.hist(gctoplot, bins=100, color='b', alpha=0.5)
plt.title("GC% histogram")
plt.xlabel("GC%")
plt.ylabel("Frequency")
plt.savefig('GCcont_hist.png')

bt35_70gc = [i for i in gctoplot if i >= 35 and i <= 70]
print len(bt35_70gc)

plt.figure(1)
plt.hist(bt35_70gc, bins=100, color='m', alpha=0.5)
plt.title('GC% = 35-70% histogram')
plt.xlabel('GC%')
plt.ylabel('Frequency')
plt.savefig('GCcont_35-70_hist.png')
