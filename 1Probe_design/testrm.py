#!/usr/bin/env python

# script shows that all repeat sequence IDs are in the file they're being removed from
# was just usng del mydat[key] but it's failing to remove all the repeats for some reason

from collections import defaultdict

RMout = '/pylon1/bi4iflp/mlim/Calypte_genome_exoncheck/RM3.out'
repeats = []
for line in open(RMout, 'r'):
	repeat_seqid = line.split()[4]
#	print repeat_seqid
	repeats.append(repeat_seqid)	
#print repeats
print 'Length of repeat seqs: ', len(repeats)

mydat = defaultdict(list)	
with open('order_above150bp_35-70GC.bed1') as f:
	for line in f:
		seqid, cdslen, gccont, seq  = line.split()
		mydat[seqid].append((int(cdslen), float(gccont), seq))

# Find matches and remove		
print len(mydat.keys())#47512

rmcol = []

sayyes1 = []
sayno1 = []
sayyes = []
sayno = []

for key in repeats:
	if key in mydat:
		sayyes1.append(key)
		del mydat[key]
	else:
		sayno1.append(key)
print 'Num yes: ', len(sayyes1)
print 'Num no: ', len(sayno1) # these are the duplicates, so the del mydat[key] got rid of the seqID in mydat (only 1 copy), but repeats has that samd ID multiple times and since it not longer has a match in mydat, it appears to not be deleted but it's actually previously deleted.
		

for key in mydat:
	if key in repeats:
		myval = 'yes'
		sayyes.append(key)
	else:
		myval = 'no'
		sayno.append(key)
	rmcol.append(myval)

print len(rmcol) #47512

print 'Num yes: ', len(sayyes)
print 'Num no: ', len(sayno)

f = open('repeats_to_remove.txt', 'w')
f.write('\n'.join(str(x) for x in sayno1))
f.close()







