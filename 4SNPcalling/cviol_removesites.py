# if sites are in combined_sites_to_remove file, remove them from sites_to_keep.keep file

# for each bed file set, you have to 
# 1) get rid of header on the sites to remove file (just manually delete)
# 2) concatenate the columns so that I can easily compare between files (do this for sites_to_remove and sites_to_keep files
#	awk '{print $1"-"$2}' < filename > newfilename
# using python script below:
# 3) compare files 
# 4) for output file, need to undo the concatenation and separate chr and site again for output file

print 'open files:'

#small_file = open('testsitestokeep.keep2','r')
#long_file = open('testcombinedsitestoremove.txt2','r')
#output_file = open('testoutput_file.txt','w')

small_file = open('/pylon5/bi4iflp/mlim/SeqCapData/Novoalign_outs/cviol_sites_to_keep.keep2','r')
long_file = open('/pylon5/bi4iflp/mlim/SeqCapData/Novoalign_outs/cviol_combined_sites_to_remove.txt2','r')
output_file = open('/pylon5/bi4iflp/mlim/SeqCapData/Novoalign_outs/cviol.keep','w')

print 'compare files:'
try:
    small_lines = small_file.readlines()
    small_lines_cleaned = [line.rstrip() for line in small_lines]
    long_file_lines = long_file.readlines()
    long_lines_cleaned = [line.rstrip() for line in long_file_lines]

    for line in small_lines_cleaned:
        if line not in long_lines_cleaned:
			line1 = line.split("-")[0]
			line2 = line.split("-")[1]
			output_file.write(line1 + '\t' + line2 + '\n')

finally:
    small_file.close()
    long_file.close()
    output_file.close()
	
print 'all done!'
