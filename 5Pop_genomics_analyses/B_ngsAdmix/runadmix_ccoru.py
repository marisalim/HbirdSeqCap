# script to run ngsAdmix
# starting with 10 iterations per K for K=1-5

import os

num_iters = [1,2,3,4,5,6,7,8,9,10]
num_K = [1,2,3,4,5]

for i in num_K:
	for j in num_iters:
		variables = dict(
		myK = str(i), myiter = str(j),
		
		# for pass 1 with all samples
		#input = '/pylon5/bi4iflp/mlim/SeqCapData/angsd_analysis/ngsAdmix/ccoru_foradmix_out.beagle.gz',
		#output = '/pylon5/bi4iflp/mlim/SeqCapData/angsd_analysis/ngsAdmix/Ccoru_admix/ccoru_K' + str(i) + '-' + str(j)) 

		# for pass 2 with 3 outliers removed
		#input = '/pylon5/bi4iflp/mlim/SeqCapData/angsd_analysis/ngsAdmix/ccoru_forAdmixn99_out.beagle.gz',
		#output = '/pylon5/bi4iflp/mlim/SeqCapData/angsd_analysis/ngsAdmix/Ccoru_admix_n99/ccoru99_K' + str(i) + '-' + str(j)) 
		
		# for pass 3 with 2 more outliers removed (may not run if PCA looks ok)
		input = '/pylon5/bi4iflp/mlim/SeqCapData/angsd_analysis/ngsAdmix/ccoru_forAdmixn97_out.beagle.gz',
		output = '/pylon5/bi4iflp/mlim/SeqCapData/angsd_analysis/ngsAdmix/Ccoru_admix_n97/ccoru97_K' + str(i) + '-' + str(j)) 
		
		
		commands = """
		echo 'Running K = {myK} and iter = {myiter}'
		
		# run ngsAdmix
		NGSadmix -likes {input} -K {myK} -P 4 -minMaf 0.05 -o {output}
		
		""".format(**variables)
		
		command_list = commands.split('\n')
		for cmd in command_list:
			os.system(cmd)
