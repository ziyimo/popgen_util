#! /usr/bin/env python

import sys

helpMsg = '''

	usage: $./msCmd2discoal.py ms_cmd.txt discoal_cmd.txt win_size

		Reads ms commandline from ms_cmd.txt, 
			e.g. "./ms nsam nreps -t θ -r ρ nsites ... "
		converts to discoal commandline and writes to discoal_cmd.txt
			e.g. "./discoal sampleSize numReplicates nSites -t θ -r ρ ..."

		win_size: # of basepairs in the simulated locus
		Neutral only for now

'''

def main(args):
	if len(args) != 4:    #3 arguments
		return helpMsg

	ms_path = args[1]
	discoal_path = args[2]
	win_size = int(args[3])

	base_cmd = ["./discoal"]
	dem_cmd = []
	#selec_cmd = 

	with open(ms_path, "r") as inF:
		ms_cmd = inF.readline().strip().split()

	ms_win_size = int(ms_cmd[7]) # i.e. nsites
	scaling = win_size/ms_win_size

	base_cmd += ms_cmd[1:3]
	base_cmd.append(str(win_size))
	base_cmd += ['-t', str(float(ms_cmd[4])*scaling), '-r', str(float(ms_cmd[6])*scaling)]

	for i in range(ms_cmd.index("-eN"), len(ms_cmd), 3):
		dem_cmd += ["-en", ms_cmd[i+1], '0', ms_cmd[i+2]]	

	discoal_cmd = base_cmd + dem_cmd

	with open(discoal_path, "w") as outF:
		outF.write(" ".join(discoal_cmd))

	return 0

sys.exit(main(sys.argv))