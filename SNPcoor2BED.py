#! /usr/bin/env python

import sys

helpMsg = '''

	usage: $./SNPcoor2BED.py chr inFile.txt outFile.bed

'''

def main(args):
	if len(args) != 4:    #3 argument
		return helpMsg

	coorF = open(args[2], 'r')
	bedF = open(args[3], 'w')

	for line in coorF:
		SNPcoor = int(line.strip())
		bedEntry = args[1] + "\t" + str(SNPcoor-1) + "\t" + str(SNPcoor) + "\n"
		bedF.write(bedEntry)

	coorF.close()
	bedF.close()

	return 0

sys.exit(main(sys.argv))

