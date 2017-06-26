#! /usr/bin/env python

import sys

helpMsg = '''

	usage: $./blast2invl.py blast_hit_table q_name s_name

	The input should be a tab-delimited hit-table (1-based) downloaded directly from NCBI with the following columns
	0.query id, 1.subject ids, 2.query acc.ver, 3.subject acc.ver, 4.% identity, 5.alignment length, 6.mismatches, 7.gap opens, 8.q. start, 9.q. end, 10.s. start, 11.s. end, 12.evalue, 13.bit score

	The output are two BED files containing the intervals in the query and the subject that correspond to all BLAST hits
	BED format:
	<chr>	<start>	<stop>
	Tab delimited, 0-based for the start coordinates

	q_name is the identifier for the query written to the output BED file
	s_name is the identifier for the subject written to the output BED file
'''

def main(args):
	if len(args) != 4:    #3 argument
		return helpMsg

	inTable = open(args[1],'r')

	qOut = open(args[1]+'_query.bed', 'w')
	sOut = open(args[1]+'_subject.bed', 'w')

	for line in inTable:
		if line[0] == '#' or line == '\n':
			continue

		cols = line.strip().split('\t')

		# starting index has to be less than the ending index in a BED file
		if int(cols[8]) < int(cols[9]):
			q_s = int(cols[8])
			q_e = int(cols[9])
		else:
			q_s = int(cols[9])
			q_e = int(cols[8])

		if int(cols[10]) < int(cols[11]):
			s_s = int(cols[10])
			s_e = int(cols[11])
		else:
			s_s = int(cols[11])
			s_e = int(cols[10])

		qOut.write(args[2]+'\t'+str(q_s-1)+'\t'+str(q_e)+'\n')	# adjust to 0-based coordinate system
		sOut.write(args[3]+'\t'+str(s_s-1)+'\t'+str(s_e)+'\n')	# ditto

	inTable.close()
	qOut.close()
	sOut.close()

	return 0


sys.exit(main(sys.argv))