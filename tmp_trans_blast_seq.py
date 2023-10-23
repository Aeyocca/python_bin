#!/bin/python
#tmp_trans_blast_seq.py

import argparse
import sys

parser = argparse.ArgumentParser(prog='PROG')
parser.add_argument('--blast', required=True, help='input blast file')
parser.add_argument('--seqID', required=True, help='input orthofinder SequenceIDs file')
parser.add_argument('--output', required=True, help='output file')
args = parser.parse_args()

#loopy loops
seq_trans = dict()

with open(args.seqID) as fh:
	for line in fh:
		la = line.strip().replace(" ","").split(":")
		seq_trans[la[0]] = la[1]

with open(args.blast) as fh:
	with open(args.output, 'w') as out:
		for line in fh:
			la = line.strip().split("\t")
			la[0] = seq_trans[la[0]]
			la[1] = seq_trans[la[1]]
			
			for col in la[:-1]:
				tmp = out.write(col + "\t")
			tmp = out.write(la[-1] + "\n")