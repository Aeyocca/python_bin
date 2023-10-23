#!/bin/python
#Alan E. Yocca
#09-04-20
#filter_fasta_length.py

import re
import argparse

parser = argparse.ArgumentParser(prog='PROG')
parser.add_argument('--input', required=True, help='input fasta file')
parser.add_argument('--length', required=False, default = 1000, type = int, help='keep entries >= to this value')
parser.add_argument('--inverse', required=False, action='store_true', help='keep entries less than this value')
parser.add_argument('--output', required=True, help='output tsv file')

args = parser.parse_args()

full_seq = dict()
header = str()
print("If there are duplicate headers, this will concatenate their sequences")
with open(args.input) as fh:
	#load fasta file? only if conditions met
	for line in fh:
		if re.match("^>", line):
			header = line.strip()
		else:
			try:
				full_seq[header] += line.strip()
			except:
				full_seq[header] = line.strip()				
	
#loopidy loopy
written = 0
with open(args.output, "w") as output:
	for seq in full_seq.keys():
		if len(full_seq[seq]) >= args.length and not args.inverse:
			output.write(seq)
			output.write("\n")
			output.write(full_seq[seq])
			output.write("\n")
			written+=1
		elif len(full_seq[seq]) < args.length and args.inverse:
			output.write(seq)
			output.write("\n")
			output.write(full_seq[seq])
			output.write("\n")
			written+=1

print("output %s seqs" % (written))