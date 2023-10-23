#!/bin/python
#Alan E. Yocca
#trim_fasta.py

#semi dummy script to get first n sequences from a fasta
#fuck it lets trim all

import argparse
import re
import sys

parser = argparse.ArgumentParser(prog='PROG')
parser.add_argument('--input', required=True, help='input fasta file')
parser.add_argument('--output', help='extracted gene file')
parser.add_argument('--trim_end', required=False, default=0, type=int, help='trim off')

args = parser.parse_args()

#load chloroplast into string

acc=""
seq=dict()
with open(args.input) as fh:
	for line in fh:
		if line.startswith(">"):
			#new fasta header
			acc = line.strip()
			seq[acc] = ""
		else:
			seq[acc] += line.strip()

with open(args.output, "w") as out_fh:
	for key in seq:
		out_fh.write(key + "\n")
		out_fh.write(seq[key][0:int(args.trim_end) + 1] + "\n")















