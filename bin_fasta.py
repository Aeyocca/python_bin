#!/bin/python
#Alan E. Yocca
#12-07-22
#bin_fasta.py
#split a fasta file into equally sized bins

import csv
import sys
import argparse
import re

parser = argparse.ArgumentParser(prog='PROG')
parser.add_argument('--input', required=True, help='input fasta file')
parser.add_argument('--out_base', required=True, help='basename for output fasta files. Will append \"_\" + args.nbins + \".fasta\"')
parser.add_argument('--nbins', required=True, type = int, help='number of bins to split into')

args = parser.parse_args()

length_dict = dict()
fasta_dict = dict()
header = []

with open(args.input) as fh:
	for line in fh:
		if line.startswith(">"):
			header = line.strip()
			length_dict[header] = 0
			fasta_dict[header]= ""
		else:
			length_dict[header] += len(line.strip())
			fasta_dict[header] += line.strip()

#partition headers
bin_assignments = dict()
bin_lengths = [0 for i in range(args.nbins)]
for i in range(args.nbins):
	bin_assignments[i] = []

#sort by length
#code from ChatGPT!
sorted_lengths = dict(sorted(length_dict.items(), key=lambda x: x[1], reverse=True))

#loop and partition into bins
for seq in sorted_lengths.keys():
	#index of smallest current bin
	i = bin_lengths.index(min(bin_lengths))
	bin_assignments[i].append(seq)
	bin_lengths[i] += sorted_lengths[seq]
	
#nice, now bin_assignments is a list of fasta headers
for i in bin_assignments.keys():
	with open(args.out_base + "_" + str(i) + ".fasta", "w") as out:
		for seq in bin_assignments[i]:
			tmp = out.write(seq + "\n")
			tmp = out.write(fasta_dict[seq] + "\n")











