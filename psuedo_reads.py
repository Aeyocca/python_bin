#!/bin/python
#psuedo_reads.py

import sys
import random
import argparse

parser = argparse.ArgumentParser(prog='PROG')
parser.add_argument('--input', required=True, help='input chloroplast assembly')
parser.add_argument('--nreads', required=False, default = 10000, help='how many reads to sim. Def 10k')
parser.add_argument('--output', required=True, help='output psuedo_reads')

args = parser.parse_args()


#just going to be super naiive by taking n substrings of length l from the chloroplast
#easier since we have a single fasta header per file and a total length of like 150kb 
#so subsetting that string 10k times is probably fast enough


read_length = 150

fasta_string = ""
with open(args.input) as fh:
	for line in fh:
		if line.startswith(">"):
			continue
		fasta_string = fasta_string + line.strip()

out_reads = []
upper_bound = len(fasta_string) - read_length - 1
for i in range(int(args.nreads)):
	slice_start = random.randint(0,upper_bound)
	loop_read = fasta_string[slice_start:(slice_start + read_length)]
	out_reads.append(loop_read)

with open(args.output, "w") as out:
	for i in range(len(out_reads)):
		tmp = out.write(">Psuedo_read_" + str(i + 1) + "\n")
		tmp = out.write(out_reads[i] + "\n")
