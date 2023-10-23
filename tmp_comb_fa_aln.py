#!/bin/python
#tmp_comb_fa_aln.py
#dummy script to combine aligned fastas with different names, do it by order
#Alan E. Yocca
#05-12-2023

import sys
import argparse

parser = argparse.ArgumentParser(prog='PROG')
parser.add_argument('--fa_list', required=True, help='list of aln fasta files')
parser.add_argument('--output', required=True, help='output bed file')

args = parser.parse_args()

out_dict = dict()
#tmp script so why not
n = 10
print("expecting %s seqs per fa" % (n))
for j in range(1,n + 1):
	out_dict[j] = ""


with open(args.fa_list) as fh:
	for line in fh:
		i = 0
		with open(line.strip()) as fa:
			for line in fa:
				if line.startswith(">"):
					i += 1
				else:
					out_dict[i] = out_dict[i] + line.strip()

with open(args.output, 'w') as out:
	for j in range(1,n+1):
		tmp = out.write(">" + str(j) + "\n")
		tmp = out.write(out_dict[j] + "\n")



