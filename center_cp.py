#!/bin/python
#center_cp.py
#script to start chloroplast fastas at the same gene
#Alan E. Yocca
#04-19-2023

import sys
import argparse

parser = argparse.ArgumentParser(prog='PROG')
parser.add_argument('--cp_fasta', required=True, help='input maf file')
parser.add_argument('--start', required=True, help='start coordinate of gene of interest, 1-indexed')
parser.add_argument('--output', required=True, help='output msa file')

args = parser.parse_args()

#add functions later automatic gene detection		

fasta = ""
header = ""

#assumes a single fasta sequence...

with open(args.cp_fasta) as fh:
	for line in fh:
		if line.startswith(">"):
			if header != "":
				sys.exit("Assumes single fasta sequence. Multiple detected.")
			header = line.strip()
		else:
			fasta += line.strip()

with open(args.output, 'w') as out:
	tmp = out.write(header + "\n")
	tmp = out.write(fasta[int(args.start) - 1:] + fasta[:int(args.start) - 1] + "\n")