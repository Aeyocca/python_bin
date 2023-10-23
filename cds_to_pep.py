#!/bin/python
#cds_to_pep.py

#use biopython

from Bio.Seq import translate
import argparse
import sys

parser = argparse.ArgumentParser(prog='PROG')
parser.add_argument('--input', required=True, help='input cds file')
parser.add_argument('--output', required=True, help='input pep file')
args = parser.parse_args()

#loopy loops
fasta_dict = dict()
header = ""
seq = ""

with open(args.input) as fh:
	for line in fh:
		if line.startswith(">"):
			if seq != "":
				fasta_dict[header] = translate(seq)
			header = line.strip()
			seq = ""
		else:
			seq += line.strip()

#straggler
fasta_dict[header] = translate(seq)

with open(args.output, 'w') as out:
	for header in fasta_dict.keys():
		tmp = out.write(header + "\n")
		tmp = out.write(fasta_dict[header] + "\n")
