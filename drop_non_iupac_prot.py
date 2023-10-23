#!/bin/python
#drop_non_iupac_prot.py

import sys
import re
import argparse

parser = argparse.ArgumentParser(prog='PROG')
parser.add_argument('--input', required=True, help='input protein fasta file')
parser.add_argument('--output', required=True, help='output protein fasta file')

args = parser.parse_args()

#load fasta into dict
header = ""
fasta_seq = dict()
with open(args.input) as fh:
	for line in fh:
		if line.startswith(">"):
			fasta_seq[line] = ""
			header = line
		else:
			fasta_seq[header] += line.strip()

#as we output check for non-iupac characters in seq string
iupac = "ACDEFGHIKLMNPQRSTVWY"
non_iupac = "BJOUZ"
matched = 0
with open(args.output, "w") as out:
	for seq in fasta_seq.keys():
		#how to check, if any of the string contains these non iupac protein codes
		if re.search("[BJOUZ]",fasta_seq[seq]):
			matched += 1
			continue
		else:
			#seq should still have line return
			tmp = out.write(seq)
			tmp = out.write(fasta_seq[seq] + "\n")
			
print("Dropped %s seqs" % (matched))