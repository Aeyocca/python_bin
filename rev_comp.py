#!/bin/python
#Alan E. Yocca
#11-17-22
#rev_comp.py

#import csv
#import sys
import argparse
#import re

parser = argparse.ArgumentParser(prog='PROG')
parser.add_argument('--input', required=True, help='input fasta file')
parser.add_argument('--output', required=True, help='output fasta file')
args = parser.parse_args()

fasta_string = ""
header = ""

#add functionality to do a specific seq but for my use case now its 1 seq per file
with open(args.input) as fh:
	header = next(fh)
	for line in fh:
		if line.startswith(">"):
			sys.exit("Multiple fasta entries in this file, one at a time please")
		fasta_string += line.strip()

def rev_comp(seq_string = ""):
	#naive reverse complement, assuming only ATCG
	#attempting to avoid IUPAC, but no garuntee
	rev_string = seq_string[::-1]
	new_string = rev_string.replace("A","Z").replace("T","Y").replace("C","P").replace("G","Q")
	return(new_string.replace("Z","T").replace("Y","A").replace("P","G").replace("Q","C"))


with open(args.output, "w") as out:
	tmp = out.write(header)
	#write unwrapped, consider wrapping later, not necessary imo
	tmp = out.write(rev_comp(fasta_string) + "\n")




