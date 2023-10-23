#!/bin/python
#Alan E. Yocca
#11-14-22
#invert_seq.py

import csv
import sys
import argparse
import re

parser = argparse.ArgumentParser(prog='PROG')
parser.add_argument('--input', required=True, help='input fasta file')
parser.add_argument('--output', required=True, help='output fasta file')
parser.add_argument('--left_bound', required=True, type = int, help='left boundary to invert (1-indexed), exclusive')
parser.add_argument('--right_bound', required=True, type = int, help='right boundary to invert (1-indexed), exclusive')

args = parser.parse_args()

#given breakpoints, break fasta, invert between breakpoints and insert 10kb Ns on either end

#3-9
#AAATTCGCTG
#AAA     TG

#not absurd to load the fasta sequence into a string, right?
#will that take a while with a wrapped fasta?
#nope, super quick

fasta_string = ""
header = ""

with open(args.input) as fh:
	header = next(fh)
	for line in fh:
		if line.startswith(">"):
			sys.exit("Multiple fasta entries in this file, one at a time please")
		fasta_string += line.strip()

#fix if left/right are switched
if args.right_bound < args.left_bound:
	tmp = args.left_bound
	args.left_bound = args.right_bound
	args.right_bound = tmp

#cut if out of bounds? no do try except when extracting roi.
#region of interest
try:
	roi = fasta_string[args.left_bound:(args.right_bound - 1)]
except IndexError:
	print(args.left_bound)
	print(args.right_bound)
	sys.exit('boundaries may exceed size of scaffold')

#left_scaff = fasta_string[:args.left_bound]
#right_scaff = fasta_string[args.right_bound:]

#do I have to reverse complement it as well? is it the 10kb break? 
#Also, are we bein inclusive with breaks?

def rev_comp(seq_string = ""):
	#naive reverse complement, assuming only ATCG
	#attempting to avoid IUPAC, but no garuntee
	rev_string = seq_string[::-1]
	new_string = rev_string.replace("A","Z").replace("T","Y").replace("C","P").replace("G","Q")
	return(new_string.replace("Z","T").replace("Y","A").replace("P","G").replace("Q","C"))

out_string = fasta_string[:args.left_bound] + "N"*10000 + \
				rev_comp(roi) + "N"*10000 + fasta_string[(args.right_bound - 1):]

with open(args.output, "w") as out:
	tmp = out.write(header)
	#write unwrapped, consider wrapping later, not necessary imo
	tmp = out.write(out_string + "\n")




