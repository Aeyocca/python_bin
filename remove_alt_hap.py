#!/bin/python
#Alan E. Yocca
#11-14-22
#remove_alt_hap.py

#import csv
#import sys
import argparse
#import re

parser = argparse.ArgumentParser(prog='PROG')
parser.add_argument('--input', required=True, help='input fasta file')
parser.add_argument('--output', required=True, help='output fasta file')
parser.add_argument('--left_bound', required=True, type = int, help='left boundary of alt hap (1-indexed), exclusive')
parser.add_argument('--center_point', required=True, type = int, help='center to break (1-indexed)')
parser.add_argument('--right_bound', required=True, type = int, help='right boundary of alt hap (1-indexed), exclusive')

args = parser.parse_args()

#given breakpoints, break fasta at midpoint, remove smaller scaffold side, pad with 10kb N

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

with open(args.output, "w") as out:
	tmp = out.write(header)
	if args.center_point < len(fasta_string) / 2:
		tmp = out.write(fasta_string[:args.center_point] + "N"*10000 + \
			fasta_string[args.right_bound:] + "\n")
	else:
		tmp = out.write(fasta_string[:args.left_bound] + "N"*10000 + \
			fasta_string[args.center_point:] + "\n")



