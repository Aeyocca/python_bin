#!/bin/python
#Alan E. Yocca
#subset_fasta.py

#extract fasta sequences. ignores whitespace

import argparse
import sys

parser = argparse.ArgumentParser(prog='PROG')
parser.add_argument('--fasta', required=True, help='input fasta file')
parser.add_argument('-i','--input', required=False, help=
	'input list of fasta headers, one per line, no carets')
parser.add_argument('-s','--string', required=False, help=
	'comma separated string of fasta headers to extract')
parser.add_argument('--output', help='output file')

args = parser.parse_args()

#specify one or the other

if ((args.input is None and args.string is None) or 
	(args.input is not None and args.string is not None)):
	sys.exit("Specify only one of --input or --string")

header_list = []
if args.input is not None:
	with open(args.input) as fh:
		for line in fh:
			header_list.append(line.strip())

if args.string is not None:
	header_list = args.string.strip().split(",")

out_fasta=dict()
header = ""
keep = 0
with open(args.fasta) as fh:
	for line in fh:
		if line.startswith(">"):
			#new fasta header
			header = line.strip().replace(">","").split(" ")[0]
			if header in header_list:
				out_fasta[header] = ""
				keep = 1
			else:
				keep = 0
		else:
			if keep:
				out_fasta[header] += line.strip()

n_seq = len(list(out_fasta.keys()))
if n_seq == 0:
	sys.exit("Unable to find any matching sequences")

print("Extracted %s sequences" % (n_seq))
with open(args.output, "w") as out_fh:
	for key in out_fasta:
		tmp = out_fh.write(">" + key + "\n")
		tmp = out_fh.write(out_fasta[key] + "\n")















