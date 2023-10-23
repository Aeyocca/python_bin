#!/bin/python
#Alan E. Yocca
#extract_trn_chloro.sh

#semi dummy script to get trnL-trnF sequences from chloroplast consensus sequences

import argparse
import re
import sys

parser = argparse.ArgumentParser(prog='PROG')
parser.add_argument('--input', required=True, help='input fasta file')
parser.add_argument('--blast', required=True, help='blast file in m8 format')
parser.add_argument('--output', help='extracted gene file')
parser.add_argument('--extend', required=False, default=0, type=int, help='extend the blast hit')

args = parser.parse_args()

#load chloroplast into string

accession=""
seq=""
with open(args.input) as fh:
	for line in fh:
		if line.startswith(">"):
			#skip? how to get accession, should I just take from here?
			#yea get here
			#funky way but works to remove carat, oh wait we can keep carrat
			accession = line.strip()
		else:
			seq += line.strip()

#get boundaries from blast file
#hmmm how to see if diff strand, well
#if one strand is diff they all should be

#just going to assume top hit is best
coords = []
with open(args.blast) as fh:
	for line in fh:
		line_array = line.strip().split("\t")
		coords = [line_array[6], line_array[7]]
		break

with open(args.output, "w") as out_fh:
	out_fh.write(accession + "_trnL_trnF\n")
	out_fh.write(seq[(int(coords[0]) - 1 - args.extend):int(coords[1]) + args.extend] + "\n")















