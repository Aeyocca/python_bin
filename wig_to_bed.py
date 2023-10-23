#!/bin/python
#wig_to_bed.py
#script to create bedgraph file from wig file
#idea is so we can use bedtools intersect to filter out CDS
#Alan E. Yocca
#04-25-2023

import sys
import argparse

parser = argparse.ArgumentParser(prog='PROG')
parser.add_argument('--wig', required=True, help='input wig file')
parser.add_argument('--chrom', required=False, help='def chromosome, for Alan\'s use')
parser.add_argument('--output', required=True, help='output bed file')

args = parser.parse_args()

out_bed = []
position = 0
#I wonder if I'll ever figure out the chromosome thing
#do they need to match for bedtools intersect?
chrom = ""

with open(args.wig) as fh:
	for line in fh:
		if line.startswith("fixedStep"):
			#new block
			#fixedStep chrom=(null) start=104833 step=1
			#check step,
			la= line.strip().split(" ")
			if la[3] != "step=1":
				sys.exit("step not equal to 1, figure something out Alan")
			
			if args.chrom != "":
				chrom = args.chrom
			else:
				chrom = la[1].replace("chrom=","")
			#minus 1 because increment at the start of else statement
			position = int(la[2].replace("start=","")) - 1
		else:
			position += 1
			out_bed.append([chrom, position, position, line.strip()])

with open(args.output, 'w') as out:
	for line in out_bed:
		for col in line[:-1]:
			tmp = out.write(str(col) + "\t")
		tmp = out.write(str(line[-1]) + "\n")

