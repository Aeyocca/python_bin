#!/bin/python
#maf_drop_self.py
 
import argparse
import sys

parser = argparse.ArgumentParser(prog='PROG')
parser.add_argument('--input', required=True, help='input maf file')
parser.add_argument('--spec', required=True, help='ref species')
parser.add_argument('--chrom', required=True, help='ref chrom, will drop all that dont match')
parser.add_argument('--output', required=True, help='output bed file')
args = parser.parse_args()

#pretty much just look for lines that start with
#s       Eriobotrya_japonica.LG12
#and have the right species but wrong chromosome and drop them
#then convert to fasta and see if its any different

out_list = []

with open(args.input) as fh:
	for line in fh:
		if line.startswith("s"):
			la = line.strip().split("\t")
			if la[1].startswith(args.spec) and args.chrom not in la[1]:
				#don't print
				continue
		#now printing out all these should capture everything else
		out_list.append(line)

with open(args.output, 'w') as out:
	for line in out_list:
		tmp = out.write(line)