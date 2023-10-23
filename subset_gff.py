#!/bin/python
#subset_gff.py
#Alan E. Yocca
#01-04-2022

import argparse
import re

parser = argparse.ArgumentParser()
parser.add_argument("--gff", required=True, help = "GFF3 formatted file")
parser.add_argument("--sub_list", required=True, help = "file where every line is a gene to subset")
parser.add_argument("--output", required=True, help = "output file")
args = parser.parse_args()

#maybe can just be a list?
sub_list = []

#stripped mRNA-1 off
with open(args.sub_list) as fh:
	for line in fh:
		sub_list.append(line.strip())

out_list = []
with open(args.gff) as fh:
	for line in fh:
		#some regex matching.. split out gene name
		#what to do if mRNA..
		#think I have to assume they are all mRNA 1
		#worst case scenario we extract them all
		
		if line.startswith("#"):
			continue
		
		line_array = line.strip().split("\t")
		
		#extract gene name
		gene = line_array[8].split("ID=")[1].split("-mRNA")[0]
		
		if gene in sub_list:
			out_list.append(line)


with open(args.output, "w") as out_fh:
	for line in out_list:
		tmp = out_fh.write(line)

#loop list of genes to subset from gff file
#assuming list of subset is small
#so fastest to get subset into dictionary