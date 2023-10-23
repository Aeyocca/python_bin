#!/bin/python
#Alan E. Yocca
#02-23-21
#rename_maker_fasta.py
#take maker generated gene id map and replace gene names in gff, trans, prot

import argparse
import re

parser = argparse.ArgumentParser(prog='PROG')
parser.add_argument('--input', required=True, help='input meta file')
parser.add_argument('--output', required=True, help='output file')
parser.add_argument('--map', required=True, help='maker id mapping file')

args = parser.parse_args()

#load map into dictionary
map_dict = {}
with open(args.map) as fh:
	for line in fh:
		line_array = line.strip().split("\t")
		map_dict[line_array[0]] = line_array[1]

#hmm loop whilst outputting? is that a thing?
with open(args.output, 'w') as out:
	with open(args.input) as fh:
		for line in fh:
			line_out = line
			if line.startswith(">"):
				gene_name = line.strip().split(">")[1].split(" ")[0]
				if gene_name in map_dict.keys():
					#replace everywhere
					line_out = line.replace(gene_name, map_dict[gene_name])
			out.write(line_out)
