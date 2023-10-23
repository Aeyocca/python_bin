#!/bin/python
#liftover_goterms.py
#dummy script to call GO terms based on orthologs to arabidopsis
#Alan E. Yocca
#05-12-2023

import sys
import argparse

parser = argparse.ArgumentParser(prog='PROG')
parser.add_argument('--trans', required=True, help='translation file of orthologs')
parser.add_argument('--output', required=True, help='output bed file')

args = parser.parse_args()

#hard code the arabidopsis go term file just assume cwd

#read in dictionary of arabidopsis GO terms

at_dict = dict()
lc = 0
with open("gene_association.tair") as fh:
	for line in fh:
		lc += 1
		if line.startswith("!"):
			continue
		la = line.strip().split("\t")
		if la[1] in at_dict.keys():
			at_dict[la[1]].append(la[4])
			
		else:
			at_dict[la[1]] = [la[4]]

ortho_dict = dict()
with open(args.trans) as fh:
	for line in fh:
		#assume each gene in left column represented once
		la = line.strip().split("\t")
		ortho_dict[la[0]] = la[1]

with open(args.output, 'w') as out:
	for gene in ortho_dict.keys():
		try:
			tmp = out.write(gene + "\t" + ",".join(at_dict[ortho_dict[gene]]) + "\n")
		except KeyError:
			#no GO term for this arabidopsis gene
			continue


