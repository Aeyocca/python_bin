#!/bin/python
#coll_to_trans.py

#read mcscanx, create translation of first column to best second column hit

import argparse
import sys

parser = argparse.ArgumentParser(prog='PROG')
parser.add_argument('--input', required=True, help='input mcscanx bed file')
parser.add_argument('--hap', required=True, help='text to match haplotype with')
parser.add_argument('--output', required=True, help='output file')
args = parser.parse_args()

#loopy loops, want translation eh just make it two dicts
score_dict = dict()
trans_dict = dict()
aln_score = ""

#think I have to back up here, make sure flip is also unique... yupp 
print("Hard coded for Apple at the moment to keep haplotypes straight")

with open(args.input) as fh:
	for line in fh:
		if "Alignment" in line:
			aln_score = line.split("score=")[1].split(" ")[0]
		else:
			if line.startswith("#"):
				continue
			la = line.strip().split("\t")
			ref_gene = la[1]
			que_gene = la[2]
			if args.hap not in ref_gene:
				tmp = ref_gene
				ref_gene = que_gene
				que_gene = tmp
			
			if ref_gene in score_dict.keys():
				if aln_score > score_dict[ref_gene]:
					#update best ortholog
					trans_dict[ref_gene] = que_gene
			else:
				#first time hitting this gene
				trans_dict[ref_gene] = que_gene
				score_dict[ref_gene] = aln_score

#lets sort these
myKeys = list(trans_dict.keys())
myKeys.sort()
sorted_dict = {i: trans_dict[i] for i in myKeys}

with open(args.output, 'w') as out:
	for ortho_a in sorted_dict.keys():
		tmp = out.write(ortho_a + "\t" + sorted_dict[ortho_a] + "\n")
