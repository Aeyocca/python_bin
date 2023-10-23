#!/bin/python
#tmp_rename_gff.py

import sys

#use fasta indexes to rename gff to match contig names in the bedgraph files

list_one = []
list_two = []
#double checked they have the same lengths in assuming the same orders
with open("Am_tri.pan.renamed.fa.fai") as fh:
	for line in fh:
		list_one.append(line.strip().split("\t")[0])

with open("Am_tri.pan.fa.fai") as fh:
	for line in fh:
		list_two.append(line.strip().split("\t")[0])

trans_dict = dict()
for i in range(len(list_one)):
	trans_dict[list_one[i]] = list_two[i]

with open("Am_tri.pan.sorted.fixed.renamed.gff") as fh:
	with open("alan.gff", "w") as out:
		for line in fh:
			if line.startswith("#"):
				continue
			la = line.strip().split("\t")
			tmp = out.write(trans_dict[la[0]] + "\t")
			for field in la[1:-1]:
				tmp = out.write(field + "\t")
			tmp = out.write(la[-1] + "\n")
