#!/bin/python
#tmp_pad_genes.py

import argparse
import sys

parser = argparse.ArgumentParser(prog='PROG')
parser.add_argument('--input', required=True, help='input mcscanx bed file')
parser.add_argument('--output', required=True, help='output file')
args = parser.parse_args()

#loopy loops
out_list = []

with open(args.input) as fh:
	for line in fh:
		la = line.strip().split("\t")
		#Maldo_hc_v1a1_chr14A_g11894_t1
		#52583 all 5 digits
		#should be: Maldo_hc_v1a1_ch10B_g00002_t1
		#is Maldo_hc_v1a1_chr11B_g3436_t1
		
		gsplit = la[3].split("g")
		under_split = gsplit[1].split("_")
		nz = 5 - len(under_split[0])
		pad = ""
		for i in range(nz):
			pad = pad + "0"
		
		gene_name = gsplit[0] + "g" + pad + gsplit[1]
		
		out_list.append([la[0],la[1],la[2],gene_name])

with open(args.output, 'w') as out:
	for line in out_list:
		for col in line[:-1]:
			tmp = out.write(col + "\t")
		tmp = out.write(line[-1] + "\n")
