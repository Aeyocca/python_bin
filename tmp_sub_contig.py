#!/bin/python
#Alan E. Yocca
#tmp_sub_contig.py

import argparse

parser = argparse.ArgumentParser(prog='PROG')
parser.add_argument('--input', required=True, help='input file')
parser.add_argument('--trans', required=True, help='translation file')
parser.add_argument('--output', required=True, help='output file')
#eventually add an argument to define a smaller set of individuals

args = parser.parse_args()

trans_dict = dict()
with open(args.trans) as fh:
	for line in fh:
		line_array = line.strip().split("\t")
		trans_dict[line_array[1]] = line_array[0]

with open(args.input) as fh:
	with open(args.output, "w") as out:
		for line in fh:
			line_array = line.strip().split("\t")
			try:
				tmp = out.write(trans_dict[line_array[0]] + "\t")
			except KeyError:
				tmp = out.write(line_array[0] + "\t")
			for i in range(1,len(line_array) - 1):
				tmp = out.write(line_array[i] + "\t")
			tmp = out.write(line_array[-1] + "\n")
