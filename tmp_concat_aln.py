#!/bin/python
#tmp_concat_aln.py
#01-31-23
#Alan E. Yocca
#temporary because some of my headers are weird for this project
#can generalize later


import sys
import re
import argparse

parser = argparse.ArgumentParser(prog='PROG')
parser.add_argument('--input', required=True, help='text file of fasta files aligned')
parser.add_argument('--output', required=True, help='output fasta alignment (not wrapped)')

args = parser.parse_args()

#hmmmm lets feed in a list of fasta files like a txt file because there are 761 of them

out_dict = {"Rosa_chinensis" : "",
			"Prunus_cerasus" : "",
			"Malus_domestica" : "",
			"Fragaria_annanassa" : ""}
#initialize headers because I'm lazy

with open(args.input) as fh_list:
	for line in fh_list:
		with open(line.strip()) as fh:
			header = ""
			for line in fh:
				if line.startswith(">"):
					header = line.strip().split(">")[1].split(" ")[0]
				else:
					out_dict[header] = out_dict[header] + line.strip()


with open(args.output,'w') as out:
	for spec in out_dict.keys():
		tmp = out.write(">" + spec + "\n")
		tmp = out.write(out_dict[spec])
		tmp = out.write("\n")	


