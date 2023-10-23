#!/bin/python
#tmp_mosdepth_to_pav.py
#Alan E. Yocca
#02-02-23

import sys
import re
import argparse

parser = argparse.ArgumentParser(prog='PROG')
parser.add_argument('--input', required=True, help='input modepth thresholds.bed file')
parser.add_argument('--output', required=True, help='output table of uniq field 4s')

args = parser.parse_args()

#dict key = seq, length, covered length

cov_gff = dict()

with open(args.input) as fh:
	##skip first
	print("skipping first line")
	first_line = fh.readline()
	#dont skip if subsetting? naw just remember to add it back in
	for line in fh:
		la = line.strip().split("\t")
		if la[3] not in cov_gff.keys():
			cov_gff[la[3]] = {"Length" : int(la[2]) - int(la[1]) + 1,
								"Covered" : int(la[4])}
		else:
			cov_gff[la[3]]["Length"] += int(la[2]) - int(la[1]) + 1
			cov_gff[la[3]]["Covered"] += int(la[4])


present = 0
total = 0
cov_threshold = 0.05

with open(args.output, 'w') as out:
	for gene in cov_gff.keys():
		total += 1
		if cov_gff[gene]["Covered"] / cov_gff[gene]["Length"] >= cov_threshold:
			present += 1
			tmp = out.write(gene + "\t1\n" )
		else:
			tmp = out.write(gene + "\t0\n" )


print("%s present out of %s" % (present, total))

#welp no I want to reproduce it with a different genome / annotation