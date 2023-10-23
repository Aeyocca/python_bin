#!/bin/bash
#Alan E. Yocca
#subset_gb.py

#"randomly" (nothing is truly random) subset a gb file so augustus training is faster

import random
import argparse

parser = argparse.ArgumentParser(prog='PROG')
parser.add_argument('--input', required=True, help='input gff file')
parser.add_argument('--ngenes', type = int, required=True, help='number of genes to subset')
parser.add_argument('--seed', type = int, required=False, default = 10, help='set see for random sampling')
parser.add_argument('--output', required=True, help='output basename for table and graph file')
args = parser.parse_args()

#get a list of random numbers from... wait.. so only way is to loop these twice, maybe can get fancy with it but these files aren't so big

#loop to count how many total genes we have

ntotal = 0
with open(args.input) as fh:
	for line in fh:
		if line.startswith("LOCUS"):
			ntotal += 1

#print("total genes = %s" % (ntotal))
		
if args.ngenes > ntotal:
	sys.exit("ngenes (%s) greater than total in genbank file (%s)" % (args.ngenes, ntotal))

#get random list
#Generate 5 random numbers between 10 and 30
random.seed(10)
randomlist = random.sample(range(0, ntotal), args.ngenes)

output = []

nloop = 0
outwrite = 0
with open(args.input) as fh:
	for line in fh:
		if line.startswith("LOCUS"):
			nloop +=1
			if nloop in randomlist:
				outwrite = 1
				output.append(line)
			else:
				outwrite = 0
		elif outwrite:
			output.append(line)

with open(args.output, "w") as out:
	for line in output:
		out.write(line)