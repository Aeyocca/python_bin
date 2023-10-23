#!/bin/python
#tmp_check_hits.py
 
import argparse
import sys

parser = argparse.ArgumentParser(prog='PROG')
parser.add_argument('--input', required=True, help='input blast file')
args = parser.parse_args()

hits = 0
with open(args.input) as fh:
	for line in fh:
		la = line.strip().split("\t")
		cns_anno = la[0].split(":")[1].split("-")
		if la[8] == cns_anno[0] and la[9] == cns_anno[1]:
			hits+=1
print("%s exact hits" % (hits))