#!/bin/python

import sys

same_same = 0
same_diff = 0
diff_diff = 0

hap1_chrom = dict()
hap2_chrom = dict()

with open("hap1_chrom_list.txt") as fh:
	for line in fh:
		la = line.strip().split("\t")
		number = la[0].replace("A","")
		hap1_chrom[la[1]] = number
with open("hap2_chrom_list.txt") as fh:
	for line in fh:
		la = line.strip().split("\t")
		number = la[0].replace("B","")
		hap2_chrom[la[1]] = number

#now to get it to output the pairs that are on the same chromosome

with open("hap1_double_pairs.txt") as fh:
	gene = ""
	loop_match = []
	for line in fh:
		la = line.strip().split("\t")
		if la[0] == gene:
			if hap1_chrom[la[0]] == hap2_chrom[la[1]]:
				loop_match.append("SAME")
			else:
				loop_match.append("DIFF")
			if loop_match == ["SAME", "SAME"]:
				same_same += 1
			elif loop_match == ["SAME", "DIFF"]:
				same_diff += 1
			elif loop_match == ["DIFF", "SAME"]:
				same_diff += 1
			elif loop_match == ["DIFF", "DIFF"]:
				diff_diff += 1
			else:
				sys.exit(loop_match)
		else:
			loop_match = []
			gene = la[0]
			if hap1_chrom[la[0]] == hap2_chrom[la[1]]:
				loop_match.append("SAME")
			else:
				loop_match.append("DIFF")

print(same_same)
print(same_diff)
print(diff_diff)

