#!/bin/python
#orthogroup_member_counts.py
#dummy_script to read in "Orthogroups_1.GeneCount.csv" and spit out bar graph data
#Alan E. Yocca
#02-17-23

import sys
out = []
tmp_out = []

with open("Orthogroups_1.GeneCount.csv") as fh:
	#get nacc from header
	header = next(fh)
	nacc = len(header.split("\t")) - 2
	out = [0 for i in range(nacc)]
	for line in fh:
		la = line.split("\t")
		loop_nacc = len([x for x in [int(j) for j in la[1:-1]] if x >= 1 ])
		out[loop_nacc - 1] += 1
		
		if loop_nacc == 2:
			tmp_out.append(line)

print(out)
with open("tmp.txt", "w") as out:
	for line in tmp_out:
		tmp = out.write(line)


