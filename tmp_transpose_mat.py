#!/bin/python
#Alan E. Yocca
#tmp_transpose_mat.py


ncol = 35
tmat = [[] for i in range(ncol)]
with open("draper_bbpan_snp_mat_maf_05.txt") as fh:
	for line in fh:
		line_array = line.strip().split("\t")
		for i in range(len(line_array)):
			tmat[i].append(line_array[i])

print("Finished loading mat")

with open("draper_bbpan_snp_mat_mat_05_trans.txt", 'w') as out:
	for line in tmat:
		for i in range(len(line) - 1):
			tmp = out.write(line[i])
			tmp = out.write("\t")
		tmp = out.write(line[-1])
		tmp = out.write("\n")