#!/bin/python

out_table = dict()
tmp_iter = 0

with open("gala_pav_raw.txt") as fh:
	for line in fh:
		la = line.strip().split("\t")
		#idx 9 onward
		core = "Auxiliary"
		if tmp_iter < 10:
			#print(sum([int(x) for x in la[8:43]]))
			tmp_iter += 1
		if sum([int(x) for x in la[8:43]]) < 33:
			core = "Core"
		if la[1].startswith("\""):
			gene_split = la[1].replace("\"","").split(",")
			for gene in gene_split:
				out_table[gene] = core
		else:
			out_table[la[1]] = core

with open("gala_pav.txt", "w") as out:
	for gene in out_table:
		tmp = out.write(gene + "\t" + out_table[gene] + "\n")




