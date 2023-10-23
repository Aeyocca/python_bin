#!/bin/python
#collect_tpm_rge.py

#dummy script to output expression table

out = dict()
acc_list = []

with open("sraid_list.txt") as fh1:
	for l1 in fh1:
		acc = l1.strip()
		acc_list.append(acc)
		out[acc] = dict()
		with open("mapped/" + acc + "_Mdom_H1.quant/abundance.tsv") as fh2:
			next(fh2)
			for l2 in fh2:
				la = l2.strip().split("\t")
				out[acc][la[0]] = la[4]

with open("SRP065793_tpm_table.txt", 'w') as ofh:
	tmp = ofh.write("Gene\t")
	for lib in acc_list[:-1]:
		tmp = ofh.write(lib + "\t")
	tmp = ofh.write(acc_list[-1] + "\n")
	for gene in out[acc_list[0]]:
		tmp = ofh.write(gene + "\t")
		for lib in acc_list[:-1]:
			tmp = ofh.write(out[lib][gene] + "\t")
		tmp = ofh.write(out[acc_list[-1]][gene] + "\n")