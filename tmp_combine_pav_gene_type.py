#!/bin/python

out_table = dict()
tmp_iter = 0

with open("gala_pav.txt") as fh:
	for line in fh:
		la = line.strip().split("\t")
		out_table[la[0]] = [la[1]]
		
no_gene_type = 0
no_pav = 0
with open("Gala_hapA.gene_type") as fh:
	for line in fh:
		la = line.strip().split("\t")
		gene = la[0].replace("A","g")
		try:
			out_table[gene].append(la[1])
		except KeyError:
			no_pav += 1

with open("gala_pav_dup_type.txt", "w") as out:
	for gene in out_table.keys():
		if len(out_table[gene]) < 2:
			no_gene_type += 1
			continue
		tmp = out.write(gene + "\t" + out_table[gene][0] + "\t" + out_table[gene][1] + "\n")

print("no pav: %s" % (no_pav))
print("no gene_type: %s" % (no_gene_type))