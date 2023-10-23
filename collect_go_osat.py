#!/bin/python
#Alan E. Yocca
#08-12-20
#collect_go_osat.py
#read rice gff and extract go terms for each gene

import re

go_dict = dict()

with open("transcripts.gff") as fh:
	for line in fh:
		if re.match("^#",line):
			continue
		line_array = line.strip().split("\t")
		if re.match("mRNA",line_array[2]):
			transcript = re.search("ID=([a-zA-Z0-9\-]*);",line_array[8]).group(1)
			#cg is to map to other values for machine learning script
			#Also to combine them all across transcripts
			cg = transcript.split("-")[0].replace("t","g")
			
			go_terms = re.findall("\((GO:[0-9]*)\)",line_array[8])
			
			try:
				#add non-repeating GO term
				#not sure why list comprehension isn't working but oh whale
				for term in go_terms:
					if term not in go_dict[cg]:
						go_dict[cg].append(term)
			except KeyError:
				go_dict[cg] = go_terms

#output
#hmm r package will compute 
with open("osat_go_list.tsv", "w") as output:
	gene_count = 1
	for gene in go_dict.keys():
		if gene_count == 1:
			#add header
			gene_count = 0
			output.write("Gene\tGO_list\n")
		output.write(gene)
		output.write("\t")
		for i in range(0,(len(go_dict[gene]) - 1)):
			output.write(go_dict[gene][i])
			output.write(",")
		output.write(go_dict[gene][(len(go_dict[gene]) - 1)])
		output.write("\n")

			
			
			