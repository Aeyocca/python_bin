#!/bin/python
#collect_tpm.py

#dummy script to output expression table


#libaries
lib_string = ["Budding_Leaves_clean", "Expanding_Leaves_clean", "Flower_Buds_clean",
				"Fruitlet_Stage_1_clean", "Fruitlet_Stage_2_clean", 
				"Open_Buds_clean", "QI_Buds_clean"]

#hap1
h1_out = dict()

for base in lib_string:
	#initialize
	h1_out[base] = dict()
	with open(base + "_hap1/abundance.tsv") as fh:
		next(fh)
		for line in fh:
			la = line.strip().split("\t")
			h1_out[base][la[0]] = la[3]
			
#hmmm this should work

with open("Pyrus_tmp_table.txt", 'w') as out:
	tmp = out.write("Gene\t")
	for lib in lib_string[:-1]:
		tmp = out.write(lib + "\t")
	tmp = out.write(lib_string[-1] + "\n")
	for gene in h1_out[lib_string[0]]:
		tmp = out.write(gene + "\t")
		for lib in lib_string[:-1]:
			tmp = out.write(h1_out[lib][gene] + "\t")
		tmp = out.write(h1_out[lib_string[-1]][gene] + "\n")