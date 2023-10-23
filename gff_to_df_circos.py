#!/bin/python
#Alan E. Yocca
#02-19-21
#gff_to_df_circos.py

import argparse

parser = argparse.ArgumentParser(prog='PROG')
parser.add_argument('--input', required=True, help='input meta file')
parser.add_argument('--output', required=True, help='output file')
parser.add_argument('--window', required=False, default = 10000, type = int, help='size of window (non-sliding)')
parser.add_argument('--all', required=False, action='store_true', help='capture all features, careful. Default is absent, so only captures \'gene\' feature')

args = parser.parse_args()

gene_bedgraph = []
with open(args.input) as fh:
	for line in fh:
		if line.startswith("#"):
			continue
		line_array = line.strip().split("\t")
		if args.all:
			gene_bedgraph.append([line_array[0], line_array[3], line_array[4]])
		elif line_array[2] == "gene":
			gene_bedgraph.append([line_array[0], line_array[3], line_array[4]])

def round_up(number = "", round_value = ""):
	#floor division should work...
	floor = number // round_value
	output = (floor + 1) * round_value
	return output


out_dict = dict()
for gene in gene_bedgraph:
	for coord in range( int(gene[1]), int(gene[2]) + 1):
		if gene[0] not in out_dict.keys():
			out_dict[gene[0]] = {}
		window = ((coord // args.window) + 1) * args.window
		try:
			out_dict[gene[0]][window] += 1
		except KeyError:
			out_dict[gene[0]][window] = 1

#window averages here??
#lets do window sums by rounding up to next highest interval of 10k
#hmm this creates an issue of going out of bounds of karyotype
#can fix that in post huh??

#hmm don't need the zeros do I??
print("Outputting")
with open(args.output, 'w') as out:
	out.write("Chromosome\t" + "chromStart\t" + "chromEnd\t" + "Data\n")
	for chr in out_dict.keys():
		for coord in out_dict[chr].keys():
			out.write(str(chr) + "\t" + str(coord - args.window + 1) + "\t" + str(coord) + 
						"\t" + str(out_dict[chr][coord]) + "\n")




