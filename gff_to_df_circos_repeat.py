#!/bin/python
#Alan E. Yocca
#gff_to_df_circos_repeat.py
#03-01-21
#add extra field of splitting by repeat class (third field of gff)


import argparse

parser = argparse.ArgumentParser(prog='PROG')
parser.add_argument('--input', required=True, help='input meta file')
parser.add_argument('--output', required=True, help='output file')
parser.add_argument('--window', required=False, default = 1000000, type = int, help='size of window (non-sliding)')
#parser.add_argument('--all', required=False, action='store_true', help='capture all features, careful. Default is absent, so only captures \'gene\' feature')
#parser.add_argument('--feature', required=False, help='specific feature, for repeat splitting')

args = parser.parse_args()

gene_bedgraph = dict()
#dictionary so we can split by repeat class

with open(args.input) as fh:
	for line in fh:
		if line.startswith("#"):
			continue
		line_array = line.strip().split("\t")
		try:
			gene_bedgraph[line_array[3]].append([line_array[0], line_array[3], line_array[4]])
		except KeyError:
			gene_bedgraph[line_array[3]] = [line_array[0], line_array[3], line_array[4]]

def round_up(number = "", round_value = ""):
	#floor division should work...
	floor = number // round_value
	output = (floor + 1) * round_value
	return output

out_dict = dict()
for repeat in gene_bedgraph.keys():
	#initialize output
	out_dict[repeat] = {}
	for bedgraph in gene_bedgraph[repeat]:
		for coord in range( int(bedgraph[1]), int(bedgraph[2]) + 1):
			#if chromosome not encountered, initialize
			if bedgraph[0] not in out_dict[repeat].keys():
				out_dict[repeat][bedgraph[0]] = {}
			window = ((coord // args.window) + 1) * args.window
			try:
				out_dict[repeat][bedgraph[0]][window] += 1
			except KeyError:
				out_dict[repeat][bedgraph[0]][window] = 1

#window averages here??
#lets do window sums by rounding up to next highest interval of 10k
#hmm this creates an issue of going out of bounds of karyotype
#can fix that in post huh??

#hmm don't need the zeros do I??
print("Outputting")
with open(args.output, 'w') as out:
	out.write("Class\t" + "Chromosome\t" + "chromStart\t" + "chromEnd\t" + "Data\n")
	for repeat in out_dict.keys():
		for chr in out_dict[repeat].keys():
			for coord in out_dict[repeat][chr].keys():
				out.write(str(repeat) + "\t" + str(chr) + "\t" + str(coord - args.window + 1) 
						+ "\t" + str(coord) + "\t" + str(out_dict[chr][coord]) + "\n")




