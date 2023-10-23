#!/bin/python
#ambiguous_base_bed.py

#dummy script to output bedgraph of ambiguous bases

import argparse

parser = argparse.ArgumentParser(prog='PROG')
parser.add_argument('--input', required=True, help='input genome fasta')
parser.add_argument('--output', required=True, help='output bed file')
args = parser.parse_args()

#dictionary of lists
out_dict = dict()

header = ""
coord = 0
prev_base = ""
group_start = 0
with open(args.input) as fh:
	for line in fh:
		if line.startswith(">"):
			#complete last interval if last fasta entry ended with "N"
			if prev_base == "N":
				out_dict[header].append([group_start,coord - 1])
			header = line.strip().replace(">","").split(" ")[0]
			out_dict[header] = []
			coord = 0
			prev_base = ""
		else:
			for base in line.strip():
				coord += 1
				if base == "N":
					if prev_base != "N":
						group_start = coord
						#else, do nothing, keep going
				elif prev_base == "N":
					#write out the group
					#print("appending")
					out_dict[header].append([group_start,coord - 1])
				prev_base = base

#get straggler
if prev_base == "N":
	out_dict[header].append([group_start,coord - 1])

with open(args.output, 'w') as out:
	for header in out_dict.keys():
		for interval in out_dict[header]:
			tmp = out.write(header + "\t" + str(interval[0]) + "\t" + 
							str(interval[1]) + "\n")

