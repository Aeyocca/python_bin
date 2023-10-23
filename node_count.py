#!/bin/python
#Alan E. Yocca
#10-21-20
#node_count.py

#count nodes between two taxa in newick file

import sys
import argparse
import re
import itertools

parser = argparse.ArgumentParser(prog='PROG')
parser.add_argument('--input', required=True, help='input newick file')
parser.add_argument('--output', required=True, help='output file')

args = parser.parse_args()


def count_nodes(tree_file = "", taxa = ""):

	#read file into string
	newick=""

	with open(tree_file) as fh:
		for line in fh:
			newick=line.strip()

	#extract string between two
	taxa_array = taxa.split(",")
	
	#add for combinations later to get 3 / 4 taxa counts
	combinations = list(itertools.combinations(range(len(taxa_array)), 2))
	node_count_list = []
	cannot_extract = 0
	for combo in combinations:
		between = ""
		between = re.search( taxa_array[combo[0]] + '(.*)' + taxa_array[combo[1]], newick)
		try:
			tmp = len(between.group(1))
		except:
			between = re.search( taxa_array[combo[1]] + '(.*)' + taxa_array[combo[0]], newick)

		#print(taxa_array[combo[0]])
		#print(taxa_array[combo[1]])
		#print(between.group(1))
		try:
			forward_parenthesis = between.group(1).count('(')
		except:
			#print(taxa)
			#print(tree_file)
			#print(combo)
			#sys.exit()
			cannot_extract = 1
		if cannot_extract:
			print("Missing a taxon in %s" % (tree_file))
			break
		reverse_parenthesis = between.group(1).count(')')
		node_count_list.append(forward_parenthesis + reverse_parenthesis + 1)
	try:
		return ([min(node_count_list), max(node_count_list)])
	except ValueError:
		return (["NA","NA"])

if __name__ == "__main__":

	print("Counting")
	output_list = []
	with open(args.input) as fh:
		for line in fh:
			line_array = line.strip().split("\t")
			gene = line_array[0].split("_")[1]
			taxa = line_array[1]
			cns = line_array[0]
			
			min_max = count_nodes(tree_file="02_trees/RAxML_bipartitionsBranchLabels."
				 + gene + "_aln", taxa = taxa)
			output_list.append([cns, gene, taxa, min_max[0], min_max[1]])

	print("outputting")
	with open(args.output, "w") as out:
		for line in output_list:
			for i in range(0,len(line) - 1):
				tmp = out.write(str(line[i]))
				tmp = out.write("\t")
			tmp = out.write(str(line[i+1]))
			tmp = out.write("\n")

