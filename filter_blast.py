#!/bin/python
#filter_blast.py

#to start, filter blast m8 format by best query hit, options for bit score and evalue

import argparse
import sys

def filter_blast(file = "", col = int(), max_value = True, keep_all = True):
	blast_dict = dict()
	#key is query, value is just the list
	#if max_value we want the hit with the highest value in the specificed column
	#so replace if difference is positive
	#if we want to keep multiple hits, have to change the way we store, 
	#value is now list of lists
	
	with open(file) as fh:
		#la for line array, parse first line to check format
		la = fh.readline().strip().split("\t")
		blast_dict[la[0]] = [la]
		if len(la) != 12:
			sys.exit("Blast file does not appear to be in m8 format")
		
		for line in fh:
			la = line.strip().split("\t")
			#if exists, overwrite if this hit is better
			#might be a faster way but these files arent huge
			if la[0] in blast_dict.keys():
				#how to change operator? check sign of difference
				#in cases of equal hits, do not overwrite
				if float(blast_dict[la[0]][0][col]) - float(la[col]) > 0 and not max_value:
					blast_dict[la[0]] = [la]
				elif float(blast_dict[la[0]][0][col]) - float(la[col]) < 0 and max_value:
					blast_dict[la[0]] = [la]
				elif float(blast_dict[la[0]][0][col]) == float(la[col]) and keep_all:
					blast_dict[la[0]].append(la)
			else:
				#if does not exist, add
				blast_dict[la[0]] = [la]
	return(blast_dict)

def write_blast(file = "", blast_dict = dict()):
	with open(file, 'w') as out:
		for query in blast_dict.keys():
			for line in blast_dict[query]:
				for col in line[:-1]:
					tmp = out.write(col + "\t")
				tmp = out.write(line[-1] + "\n")
			

if __name__ == "__main__":
	parser = argparse.ArgumentParser(prog='PROG')
	parser.add_argument('--input', required=True, help='input fasta file')
	parser.add_argument('--column', required=False, default = 11, help='column of blast m8 format, 1-indexed. bitscore == 12, e-value == 11')
	parser.add_argument('--max_value', required=False, action='store_true', help='whether to keep the largest or smallest value. if smallest desired (eg evalue) then do not specify')
	parser.add_argument('--keep_all', required=False, action='store_true', help='keep all of the best hits if equal')
	parser.add_argument('--output', required=True, help='output file')
	args = parser.parse_args()
	
	if args.max_value:
		print("Keeping maximun value in column %s" % (args.column))
	else:
		print("Keeping minimum value in column %s" % (args.column))

	blast_dict = filter_blast(file = args.input, col = int(args.column) - 1, 
							  max_value = args.max_value, keep_all = args.keep_all)
	write_blast(file = args.output, blast_dict = blast_dict)
