#!/bin/python
#kmer_diff.py

#take two jellyfish tables and output kmer differences

import argparse
import sys
import pandas as pd

def read_jf(file = ""):
	"""read in column formatted jellyfish dump file and return list of kmers"""
	#function list, but its fun too, no?
	fun_list = []
	with open(file) as fh:
		for line in fh:
			fun_list.append(kmer)
	
	return(fun_list)

def write_kmer_diff(file = "", kmer_list = []):
	"""take a list of kmers and write to file, straightforward"""
	with open(file, 'w') as out:
		for kmer in kmer_list:
			tmp = out.write(kmer + "\n")

if __name__ == "__main__":
	parser = argparse.ArgumentParser(prog='PROG')
	parser.add_argument('--file_one', required=True, help='One jellyfish file, created from jellyfish dump with -c argument. Order does not matter')
	parser.add_argument('--file_two', required=True, help='Other jellyfish file, created from jellyfish dump with -c argument. Order does not matter')
	parser.add_argument('--out_one', required=True, help='file_one specific kmers')
	parser.add_argument('--out_two', required=True, help='file_two specific kmers')
	args = parser.parse_args()
	
	list_one = read_jf(args.file_one)
	list_two = read_jf(args.file_two)
	
	#this cannot be the fastest way.. especially since looping ~0.5Gb twice.. we'll see
	#spec for specific
	list_one_spec = [x for x in list_one if x not in list_two]
	list_two_spec = [x for x in list_two if x not in list_one]
	
	#hmm I guess we have kmers specific to each
	write_kmer_diff(file = args.out_one, kmer_list = list_one_spec)
	write_kmer_diff(file = args.out_two, kmer_list = list_two_spec)
