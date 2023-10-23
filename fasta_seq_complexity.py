#!/bin/python
#fasta_seq_complexity.py

#use biopython
#have the same script for dropping "low" complexity seqs and for spitting them all out
#working with short seqs and on the order of 10^5 seqs so just load everything into memory then process

from Bio.SeqUtils import lcc
import argparse
import sys

def load_fasta(file):
	fasta_dict = dict()
	header = ""
	seq = ""

	with open(args.input) as fh:
		for line in fh:
			if line.startswith(">"):
				#skips first sequence
				if seq != "":
					fasta_dict[header] = seq
				header = line.strip()
				seq = ""
			else:
				seq += line.strip()
	#grab last sequence
	fasta_dict[header] = seq
	return(fasta_dict)

def calc_complexity(fasta_dict):
	complexity_dict = dict()
	for gene in fasta_dict.keys():
		complexity_dict[gene] = lcc.lcc_simp(fasta_dict[gene])
	return(complexity_dict)

def write_complexity_tsv(outfile = str(), fasta_dict = dict(),
							 complexity_dict = dict()):
	with open(outfile, 'w') as out:
		for gene in fasta_dict.keys():
			tmp = out.write(gene + "\t" + str(complexity_dict[gene]) + "\n")

def drop_low_seqs(threshold = str(), fasta_dict = dict(),
							 complexity_dict = dict()):
	dropped = 0
	for gene in fasta_dict.keys():
		if complexity_dict[gene] <= float(threshold):
			del fasta_dict[gene]
			dropped+=1
	print("Removed %s sequences below a threshold of %s" % (dropped,threshold))
	return(fasta_dict)

def write_fasta(outfile = str(), fasta_dict = dict()):
	with open(outfile, 'w') as out:
		for gene in fasta_dict.keys():
			tmp = out.write(">" + gene + "\n" + str(fasta_dict[gene]) + "\n")

if __name__ == "__main__":
	parser = argparse.ArgumentParser(prog='PROG')
	parser.add_argument('--input', required=True, help='input fasta file')
	parser.add_argument('--threshold', required=False, help='if specified, remove seqs with complexity below this threshold. otherwise writes out a tsv of seq complexity')
	parser.add_argument('--output', required=True, help='output file')
	args = parser.parse_args()
	
	fasta_dict = load_fasta(args.input)
	complexity_dict = calc_complexity(fasta_dict)
	if args.threshold != "":
		write_complexity_tsv(outfile = args.output, fasta_dict = fasta_dict,
							 complexity_dict = complexity_dict)
	else:
		fasta_dict = drop_low_seqs(threshold = args.threshold, fasta_dict = fasta_dict,
								   complexity_dict = complexity_dict)
		write_fasta(outfile = args.output, fasta_dict = fasta_dict)
	
	
	
	
	
	
	
	
	
	