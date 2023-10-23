#!/bin/python
#Alan E. Yocca
#drop_short_seqs.py

import argparse

parser = argparse.ArgumentParser(prog='PROG')
parser.add_argument('--input', required=True, help='input fasta file')
parser.add_argument('--length', default = 50000, type = int, required=False, help='length of seq to drop inclusive')
parser.add_argument('--trashbin', required=False, help='if specified, where to put filtered out contigs')
parser.add_argument('--output', required=True, help='output fasta file')
args = parser.parse_args()


length = 0
seq = []

with open(args.input) as fh:
	with open(args.output, "w") as out:
		#so we don't have to check if its the first seq 
		header = next(fh)
		for line in fh:
			if line.strip().startswith(">"):
				#calc length
				seq_lens = [len(x) for x in seq]
				sum_len = sum(seq_lens)
				if sum_len > args.length:
					#write out
					tmp = out.write(header)
					for entry in seq:
						tmp = out.write(entry + "\n")
				else:
					if args.trashbin is not None:
						with open(args.trashbin, "a") as trash:
							tmp = trash.write(header)
							for entry in seq:
								tmp = trash.write(entry + "\n")
				#reset
				header = line
				seq = []
			else:
				seq.append(line.strip())
