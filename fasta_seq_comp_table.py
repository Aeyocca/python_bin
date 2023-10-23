#!/bin/python
#fasta_seq_comp_table.py

#take in fasta, output 
import argparse
import sys

parser = argparse.ArgumentParser(prog='PROG')
parser.add_argument('--fasta', required=True, help='input gene bed file')
parser.add_argument('--output', required=True, help='output bed file')
args = parser.parse_args()

#whats our output look like?
#chrom A T C G GC_per AT_per A_per T_per C_per G_per Seq_length
#yea do all these here


out_dict = dict()
fasta_dict = dict()
header = ""
seq = ""

def add_seq(out_dict, seq, header):
	out_dict[header] = dict()
	
	a_count = seq.count("A")
	c_count = seq.count("C")
	t_count = seq.count("T")
	g_count = seq.count("G")
	seq_length = len(seq)
	
	out_dict[header]["A"] = a_count
	out_dict[header]["T"] = t_count
	out_dict[header]["C"] = c_count
	out_dict[header]["G"] = g_count
	out_dict[header]["GC_per"] = (c_count + g_count) / seq_length
	out_dict[header]["A_per"] = a_count / seq_length
	out_dict[header]["T_per"] = t_count / seq_length
	out_dict[header]["C_per"] = c_count / seq_length
	out_dict[header]["G_per"] = g_count / seq_length
	out_dict[header]["Length"] = seq_length
	
	return(out_dict)
	

with open(args.fasta) as fh:
	for line in fh:
		if line.startswith(">"):
			if seq != "":
				out_dict = add_seq(out_dict, seq, header)
			header = line.strip()
			seq = ""
		else:
			seq += line.strip()

#output
out_header = ["Header","A","T","C","G","GC_per","A_per","T_per","C_per","G_per","Length"]
with open(args.output, 'w') as out:
	for line in out_header[:-1]:
		tmp = out.write(line + "\t")
	tmp = out.write(out_header[-1] + "\n")
	for header in out_dict.keys():
		tmp = out.write(header + "\t" + str(out_dict[header]["A"]) + "\t"
						+ str(out_dict[header]["A"]) + "\t" + str(out_dict[header]["T"]) + "\t"
						+ str(out_dict[header]["C"]) + "\t" + str(out_dict[header]["G"]) + "\t"
						+ str(out_dict[header]["GC_per"]) + "\t" + str(out_dict[header]["A_per"]) + "\t"
						+ str(out_dict[header]["T_per"]) + "\t" + str(out_dict[header]["C_per"]) + "\t"
						+ str(out_dict[header]["G_per"]) + "\t" + str(out_dict[header]["Length"]) + "\n")



