#!/bin/python
#samtools_faidx_from_bed.py

#takes in gene_bed, exon_bed and generates an intron_bed
#couldn't exactly figure out how to get bedtools intersect to do this for me

import argparse
import sys

parser = argparse.ArgumentParser(prog='PROG')
parser.add_argument('--bed', required=True, help='input bed file')
parser.add_argument('--fasta', required=True, help='input fasta file')
parser.add_argument('--output', required=True, help='output bed file')
args = parser.parse_args()

#loop list of introns and? No use LOL from bed as indicies to subset a large fasta string, that should be fastest huh

fasta_dict = dict()
header = []

#what the hell have I been testing on, why is this soooo slow all of the suddent
#you know what? It might be that its trying to read wrapped files so the append operation
#is ocurring many times...
#yupp thats probably it!
with open(args.fasta) as fh:
	for line in fh:
		if line.startswith(">"):
			header = line.strip().replace(">","")
			fasta_dict[header]= ""
		else:
			fasta_dict[header] += line.strip()

bed_dict = dict()
with open(args.bed) as fh:
	for line in fh:
		la = line.strip().split("\t")
		if la[0] not in bed_dict:
			#initialize
			bed_dict[la[0]] = []
		#subtracting 1 from everything because will use the 1-based integers to subset the fasta file
		bed_dict[la[0]].append([int(la[1]) - 1,int(la[2]) - 1])


out_dict = dict()
for chrom in fasta_dict.keys():
	#alrighty, should work?
	print("Subsetting %s" % (chrom))
	#hmm this gets them all as one concatenated string, thats fine right?
	#ugh, no probably want them separate, great also if genes were listed in the fasta header
	#this is of pretty limited utility huh? All for motif enrichment
	#check how fast this runs
	#fast, nice, but probably not using it anymore
	#output_list = [fasta_dict[chrom][x] for x in bed_dict[chrom]]
	#out_dict[chrom] = ''.join(output_list)

	for interval in bed_dict[chrom]:
		loop_list = [fasta_dict[chrom][x] for x in range(interval[0],interval[1] + 1)]
		if chrom not in out_dict:
			#initialize
			out_dict[chrom] = [[''.join(loop_list)]]
		out_dict[chrom].append([''.join(loop_list)])

i = 0
with open(args.output, 'w') as out:
	for chrom in out_dict.keys():
		for fasta_string in out_dict[chrom]:
			i += 1
			tmp = out.write(">" + chrom + "_bed_line_" + str(i) + "\n" + 
							''.join(fasta_string) + "\n")

