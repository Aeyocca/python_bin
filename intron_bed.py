#!/bin/python
#intron_bed.py

#takes in gene_bed, exon_bed and generates an intron_bed
#couldn't exactly figure out how to get bedtools intersect to do this for me

import argparse
import sys

parser = argparse.ArgumentParser(prog='PROG')
parser.add_argument('--gene_bed', required=True, help='input gene bed file')
parser.add_argument('--exon_bed', required=True, help='input exon bed file')
parser.add_argument('--output', required=True, help='output bed file')
args = parser.parse_args()

#dictionary of lists
#gene_bed_dict = dict()
#exon_bed_dict = dict()
#out_dict = dict()

#using gene_bed as a basis, so if gene_bed doesn't have a certain chromosome, then it won't be reported in the output

def load_bed(file = ""):
	out_dict = dict()
	with open(file) as fh:
		for line in fh:
			la = line.strip().split("\t")
			#check if need to initialize chrom
			if la[0] not in out_dict:
				#initialize
				out_dict[la[0]] = []
			#these are all consecutive...
			out_dict[la[0]] += list(range(int(la[1]),int(la[2]) + 1))
	return(out_dict)

def drop_exon(list_one = [], list_two = []):
	#list_one = gd
	#one line, why make a function
	#one for checks of chromosome matching?
	#[x for x in list_one if x not in list_two]
	
	print("Assuming sorted inputs")
	
	#Hmm something to speed it up like a sliding window assuming everything is sorted
	#maybe one single list instead of list of lists
	out_list = []
	j = 0
	i = 0
	#also what happens when j > len(list_two) ? 
	while i < len(list_one):
		if list_one[i] > list_two[j]:
			#haven't gotten to the exon yet
			#increment and next, wait don't want to increment i though
			j += 1
		elif i < list_two[j]:
			#intron!!
			out_list += i
			i += 1
		elif i == list_two[j]:
			#exon sequence, increment and skip
			i += 1
	
	return out_list

def write_out(file = "",out_dict = dict()):
	#outputting.. every time we get non-consecutive values, add linebreak
	#how to...
	#for first value to last value, max of a few million so shouldn't take that long
	with open(file,'w') as out:
		for chrom in out_dict.keys():
			interval = []
			#print("chrom: %s" % (chrom))
			#print("length introns: %s" % (len(out_dict[chrom])))
			#sort
			out_dict[chrom].sort()
			#print("Finished sorting %s" % (chrom))
			#testing length of interval is incredibly slow for some reason?
			#or maybe this loop is just crazy slow
			
			#sys.exit(out_dict[chrom][0:10])
			
			#how to make this faster, want intervals of continuous sequence
			min_coord = out_dict[chrom][0]
			last_coord = out_dict[chrom][0]
			
			for i in out_dict[chrom][1:]:
				if i - last_coord > 1:
					#write out
					tmp = out.write(chrom + "\t" + str(min_coord) + "\t" + 
									str(last_coord) + "\n")
					min_coord = i
				last_coord = i
			#straggler
			tmp = out.write(chrom + "\t" + str(min_coord) + "\t" + 
									str(last_coord) + "\n")
			

if __name__ == "__main__":
	gbd = load_bed(file = args.gene_bed)
	ebd = load_bed(file = args.exon_bed)
	
	out_dict = dict()
	for chrom in gbd.keys():
		#print("Subsetting %s" % (chrom))
		#out_dict[chrom] = [x for x in gbd[chrom] if x not in ebd[chrom]]
		#out_dict[chrom] = drop_exon(list_one = gbd[chrom], list_two = ebd[chrom])
		#def a better way to do this, too slow it seems
		out_dict[chrom] = list(set(gbd[chrom]) - set(ebd[chrom]))
		
		#hmm any concern this might leave "UTR" seqs as "introns" still based on gff naming conventions? grepped for "exon", but maybe people labeled UTR instead of "exon"
	write_out(file = args.output, out_dict = out_dict)