#!/bin/python
#Alan E. Yocca
#04-15-2021
#first_intron_coords.py

import csv
import sys
import argparse
import re

parser = argparse.ArgumentParser(prog='PROG')
parser.add_argument('--gff_file', required=True, help='osat or bdis gff file')
parser.add_argument('--species', required=True, help='osat or bdis for controlling regex')
args = parser.parse_args()

#tragic.. this gets the longest first intron, not the longest transcript.........
#AHHHHH
#if we make the fourth array entry the gene length we can get this but still fucks with the intronless 
#gene calc

#read in gff, output bed file? Or should I just do the extraction here too?
#no think I can use bedtools getfasta for that

#Key Gene, Value with chrom, start, stop, length
#if zero length... grrrr well just put mRNA coords
#need to filter this bed file later for which genes have no introns
#so can get concrete value and not that its just missing
#how do we know zero length introns then????
#ugh. another array element
gff_dict = dict()
#initialize some variables
cg_length = ""
cg = ""
#current_transcript_intron_start = ""
ct_int_start = ""
ct_int_stop = ""
ct_chrom = ""
#yay series of logic switches!
#well I can't think of any other distinguishable characteristics
first_exon = True
second_exon = False

exon_anno = ""
gff_file = ""
if args.species == "osat":
	exon_anno = "exon"
	gff_file = "/mnt/research/edgerpat_lab/AlanY/00_collab/02_core_pred/06_rice/IRGSP-1.0_representative/transcripts_exon.gff"
elif args.species == "bdis":
	exon_anno = "CDS"
	gff_file = "Bdistachyon_283_v2.1.gene.gff3"
else:
	sys.exit("Please specify one of osat or bdis for species")

with open(args.gff_file) as fh:
	for line in fh:
		if re.match("^#",line):
			continue
		line_array = line.strip().split("\t")
		if re.match("mRNA",line_array[2]):
			cg_length = int(line_array[4]) - int(line_array[3])
			if second_exon:
				#second exon doesn't exist, length to zero
				try:
					if gff_dict[4] < cg_length:
						gff_dict[cg] = (ct_chrom, line_array[3], line_array[4], 0, cg_length)
				except KeyError:
					gff_dict[cg] =  (ct_chrom, line_array[3], line_array[4], 0, cg_length)
		
			#next gene, next cds annotation will be first exon
			first_exon = True
			second_exon = False
			
			if ct_int_start == "":
				#print("Passing first")
				pass
	
			if args.species == "osat":
				cg = re.search("ID=([a-zA-Z0-9\-]*);",line_array[8]).group(1)
				cg = cg.split("-")[0].replace("t","g")
			else:
				cg = re.search("ID=([a-zA-Z0-9\-]*)\.",line_array[8]).group(1)
			
		#manual checking of gff file show first "CDS" annotation is first exon regardless of strand
		#going to measure between first and second "CDS" annotation
		#hmm should they all be +/- 1 nucleotide?
		#need to double check bedtools getfasta, I imagine I should
		if re.match(exon_anno,line_array[2]):
			if first_exon:
				first_exon = False
				#hmm only overwrite if will be longer
				ct_int_start = int(line_array[4])
				#next cds will be our second exon
				#WAIT what if it doesn't exist
				second_exon = True
				ct_chrom = line_array[0]
			elif second_exon:
				second_exon = False
				ct_int_stop = int(line_array[3])
				#if longest replace so yeah this has to be a dictionary
				if line_array[6] == "-":
					#switch start / stop
					tmp = ct_int_start
					ct_int_start = ct_int_stop
					ct_int_stop = tmp
				
				ct_length = (ct_int_stop - ct_int_start - 2)
				try:
					#minus two here because not including last base of first/second exon
					if gff_dict[cg][3] < cg_length:
						gff_dict[cg] = (line_array[0], ct_int_start + 1, ct_int_stop - 1, ct_length, cg_length)
				except KeyError:
					gff_dict[cg] = (line_array[0], ct_int_start + 1, ct_int_stop - 1, ct_length, cg_length)

with open(args.gff_file + "_first_intron.bed", "w") as out:
	for gene in gff_dict.keys():
		#a little funky
		intron_str = [str(x) for x in gff_dict[gene]]
		tmp = out.write(intron_str[0] + "\t" + intron_str[1] + "\t"
						+ intron_str[2] + "\t" + gene + "\t" + intron_str[3] + "\n")
