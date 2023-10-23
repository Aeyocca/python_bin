#!/bin/python
#Alan E. Yocca
#multi_vcf_to_snp_matrix.py
#convert multi vcf file to SNP matrix
#SNP PCA is my goal
#which is rows and which is columns?

import csv
import sys
import argparse
import re

parser = argparse.ArgumentParser(prog='PROG')
parser.add_argument('--input', required=True, help='input vcf file')
parser.add_argument('--maf', required=False, default = "0", help='Drop sites where MAF below this')
parser.add_argument('--output', required=True, help='output txt file')
#eventually add an argument to define a smaller set of individuals

args = parser.parse_args()

def binarize_vcf(line = "", maf = ""):
	end_fields = line.strip().split("\t")[9:]
	
	#MAF filter... how to handle >2 alleles
	#"Minor allele frequency (MAF) is the frequency at which the second most common allele occurs in a given population"
	
	allele_dict = dict()
	#hmm have set amount of accessions so can just add fractions...
	fractional_amount = 1/len(end_fields)

	out_line = []
	for acc in end_fields:
		genotype = acc.split(":")[0]
		if genotype == ".":
			out_line.append("0")
			try:
				allele_dict["0"] += fractional_amount
			except KeyError:
				allele_dict["0"] = fractional_amount
		else:
			out_line.append(genotype)
			try:
				allele_dict[genotype] += fractional_amount
			except KeyError:
				allele_dict[genotype] = fractional_amount
	
	#get value of second highest frequency
	maf_locus = sorted(list(allele_dict.values()), reverse = True)[1]
	if maf_locus < float(maf):
		out_line = []

	return out_line

with open(args.input) as fh:
	with open(args.output, "w") as out:
		for line in fh:
			if line.startswith("##"):
				continue
			elif line.startswith("#C"):
				acc_list = line.strip().split("\t")[9:]
				tmp = out.write("Marker" + "\t")
				for i in range(len(acc_list) - 1):
					tmp = out.write(acc_list[i] + "\t")
				tmp = out.write(acc_list[-1] + "\n")
			else:
				marker = line.strip().split("\t")[0] + "_" + line.strip().split("\t")[1]
				
				out_line = binarize_vcf(line = line, maf = args.maf)
				if len(out_line) > 0:
					tmp = out.write(marker + "\t")
					for i in range(len(out_line) - 1):
						out.write(out_line[i] + "\t")
					
					out.write(out_line[-1] + "\n")

