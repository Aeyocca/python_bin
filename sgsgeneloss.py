#!/bin/python
#sgsgeneloss.py

import sys
import re
import argparse

parser = argparse.ArgumentParser(prog='PROG')
parser.add_argument('--input', default = "none", required=False, help='input bedgraph files, comma separated')
parser.add_argument('--input_list', required=False, default = "none", help='file of input bedgraph files, one per line separated')
parser.add_argument('--gff', required=True, help='reference annotation')
parser.add_argument('--depth_threshold', required=False, default = 1, help='depth threshold for presence, default 1 read')
parser.add_argument('--cov_threshold', required=False, default = 0.8, help='horizontal coverage threshold for presence. Default 0.8 of exon bases passing depth_threshold')
parser.add_argument('--output', required=True, help='output PAV table')

args = parser.parse_args()

def load_gff(gff = ""):
	#dictionary, key = gene, value = list of exon boundaries
	#two objects returned, a gff_bg to index bases, and a gff_cov to collect length
	#eh, can be one object but lets make it two
	gene_id = ""
	ref_gff = dict()
	cov_gff = dict()
	with open(gff) as fh:
		for line in fh:
			if line.startswith("#"):
				continue
			la = line.strip().split("\t")
			###sooo 5' utr / 3' were messing with it
			#if la[2] == "exon" or la[2] == "five_prime_UTR" or la[2] == "three_prime_UTR":
			if la[2] == "exon":
				trans_id_tmp = la[8].split("Parent=")[1]
				trans_id = trans_id_tmp.split(";")[0]
				#should be initialized at the same time so just check one
				if trans_id in cov_gff:
					for bp in range(int(la[3]),int(la[4]) + 1):
						cov_gff[trans_id]["Length"] += 1
						#if overlapping gene
						if str(bp) in ref_gff[la[0]].keys():
							ref_gff[la[0]][str(bp)].append(trans_id)
						else:
							ref_gff[la[0]][str(bp)] = [trans_id]
				else:
					cov_gff[trans_id] = {"Length" : int(la[4]) - int(la[3]) + 1, "Cov" : 0, "Cov_bases" : 0}
					if la[0] not in list(ref_gff.keys()):
						ref_gff[la[0]] = dict()
					for bp in range(int(la[3]),int(la[4]) + 1):
						if str(bp) in ref_gff[la[0]].keys():
							ref_gff[la[0]][str(bp)].append(trans_id)
						else:
							ref_gff[la[0]][str(bp)] = [trans_id]
	return(ref_gff,cov_gff)

#that should be it, gene set to a list of exon boundaries
def loop_bg(bg_file = "", ref_gff = dict(), cov_gff = dict()):
	#return an updated coverage gff
	with open(bg_file) as fh:
		for line in fh:
			la = line.strip().split("\t")
			for base in range(int(la[1]),int(la[2])):
				try:
					if int(la[3]) >= 1:
						for gene in ref_gff[la[0]][str(base)]:
							cov_gff[gene]["Cov"] += int(la[3])
							cov_gff[gene]["Cov_bases"] += 1
				except KeyError:
					#no gene overlaps this base
					continue
	return(cov_gff)

def calc_pav(cov_gff = dict(), cov_threshold = float(), depth_threshold = int()):
	#simple division here
	out_dict = dict()
	for gene in cov_gff.keys():
		if cov_gff[gene]["Cov"] / cov_gff[gene]["Length"] >= float(depth_threshold) and cov_gff[gene]["Cov_bases"] / cov_gff[gene]["Length"] >= float(cov_threshold):
			out_dict[gene] = 1
		else:
			out_dict[gene] = 0
	return(out_dict)

def transpose_dict(in_dict = dict()):
	out_dict = dict()
	i = 0
	for tag in in_dict.keys():
		i = i + 1
		if i == 1:
			#first loop so initializing
			out_dict["File"] = [tag]
			for gene in in_dict[tag].keys():
				out_dict[gene] = [in_dict[tag][gene]]
		else:
			out_dict["File"].append(tag)
			for gene in in_dict[tag].keys():
				out_dict[gene].append(in_dict[tag][gene])
	return(out_dict)

def main():
	#print("WARNING: ensure only primary transcripts are in gff... ? wait should that matter?")
	#nope that actually shouldn't affect the calculation for the primary, can filter after this script
	output = dict()
	#print("%s" % (args.depth_threshold))
	#ref_gff, cov_gff = load_gff(gff = args.gff)

	in_list = []
	if args.input_list == "none" and args.input == "none":
		sys.exit("Must specify one of input or input_list")
	elif args.input == "none":
		with open(args.input_list) as fh:
			for line in fh:
				in_list.append(line.strip())
	else:
		in_list = args.input.split(",")
	
	for bg_file in in_list:
		tag = bg_file.replace(".bed","")
		output[tag] = dict()
		#reset every iteration to be sure
		ref_gff, cov_gff = load_gff(gff = args.gff)
		cov_gff = loop_bg(bg_file = bg_file, ref_gff = ref_gff, cov_gff = cov_gff)
		output[tag] = calc_pav(cov_gff = cov_gff, cov_threshold = args.cov_threshold, depth_threshold = args.depth_threshold)		
		
	out_t = transpose_dict(output)
	
	with open(args.output, 'w') as out:
		for gene in out_t.keys():
			tmp = out.write(gene + "\t")
			for col in out_t[gene][:-1]:
				tmp = out.write(str(col) + "\t")
			tmp = out.write(str(out_t[gene][-1]) + "\n")



if __name__ == "__main__":
	main()












