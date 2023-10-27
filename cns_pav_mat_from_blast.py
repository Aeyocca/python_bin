#!/bin/python
#cns_pav_mat_from_blast.py

#take list of accessions and hard code blast strings to read
#add a flag for only retaining single-copy hits

import argparse
import sys
import pandas as pd

def init_matrix(file = ""):
	"""
	intialize CNS matrix, first column reference CNS, each additional is 
	binary PAV per accession. Returns empty dictionary
	"""
	
	cns_list = []
	#make reference all present
	with open(file) as fh:
		for line in fh:
			cns_list.append(line.strip())
			
	return(pd.DataFrame(cns_list, columns = ['Ref_CNS']))

def add_accession(file = "", acc = "", df = pd.DataFrame()):
	"""
	Reads filtered blast file and adds binary column for present and absent
	optionally only keep single hits.. eh issue with that is multi hits will 
	be indistinguishable from non-hits.. think of a different way to notate
	other than simply dropping...
	Wow I'm dumb, just list the number of times present
	changes my list comprehension though
	"""
	
	print("Loading %s" % (acc))
	acc_list = []
	with open(file) as fh:
		for line in fh:
			acc_list.append(line.strip().split("\t")[0])
	
	#going to be a little slow I think?
	acc_column = [acc_list.count(x) for x in df['Ref_CNS'].tolist()]
	
	#order should be maintained here
	df[acc] = acc_column
	return(df)
	
if __name__ == "__main__":
	parser = argparse.ArgumentParser(prog='PROG')
	parser.add_argument('--accession_list', required=True, help='list of accessions')
	parser.add_argument('--ref_cns_list', required=True, help='list of reference CNS to search for. One line per CNS')
	parser.add_argument('--query_base', required=True, help='part of blast file name')
	parser.add_argument('--output', required=True, help='output file')
	args = parser.parse_args()
	
	"""
	whats our data structure here? dictionary of CNS - accession - value or 
	accession - CNS - value? What if we just did a pandas dataframe and some 
	sort of %in% operator like in perl???
	"""
	
	cns_df = init_matrix(args.ref_cns_list)
	
	with open(args.accession_list) as fh:
		for acc in fh:
			#hard code for now
			blast_dir = "02_blast/"
			blast_string = (
				blast_dir + acc.strip() + "." + args.query_base + ".ws7.filter.blast")
			cns_df = add_accession(file = blast_string, acc = acc.strip(), df = cns_df)
	
	cns_df.to_csv(args.output, sep = "\t", index = False, header = True)