#!/bin/python
#CHATGPT GENERATED!
#this is exciting
#bed_sliding_window.py
#hmm needs a few tweaks but should help

import pandas as pd
import argparse

parser = argparse.ArgumentParser(prog='PROG')
parser.add_argument('--in_bed', required=True, help='input bed feature file')
parser.add_argument('--faidx', required=True, help='fasta index file')
parser.add_argument('--window_size', required=False, default="50000", help='window_size')
parser.add_argument('--out_bed', required=True, help='output sliding window file')

args = parser.parse_args()

#going to read bed file instead
in_bed = pd.read_table(args.in_bed,header=None)
in_bed.columns = ["chrom", "start", "end", "gene"]

# Define sliding window parameters
#have to be non-overlapping because... easier to calculate overlap tbh
window_size = int(args.window_size)

# Create empty DataFrame to store depth values
depth_df = pd.DataFrame(columns=["chrom", "start", "end", "depth"])

#fill in sliding windows using fasta index, initialize all windows with 0 depth
with open(args.faidx) as fh:
	for line in fh:
		la = line.strip().split("\t")
		for i in range(0,int(la[1]),window_size):
			depth_df = pd.concat([depth_df, pd.DataFrame({'chrom': [la[0]], 'start': [i], 'end': [i + window_size - 1],"depth": [0]})], ignore_index = True)

#convert to numeric, just have to do start and depth
depth_df["start"] = pd.to_numeric(depth_df["start"])
depth_df["depth"] = pd.to_numeric(depth_df["depth"])

# Iterate over chromosomes
for chrom in in_bed["chrom"].unique():
    # Filter gene coordinates for current chromosome
    chrom_genes = in_bed[in_bed["chrom"] == chrom]
    loop_depth = depth_df[depth_df["chrom"] == chrom]
    chrom_genes = chrom_genes.sort_values(by=["start"])
    
    for i, gene in chrom_genes.iterrows():
     	#subtract gene start from all windows but.. ahh
    	#smallest positive difference
    	subs = gene["start"] - loop_depth["start"]
    	left_window = loop_depth["start"][subs == min([i for i in subs if i >= 0])].iloc[0]
    	
    	subs = gene["end"] - loop_depth["start"]
    	right_window = loop_depth["start"][subs == min([i for i in subs if i >= 0])].iloc[0]
    	
    	#add one to all windows between these
    	for i in range(left_window,right_window,window_size):
    		depth_df.loc[(depth_df.chrom == chrom) & (depth_df.start == i),"depth"] += 1
    	
    	if left_window == right_window:
    		depth_df.loc[(depth_df.chrom == chrom) & (depth_df.start == left_window),"depth"] += 1
    	else:
    		#need to capture the right window
    		depth_df.loc[(depth_df.chrom == chrom) & (depth_df.start == right_window),"depth"] += 1

# Save depth values to BED file
depth_df.to_csv(args.out_bed, sep="\t", header=False, index=False)

#print(depth_df)