#!/bin/python
#filter_exon_bed_intersect.py

#take our bedtools intersect file of exons in one haplotypes with lifted over exons from the other
#check how many retain coverage across 50% of the entire gene model
#check how many are overlapped by the same gene model

import argparse
import sys
import pandas as pd



if __name__ == "__main__":
	parser = argparse.ArgumentParser(prog='PROG')
	parser.add_argument('--input', required=True, help='input bedtools intersect output file')
	parser.add_argument('--a_bed', required=True, help='original -a bed file')
	parser.add_argument('--threshold', required=False, default = "0.5", help='proportion overlap, default = 0.5')	
	parser.add_argument('--output', required=True, help='output file')
	args = parser.parse_args()
	
	#might be able to do all in pandas somehow
	#not very large files so looping them might be fine
						 
	#load in original a_bed file, may be unreported exons in bedtools intersect/
	#its a bit messy so consult the original
	a_bed = pd.read_csv(args.a_bed, sep = "\t")
	a_bed.columns = ['a_chrom','a_start','a_end','a_feature']
	a_bed['a_exon_length'] = a_bed['a_end'] - a_bed['a_start'] + 1
	a_feature_grouped = a_bed.groupby(['a_feature'])['a_exon_length'].sum()
	
	#that is a series of a_feature and its length, use later
	
	#next, we need to load in bedtools intersect file
	bed_table = pd.read_csv(args.input, sep = "\t")
	bed_table.columns = ['a_chrom','a_start','a_end','a_feature',
						 'b_chrom','b_start','b_end','b_feature',
						 'bp_overlap']
	
	#drop these lines where I assume there are just no b_bed overlaps
	bed_table = bed_table.drop(bed_table[bed_table.b_start == -1].index)
	
	#for each unique combination of a_feature and b_feature,
	#sum the bp_overlap
	bed_table_grouped = bed_table.groupby(['a_feature','b_feature'])['bp_overlap'].sum()

	#Lets convert them both to dataframes, then merge based on 'a_feature'	
	bed_group_df = pd.DataFrame(bed_table_grouped)
	a_group_df = pd.DataFrame(a_feature_grouped)
	
	#change rownames to their own column
	a_group_df['a_feature'] = a_group_df.index.values
	
	#since grouped on 2 variables, have MultiIndex so parse twice
	#index.get_level_values(0)
	bed_group_df['a_feature'] = bed_group_df.index.get_level_values(0)
	bed_group_df['b_feature'] = bed_group_df.index.get_level_values(1)
	
	#reset idx, because apparently 'a_feature' is the name of one of the indicies
	a_group_df.reset_index(drop = True, inplace = True)
	bed_group_df.reset_index(drop = True, inplace = True)
	
	#integrate a_exon_length
	merged_df = pd.merge(bed_group_df,a_group_df)
	
	#filter for combos where coverage >0 50% of total exon length
	over_fifty = merged_df[merged_df['bp_overlap'] / merged_df['a_exon_length'] >= float(args.threshold)]
	
	#write out
	over_fifty.to_csv(args.output, index = False, sep = "\t")
	
	
