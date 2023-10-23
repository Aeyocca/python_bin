#!/bin/python
#Alan E. Yocca
#11-14-22
#combine_EDTA.py

#import csv
#import sys
#import argparse
#import re

#parser = argparse.ArgumentParser(prog='PROG')
#parser.add_argument('--input', required=True, help='input fasta file')
#parser.add_argument('--output', required=True, help='output fasta file')

#args = parser.parse_args()



#dont think I'm even going to argparse

def add_rep(out_df = dict(), line = "", rep = "", unknown = ""):
	if rep in line:
		if rep == "unknown":
			rep = unknown + "_unknown"
		out_df[rep]["Count"] += int(line.strip().split()[1])
		out_df[rep]["bpMasked"] += int(line.strip().split()[2])
	return out_df


def read_summary(out_df = dict(), filename = ""):

	rep_class = ""
	sub_class = ""
	#no we can just do subclasses, ugh, two unknowns
	#are any substrings of the other?
	with open(filename) as fh:
		for line in fh:
			if line.startswith("LTR"):
				rep_class = "LTR"
				continue
			elif line.startswith("nonLTR"):
				rep_class = "nonLTR"
				continue
			#start the nonsensical if/else matching
			rep_list = ["Copia","Gypsy","unknown","CACTA","Mutator",
						"PIF_Harbinger","Tc1_Mariner","hAT","LINE_element","helitron",
						"repeat_region"]
			
			for rep in rep_list:
				out_df = add_rep(out_df = out_df, line = line, 
				rep = rep, unknown = rep_class)
	return out_df

rep_list = ["Copia","Gypsy","LTR_unknown","CACTA","Mutator",
			"PIF_Harbinger","Tc1_Mariner","hAT","LINE_element","helitron",
			"repeat_region","nonLTR_unknown"]
#initialize
out_df = dict()
for rep in rep_list:
	out_df[rep] = {"Count" : 0, "bpMasked" : 0}

for i in range(1,18):
	filename = "chr" + str(i) + "A/chr" + str(i) + "A.fasta.mod.EDTA.TEanno.sum"
	out_df = read_summary(out_df = out_df, filename = filename)

with open("hap1_edta_sum.txt", "w") as out:
	for rep in out_df.keys():
		tmp = out.write(rep + "\t" + str(out_df[rep]["Count"]) + 
				"\t" + str(out_df[rep]["bpMasked"]) + "\n")

out_df = dict()
for rep in rep_list:
	out_df[rep] = {"Count" : 0, "bpMasked" : 0}

for i in range(1,18):
	filename = "chr" + str(i) + "B/chr" + str(i) + "B.fasta.mod.EDTA.TEanno.sum"
	out_df = read_summary(out_df = out_df, filename = filename)

with open("hap2_edta_sum.txt", "w") as out:
	for rep in out_df.keys():
		tmp = out.write(rep + "\t" + str(out_df[rep]["Count"]) + 
				"\t" + str(out_df[rep]["bpMasked"]) + "\n")

"""
#LINE_element    0       0
#repeat_region
#nonLTR_unknown

LTR                    --           --           --   
    Copia              2985         2840585      7.91% 
    Gypsy              4377         5473133      15.24% 
    unknown            4467         4260100      11.86% 
TIR                    --           --           --   
    CACTA              2125         811547       2.26% 
    Mutator            5859         1433060      3.99% 
    PIF_Harbinger      2313         812208       2.26% 
    Tc1_Mariner        22           8166         0.02% 
    hAT                2216         766035       2.13% 
nonLTR                 --           --           --   
    LINE_element       99           70656        0.20% 
    unknown            63           61325        0.17% 
nonTIR                 --           --           --   
    helitron           502          198842       0.55% 
repeat_region          4687         1084594      3.02% 
                      ---------------------------------
    total interspersed 29715        17820251     49.63%
    
    
    
Class                  Count        bpMasked    %masked
=====                  =====        ========     =======
LTR                    --           --           --   
    Copia              3082         3276528      9.03% 
    Gypsy              4370         5221611      14.39% 
    unknown            4766         4087867      11.27% 
TIR                    --           --           --   
    CACTA              1527         629424       1.73% 
    Mutator            5676         1387037      3.82% 
    PIF_Harbinger      2129         748015       2.06% 
    Tc1_Mariner        27           8791         0.02% 
    hAT                2537         832384       2.29% 
nonTIR                 --           --           --   
    helitron           585          174335       0.48% 
                      ---------------------------------
    total interspersed 24699        16365992     45.11%

---------------------------------------------------------
Total                  24699        16365992     45.11%
"""



