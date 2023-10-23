#!/bin/python
#collect_edta.py

#dummy script to output edta summary 
#since we ran edta on scaffolds separately

"""
Repeat Classes
==============
Total Sequences: 1
Total Length: 25033000 bp
Class                  Count        bpMasked    %masked
=====                  =====        ========     =======
LTR                    --           --           --   
    Copia              1449         1369252      5.47% 
    Gypsy              2585         3545548      14.16% 
    unknown            3339         2464697      9.85% 
TIR                    --           --           --   
    CACTA              1368         485135       1.94% 
    Mutator            4036         893940       3.57% 
    PIF_Harbinger      1057         316493       1.26% 
    Tc1_Mariner        266          82445        0.33% 
    hAT                614          162691       0.65% 
nonLTR                 --           --           --   
    LINE_element       53           24778        0.10% 
    unknown            33           39824        0.16% 
nonTIR                 --           --           --   
    helitron           1136         341972       1.37% 
repeat_region          4350         1143306      4.57% 
                      ---------------------------------
    total interspersed 20286        10870081     43.42%

---------------------------------------------------------
Total                  20286        10870081     43.42%
"""

#welp.. if everything is in the same order....

#eh, series of if/else statements incase spacing changes

TE_dict = {"LTR_Copia" : {"Count" : 0, "bpMasked" : 0},
		   "LTR_Gypsy" : {"Count" : 0, "bpMasked" : 0},
		   "LTR_unknown" : {"Count" : 0, "bpMasked" : 0},
		   "TIR_CACTA" : {"Count" : 0, "bpMasked" : 0},
		   "TIR_Mutator" : {"Count" : 0, "bpMasked" : 0},
		   "TIR_PIF_Harbinger" : {"Count" : 0, "bpMasked" : 0},
		   "TIR_Tc1_Mariner" : {"Count" : 0, "bpMasked" : 0},
		   "TIR_hAT" : {"Count" : 0, "bpMasked" : 0},
		   "nonLTR_LINE" : {"Count" : 0, "bpMasked" : 0},
		   "nonLTR_unknown" : {"Count" : 0, "bpMasked" : 0},
		   "nonTIR_helitron" : {"Count" : 0, "bpMasked" : 0},
		   "repeat_region" : {"Count" : 0, "bpMasked" : 0}}

def add_count_seq(TE_dict = dict(), ls = "", tag = ""):
	la = ls.split(" ")
	#remove empty fields
	la = [x for x in la if x != ""]
	TE_dict[tag]["Count"] = TE_dict[tag]["Count"] + int(la[1])
	TE_dict[tag]["bpMasked"] = TE_dict[tag]["bpMasked"] + int(la[2])
	return TE_dict
	
tot_length = 0

hap = "A"
#hap = "B"
	
for i in range(1,18):
	with open("chr" + str(i) + hap + "/chr" + str(i) + hap + ".fasta.mod.EDTA.TEanno.sum") as fh:
		nonltr = False
		for line in fh:
			ls = line.strip()
			if ls.startswith("Total Length"):
				tot_length += int(ls.split("Total Length: ")[1].replace(" bp",""))
			elif ls.startswith("Copia"):
				TE_dict = add_count_seq(tag = "LTR_Copia", TE_dict = TE_dict, ls = ls)
			elif ls.startswith("Gypsy"):
				TE_dict = add_count_seq(tag = "LTR_Gypsy", TE_dict = TE_dict, ls = ls)
			elif ls.startswith("unknown"):
				if nonltr:
					TE_dict = add_count_seq(tag = "nonLTR_unknown", TE_dict = TE_dict, ls = ls)
				else:
					TE_dict = add_count_seq(tag = "LTR_unknown", TE_dict = TE_dict, ls = ls)
			elif ls.startswith("CACTA"):
				TE_dict = add_count_seq(tag = "TIR_CACTA", TE_dict = TE_dict, ls = ls)
			elif ls.startswith("Mutator"):
				TE_dict = add_count_seq(tag = "TIR_Mutator", TE_dict = TE_dict, ls = ls)
			elif ls.startswith("PIF_Harbinger"):
				TE_dict = add_count_seq(tag = "TIR_PIF_Harbinger", TE_dict = TE_dict, ls = ls)
			elif ls.startswith("Tc1_Mariner"):
				TE_dict = add_count_seq(tag = "TIR_Tc1_Mariner", TE_dict = TE_dict, ls = ls)
			elif ls.startswith("hAT"):
				TE_dict = add_count_seq(tag = "TIR_hAT", TE_dict = TE_dict, ls = ls)
			elif ls.startswith("LINE_element"):
				TE_dict = add_count_seq(tag = "nonLTR_LINE", TE_dict = TE_dict, ls = ls)
			elif ls.startswith("helitron"):
				TE_dict = add_count_seq(tag = "nonTIR_helitron", TE_dict = TE_dict, ls = ls)
			elif ls.startswith("nonLTR"):
				nonltr = True
			elif ls.startswith("repeat_region"):
				TE_dict = add_count_seq(tag = "repeat_region", TE_dict = TE_dict, ls = ls)	
				break

#print(str(tot_length) + "\n")
with open("hap" + hap + "_edta_summary.txt", 'w') as out:
	tmp = out.write("Repeat\tCount\tbpMasked\t% masked\n")
	for repeat in TE_dict.keys():
		tmp = out.write(repeat + "\t" + str(TE_dict[repeat]["Count"]) + "\t" + 
						str(TE_dict[repeat]["bpMasked"]) + "\t" + 
						str(round(TE_dict[repeat]["bpMasked"] / tot_length, 4)) + "\n")

