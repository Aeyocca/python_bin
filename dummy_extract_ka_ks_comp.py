#!/bin/python
#dummy_extract_ka_ks_comp.py
#dummy script to read some output
#get data tidyish to plot

import sys
import re
#model	accuracy	class	feat_imp	arb_value Bal_Type	Spec_type

"""
{99 : 
	{"train_test" : {
		"Train Osat Test Osat" : {
			"SVC" : num(),
			"GNB" : num(),
			"RFT" : num(),
			"Class" : {"Ka" : num()}
		}}}}
"""
#nested dictionary? eh
#write out while reading? naw
#model accuracy	arb_value
#model	class	feat_imp arb_value

tidy_dict = dict()
train_test = ""
arb_test_list = [0,0.25,0.5,0.75,1,1.5,2,3,5,10,20,50,99]
wkdir = "/mnt/research/edgerpat_lab/AlanY/Error_files/"

for i in range(len(arb_test_list)):
	slurm_id = i + 1
	tidy_dict[arb_test_list[i]] = dict()
	with open(wkdir + "bdis_model_test-12790243_" + str(slurm_id) + ".SLURMout") as fh:
		for line in fh:
			if line.startswith("Training"):
				line_array = line.strip().split(" ")
				train_test = ' '.join([str(elem) for elem in line_array[2:6]]) 
				#initiate dict
				tidy_dict[arb_test_list[i]][train_test] = {"Class" : {}}
			if line.startswith("SVC"):
				line_array = line.strip().split(" ")
				accuracy = line.strip().split(" ")
				tidy_dict[arb_test_list[i]][train_test]["SVC"] = line_array[-1]
			if line.startswith("GNB"):
				line_array = line.strip().split(" ")
				accuracy = line.strip().split(" ")
				tidy_dict[arb_test_list[i]][train_test]["GNB"] = line_array[-1]
			if line.startswith("RFT"):
				line_array = line.strip().split(" ")
				accuracy = line.strip().split(" ")
				tidy_dict[arb_test_list[i]][train_test]["RFT"] = line_array[-1]
			if re.match(r"[0-9]+", line):
				line_array = re.split('\s+', line)
				tidy_dict[arb_test_list[i]][train_test]["Class"][line_array[1]] = line_array[2]


with open("osat_bdis_ka_ks_arb_change_tidy.txt", "w") as out:
	out.write("Arb_Value\tTrain_Test\tModel\tAccuracy\tClass\tFeat_Imp")
	for arb in tidy_dict:
		#ttt == train_test_type
		for ttt in tidy_dict[arb]:
			for key in tidy_dict[arb][ttt]:
				if key == "Class":
					for ka_ks_class in tidy_dict[arb][ttt][key]:
						out.write("\n" + str(arb) + "\t" + str(ttt) + "\tSVC\t" + 
								str(tidy_dict[arb][ttt]["SVC"] + "\t" + str(ka_ks_class) +
								"\t" + str(tidy_dict[arb][ttt]["Class"][ka_ks_class])))
						out.write("\n" + str(arb) + "\t" + str(ttt) + "\tGNB\t" + 
								str(tidy_dict[arb][ttt]["GNB"] + "\t" + str(ka_ks_class) +
								"\t" + str(tidy_dict[arb][ttt]["Class"][ka_ks_class])))
						out.write("\n" + str(arb) + "\t" + str(ttt) + "\tRFT\t" + 
								str(tidy_dict[arb][ttt]["RFT"] + "\t" + str(ka_ks_class) +
								"\t" + str(tidy_dict[arb][ttt]["Class"][ka_ks_class])))
