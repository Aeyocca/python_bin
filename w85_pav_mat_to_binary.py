#!/bin/python
#w85_pav_mat_to_binary.py
#10-22-22
#Alan E. Yocca

import sys

#nhb_idx = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]
#dropping Draper
nhb_idx = [1,2,3,4,5,6,7,8,9,10,11,12,13,14]

shb_idx = [16,17,18,19,20,21,22,23,24,25,26,27,28,29]
#nice they're ordered

header = []
nhb_out = []
shb_out = []
i = 0		

#who am I kidding this will all be hardcoded
with open("Pangenome_Core_Summary_BBonly_4Alan.txt") as fh:
	#header line
	header = next(fh).strip().split("\t")
	nhb_header = [header[x] for x in nhb_idx]
	shb_header = [header[x] for x in shb_idx]
	nhb_header.insert(0,header[0])
	shb_header.insert(0,header[0])
	nhb_out.append(nhb_header)
	shb_out.append(shb_header)
	
	for line in fh:
		i+=1
		#don't strip or will remove empty fields
		#since subsetting to nhb / shb don't care about stripping anyway
		line_array = line.split("\t")
		bin_output = [0 if x == "" else 1 for x in line_array]
		nhb_sub = [bin_output[x] for x in nhb_idx]
		shb_sub = [bin_output[x] for x in shb_idx]
		
		nhb_sub.insert(0,line_array[0])
		shb_sub.insert(0,line_array[0])
		
		#if i == 1:
		#	nhb_out = [nhb_sub]
		#	shb_out = [shb_sub]
		#else:
		#	nhb_out.append(nhb_sub)
		#	shb_out.append(shb_sub)
		nhb_out.append(nhb_sub)
		shb_out.append(shb_sub)

i = 0		
with open("w85_ref_NHB_pan_mat.txt", "w") as out:
	for line in nhb_out:
		i+=1
		try:
			for item in line[:-1]:
				tmp = out.write(str(item) + "\t")
			tmp = out.write(str(line[-1]) + "\n")
		except TypeError:
			print(i)
			print(len(nhb_out))
			sys.exit()

with open("w85_ref_SHB_pan_mat.txt", "w") as out:
	for line in shb_out:
		for item in line[:-1]:
			tmp = out.write(str(item) + "\t")
		tmp = out.write(str(line[-1]) + "\n")

