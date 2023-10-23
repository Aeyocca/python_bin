#!/bin/python
#Alan E. Yocca
#bbpan_snp_mat_sub.py

#import pandas as pd
#import numpy as np
#import matplotlib.pyplot as plt
#from sklearn.decomposition import PCA
from random import sample
import argparse

parser = argparse.ArgumentParser(prog='PROG')
parser.add_argument('--input', required=True, help='input vcf file')
parser.add_argument('--sub_size', required=False, default = 1000000, help='sites to sub')
parser.add_argument('--line_count', required=False, default = 61846878, help='total number of sites')
parser.add_argument('--output', required=True, help='output txt file')

args = parser.parse_args()


sub_mat = []
header = []
random_sample = sample(range(0,int(args.line_count)), int(args.sub_size))
#sort
sorted_random = sorted(random_sample,key = int)
#a little funky, but to avoid list index out of range
#add a final item that won't match so j stops being incremented
sorted_random.append(0)
j=0
i=0
with open(args.input) as fh:
	header = next(fh)
	#length of file, any way to not hard code this?
	for line in fh:
		if i == sorted_random[j]:
			#print
			line_array = line.strip().split("\t")
			sub_mat.append(line_array)
			j+=1
			
		if (i % 10000000) == 0:
			print("through line " + str(i))
		i+=1

print("Outputting")

with open(args.output, 'w') as out:
	tmp = out.write(header)
	for line in sub_mat:
		for i in range(0,len(line) -1):
			tmp = out.write(line[i] + "\t")
		tmp = out.write(line[i] + "\n")

