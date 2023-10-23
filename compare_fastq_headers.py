#!/bin/bash
#Alan E. Yocca
#compare_fastq_headers.py

import sys
import argparse

parser = argparse.ArgumentParser(prog='PROG')
parser.add_argument('--read_one', required=True, help='read_one file')
parser.add_argument('--read_two', required=True, help='read_two file')
args = parser.parse_args()



file_one_list=[]
file_two_list=[]

file_one_headers = []
file_two_headers = []


#up to first space should match
#right now just check they match through the whole file

i=0
with open(args.read_one) as f1, open(args.read_two) as f2:
	for f1_line in f1:
		i+=1
		f2_line = f2.readline()
		if i+3 % 4 == 0:
			#header sequence
			f1_string = f1_line.split(" ")[0]
			f2_string = f2_line.split(" ")[0]
			
			if f1_string != f2_string:
				sys.exit("%s not equal %s" % (f1_string,f2_string))

#i=0
#j=0
#with open(args.read_two) as fh:
#	for line in fh:
#		i+=1
#		if i % 4 == 0:
#			#header sequence
#			match_string = line.split(" ")[0]
#			if match_string != file_one_headers[j]:
#				sys.exit("%s not equal %s" % (match_string,file_one_headers[j]))
#			j +=1
