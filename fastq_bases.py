#!/bin/python
#fastq_bases.py

#dummy script to count bases in fastq file

#only reads the line after the one starting with an "@"
#ugh.. that doesnt work it seems since "@" used as a quality value...
#if we do see it though, don't count double lines
#wait! this will still work becase it should pass until it gets a line without an @
#if we hit a quality line starting with @, the next line will be the header of the next read

import sys
import argparse
import statistics

parser = argparse.ArgumentParser(prog='PROG')
parser.add_argument('--input', required=True, help='input fastq file')
parser.add_argument('--median', required=False, default = 0, help='if set to 1, report median also. Assumes not wrapped')
args = parser.parse_args()

base_count = 0
base_list = []
switch = 0

with open(args.input) as fh:
	for line in fh:
		if line.startswith("@"):
			switch = 1
		else:
			if switch:
				read_length = len(line.strip())
				base_count += read_length
				base_list.append(read_length)
				switch = 0

print(str(base_count))
if args.median:
	print(statistics.median(base_list))