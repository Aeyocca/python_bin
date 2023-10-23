#!/bin/python
#wig_bed_to_cns.py
#script to create bed of "CNS" from wig file
#Alan E. Yocca
#04-25-2023

import sys
import argparse
import numpy as np

parser = argparse.ArgumentParser(prog='PROG')
parser.add_argument('--wig_bed', required=True, help='Alan\'s bed created from wig')
parser.add_argument('--output', required=True, help='output bed file')

args = parser.parse_args()

high_score = 0.75
low_score = 0.5
min_len = 7
low_len = 12

#regions of at least min_len with avg >= high_score
#with no run of low_len avg <= low_score
#not sure if avg for either of these runs but seems silly if they aren't avgs

position = 0
loop_scores = []
cns_bed = []

with open(args.wig_bed) as fh:
	for line in fh:
		la = line.strip().split("\t")
		
		if int(la[1]) - position > 1:
			#non consecutive bases
			loop_scores = [float(la[3])]
			position = int(la[1])
			continue
		
		position = int(la[1])
		loop_scores.append(float(la[3]))
		#now the series of checks!
		#is it long enough?
		if len(loop_scores) < min_len:
			#not long enough, keep going
			continue
		
		#is the average above threshold...
		#when to cut it off though??
		if np.mean(loop_scores) >= high_score:
			if np.mean(loop_scores[-low_len:]) <= low_score:
				#cut it here
				cns_bed.append([position - len(loop_scores[:-low_len]), position])
				#reset loop_scores
				loop_scores = []
			
			continue
		
		#so its long, but the average isn't above threshold
		#ummmmm, if it was a CNS, the above if statement would have caught it
		#soo the assumption here is we haven't hit the threshold yet,
		#so remove the first element and move on
		loop_scores = loop_scores[1:]

with open(args.output, 'w') as out:
	for cns in cns_bed:
		tmp = out.write(str(cns[0]) + "\t" + str(cns[1]) + "\n")

