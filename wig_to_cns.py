#!/bin/python
#wig_to_cns.py
#script to create bed of "CNS" from wig file
#Alan E. Yocca
#04-25-2023

import sys
import argparse
import numpy as np

parser = argparse.ArgumentParser(prog='PROG')
parser.add_argument('--wig', required=True, help='input wig file')
parser.add_argument('--exon_bed', required=False, help='bed file of exon boundaries')
parser.add_argument('--output', required=True, help='output bed file')

args = parser.parse_args()

#load in exon boundaries
#hmm how to efficiently parse these
#if we hit an exon, if we are in an exon
#sooo just one giant list of coords, how much time does that take?!?!
#lets try it
exon_list = []
if args.exon_bed != "":
	with open(args.exon_bed) as fh:
		for line in fh:
			#assume 2nd and 3rd column
			la = line.strip().split("\t")
			for coord in list(range(int(la[1]),int(la[2]) + 1)):
				exon_list.append(coord)
	print("Finished loading exons")

#hmmmm our wig file only has fixedStep and lists single bases at a time, so a bit of a dummy script

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

with open(args.wig) as fh:
	for line in fh:
		if line.startswith("fixedStep"):
			#new block
			#fixedStep chrom=(null) start=104833 step=1
			#check step,
			la= line.strip().split(" ")
			if la[3] != "step=1":
				sys.exit("step not equal to 1, figure something out Alan")
			step_start = int(la[2].replace("start=",""))
			if step_start - position > 1:
				#reset loop_scores
				loop_scores = []
				#set position, -1 bc adding at the start of next loop
				position = step_start - 1

		else:
			position += 1
			
			if position in exon_list:
				#reeeeeset
				loop_scores = []
				continue
			#add score to current loop_scores
			loop_scores.append(float(line.strip()))
			
			#now the series of checks!
			#is it long enough?
			if len(loop_scores) < min_len:
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

