#!/bin/python
#Alan E. Yocca
#04-04-20

#subset_fastq.py
import argparse
import subprocess
import random

parser = argparse.ArgumentParser(prog='PROG')
parser.add_argument('--input', required=True, help='input fastq file')
parser.add_argument('--nreads', type = int, required=True, help='number of reads to subset')
parser.add_argument('--output', required=True, help='output fastq file')
args = parser.parse_args()

#ugh, need to know the total size huh... if doing line by line then yes
#could just take the first X lines...
#could sys(wc -l), that would add a few minutes but memory usage would still be low I believe
#lets do that

wc_string = subprocess.run(("wc -l %s" % (args.input)), stdout=subprocess.PIPE, universal_newlines=True, shell=True)

print("wc string: %s" % (wc_string.stdout))

wc_lines = [int(i) for i in wc_string.stdout.split() if i.isdigit()]
total_reads = wc_lines[0] / 4


print("total reads: %s" % (total_reads))

#get random list of numbers in set 0 to total_reads
if args.nreads > total_reads:
	sys.exit("number of reads (%s) must be less than total reads (%s)" % (args.nreads, total_reads))

sample_space = range(1,int(1+total_reads))
subset = random.sample(sample_space, k = int(args.nreads))

print("len sub: %s" % (len(subset)))

with open(args.input) as input:
	with open(args.output, "w") as output:
	
		for line in range(1,int(total_reads + 1)):
			line1 = input.readline()
			line2 = input.readline()
			line3 = input.readline()
			line4 = input.readline()
		
			if line in subset:
				output.write(line1)
				output.write(line2)
				output.write(line3)
				output.write(line4)
