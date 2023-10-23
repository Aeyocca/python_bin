#!/bin/python
#05-15-20
#fst_calc_all.py

import argparse
import dendropy
from datetime import datetime

parser = argparse.ArgumentParser(prog='PROG')
parser.add_argument('--input', required=True, type = int, help='how many accessions')
#parser.add_argument('--nreads', type = int, required=True, help='number of reads to subset')
parser.add_argument('--output', required=True, help='output fst file')
args = parser.parse_args()

rounded_up = -(- args.input // 2)
#doing only half got max fst of 0.5... oddd... try the full gamet
print("Running all instead of half")
rounded_up = 1136
output = []
for i in range(0,rounded_up):
  if (i % 100) == 0:
    print("Starting position %s" % (i))
  maf = i / args.input
  phylip_string= str(args.input) + " 1\n"
  simulated_snp_list = [1] * i
  simulated_snp_list = simulated_snp_list + [0] * (args.input - i)
  for taxa in range(0,args.input):
    ecotype_string = str(taxa) + " " * (10 - len(str(taxa)))
    phylip_string += ecotype_string + str(simulated_snp_list[taxa]) + "\n"
    loop_snps = dendropy.StandardCharacterMatrix.get(
    data=phylip_string,
    schema="phylip")
  loop_fst = dendropy.calculate.popgenstat.nucleotide_diversity(loop_snps)
  output.append([maf,loop_fst])

with open(args.output, 'w', newline="") as f:
  for sublist in output:
    #assign return value to variable
    tmp=f.write(str(sublist[0]))
    tmp=f.write('\t')
    tmp=f.write(str(sublist[1]))
    tmp=f.write('\n')


