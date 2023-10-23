#!/bin/python
#05-15-20
#loop_snp_files.py

import argparse
import dendropy
from datetime import datetime

parser = argparse.ArgumentParser(prog='PROG')
parser.add_argument('--input', required=True, help='input snp split file')
#parser.add_argument('--nreads', type = int, required=True, help='number of reads to subset')
parser.add_argument('--output', required=True, help='output snp split file')
args = parser.parse_args()

def read_list(filename = ""):
  lines=[]
  print("Loading snps")
  with open(filename) as fh:
    lines = [ [split.strip() for split in line.split('\t')] for line in fh]
  #split snp array
  #remove that last element that was produced accidentally
  #print("Splitting line")
  #now = datetime.now()
  #current_time = now.strftime("%H:%M:%S")
  #print("Current Time =", current_time)
  split_line = [[x[0], x[1], x[2].split(',')[0:1135]] for x in lines]
  #print("Finished splitting")
  #now = datetime.now()
  #current_time = now.strftime("%H:%M:%S")
  #print("Current Time =", current_time)
  return split_line

def calc_fst(snp_list = [], add_maf = True):
  fst_list = []
  ntaxa = str(len(snp_list[0][2]))
  pos=0
  for scp in snp_list:
    if (pos % 1000) == 0:
      print("Starting position %s" % (pos))
    pos+=1
    scp[2] = [int(x) for x in scp[2]]
    maf = sum(scp[2]) / len(scp[2])
    if maf > 0.5:
      maf = 1 - maf
    if maf == 0:
      #print("MAF == 0, therefore Fst == 0. Not calculating to save time")
      loop_fst = 0
      if add_maf:
        fst_list.append([scp[0],scp[1],loop_fst,maf])
      else:
        fst_list.append([scp[0],scp[1],loop_fst])
      continue
    phylip_string= ntaxa + " 1\n"
    for taxa in range(0,int(ntaxa)):
      ecotype_string = str(taxa) + " " * (10 - len(str(taxa)))
      phylip_string += ecotype_string + str(scp[2][taxa]) + "\n"
      loop_snps = dendropy.StandardCharacterMatrix.get(
      data=phylip_string,
      schema="phylip")
    loop_fst = dendropy.calculate.popgenstat.nucleotide_diversity(loop_snps)
    if add_maf:
      fst_list.append([scp[0],scp[1],loop_fst,maf])
    else:
      fst_list.append([scp[0],scp[1],loop_fst])
  return fst_list


if __name__ == "__main__":
  snp_list = read_list(filename = args.input)
  add_maf = True
  fst_list = calc_fst(snp_list = snp_list, add_maf = add_maf)
  with open(args.output, 'w', newline="") as f:
    for sublist in fst_list:
      #assign return value to variable
      tmp=f.write(str(sublist[0]))
      tmp=f.write('\t')
      tmp=f.write(str(sublist[1]))
      tmp=f.write('\t')
      tmp=f.write(str(sublist[2]))
      if add_maf:
        tmp=f.write('\t')
        tmp=f.write(str(sublist[3]))
      tmp=f.write('\n')
  print("Finished snp list %s" % (args.input))
