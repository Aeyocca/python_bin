#!/bin/python
#05-15-20
#fst_calc_from_table.py

import sys
import csv
import numpy as np
import h5py
import pandas as pd
#import multiprocessing

print("top of hdf5")
f = h5py.File('imputed_snps_binary.hdf5','r')
positions = f['positions'][:]
snps = f['snps'][:]
accessions = f['accessions'][:]
chr_regions = f['positions'].attrs['chr_regions']

"""
  #get index of accessions of interest:
  id_ext_list = list()
  with open('eco_id_list.txt') as eco_fh:
    for line in eco_fh:
        #print("goal ID: %s" % (str(line.strip())))
        id_ext_list.append(str(line.strip()))
  ntaxa = str(len(id_ext_list))

  iter_index = 0
  indicies_of_interest = list()
  for acc in accessions:
    acc_strip = str(acc.strip()).replace("b","")
    acc_strip = acc_strip.replace("'","")
    if acc_strip in id_ext_list:
      indicies_of_interest.append(iter_index)
    iter_index += 1
"""
#capture all 1135 accessions
id_ext_list = list()
iter_index = 0
indicies_of_interest = list()
for acc in accessions:
  acc_strip = str(acc.strip()).replace("b","")
  acc_strip = acc_strip.replace("'","")
  indicies_of_interest.append(iter_index)
  id_ext_list.append(acc_strip)
  
ntaxa = str(len(id_ext_list))
#use range so can get numerical chromosome for output
#need external variable to keep snp indexing incrementing as we jump chromosomes
snp_index = 0
all_fst = []
progress = 0
  
#loop through snp structure, create new 2d list with snps, chr_pos?
#could do dictionary chr_pos = snps
#no do list because splitting that
snp_chr_pos = []
  
for i in range(0,len(chr_regions)):
  chromosome = "Chr" + str((i + 1))
  for chr_pos in positions[chr_regions[i][0]:chr_regions[i][1]]:
    #extract SNPs at this location, get into phylip string
    snp_subset = snps[snp_index][indicies_of_interest]
    snp_index+=1
    snp_chr_pos.append([chromosome,chr_pos,snp_subset])

print("Finished creating snp_chr_pos list")
#print("splitting list")
#chunk_size = 100000
#snp_splits = [snp_chr_pos[i:i + chunk_size] for i in range(0, len(snp_chr_pos), chunk_size)]

#see if we can do this without splitting
#fst table is MAF<tab>Fst
#since all sites should have 1135 values, all mafs will be represented
#load into dictionary
fst_dict = {}
with open("1135_fst_calc_full.txt") as fh:
  #lines = [ [split.strip() for split in line.split('\t')] for line in fh]
  #no list comprehension for this one
  for line in fh:
    split = line.split('\t')
    fst_dict[split[0]] = split[1].strip()

all_fst = []
progress = 0
for scp in snp_chr_pos:
  progress += 1
  if (progress % 100000) == 0:
    print("Position: %s" % progress)
  #calc maf
  maf = float(sum(scp[2])) / float(1135)
  loop_fst = fst_dict[str(maf)]
  if (maf != 1 and
  maf != 0):
    sys.exit("Found non-zero fst!: %s %s %s %s" % (scp[0],scp[1],loop_fst, maf))
  all_fst.append([scp[0],scp[1],loop_fst])

    
#use this to multiprocess
def calc_fst(snp_chr_pos = []):
  sub_process_fst = []
  for scp in snp_chr_pos:
    phylip_string= ntaxa + " 1\n"
    for eco_idx in range(0,len(id_ext_list)):
      ecotype_string = id_ext_list[eco_idx] + " " * (10 - len(id_ext_list[eco_idx]))
      phylip_string += ecotype_string + str(scp[2][eco_idx]) + "\n"
      loop_snps = dendropy.StandardCharacterMatrix.get(
      data=phylip_string,
      schema="phylip")
    loop_fst = dendropy.calculate.popgenstat.nucleotide_diversity(loop_snps)
    sub_process_fst.append([scp[0],scp[1],loop_fst])
  print("Finished a chunk of snps")
  return sub_process_fst
    
#now feeding calc_fst this splitup list
#pool = multiprocessing.Pool(processes = 20)
#all_fst = pool.map(calc_fst,snp_splits)
  
#print("Finished creating all_fst")

"""  
  for i in range(0,len(chr_regions)):
    chromosome = "Chr" + str((i + 1))
    for chr_pos in positions[chr_regions[i][0]:chr_regions[i][1]]:
      #extract SNPs at this location, get into phylip string
      snp_subset = snps[snp_index][indicies_of_interest]
      snp_index+=1
      #length of snp subset should be ntaxa
      #get ntaxa outside this block so don't have to compute 27 million times
      phylip_string= ntaxa + " 1\n"
      for eco_idx in range(0,len(id_ext_list)):
        ecotype_string = id_ext_list[eco_idx] + " " * (10 - len(id_ext_list[eco_idx]))
        phylip_string += ecotype_string + str(snp_subset[eco_idx]) + "\n"
      
      loop_snps = dendropy.StandardCharacterMatrix.get(
        data=phylip_string,
        schema="phylip")
      loop_fst = dendropy.calculate.popgenstat.nucleotide_diversity(loop_snps)
      #future me problem to figure out how to convert to genomic location...
      all_fst.append([chromosome,chr_pos,loop_fst])
      
      progress += 1
      if (progress % 100000) == 0:
        print("Position: %s" % progress)
"""
print("Starting to write final output")
#  with open("fst_sub_srr.csv", "w", newline="") as fh:
#    writer = csv.writer(fh)
#    writer.writerows(all_fst)

with open("fst_all_acc.csv", "w", newline="") as fh:
  writer = csv.writer(fh)
  writer.writerows(all_fst)
"""
if __name__ == "__main__":
  print("In main block")
  loop_hdf5()
  print("post hdf5 call")
"""


"""
#lets see if I can make a phylip with a string
phylip_string="2 3\n430       000\n500       913"
test_snps = dendropy.StandardCharacterMatrix.get(
        data=phylip_string,
        schema="phylip")

print(dendropy.calculate.popgenstat.nucleotide_diversity(test_snps))
"""