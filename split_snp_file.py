#!/bin/python
#05-15-20
#split_snp_file.py

import dendropy
import sys
import csv
import numpy as np
import h5py
import pandas as pd
import multiprocessing

np.set_printoptions(threshold=np.inf)

#split snps into 100 separate files to run fst calculation in parallel


print("top of hdf5")
f = h5py.File('imputed_snps_binary.hdf5','r')
positions = f['positions'][:]
snps = f['snps'][:]
accessions = f['accessions'][:]
chr_regions = f['positions'].attrs['chr_regions']

#capture all 1135 accessions
id_ext_list = list()
iter_index = 0
indicies_of_interest = list()
for acc in accessions:
  acc_strip = str(acc.strip()).replace("b","")
  acc_strip = acc_strip.replace("'","")
  indicies_of_interest.append(iter_index)
  id_ext_list.append(acc_strip)
  iter_index += 1
  
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
print("splitting list")
chunk_size = 20000
snp_splits = [snp_chr_pos[i:i + chunk_size] for i in range(0, len(snp_chr_pos), chunk_size)]

for i in range(0,len(snp_splits)):
  with open('snp_split_%s.txt' % i,'w', newline="") as f:
    for sublist in snp_splits[i]:
      #assign return value to variable
      tmp=f.write(str(sublist[0]))
      tmp=f.write('\t')
      tmp=f.write(str(sublist[1]))
      tmp=f.write('\t')
      #add iterator to avoid comma at the end of the line
      snp_iter = 1
      sub_list_len = len(sublist[2])
      for snp in sublist[2]:
        tmp=f.write(str(snp))
        if snp_iter < sub_list_len:
          tmp=f.write(',')
        snp_iter+=1
      tmp=f.write('\n')
  print("Finished subset %s" % (i))

"""
for i in range(0,len(snp_splits)):
  with open('snp_split_%s.txt' % i,'w', newline="") as f:
    for sublist in snp_splits[i]:
      #assign return value to variable
      tmp=f.write(str(sublist[0]))
      tmp=f.write('\t')
      tmp=f.write(str(sublist[1]))
      tmp=f.write('\t')
      for snp in sublist[2]:
        tmp=f.write(str(snp))
        tmp=f.write(',')
      tmp=f.write('\n')
  print("Finished subset %s" % (i))
"""
"""
for i in range(0,len(snp_splits)):
  with open('snp_split_%s.txt' % i,'w', newline="") as f:
    for sublist in snp_splits[i]:
      for snp in sublist[2]:
        f.write(str(snp))
        f.write(',')
        break
      break
  break
"""