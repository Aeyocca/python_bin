#!/bin/python
#05-15-20
#fst_calc.py

import dendropy

#hmmm sliding window. non-overlapping?

#lets just do each site
#whats the structure of the character matrix
"""
all_snps = dendropy.StandardCharacterMatrix.get(
        path="srr_sub_snp.phylip",
        schema="phylip")

all_fst = []
for i in range(0,length(all_snps)):
  loop_snp = all_snps.export_character_subset(i)
  loop_fst = dendropy.calculate.popgenstat.nucleotide_diversity(loop_snp)
  all_fst.append[(i,loop_fst)]
"""  

phylip_dict = {}

with open("srr_sub_snp.phylip") as fh:
  next(fh)
  for line in fh:
    #split by whitespace, should result in two elements
    tmp_array = line.split()
    #ecotype now set to the snps as a string
    phylip_dict[tmp_array[0]] = tmp_array[1]

print("length phylip dictionary: %s" % (len(phylip_dict)))
print("length ecotype 430: %s" % (len(phylip_dict[430])))