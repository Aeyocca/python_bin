#!/bin/python
#Alan E. Yocca
#07-20-20
#load_ka_ks_osat.py

import csv
import sys
import re

def load_pav_file(filename = "", acc_idx_array = [i for i in range(1,454)], gene_idx = 0):
  pav_dict = dict()
  with open(filename) as fh:
    next(fh)
    for line in fh:
      line_array = line.strip().split("\t")
      #select specific accessions from header
      #define this gene as core or dispenable      
      if set([line_array[i] for i in acc_idx_array]) == {"1"}:
        pav_dict[line_array[gene_idx]] = {"Membership" : "Core"}
      else:
        pav_dict[line_array[gene_idx]] = {"Membership" : "Dispensable"}
  return pav_dict


def add_ka_ks(filename = "", pav_dict = dict(), comparator = ""):
  print("Reading ka/ks file: %s" % (filename))
  #4 is ka/ks 5 is bna gene
  ka_ks_dict = dict()
  ka_dict = dict()
  ks_dict = dict()
  with open(filename) as fh:
    for line in fh:
      next(fh)
      for line in fh:
        line_array = line.strip().split("\t")
        gene_no_trans = line_array[5].split('-')[0].replace("t","g")
        ka_ks_dict[gene_no_trans] = line_array[4]
        ka_dict[gene_no_trans] = line_array[2]
        ks_dict[gene_no_trans] = line_array[3]
  no_ortho = 0
  for gene in pav_dict.keys():
    #a little difference in the string of gene names, edit here
    #gene_no_trans = gene.split('-')[0].replace("t","g")
    #pretty sure all transcript identifiers just 0.1
    #will this initialize properly?
    try:
      pav_dict[gene][comparator] = {"Ka_Ks" : ka_ks_dict[gene],
      								"Ka" : ka_dict[gene],
      								"Ks" : ks_dict[gene]}
    except KeyError:
      #delete gene?
      #No, we can handle these when we output, just skip genes without a second level?
      no_ortho += 1
  print("No ka/ks for %s genes for %s" % (no_ortho, comparator))
  return pav_dict

if __name__ == "__main__":
  #load feature files
  pav_dict = dict()
  wkdir="/mnt/research/edgerpat_lab/AlanY/00_collab/02_core_pred/06_rice/"
 
  pav_dict = load_pav_file(filename = wkdir + "/pav_info_list.txt")
  pav_dict = add_ka_ks(filename = wkdir + "../02_ka_ks/07_osat_brachy/03_codeml/osat.osat.ka.ks.txt", 
  						comparator = "osat", pav_dict = pav_dict)
  pav_dict = add_ka_ks(filename = wkdir + "../02_ka_ks/07_osat_brachy/03_codeml/osat.oruf.ka.ks.txt",
  						comparator = "oruf", pav_dict = pav_dict)
  pav_dict = add_ka_ks(filename = wkdir + "../02_ka_ks/07_osat_brachy/03_codeml/osat.brachy.ka.ks.txt", 
  						comparator = "bdis", pav_dict = pav_dict)
  pav_dict = add_ka_ks(filename = wkdir + "../02_ka_ks/07_osat_brachy/03_codeml/osat.zea.ka.ks.txt",
  						comparator = "zmay", pav_dict = pav_dict)

   
  with open(wkdir + "/ka_ks_all.txt", "w") as output:
    key_count = 1
    output.write("Gene")
    output.write("\t")
    output.write("Membership")
    output.write("\t")
    output.write("Comparator")
    output.write("\t")
    output.write("Ka_Ks")
    output.write("\t")
    output.write("Ka")
    output.write("\t")
    output.write("Ks")
    output.write("\n")
    for gene in pav_dict.keys():
      if len(pav_dict[gene].keys()) == 1:
        continue
      for comparator in pav_dict[gene].keys():
        if comparator == "Membership":
          continue
        output.write(gene)
        output.write("\t")
        output.write(pav_dict[gene]["Membership"])
        output.write("\t")
        output.write(comparator)
        output.write("\t")
        output.write(pav_dict[gene][comparator]["Ka_Ks"])
        output.write("\t")
        output.write(pav_dict[gene][comparator]["Ka"])
        output.write("\t")
        output.write(pav_dict[gene][comparator]["Ks"])
        output.write("\n")
  #...that should do it..
  
  print("Finished")






