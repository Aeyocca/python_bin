#!/bin/python
#Alan E. Yocca
#07-25-20
#collect_ka_ks_table.py

import csv
import sys
import argparse
import re


parser = argparse.ArgumentParser(prog='PROG')
#parser.add_argument('--input', required=True, help='input vcf file')
parser.add_argument('--arb', required=False, default = 99, help='arbitrary ka/ks value')
parser.add_argument('--output', required=True, help='output tsv file')

args = parser.parse_args()

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


def add_ka_ks(filename = "", pav_dict = dict(), arb_value = 99, comparator = ""):
  print("Adding arbitrary value for genes without ortholog: %s" % (arb_value))
  print("Reading ka/ks file: %s" % (filename))
  print("For comparator: %s" % (comparator))
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
  for gene in pav_dict.keys():
    #a little difference in the string of gene names, edit here
    #gene_no_trans = gene.split('-')[0].replace("t","g")
    #pretty sure all transcript identifiers just 0.1
    try:
      pav_dict[gene][comparator]  = {"Ka_Ks" : ka_ks_dict[gene]}
      pav_dict[gene][comparator]["Ka"] = ka_dict[gene]
      pav_dict[gene][comparator]["Ks"] = ks_dict[gene]
    except:
      pav_dict[gene][comparator]  = {"Ka_Ks" : arb_value}
      pav_dict[gene][comparator]["Ka"] = arb_value
      pav_dict[gene][comparator]["Ks"] = arb_value
  return pav_dict


if __name__ == "__main__":
  #load feature files
  pav_dict = dict()
  wkdir="/mnt/research/edgerpat_lab/AlanY/00_collab/02_core_pred/06_rice/"
  pav_dict = load_pav_file(filename = wkdir + "/pav_info_list.txt")
  pav_dict = add_ka_ks(filename = wkdir + "../02_ka_ks/07_osat_brachy/03_codeml/osat.osat.ka.ks.txt", 
                       pav_dict = pav_dict, arb_value = args.arb, comparator = "osat")

  pav_dict = add_ka_ks(filename = wkdir + "../02_ka_ks/07_osat_brachy/03_codeml/osat.oruf.ka.ks.txt", 
                       pav_dict = pav_dict, arb_value = args.arb, comparator = "oruf")
                       
  pav_dict = add_ka_ks(filename = wkdir + "../02_ka_ks/07_osat_brachy/03_codeml/osat.brachy.ka.ks.txt", 
                       pav_dict = pav_dict, arb_value = args.arb, comparator = "bdis")


  pav_dict = add_ka_ks(filename = wkdir + "../02_ka_ks/07_osat_brachy/03_codeml/osat.zea.ka.ks.txt", 
                       pav_dict = pav_dict, arb_value = args.arb, comparator = "zmay")

  """
pav_dict = add_ka_ks(filename = wkdir + "../02_ka_ks/07_osat_brachy/03_codeml/osat.osat.ka.ks.txt", pav_dict = pav_dict, arb_value = 99, comparator = "osat")

pav_dict = add_ka_ks(filename = wkdir + "../02_ka_ks/07_osat_brachy/03_codeml/osat.oruf.ka.ks.txt", pav_dict = pav_dict, arb_value = 99, comparator = "oruf")

pav_dict = add_ka_ks(filename = wkdir + "../02_ka_ks/07_osat_brachy/03_codeml/osat.brachy.ka.ks.txt", pav_dict = pav_dict, arb_value = 99, comparator = "brachy")


with open(args.output, "w") as output:
  key_count = 1
  for key in pav_dict.keys():
    if key_count == 1:
      #add header
      key_count = 0
      output.write("Gene")
      output.write("\t")
      output.write("Comparator")
      for colname in list(pav_dict[key].keys()):
        output.write("\t")
        if colname == "Membership":
          output.write(colname)
        else:
          for data_type in pav_dict[key][colname].keys():
            output.write("\t")
            output.write(colname + "_" + data_type)
      output.write("\n")
    output.write(key)
    for sub_key in pav_dict[key].keys():
      output.write("\t")
      if sub_key == "Membership":
        output.write(sub_key)
      else:
        for data_type in pav_dict[key][sub_key].keys():
          output.write("\t")
          output.write(pav_dict[key][sub_key][data_type])
    output.write("\n")
  """                 

  #sucks to suck if you want to come back to this, but works for now      
  with open(args.output, "w") as output:
    key_count = 1
    for key in pav_dict.keys():
      if key_count == 1:
        #add header
        key_count = 0
        output.write("Gene")
        output.write("\t")
        output.write("Comparator")
        for colname in list(pav_dict[key].keys())[0:2]:
          if colname == "Membership":
            output.write("\t")
            output.write(colname)
          else:
            for data_type in pav_dict[key][colname].keys():
              output.write("\t")
              output.write(data_type)
        output.write("\n")
      
      #for comparator, print
      #gene \t membership \t comparator \t ka \t ks, etc
      #output.write(key)      
      for sub_key in list(pav_dict[key].keys())[1:]:
        output.write(key)
        output.write("\t")
        output.write(sub_key)
        output.write("\t")
        output.write(pav_dict[key]["Membership"])
        for data_type in pav_dict[key][sub_key].keys():
          output.write("\t")
          output.write(str(pav_dict[key][sub_key][data_type]))
        output.write("\n")

  print("Finished")
















