#!/bin/python
#Alan E. Yocca
#07-06-20
#parse_vac_gff.py

import csv
import sys
import argparse
import re


parser = argparse.ArgumentParser(prog='PROG')
parser.add_argument('--input', required=True, help='input gff file')
#parser.add_argument('--arb', required=False, default = 99, help='arbitrary ka/ks value')
#parser.add_argument('--ka_ks_file', required=True, help='input ka ks file')
#parser.add_argument('--window', type = int, default = 0, required=False, help='increase size of features by n bp')
parser.add_argument('--output', required=True, help='output tsv file')

args = parser.parse_args()


def calc_ec_el_il(cg_breaks = []):
  if len(cg_breaks) == 0:
    print("Failed in calc_ec_el_il")
    raise
  exon_length = 0
  intron_length = 0
  exon_count = 0
  #initialize
  previous_stop = 0
  for interval in cg_breaks:
    exon_count +=1
    exon_length += (int(interval[1]) - int(interval[0]) + 1)
    if exon_count > 1:
      intron_length += (int(interval[0]) - int(previous_stop))
    previous_stop = interval[1]
  return (exon_count, exon_length, intron_length)
    


def parse_gff(filename = "", pav_dict = dict()):
  no_pav_info = 0
  no_pav_info_genes = []
  test = 1
  cg_breaks = []
  cg = ""
  cg_aed = ""
  cg_length = ""
  with open(filename) as fh:
    for line in fh:
      if re.match("^#",line):
        continue
      line_array = line.strip().split("\t")
      if re.match("gene",line_array[2]):
        #ugh, load in things, when the variables made?
          
        if len(cg_breaks) == 0:
          print("Passing first")
          pass
        else:
          (cg_exon_count, cg_exon_length, cg_intron_length) = calc_ec_el_il(cg_breaks = cg_breaks)
          try:
            pav_dict[cg] = (cg_length, cg_exon_count, cg_exon_length,
            				cg_intron_length)
          except:
            #print("Can't find gene: %s" % (cg))
            no_pav_info_genes.append(cg)
            no_pav_info += 1
            #
        
        #make new 2d array
        cg = re.search("ID=([a-zA-Z0-9_\-\.]*);",line_array[8]).group(1)
        cg_breaks = []
        cg_length = int(line_array[4]) - int(line_array[3])
        if cg_length == 0:
          sys.exit("Zero length gene huh?? %s %s" % (cg, cg_length))

      if re.match("exon",line_array[2]):
        cg_breaks.append([line_array[3],line_array[4]])
  
  print("Finished parsing gff. Failed to find pav info for %s genes" % (no_pav_info))
  return pav_dict

def add_go_terms(filename = "", pav_dict = dict()):
  key_error = 0
  with open(filename) as fh:
    for line in fh:
      line_array = line.strip().split("\t")
      #add transcript identifier so we can find it
     # gene_trans = line_array[0].replace("D2",".1D2")
      try:
        pav_dict[line_array[0]].extend([line_array[1]])
      except KeyError:
        key_error+=1

  print("Key errors adding go terms: %s" % (key_error))
  return pav_dict

if __name__ == "__main__":
  pav_dict = dict()
  pav_dict = parse_gff(filename = args.input, pav_dict = pav_dict)

  with open(args.output, "w") as output:
    for key in pav_dict.keys():
      output.write(key)
      for item in pav_dict[key]:
        output.write("\t")
        output.write(str(item))
      output.write("\n")

  print("Finished")










