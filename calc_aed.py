#!/bin/python
#Alan E. Yocca
#06-11-20
#calc_aed.py

#read in gff, output cumulative aed table and graph
#going to use some regex to get those values so this is probably going to be fragile

import re
import argparse
#from datetime import datetime
from matplotlib import pyplot as plt

parser = argparse.ArgumentParser(prog='PROG')
parser.add_argument('--input', required=True, help='input gff file')
#parser.add_argument('--nreads', type = int, required=True, help='number of reads to subset')
parser.add_argument('--output', required=True, help='output basename for table and graph file')
args = parser.parse_args()

def read_gff_aed(filename = ""):
  #get it into list?
  #need every mRNA annotation, 
  #just need to keep the 9th field
  lines_aed = []
  with open(filename) as fh:
    #lines = [ [split.strip() for split in line.split('\t')] for line in fh]
    #no list comprehension for this one
    for line in fh:
      split = line.split('\t')
      if re.match('^#',split[0]) is not None:
        continue
      if re.match('^>',split[0]) is not None:
        print("Looks like we hit fasta sequence, exiting read_gff_aed")
        break
      if split[2] == "mRNA":
        lines_aed.append(split[8])
    
  return lines_aed

def output_table(list = [], tag = ""):
  list_length = len(list)
  output = []
  for sub in range(0,11,1):
    sub = sub / 10
    #percentage of gene models with aed this value or less
    fraction = len([i for i in list if i <= sub]) / list_length
    output.append([sub,fraction])
  
  with open(tag + '_table.txt', 'w', newline="") as f:
    for line in output:
      #assign return value to variable
      tmp=f.write(str(line[0]))
      tmp=f.write('\t')
      tmp=f.write(str(line[1]))
      tmp=f.write('\n')

def output_graph(list = [], tag = ''):
  list_length = len(list)
  #get table at finer subdivision
  x_list = []
  y_list = []
  for sub in [x * 0.01 for x in range(0, 101)]:
    #sub = sub / 100
    fraction = len([i for i in list if i <= sub]) / list_length
    x_list.append(sub)
    y_list.append(fraction)
  
  plt.plot(x_list,y_list)
  plt.ylabel('Cumulative Fraction of Annotations')
  plt.xlabel('(e)AED')
  plt.savefig(tag + '_graph.png')


if __name__ == '__main__':
  lines_aed = read_gff_aed(args.input)
  aed_list = [float(re.search(';_AED=([0-9\.]+);',string).group(1)) for string in lines_aed ]
  no_aed = []
  exon_aed_list = []
  for string in lines_aed:
    try:
      exon_aed_list.append(float(re.search(';_eAED=([0-9\.]+);',string).group(1)))
    except:
      no_aed.append(re.search(';_eAED=([0-9\.]+);',string))
  #exon_aed_list = [float(re.search(';_eAED=([0-9\.]+);',string).group(1)) for string in lines_aed ]
  print("Gene models without eAED: %s" % (len(no_aed)))

  output_table(list = aed_list, tag = args.output + "_AED")
  output_table(list = exon_aed_list, tag = args.output + "_eAED")
  output_graph(list = aed_list, tag = args.output + "_AED")
  output_graph(list = exon_aed_list, tag = args.output + "_eAED")
