#!/bin/python
#Alan E. Yocca
#gff_exon_distribution.py
#get distribution of exon counts in a gff file
#get some dummy flags to get the values I want (how many single exon gene models are there)
#don't need to rely on regex to get gene names right?
#could trust the gff file

import csv
import sys
import argparse
import re

parser = argparse.ArgumentParser(prog='PROG')
parser.add_argument('--input', required=True, help='input gff file')
parser.add_argument('--transcript_tag', required=False, default = "mRNA", help='gff file third column tag for transcripts')
parser.add_argument('--exon_tag', required=False, default = "CDS", help='gff file third column tag for exons')
parser.add_argument('--trans_split', required=False, default = ".", help='split character for transcript to gene')
parser.add_argument('--gene_regex', required=False, default = "ID=([a-zA-Z0-9]*)\.", 
					help='regex to get the base gene name (drop transcript identifier)')
parser.add_argument('--output', required=True, help='output tsv file')
#eventually add an argument to define a smaller set of individuals

args = parser.parse_args()


def calc_ec_el_il(cg_breaks = [], cg_strand = "", dummy_transcript = "poop"):
  #if minus strand, reverse order of list
  if cg_strand == "-":
    cg_breaks = cg_breaks[::-1]

  if len(cg_breaks) == 0:
    print("Failed in calc_ec_el_il")
    print(cg_breaks)
    print(dummy_transcript)
    sys.exit()
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
    


def load_gff(filename = "", trans_tag = "args.transcript_tag", 
				exon_tag = "args.exon_tag", gene_regex = "args.gene_regex"):
  loop_dict = dict()
  cg_breaks = []
  cg = ""
  transcript = ""
  cg_aed = ""
  cg_length = ""
  cg_strand = ""
  with open(filename) as fh:
    for line in fh:
      if re.match("^#",line):
        continue
      line_array = line.strip().split("\t")
      if re.match(trans_tag,line_array[2]):  
        #initiate the transcript in dict
        #first element will be gene length
        try:
          transcript = re.search(gene_regex,line_array[8]).group(1)
        except AttributeError:
          print("Gene regex: %s" % (gene_regex))
          sys.exit("String: %s" % (line_array[8]))
        loop_dict[transcript] = [int(line_array[4]) - int(line_array[3])]
      elif re.match(exon_tag,line_array[2]):
        #print("Appending")
        loop_dict[transcript].append([line_array[3],line_array[4]])
  return loop_dict

def calc_gene_chars(gff_dict = dict(), trans_split = "trans_split"):
  #also filters for the longest transcript
  gene_char_dict = dict()
  #key will be basename gene (no transcript identifier)
  for transcript in gff_dict.keys():
    gene = transcript.split(trans_split)[0]
    #initiate if undefined
    tmp = gene_char_dict.setdefault(gene, [0,1])
    if gene_char_dict[gene][1] < gff_dict[transcript][0]:
      #overwrite
      (cg_exon_count, cg_exon_length, cg_intron_length) = calc_ec_el_il(cg_breaks = gff_dict[transcript][1:], dummy_transcript = transcript)
      gene_char_dict[gene] = (transcript, gff_dict[transcript][0], cg_exon_count, 
      							cg_exon_length,cg_intron_length)
    else:
      #initiate
      (cg_exon_count, cg_exon_length, cg_intron_length) = calc_ec_el_il(cg_breaks = gff_dict[transcript][1:], dummy_transcript = transcript)
      gene_char_dict[gene] = (transcript, gff_dict[transcript][0], cg_exon_count, 
      							cg_exon_length,cg_intron_length)
  return gene_char_dict


def filter_longest(dict = dict()):
  loop_dict = dict()
  for transcript in dict.keys():
    gene = transcript.split(".")[0]
    try:
      if loop_dict[gene][1] < dict[transcript][1]:
      	#overwrite
      	loop_dict[gene] = dict[transcript]
    except KeyError:
      #initiate
      loop_dict[gene] = dict[transcript]
  return loop_dict

def echo_single_exon(dict = dict()):
  sec = 0
  total = len(dict.keys())
  for transcript in dict.keys():
    if dict[transcript][2] == 1:
      sec += 1
  print("N single exon gene models: %s" % (sec))
  print("Prop single exon gene models: %s" % (sec/total))
  
def output_dict(out_filename = "", out_dict = dict()):
  with open(out_filename, "w") as out_fh:
    key_count = 1
    for key in out_dict.keys():
      if key_count == 1:
        #add header
        key_count = 0
        tmp = out_fh.write("Gene\tLength\tExon_Count\tExon_Length\tIntron_Length\n")
        tmp = out_fh.write(str(out_dict[key][0]))
        tmp = out_fh.write("\t")
        tmp = out_fh.write(str(out_dict[key][1]))
        tmp = out_fh.write("\t")
        tmp = out_fh.write(str(out_dict[key][2]))
        tmp = out_fh.write("\t")
        tmp = out_fh.write(str(out_dict[key][3]))
        tmp = out_fh.write("\t")
        tmp = out_fh.write(str(out_dict[key][4]))
        tmp = out_fh.write("\n")

if __name__ == "__main__":
  gff_dict = load_gff(filename = args.input, trans_tag = args.transcript_tag, 
  						exon_tag = args.exon_tag, gene_regex = args.gene_regex)

  out_dict = calc_gene_chars(gff_dict =gff_dict, trans_split = args.trans_split)
  print("len dict: %s" % (len(out_dict.keys())))
  
  echo_single_exon(out_dict)
  
  output_dict(out_filename = args.output, out_dict = out_dict)



