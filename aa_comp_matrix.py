#!/bin/python
#Alan E. Yocca
#07-13-20
#aa_comp_matrix.py
#read in pep file and output a gene by amino acid matrix
#where each cell is the proportion of an amino acid in that gene
#hmmm should be able to do this for cds also huh, don't feel like loading in
#a translation table


import csv
import sys
import argparse
import re

parser = argparse.ArgumentParser(prog='PROG')
parser.add_argument('--input', required=True, help='input csv file')
parser.add_argument('--output', required=True, help='output tsv file')

args = parser.parse_args()

def load_fasta(filename = ""):
  print("Default stripping header after first space")
  header = ""
  pep_seq = dict()
  with open(filename) as fh:
    for line in fh:
      if re.match("^>",line):
        line_array = line.strip().split(" ")
        header = line_array[0].strip(">")
      else:
        #also need to remove stop codons
        line = line.replace("*","")
        try:
          pep_seq[header] += line.strip()
        except KeyError:
          pep_seq[header] = line.strip()
  #print("Length pep_seq %s" % (len(pep_seq.keys())))          
  print("Loaded %s genes from %s" % (len(pep_seq.keys()), filename))
  return pep_seq

def create_prop_mat(pep_seq = dict()):
  #hmm could just get the unique amino acids in the genome, nah
  #X is last because unspecified AA
  aa_array = ["A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y","X"]
  aa_mat = dict()
  for gene in pep_seq.keys():
    prop_array = [float(pep_seq[gene].count(x))/float(len(pep_seq[gene])) for x in aa_array]
    if sum(prop_array) < 0.95:
      print("Gene length: %s" % (float(len(pep_seq[gene]))))
      print("Gene: %s " % (gene))
      print("Seq: %s" % (pep_seq[gene]))
      print("%s" % (prop_array))
      sys.exit("Uh oh!: %s" % (sum(prop_array)))
    aa_mat[gene] = prop_array
  return aa_mat



if __name__ == "__main__":
  #pep_seq is dictionary with gene as key, sequence as value
  pep_seq = load_fasta(filename = args.input)

  #hmmm what format for the matrix. pandas?
  aa_mat = create_prop_mat(pep_seq = pep_seq)

  #dictionary will be simplest to work with, wait we need the amino acids to align
  #can handle that, they should if we make the values and array
  
  with open(args.output, "w") as output:
    for gene in aa_mat.keys():
      output.write(gene)
      for aa in aa_mat[gene]:
        output.write("\t")
        output.write(str(aa))
      output.write("\n")






