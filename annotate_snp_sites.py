#!/bin/python
#Alan E. Yocca
#05-19-20
#annotate_snp_sites.py
#load in a whole host of bed files and annotate the fst_sub_srr.csv file
#convert that to bed also... whats it called when bed lists every site?

import csv
import sys
import argparse

parser = argparse.ArgumentParser(prog='PROG')
parser.add_argument('--input', required=True, help='input tsv file')
parser.add_argument('--window', type = int, default = 0, required=False, help='increase size of features by n bp')
parser.add_argument('--output', required=True, help='output tsv file')

args = parser.parse_args()


def load_annotation(bed_file = str(), feature = str(), feature_dict = dict()):
  #loop through bed file
  #load annotation at each site into the feature dictionary
  print("Loading bed file: %s" % (bed_file))
  with open(bed_file) as fh:
    for line in fh:
      line_array = line.split()
      #alter window size
      line_array[1] = int(line_array[1]) - args.window
      if (line_array[1] < 0):
        line_array[1] = 0
      line_array[2] = int(line_array[2]) + args.window
      for position in range(int(line_array[1]),(int(line_array[2]) + 1)):
        chrom_pos = line_array[0] + "_" + str(position)
        try:
          feature_dict[chrom_pos] += "," + feature
        except:
          feature_dict[chrom_pos] = feature
  return feature_dict


def append_snp_annotation(snp_file = str(), feature_dict = dict(), output = str()):
  #loop through snp file
  print("Appending snp annotation: %s" % (snp_file))
  out_array = []
  with open(snp_file) as snp_fh:
    for line in snp_fh:
      line = line.strip()
      line_array = line.split(',')
      chrom_pos = line_array[0] + "_" + line_array[1]
      try:
        for feature in feature_dict[chrom_pos].split(","):
          out_array.append(line_array + [feature])
      except:
        out_array.append(line_array + ["Other"])

  print("Writing final output: %s" % (output))
  with open(output, "w", newline="") as fh:
    writer = csv.writer(fh, delimiter = "\t")
    writer.writerows(out_array)

if __name__ == "__main__":
  #load feature files
  fd = dict()
  fd = load_annotation(bed_file = "tair10_cns_gene.bed", feature = "CNS_gene", feature_dict = fd)
  fd = load_annotation(bed_file = "pav_cns.txt", feature = "PAV_CNS", feature_dict = fd)
  fd = load_annotation(bed_file = "posv_cns.txt", feature = "PosV_CNS", feature_dict = fd)
  fd = load_annotation(bed_file = "tair10_pav_cns_gene.bed", feature = "PAV_CNS_gene", feature_dict = fd)
  fd = load_annotation(bed_file = "tair10_posv_cns_gene_original.bed", feature = "PosV_CNS_gene_og", feature_dict = fd)
  fd = load_annotation(bed_file = "tair10_posv_cns_gene_new_closest.bed", feature = "PosV_CNS_gene_new", feature_dict = fd)
  fd = load_annotation(bed_file = "tair10_non_cns_gene.bed", feature = "non_CNS_gene", feature_dict = fd)
  fd = load_annotation(bed_file = "coll_cns.bed", feature = "Coll_CNS", feature_dict = fd)
  fd = load_annotation(bed_file = "tair10_coll_cns_gene.bed", feature = "Coll_CNS_gene", feature_dict = fd)

  fd = load_annotation(bed_file = "tair10_cns_gain_gene.bed", feature = "CNS_Gain", feature_dict = fd)
  fd = load_annotation(bed_file = "tair10_cns_loss_gene.bed", feature = "CNS_Loss", feature_dict = fd)
  fd = load_annotation(bed_file = "tair10_cns_nc_gene.bed", feature = "CNS_NC", feature_dict = fd)


  #annotate SNP file
  #sub 30 accessions
  append_snp_annotation(snp_file = args.input, feature_dict = fd, output = args.output)
  #all 1135 accessions
  #append_snp_annotation(snp_file = "fst_all_acc.csv", feature_dict = fd, output = "fst_all_acc_annotated.tsv")
  
  






