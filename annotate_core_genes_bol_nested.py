#!/bin/python
#Alan E. Yocca
#07-20-20
#annotate_core_genes_bol_nested.py

import csv
import sys
import argparse
import re


parser = argparse.ArgumentParser(prog='PROG')
#parser.add_argument('--input', required=True, help='input vcf file')
parser.add_argument('--arb', required=False, default = 99, help='arbitrary ka/ks value')
parser.add_argument('--ka_ks', required=False, default = "/mnt/research/edgerpat_lab/AlanY/00_collab/02_core_pred/02_ka_ks/07_osat_brachy/03_codeml/tmp_ka_ks.txt", help = 'Ka/ks file, relative to wkdir')
#parser.add_argument('--ka_ks_file', required=True, help='input ka ks file')
#parser.add_argument('--window', type = int, default = 0, required=False, help='increase size of features by n bp')
parser.add_argument('--output', required=True, help='output tsv file')

#eventually add an argument to define a smaller set of individuals

args = parser.parse_args()

def load_vcf_file(filename = "", acc_idx_array = [i for i in range(9,19)], gene_idx = 2):
  pav_dict = dict()
  
  with open(filename) as fh:
    next(fh)
    for line in fh:
      line_array = line.strip().split("\t")
      #select specific accessions from header
      #define this gene as core or dispenable      
      if set([line_array[i] for i in acc_idx_array]) == {"1/1"}:
        pav_dict[line_array[gene_idx]] = {"Membership" : "Core"}
      else:
        pav_dict[line_array[gene_idx]] = {"Membership" : "Dispensable"}

  return pav_dict

  
def add_ka_ks(filename = "", pav_dict = dict(), arb_value = 99):
  print("Adding arbitrary value for genes without arabidopsis ortholog: %s" % (arb_value))
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
        ka_ks_dict[line_array[5]] = line_array[4]
        ka_dict[line_array[5]] = line_array[2]
        ks_dict[line_array[5]] = line_array[3]
    for gene in pav_dict.keys():
      #a little difference in the string of gene names, edit here
      #gene_no_trans = gene.split('-')[0].replace("t","g")
      #pretty sure all transcript identifiers just 0.1
      try:
        pav_dict[gene]["Ka_Ks"] = ka_ks_dict[gene]
        pav_dict[gene]["Ka"] = ka_dict[gene]
        pav_dict[gene]["Ks"] = ks_dict[gene]
      except:
        pav_dict[gene]["Ka_Ks"] = str(arb_value)
        pav_dict[gene]["Ka"] = str(arb_value)
        pav_dict[gene]["Ks"] = str(arb_value)
  return pav_dict




def calc_ec_el_il(cg_breaks = [], cg_strand = ""):
  #if minus strand, reverse order of list
  if cg_strand == "-":
    cg_breaks = cg_breaks[::-1]

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
      if re.match("mRNA",line_array[2]):
        #ugh, load in things, when the variables made?
          
        if len(cg_breaks) == 0:
          print("Passing first")
          pass
        else:
          (cg_exon_count, cg_exon_length, cg_intron_length) = calc_ec_el_il(cg_breaks = cg_breaks, cg_strand = cg_strand)
          
          try:
            if loop_dict[cg][1] < cg_length:
              loop_dict[cg] = (transcript, cg_length, cg_exon_count, cg_exon_length,
            				cg_intron_length)
          except KeyError:
            loop_dict[cg] = (transcript, cg_length, cg_exon_count, cg_exon_length,
            				cg_intron_length)

          """
          except:
            print("First gene")
            if test > 1:
              sys.exit("something else caused appending to fail")
            test +=1
          """
      
      #run calculations and append to dictionary
      #_QI=0|1|0|1|1|1|2|0|358;Parent=BnaA01g00060D2;_AED=0.08;ID=BnaA01g00060.1D2;Name=BnaA01g00060.1D2;_eAED=0.08
      #Name=BnaA03g40410.1D2;_eAED=0.19;_AED=0.19;ID=BnaA03g40410.1D2;Parent=BnaA03g40410D2;_QI=0|0.6|0.5|0.83|1|1|6|0|395
      
      #make new 2d array
        transcript = re.search("ID=([a-zA-Z0-9\-]*);",line_array[8]).group(1)
        cg = transcript.split("-")[0].replace("t","g")
        #cg = cg.replace("t","g")
        cg_breaks = []
        cg_length = int(line_array[4]) - int(line_array[3])
        cg_strand = line_array[6]
        if cg_length == 0:
          sys.exit("Zero length gene huh?? %s %s" % (cg, cg_length))

      if re.match("exon",line_array[2]):
        cg_breaks.append([line_array[3],line_array[4]])
  """
  with open("tmp_no_pav_info.txt", "w") as output:
    for item in no_pav_info_genes:
      output.write(item)
      output.write("\n")
  """
  #A Little clunky, but easiest way to replace info if found a longer isoform
  no_gene_info = 0
  longest_transcript_list = []
  """
  for gene in loop_dict.keys():
    longest_transcript_list.append(loop_dict[gene][0])
    #gene = transcript.split("-")[0].replace("t","g")
    try:
      pav_dict[gene].extend((loop_dict[gene][1:5]))
    except KeyError:
      no_gene_info += 1
  """
  membership = []
  remove = []
  for gene in pav_dict.keys():
    try:
      #(transcript, cg_length, cg_exon_count, cg_exon_length,cg_intron_length)
      pav_dict[gene]["Length"] = loop_dict[gene][1]
      pav_dict[gene]["Exon_Count"] = loop_dict[gene][2]
      pav_dict[gene]["Exon_Length"] = loop_dict[gene][3]
      pav_dict[gene]["Intron_Length"] = loop_dict[gene][4]

    except:
      no_gene_info += 1
      membership.append(pav_dict[gene]["Membership"])
      remove.append(gene)
  for gene in remove:
    del pav_dict[gene]
  print("No pav info for %s genes in gff file, removing" % (no_gene_info))
  print("%s core, %s dispensable" % (len([x for x in membership if x == "Core"]),
                                     len([x for x in membership if x != "Core"])))
  """
  with open("osat_longest_transcript.txt", "w") as output:
    for item in longest_transcript_list:
      output.write(item)
      output.write("\n")
  """

  return pav_dict

def parse_gff(filename = "", pav_dict = dict()):
  loop_dict = dict()
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
            loop_dict[cg] = (transcript, cg_length, cg_exon_count, cg_exon_length,
            				cg_intron_length)
          except KeyError:
            #print("Can't find gene: %s" % (cg))
            no_pav_info_genes.append(cg)
            no_pav_info += 1
      #make new 2d array
        cg = re.search("ID=gene:([a-zA-Z0-9]*);",line_array[8]).group(1)
        cg_breaks = []
        cg_length = int(line_array[4]) - int(line_array[3])
        if cg_length == 0:
          sys.exit("Zero length gene huh?? %s %s" % (cg, cg_length))

      if re.match("exon",line_array[2]):
        cg_breaks.append([line_array[3],line_array[4]])
  
  membership = []
  remove = []
  for gene in pav_dict.keys():
    try:
      #(transcript, cg_length, cg_exon_count, cg_exon_length,cg_intron_length)
      pav_dict[gene]["Length"] = loop_dict[gene][1]
      pav_dict[gene]["Exon_Count"] = loop_dict[gene][2]
      pav_dict[gene]["Exon_Length"] = loop_dict[gene][3]
      pav_dict[gene]["Intron_Length"] = loop_dict[gene][4]

    except:
      no_gene_info += 1
      membership.append(pav_dict[gene]["Membership"])
      remove.append(gene)
  for gene in remove:
    del pav_dict[gene]
  print("No pav info for %s genes in gff file, removing" % (no_gene_info))
  print("%s core, %s dispensable" % (len([x for x in membership if x == "Core"]),
                                     len([x for x in membership if x != "Core"])))

  return pav_dict

def calc_gc_per(filename = "", pav_dict = dict()):
  header = ""
  cds_seq = dict()
  with open(filename) as fh:
    for line in fh:
      if re.match("^>",line):
        line_array = line.strip().split(" ")
        header = line_array[0].strip(">").split("-")[0].replace("t","g")
      else:
        #also need to remove stop codons
        line = line.replace("*","")
        try:
          cds_seq[header] += line.strip()
        except KeyError:
          cds_seq[header] = line.strip()
  failed = 0
  membership = []
  for gene in pav_dict.keys():
    gc_per = ""
    try:
      gc_per = float(cds_seq[gene].count("G") + cds_seq[gene].count("C")) / float(len(cds_seq[gene]))
    except KeyError:
      failed += 1
      membership.append(pav_dict[gene]["Membership"])
    #ugh... going to just take the first transcript
    try:
      pav_dict[gene]["GC_Per"]
    except:
      #set if undefined
      pav_dict[gene]["GC_Per"] = gc_per
      
  print("No GC info for %s genes" % (failed))
  print("%s core, %s dispensable" % (len([x for x in membership if x == "Core"]),
                                     len([x for x in membership if x != "Core"])))

  return pav_dict
  
def calc_dinucleotide_per(filename = "", pav_dict = dict()):
  header = ""
  cds_seq = dict()
  dinuc_pair_list = ["AA","AC","AG","AT",
                     "CA","CC","CG","CT",
                     "GA","GC","GG","GT",
                     "TA","TC","TG","TT"]

  with open(filename) as fh:
    for line in fh:
      if re.match("^>",line):
        line_array = line.strip().split(" ")
        header = line_array[0].strip(">").split("-")[0].replace("t","g")
      else:
        #also need to remove stop codons
        line = line.replace("*","")
        try:
          cds_seq[header] += line.strip()
        except KeyError:
          cds_seq[header] = line.strip()
  failed = 0
  membership = []
  for gene in pav_dict.keys():
    try:
      gene_length = float(len(cds_seq[gene]))
      for dinuc_pair in dinuc_pair_list:
        #ugh... going to just take the first transcript 
        try:
          pav_dict[gene][dinuc_pair]
        except:
          #set if undefined
          pav_dict[gene][dinuc_pair] = float(cds_seq[gene].count(dinuc_pair)) / gene_length
    except KeyError:
      failed += 1
      membership.append(pav_dict[gene]["Membership"])  
  print("No dinucleotide info for %s genes" % (failed))
  print("%s core, %s dispensable" % (len([x for x in membership if x == "Core"]),
                                     len([x for x in membership if x != "Core"])))

  return pav_dict

if __name__ == "__main__":
  #load feature files
  pav_dict = dict()
  wkdir="/mnt/research/edgerpat_lab/AlanY/00_collab/02_core_pred/"
  pav_dict = load_vcf_file(filename = wkdir + "/01_reference/BOLEPan.pav.13062016.vcf")
  
  print(len(pav_dict.keys()))
  
  #pav_dict = calc_gc_per(filename = wkdir + "/IRGSP-1.0_cds_2020-06-03.fasta", pav_dict = pav_dict)
  #pav_dict = calc_dinucleotide_per(filename = wkdir + "/IRGSP-1.0_cds_2020-06-03.fasta", pav_dict = pav_dict)
  pav_dict = add_ka_ks(filename = args.ka_ks, 
                       pav_dict = pav_dict, arb_value = args.arb)
  #pav_dict = parse_gff(filename = wkdir + "/IRGSP-1.0_representative/transcripts_exon.gff", pav_dict = pav_dict)
  #pav_dict = add_exp(filename = wkdir + "/rice_exp_table.txt", pav_dict = pav_dict)              
  #pav_dict = add_aa_table(filename = wkdir + "/osat_longest_trans_aa_table.txt", pav_dict = pav_dict)
    
  with open(args.output, "w") as output:
    key_count = 1
    for key in pav_dict.keys():
      if key_count == 1:
        #add header
        key_count = 0
        output.write("Gene")
        for colname in pav_dict[key].keys():
          output.write("\t")
          output.write(colname)
        output.write("\n")
      output.write(key)
      for sub_key in pav_dict[key].keys():
        output.write("\t")
        output.write(str(pav_dict[key][sub_key]))
      output.write("\n")

  print("Finished")






