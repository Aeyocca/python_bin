#!/bin/python
#Alan E. Yocca
#07-20-20
#annotate_core_genes_osat.py

import csv
import sys
import argparse
import re

sys.exit("Use annotate_core_genes_osat_nested.py instead")

parser = argparse.ArgumentParser(prog='PROG')
#parser.add_argument('--input', required=True, help='input vcf file')
parser.add_argument('--arb', required=False, default = 99, help='arbitrary ka/ks value')
parser.add_argument('--ka_ks_file', required=False, default = "/mnt/research/edgerpat_lab/AlanY/00_collab/02_core_pred/02_ka_ks/07_osat_brachy/03_codeml/tmp_ka_ks.txt", help='input ka ks file')
#parser.add_argument('--window', type = int, default = 0, required=False, help='increase size of features by n bp')
parser.add_argument('--output', required=True, help='output tsv file')

#eventually add an argument to define a smaller set of individuals

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
        pav_dict[line_array[gene_idx]] = ["Core"]
      else:
        pav_dict[line_array[gene_idx]] = ["Dispensable"]
   
  return pav_dict
  
def add_ka_ks(filename = "", pav_dict = dict(), arb_value = 99):
  print("Adding arbitrary value for genes without ortholog: %s" % (arb_value))
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
  for gene in pav_dict.keys():
    #a little difference in the string of gene names, edit here
    #gene_no_trans = gene.split('-')[0].replace("t","g")
    #pretty sure all transcript identifiers just 0.1
    try:
      pav_dict[gene].extend([ka_ks_dict[gene],ka_dict[gene],ks_dict[gene]])
    except:
      pav_dict[gene].extend([arb_value,arb_value,arb_value])
      

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
      pav_dict[gene].extend((loop_dict[gene][1:5]))
    except:
      no_gene_info += 1
      membership.append(pav_dict[gene][1])
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



"""    
    #read through each gene, use some structure of each gene to calculate
    #certain stats with different functions
    #current_gene_start
    #current_gene_stop
    #cg_exon_boundaries
    #just get a 2d array?


    #how to get a single gene into 2d array?
    #need to extract the gene ID    
    
    #gene length
    #exon length
    #exon count
    #intron length
    #thats it right?
    #lets extract AED whilst we are at it
"""

def add_aa_table(filename = "", pav_dict = dict()):
  failed = 0
  aa_dict = dict()
  with open(filename) as fh:
    for line in fh:
      line_array = line.strip().split("\t")
      for i in range(1,len(line_array)):
        try:
          aa_dict[line_array[0]].extend([line_array[i]])
        except KeyError:
          aa_dict[line_array[0]] = [line_array[i]]
  membership = []
  remove = []
  for gene in pav_dict.keys():
    try:
      for item in aa_dict[gene]:
        pav_dict[gene].extend([item])
    except KeyError:
      failed += 1
      membership.append(pav_dict[gene][1])
      remove.append(gene)
  for gene in remove:
    del pav_dict[gene]
  print("No aa info added for %s genes" % (failed))
  print("%s core, %s dispensable" % (len([x for x in membership if x == "Core"]),
                                     len([x for x in membership if x != "Core"])))
  return pav_dict
	
def add_exp(filename = "", pav_dict = dict()):
  failed = 0
  exp_dict = dict()
  with open(filename) as fh:
    for line in fh:
      line_array = line.strip().split("\t")
      #strip_gene
      gene = line_array[0].split("-")[0].replace("t","g")
      exp_dict[gene] = line_array[3]
  membership = []
  remove = []
  for gene in pav_dict.keys():
    try:
      pav_dict[gene].extend([exp_dict[gene]])
    except KeyError:
      failed += 1
      membership.append(pav_dict[gene][1])
      remove.append(gene)
  for gene in remove:
    del pav_dict[gene]
  print("No expression info for %s genes" % (failed))
  print("%s core, %s dispensable" % (len([x for x in membership if x == "Core"]),
                                     len([x for x in membership if x != "Core"])))
  return pav_dict

if __name__ == "__main__":
  #load feature files
  pav_dict = dict()
  wkdir="/mnt/research/edgerpat_lab/AlanY/00_collab/02_core_pred/06_rice/"
  pav_dict = load_pav_file(filename = wkdir + "/pav_info_list.txt")
  pav_dict = add_ka_ks(filename = args.ka_ks_file,
                       pav_dict = pav_dict, arb_value = args.arb)
  pav_dict = parse_gff(filename = wkdir + "/IRGSP-1.0_representative/transcripts_exon.gff", pav_dict = pav_dict)
  pav_dict = add_exp(filename = wkdir + "/rice_exp_table.txt", pav_dict = pav_dict)              
  pav_dict = add_aa_table(filename = wkdir + "/osat_longest_trans_aa_table.txt", pav_dict = pav_dict)

  """
  pav_dict = add_ka_ks(filename = wkdir + "/02_ka_ks/05_bol_at/03_codeml/bol.at.ka.ks.txt", pav_dict = pav_dict, arb_value = args.arb)
  #think of add_feature function
  #perhaps some gff_feature_extraction

  """
  """
  #what do I want output to be? should I have headers? maybe eventually, not now
  with open(args.output, 'w') as csv_file:  
    writer = csv.writer(csv_file, delimiter = "\t")
    for key in pav_dict.keys():
      for 
      writer.writerow([key, value])
  """

  """  
  pav_dict = parse_gff(filename = wkdir + "/01_reference/Brassica_oleracea.BOL.42.gff3", pav_dict = pav_dict)
 
  pav_dict = add_aa_table(filename = wkdir + "/01_reference/bol_aa_table.txt", pav_dict = pav_dict)

  pav_dict = add_go_terms(filename = wkdir + "/02_ka_ks/05_bol_at/01_data/bol.at.go.tsv", pav_dict = pav_dict)

  """
  
  with open(args.output, "w") as output:
    for key in pav_dict.keys():
      output.write(key)
      for item in pav_dict[key]:
        output.write("\t")
        output.write(str(item))
      output.write("\n")

  print("Finished")






