#!/bin/python
#Alan E. Yocca
#11-11-22
#annotate_core_genes_w85.py

import csv
import sys
import argparse
import re

#make every file a cmd line arg so can test different accessions
#defaults will be for the reference
parser = argparse.ArgumentParser(prog='PROG')
parser.add_argument('--arb', required=False, default = 99, help='arbitrary ka/ks value')
parser.add_argument('--pav_file', required=True, 
  help='SPACE SEPARATED! PAV matrix, first column is gene names, additional cols are binary calls for each accession')
parser.add_argument('--core_threshold', type = float, required=False, default = 1, help='proportion of individuals to call gene Core')
parser.add_argument('--softcore_threshold', type = float, required=False, default = 1, help='proportion of individuals to call gene soft core, default ignore')
parser.add_argument('--ortho_ka_ks', required=False, 
  help = 'Ka/ks file, absolute path, outgroup as reference')
parser.add_argument('--cds_file', required=False, 
  help = 'cds file, absolute path')
parser.add_argument('--para_ka_ks', required=False, 
  help = 'Ka/Ks file, self vs. self, absolute path')
parser.add_argument('--gff', required=False, 
  help = 'gff file, absolute path')
parser.add_argument('--gene_type', required=False, 
  help = 'gene_type file from MCScanX, absolute path')
parser.add_argument('--intron_fasta', required=False,
  help = 'fasta file of first introns')
parser.add_argument('--gtf', required=False,
  help = 'gtf file of expression. Could be comma separated list')
parser.add_argument('--orthogroups', required=False,
  help = 'orthofinder Orthogroups.tsv output file for phylostrat assignment')
parser.add_argument('--aa_mat', required=False,
  help = 'amino acid matrix')
parser.add_argument('--add_n_present', required=False, default = 0,
  help = 'if specified as 1, add column for n accessions gene is present in')

parser.add_argument('--output', required=True, help='output tsv file')

#eventually add an argument to define a smaller set of individuals

args = parser.parse_args()

def load_pav_file(filename = "", core_threshold = 1, softcore_threshold = 1, add_n_present = 0):
  pav_dict = dict()
  with open(filename) as fh:
    next(fh)
    for line in fh:
      line_array = line.replace("\n","").split(" ")
      #number of accessions will be ncol minus first one
      #nacc=len(line_array) - 1
      npresent=sum([int(x) for x in line_array[1:]])
      try:
        prop_present=npresent/(len(line_array) - 1)
      except ZeroDivisionError:
        sys.exit(line_array)
      if prop_present >= core_threshold:
      	pav_dict[line_array[0]] = {"Membership" : "Core"}
      elif prop_present >= softcore_threshold:
      	pav_dict[line_array[0]] = {"Membership" : "Softcore"}
      else:
        pav_dict[line_array[0]] = {"Membership" : "Variable"}
      if add_n_present:
        pav_dict[line_array[0]]["N_present"] = npresent
  return pav_dict
  
def add_ka_ks(filename = "", pav_dict = dict(), arb_value = 99, tag = ""):
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
        gene_no_trans = line_array[5].split('.')[0]
        ka_ks_dict[gene_no_trans] = line_array[4]
        ka_dict[gene_no_trans] = line_array[2]
        ks_dict[gene_no_trans] = line_array[3]
  for gene in pav_dict.keys():
    #a little difference in the string of gene names, edit here
    #gene_no_trans = gene.split('-')[0].replace("t","g")
    #pretty sure all transcript identifiers just 0.1
    try:
      pav_dict[gene][tag + "Ka_Ks"] = ka_ks_dict[gene]
      pav_dict[gene][tag + "Ka"] = ka_dict[gene]
      pav_dict[gene][tag + "Ks"] = ks_dict[gene]
    except:
      pav_dict[gene][tag + "Ka_Ks"] = arb_value
      pav_dict[gene][tag + "Ka"] = arb_value
      pav_dict[gene][tag + "Ks"] = arb_value
  
  return pav_dict

def calc_ec_el_il(cg_breaks = []):
  
  #print(cg_breaks)
  #if minus strand, reverse order of list
  #COLUMN INDICIES STAY THE SAME JSUT THE ORDER THEY APPEAR IN GFF diff for - strand
  
  #screw strand, just sort
  cg_breaks = [[int(x) for x in y] for y in cg_breaks]
  cg_breaks = sorted(cg_breaks, key=lambda x: (x[1], x[0]))
  
  if len(cg_breaks) == 0:
    print("Failed in calc_ec_el_il")
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

#think I will just redo this whole function...
def parse_gff(filename = "", pav_dict = dict()):
  #calc_el_ec_il works, just need to load those in accurately
  #GFF3 IS A HIGHLY IRREGULAR FORMAT, DOIN MY BEST WITH WHAT I GOT
  #make a gff datastructure
  #Assuming single mRNA per gene
  #hmm for w85 maybe we want to keep the mRNA part? doesnt matter, one per
  #lets keep it
  #strand also seems to not matter here
  gff_dict = dict()
  #key will be gene,
  #value will be cg breaks...?
  #just need exon breaks, TSS - TES, and CDS breaks
  #exon / CDS yes
  with open(filename) as fh:
    for line in fh:
      if re.match("^#",line):
        continue
      #collect exon lengths and cds length separately
      line_array = line.strip().split("\t")
      if re.match("mRNA",line_array[2]):
        #initialize gene
        gene = re.search("ID=([a-zA-Z0-9\-_\.]*);",line_array[8]).group(1)
        gff_dict[gene] = {"Exon_breaks" : [], "CDS_breaks" : [],
                          "Exon_count" : "", "Exon_length" : "",
                          "Intron_length" : "", "CDS_Length" : "",
                          "Gene_length" : int(line_array[4]) - int(line_array[3]) + 1,
                          "Strand" : line_array[6]}
      if re.match("exon",line_array[2]):
        gene = re.search("Parent=([a-zA-Z0-9\-_\.]*);",line_array[8]).group(1)
        gff_dict[gene]["Exon_breaks"].append([line_array[3],line_array[4]])
      if re.match("CDS",line_array[2]):
        gene = re.search("Parent=([a-zA-Z0-9\-_\.]*);",line_array[8]).group(1)
        gff_dict[gene]["CDS_breaks"].append([line_array[3],line_array[4]])
  
  skipped = 0
  for gene in gff_dict.keys():
    #(exn_count, exn_intra_length, exn_inter_length) = calc_ec_el_il(cg_breaks = gff_dict[gene]["Exon_breaks"])
    
    #keep only the first isoform
    #will end in \.[0-9]
    if bool(re.search('\.[2-9]$', gene)):
      continue
    
    #strip off the \.1
    strip_gene = re.sub("\.[1-9]$","",gene)

    (cds_count, cds_intra_length, cds_inter_length) = calc_ec_el_il(cg_breaks = gff_dict[gene]["CDS_breaks"])
    
    try:
      pav_dict[strip_gene]["Length"] = gff_dict[gene]["Gene_length"]
    except KeyError:
      #no pav info, skip
      skipped +=1
      continue
    #pav_dict[gene]["Exon_Count"] = exn_count
    #pav_dict[gene]["Exon_Length"] = exn_intra_length
    pav_dict[strip_gene]["Intron_Length"] = cds_inter_length
    pav_dict[strip_gene]["CDS_Length"] = cds_intra_length
    pav_dict[strip_gene]["CDS_Count"] = cds_count

  
  print("Skipped %s genes with no pav info, but had a gff entry" % (skipped))
  return pav_dict

def add_go_terms(filename = "", pav_dict = dict()):
  sys.exit("You don't want to add go terms, need to make nested version later if you do")
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

def add_aa_table(filename = "", pav_dict = dict()):
  failed = 0
  aa_dict = dict()
  aa_list = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", "X"]
  with open(filename) as fh:
    for line in fh:
      line_array = line.strip().split("\t")
      #gene = line_array[0].split(".")[0]
      gene = line_array[0]
      if gene in aa_dict:
      	continue
      for i in range(1,len(line_array)):
        try:
          aa_dict[gene].extend([line_array[i]])
        except KeyError:
          aa_dict[gene] = [line_array[i]]
  membership = []
  remove = []
  for gene in pav_dict.keys():
    try:
      for i in range(0,len(aa_dict[gene])):
        pav_dict[gene][aa_list[i]] = aa_dict[gene][i]
    except KeyError:
      failed += 1
      membership.append(pav_dict[gene]["Membership"])
      remove.append(gene)
  for gene in remove:
    del pav_dict[gene]
  print("No aa info added for %s genes" % (failed))
  print("%s core, %s soft/variable" % (len([x for x in membership if x == "Core"]),
                                     len([x for x in membership if x != "Core"])))
  return pav_dict

def calc_gc_per(filename = "", pav_dict = dict()):
  header = ""
  cds_seq = dict()
  with open(filename) as fh:
    for line in fh:
      if re.match("^>",line):
        line_array = line.strip().split(" ")
        header = line_array[0].strip(">").split(".")[0]
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
       gc_per = float(cds_seq[gene].count("G") + cds_seq[gene].count("C") +
      				 cds_seq[gene].count("g") + cds_seq[gene].count("c")) / float(len(cds_seq[gene]))
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
  print("%s core, %s soft/variable" % (len([x for x in membership if x == "Core"]),
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
        header = line_array[0].strip(">").split(".")[0]
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
          #why do I have this try except here?
          pav_dict[gene][dinuc_pair] = float(cds_seq[gene].upper().count(dinuc_pair)) / gene_length
    except KeyError:
      failed += 1
      membership.append(pav_dict[gene]["Membership"])  
  print("No dinucleotide info for %s genes" % (failed))
  print("%s core, %s soft/variable" % (len([x for x in membership if x == "Core"]),
                                     len([x for x in membership if x != "Core"])))
  
  return pav_dict

def add_dup_type(filename = "", pav_dict = dict()):
  failed = 0
  dup_dict = dict()
  with open(filename) as fh:
    for line in fh:
      line_array = line.strip().split("\t")
      dup_dict[line_array[0]] = line_array[1]
  membership = []
  remove = []
  for gene in pav_dict.keys():
    try:
      pav_dict[gene]["Dup_Type"] = dup_dict[gene]
    except KeyError:
      failed += 1
      membership.append(pav_dict[gene]["Membership"])
      remove.append(gene)
  for gene in remove:
    del pav_dict[gene]
  print("No Duplication info for %s genes" % (failed))
  print("%s core, %s dispensable" % (len([x for x in membership if x == "Core"]),
                                     len([x for x in membership if x != "Core"])))
  return pav_dict

def add_dup_type_one_hot_encoding(filename = "", pav_dict = dict()):
  #per reviewer request, try one hot encoding instead of single "dup_type class"
  #Type of dup	Code	Number
  #Singleton	0	
  #Dispersed	1	
  #Proximal	2	
  #Tandem	3	
  #WGD or segmental	4	
  trans_dict = {"0" : "Singleton",
  				"1" : "Dispersed",
  				"2" : "Proximal",
  				"3" : "Tandem",
  				"4" : "WGD"}
  failed = 0
  dup_dict = dict()
  with open(filename) as fh:
    for line in fh:
      line_array = line.strip().split("\t")
      dup_dict[line_array[0]] = line_array[1]
  membership = []
  remove = []
  for gene in pav_dict.keys():
    for dup_type in trans_dict.keys():
      try:
        if dup_dict[gene] == dup_type:
          pav_dict[gene][trans_dict[dup_type]] = 1
        else:
          pav_dict[gene][trans_dict[dup_type]] = 0
      except KeyError:
        failed += 1
        membership.append(pav_dict[gene]["Membership"])
        remove.append(gene)
  for gene in remove:
    del pav_dict[gene]
  print("No Duplication info for %s genes" % (failed))
  print("%s core, %s soft/variable" % (len([x for x in membership if x == "Core"]),
                                     len([x for x in membership if x != "Core"])))
  return pav_dict

def calc_intron_one_stats(filename = "", pav_dict = dict()):
  header = ""
  intron_seq = dict()
  with open(filename) as fh:
    for line in fh:
      if re.match("^>",line):
        line_array = line.strip().split(" ")
        header = line_array[0].strip(">").split(".")[0]
      else:
        #also need to remove stop codons
        line = line.replace("*","")
        try:
          intron_seq[header] += line.strip()
        except KeyError:
          intron_seq[header] = line.strip()
  failed = 0
  membership = []
  
  #print("Length pav_dict: %s" % (len(pav_dict.keys())))
  #print("Length intron_seq: %s" % (len(intron_seq.keys())))
  
  for gene in pav_dict.keys():
    gc_per = ""
    try:
      if len(intron_seq[gene]) == 0:
        gc_per = 0
      else:
        gc_per = float(intron_seq[gene].count("G") + intron_seq[gene].count("C") +
      				 intron_seq[gene].count("g") + intron_seq[gene].count("c")) / float(len(intron_seq[gene]))
    except KeyError:
      failed += 1
      membership.append(pav_dict[gene]["Membership"])
    #ugh... going to just take the first transcript
    try:
      pav_dict[gene]["Intron_one_GC"]
    except:
      #set if undefined
      pav_dict[gene]["Intron_one_GC"] = gc_per
      pav_dict[gene]["Intron_one_Length"] = len(intron_seq[gene])
  
  print("No First intron info for %s genes" % (failed))
  print("%s core, %s soft/variable" % (len([x for x in membership if x == "Core"]),
                                     len([x for x in membership if x != "Core"])))
  
  return pav_dict

def add_exp(filename = "", pav_dict = dict(), tag = ""):
  #Assuming stringtie output
  #Adding TPM and FPKM
  failed = 0
  new_gene = 0
  exp_dict = dict()
  with open(filename) as fh:
    for line in fh:
      if line.startswith("#"):
        continue
      line_array = line.strip().split("\t")
      if line_array[2] == "transcript":
        gtf_array = line_array[8].split("; ")
        #0 - gene_id
        #1 - transcript_id
        #2 - reference_id
        #3 - ref_gene
        #4 - cov
        #5 - FPKM
        #6 - TPM
        trans_id = [x for x in gtf_array if "transcript_id" in x]
        FPKM = [x for x in gtf_array if "FPKM" in x]
        TPM = [x for x in gtf_array if "TPM" in x]
        
        gene = re.search('"([a-zA-Z0-9\-_]*)"',str(trans_id)).group(1)
        fpkm_value = re.search('"([\.0-9]*)"',str(FPKM)).group(1)
        tpm_value = re.search('"([\.0-9]*)"',str(TPM)).group(1)
        exp_dict[gene] = {"FPKM_" + tag : fpkm_value, "TPM_" + tag : tpm_value}
        
  membership = []
  remove = []
  for gene in pav_dict.keys():
    try:
      pav_dict[gene]["FPKM_" + tag] = exp_dict[gene]["FPKM_" + tag]
      pav_dict[gene]["TPM_" + tag] = exp_dict[gene]["TPM_" + tag]
    except KeyError:
      failed += 1
      #these aren't failed, just no expression listed, list a zero
      #keeping variable names though
      membership.append(pav_dict[gene]["Membership"])
      remove.append(gene)
  for gene in remove:
    pav_dict[gene]["FPKM_" + tag] = 0
    pav_dict[gene]["TPM_" + tag] = 0
  print("No expression info for %s genes, setting to zero" % (failed))
  print("%s core, %s soft/variable" % (len([x for x in membership if x == "Core"]),
                                     len([x for x in membership if x != "Core"])))
  return pav_dict

def phylostrat(filename = "", pav_dict = dict()):

  #probably just a temporary function so going to hard code a few things


  with open(filename) as fh:
    #skip first line
    header = next(fh).replace("\n","").split("\t")[1:]
    for line in fh:
      #cant strip because of empty columns
      line_array=line.replace("\n","").split("\t")[1:]
      #assuming reference column
      ref_idx = 4
      ref_genes = line_array[ref_idx].split(", ")
      #strata_order = ["Atrich","Ppat","Osat","Ppers","vmacro"]
      strat_ages =  [3, 2, 4, 1, 0]
      #whats the oldest strata for thi orthogroup
      present_strats = [strat_ages[i] for i in range(len(line_array)) if len(line_array[i]) > 0 ]
      oldest_strat = max(present_strats)
      if len(ref_genes[0]) == 0:
        continue
      
      for gene in ref_genes:
        try:
          if pav_dict[gene]["Phylostrata"] < oldest_strat:
            pav_dict[gene]["Phylostrata"] = oldest_strat
        except KeyError:
          pav_dict[gene]["Phylostrata"] = oldest_strat
  
  #if not in the orthogroups file, phylostrata 0
  for gene in pav_dict.keys():
    if "Phylostrata" not in pav_dict[gene]:
      pav_dict[gene]["Phylostrata"] = 0  
  return pav_dict


if __name__ == "__main__":

  #load feature files
  pav_dict = dict()
  
  softcore = args.softcore_threshold
  if softcore == "":
    softcore = args.core_threshold
  #if softcore > args.core_theshold:
  #  print("softcore > core thresholds, all softcore will be core")
  pav_dict = load_pav_file(filename = args.pav_file, core_threshold = args.core_threshold, softcore_threshold = softcore, add_n_present = args.add_n_present)
  
  if args.cds_file is not None:
    pav_dict = calc_gc_per(filename = args.cds_file, pav_dict = pav_dict)
  
  
  if args.cds_file is not None:
    pav_dict = calc_dinucleotide_per(filename = args.cds_file, pav_dict = pav_dict)
  
  if args.ortho_ka_ks is not None:
    pav_dict = add_ka_ks(filename = args.ortho_ka_ks, tag = "Ortho_",
                       pav_dict = pav_dict, arb_value = args.arb)
  if args.para_ka_ks is not None:
    pav_dict = add_ka_ks(filename = args.para_ka_ks, tag = "Para_",
                       pav_dict = pav_dict, arb_value = args.arb)
  if args.gff is not None:
    pav_dict = parse_gff(filename = args.gff, pav_dict = pav_dict)
  
  if args.gtf is not None:
    gtf_list = args.gtf.split(",")
    for file in gtf_list:
      tag = file.split("/")[-1]
      print(tag)
      pav_dict = add_exp(filename = file, pav_dict = pav_dict, tag = tag)
  
  if args.orthogroups is not None:
    pav_dict = phylostrat(filename = args.orthogroups, pav_dict = pav_dict)
    
  #print(list(pav_dict.keys())[0])
  
  #pav_dict = add_exp(filename = wkdir + "/rice_exp_table.txt", pav_dict = pav_dict)     
  if args.aa_mat is not None:  
    pav_dict = add_aa_table(filename = args.aa_mat, pav_dict = pav_dict)
  #pav_dict = add_dup_type(pav_dict = pav_dict, filename = args.gene_type)
  
  if args.gene_type is not None:
    pav_dict = add_dup_type_one_hot_encoding(pav_dict = pav_dict, filename = args.gene_type)
  
  if args.intron_fasta is not None:
    pav_dict = calc_intron_one_stats(pav_dict = pav_dict, filename = args.intron_fasta)

  print(len(list(pav_dict.keys())))
  
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


