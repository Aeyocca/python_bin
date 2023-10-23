#!/bin/python
#tmp_rename_hc_chrom.py

#dummy script to rename honeycrisp genes

#Create dictionary for gene : chromosome

gene_to_chrom = dict()
#uh ohhh! gene names match across haplotypes!
#can we get species info from orthofinder?

with open("Honeycrisp_HAP1_braker1+2_combined_fullSupport_renamed.gff3") as fh:
	for line in fh:
		if line.startswith("#"):
			continue
		la = line.strip().split("\t")
		if la[2] == "mRNA":
			gene = la[8].replace("ID=","").split(";")[0]
			gene_to_chrom[gene] = la[0] + "_" + gene

out_dict = dict()
with open("Honeycrisp_HAP1_braker1+2_combined_fullSupport_renamed_filtered.pep.fa") as fh:
	header = ""
	for line in fh:
		if line.startswith(">"):
			header = line.strip().replace(">","")
		else:
			try:
				out_dict[header] = out_dict[header] + line.strip()
			except KeyError:
				out_dict[header] = line.strip()

with open("Mdom_HC_hap1.pep",'w') as out:
	for old_gene in out_dict.keys():
		tmp = out.write(">" + gene_to_chrom[old_gene] + "\n")
		tmp = out.write(out_dict[old_gene] + "\n")

#reset
gene_to_chrom = dict()
with open("Honeycrisp_HAP2_braker1+2_combined_fullSupport_renamed.gff3") as fh:
	for line in fh:
		if line.startswith("#"):
			continue
		la = line.strip().split("\t")
		if la[2] == "mRNA":
			gene = la[8].replace("ID=","").split(";")[0]
			gene_to_chrom[gene] = la[0] + "_" + gene

#reset
out_dict = dict()
with open("Honeycrisp_HAP2_braker1+2_combined_fullSupport_renamed_filtered.pep.fa") as fh:
	header = ""
	for line in fh:
		if line.startswith(">"):
			header = line.strip().replace(">","")
		else:
			try:
				out_dict[header] = out_dict[header] + line.strip()
			except KeyError:
				out_dict[header] = line.strip()

with open("Mdom_HC_hap2.pep",'w') as out:
	for old_gene in out_dict.keys():
		tmp = out.write(">" + gene_to_chrom[old_gene] + "\n")
		tmp = out.write(out_dict[old_gene] + "\n")
