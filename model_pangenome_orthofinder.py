#!/bin/python
#model_pangenome_orthofinder.py
#Alan E. Yocca
#03-24-2023

#actually might be useful to others
#... does GENESPACE have a script for this?

import numpy as np
import statistics
import itertools
import multiprocessing
from multiprocessing import Process, Queue
import argparse

parser = argparse.ArgumentParser(prog='PROG')
parser.add_argument('--orthogroups', required=False, default = "Orthogroups.csv", help='Orthogroups.csv file')
parser.add_argument('--unassigned', required=False, default = "Orthogroups_UnassignedGenes.csv", help='Orthogroups_UnassignedGenes.csv file')
parser.add_argument('--output', required=False, default = "bbpan_model.txt", help='output')

args = parser.parse_args()


#logic:

#for a combination of n accessions, calculate efficiently (lets hope)
#the number of core and number of variable genes in this subset
#reducing the number of variable genes is the toughest part, how to
#count number of variable genes????


#can we do it by building out panmats?

#we can do this by averaging nvar per accession for each split
#that I think is how we did it with cranberry

#so build out panmaps that are subsettable by index
#need the column order to be consistent

#loop panmats, sub them on a certain index, and just check the number
#present in all, each index nvar/ncore will be averages across individuals

#hmm. isn't that what the other script does? but it doesn't count variable, huh

#loop the orthogroup file to build out pan mats
#but we need genes specific to each individual also, should be provided by orthofinder

#unassigned genes
#but isn't this count across accessions misleading? Its not total genes, but linking 
#variable genes is tough

#this problem goes away if genes are only represented in a single orthogroup
#to efficiently screen...
#load in all genes, and see if unique == ntotal

def check_uniq_assignment(orthogroup_file = ""):
	gene_list = []
	with open(orthogroup_file) as fh:
		#skip header
		next(fh)
		for line in fh:
			#some columns are empty, so strip is a bad idea
			la = line.replace("\n","").split("\t")[1:]
			for col in la:
				if col == "":
					continue
				for gene in col.split(", "):
					gene_list.append(gene)
	if len(gene_list) != len(list(set(gene_list))):
		print("WARNING: some genes appear in >1 orthogroup, this will inflate variable gene counts")
	print("%s %s" % (len(gene_list), len(list(set(gene_list)))))

def load_orthogroups(orthogroup_file = "", unassigned_file = ""):
	pav_mats = []
	header = ""
	#
	with open(orthogroup_file) as fh:
		#skip header
		header = next(fh)
		
		#initialize pav_mats object
		pav_mats = [[] for x in range(len(header.split("\t")) - 1)]
		
		for line in fh:
			la = line.replace("\n","").split("\t")[1:]
			bin_list = [1 if x != "" else 0 for x in la]
			for i in range(len(la)):
				for gene in la[i].split(", "):
					pav_mats[i].append(bin_list)
	#
	#wait its the same exact code, but duplicate
	#because checking headers match
	with open(unassigned_file) as fh:
		un_header = next(fh)
		if un_header != header:
			print("WARNING: header of %s and %s do not match" % 
								(orthogroup_file, unassigned_file))
		
		for line in fh:
			la = line.replace("\n","").split("\t")[1:]
			bin_list = [1 if x != "" else 0 for x in la]
			#this should only return a single presence...
			if sum(bin_list) >1:
				print("WARNING: multiple accessions present for this supposed unassigned gene... %s" 
						% (bin_list))
			for i in range(len(la)):
				for gene in la[i].split(", "):
					pav_mats[i].append(bin_list)
	np_mats = np.array(pav_mats)
	return(np_mats)

#def count_core_var(np_mats = [], combo = [], queue = ""):
def count_core_var(np_mats, combo, queue):
	#np_mats = np.array(pav_mats)
	
	ncore_array = []
	nvar_array = []
	ncol = len(combo)
	
	#need to multithread this loop...
	
	for mat in np_mats[combo]:
		#convert back to numpy array
		#need to eliminate this step
		#mat = np.array(mat)
		
		#does it have to be a numpy array?
		#can we listcomp from here?
		ncore = 0
		nvar = 0
		#print("looping all genes")
		
		for gene in mat:
			bin_list = [gene[i] for i in combo]
			if sum(bin_list) == ncol:
				ncore = ncore + 1
			else:
				nvar = nvar + 1
		
		#sub_mat = [mat[x][i] for i in combo for x in len(mat[0])]
		
		#sub_mat = mat[:,combo]
		#colsums = sub_mat.sum(axis = 1)
		
		
		ncore_array.append(ncore)
		nvar_array.append(nvar)
	
	#I think thats slightly faster
	#across 100 procs estimated 5 hours
	
	#lets hope that efficient enough
	#return for a given combination, the average of these arrays
	#why not return them all?
	#print(str(ncore))
	#print(str(nvar))
	#return(ncore_array, nvar_array)
	queue.put([ncore_array, nvar_array])

def loop_combos(nacc = int(), np_mats = []):
	core_means = []
	core_stds = []
	var_means = []
	var_stds = []
	accs = list(range(nacc))
	for i in accs:
		loop_core = []
		loop_var = []
		combos_iter = itertools.combinations(accs, i + 1)
		combos = list(combos_iter)
		print("%s accs; %s combos" % (i + 1, len(list(combos))))
		
		#should figure how to multithread this loop
		#could probably define outside this loop but lets keep it here
		
		#ok so spawning 2 million processes isn't working or I'm not doing it right
		#lets collect processes after every 1k jobs
		q = Queue()
		processes = []
		j = 0
		for combo in list(combos):
			j = j + 1
			if j % 1000:
				#collect processes and reset queue / procs
				for p in processes:
					lol = q.get()
					loop_core.append(lol[0])
					loop_var.append(lol[1])
				for p in processes:
					p.join()
				
				q = Queue()
				processes = []
			
			#print("combo %s" % (j))
			p = Process(target=count_core_var, args=(np_mats,list(combo),q))
			#(ncore,nvar) = count_core_var(np_mats = np_mats, 
			#										combo = list(combo))
			processes.append(p)
			p.start()
			#loop_core.append(ncore)
			#loop_var.append(nvar)
		
		#collect stragglers
		for p in processes:
			lol = q.get()
			loop_core.append(lol[0])
			loop_var.append(lol[1])
		
		for p in processes:
			p.join()
		
		#think we have to flatten them?
		loop_core = [item for sublist in loop_core for item in sublist]
		loop_var = [item for sublist in loop_var for item in sublist]
		
		core_means.append(sum(loop_core)/len(loop_core))
		core_stds.append(statistics.stdev(loop_core))
		var_means.append(sum(loop_var)/len(loop_var))
		var_stds.append(statistics.stdev(loop_var))
	
	#just return the output matrix
	header = ["Nacc","Mean_core","Stdev_core","Mean_var","Stdev_var"]
	out_mat = [header]
	out_data = [list(range(1,len(accs) + 1)), core_means, core_stds, var_means, var_stds]
	out_data_t = list(map(list, zip(*out_data)))
	out_data_t.insert(0,header)
	return(out_data_t)

if __name__ == "__main__":
	#check_uniq_assignment(orthogroup_file = "Orthogroups.csv")
	#PASSED! nice, much less complicated, no?
	
	np_mats = load_orthogroups(orthogroup_file = args.orthogroups, 
								unassigned_file = args.unassigned)
		
	out_mat = loop_combos(nacc = len(np_mats), np_mats = np_mats)
	
	with open(args.output, 'w') as out:
		for line in out_mat:
			for val in line[:-1]:
				tmp = out.write(str(val) + "\t")
			tmp = out.write(str(line[-1]) + "\n")