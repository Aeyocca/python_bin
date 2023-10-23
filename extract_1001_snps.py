#!/bin/python
#Alan E. Yocca
#02-06-20
#extract_1001_snps.py
#grab SNP matrix of the accessions you want
#feeding into PCA, so accessions are column names
#positions are row names
#02-07-20
#too mem intensive to save to objects, write to file and parse later

import numpy as np
import h5py
import pandas as pd

#still no good with argparse so just hard coding in files


#this was done earlier, now need to load in files, loop each line
#np.savetxt('1001_accessions.csv', h5py.File('imputed_snps_binary.hdf5')['accessions'], '%s', ',')
#print("Finished loading accessions")

#np.savetxt('1001_positions.csv', h5py.File('imputed_snps_binary.hdf5')['positions'], '%s', ',')
#print("Finished loading positions")

#np.savetxt('1001_snps.csv', h5py.File('imputed_snps_binary.hdf5')['snps'], '%s', ',')
#print("Finished loading snps")

#load files in
acc_list = list()
pos_list = list()
snp_list = list()

with open('1001_accessions.csv') as fh:
    for line in fh:
        line = str(line.strip()).replace("b","")
        line = line.replace("'","")
        #strip b'' whatever the hell that is...
        acc_list.append(line)

print("Finished loaded accessions")
	
with open('1001_positions.csv') as positions:
    for line in positions:
        pos_list.append(line.strip())

print("Finished	loaded positions")

with open('1001_snps.csv') as snps:
    for line in snps:
        snp_list.append(line.strip())

print("Finished	loaded SNPs")

#print("Finished loading in all lists")
#acc_list = list(h5py.File('imputed_snps_binary.hdf5','r')['accessions'])
#pos_list = list(h5py.File('imputed_snps_binary.hdf5','r')['positions'])
#snp_list = list(h5py.File('imputed_snps_binary.hdf5','r')['snps'])

#read in desired list of accessions
id_ext_list = list()
with open('eco_id_list.txt') as eco_fh:
    for line in eco_fh:
        #print("goal ID: %s" % (str(line.strip())))
        id_ext_list.append(str(line.strip()))


#get row numbers of accessions we want to extract
iter_index = 0
indicies_of_interest = list()
for acc in acc_list:
    if acc in id_ext_list:
        indicies_of_interest.append(iter_index)
    iter_index += 1

#print("indicies of interest: %s" % (indicies_of_interest))
#subset = [match for match in list_one if match in list_two]

#loop through position
#	for each position, loop through SNP row for that position
#	only keep values at positions matching row numbers we want to extract
#	careful index base is same
#	actually don't have to loop positions yet, can just loop SNPs
#	add positions as row names afterward

out_mat = list()
iterator = 0
#print("Beginning to loop snps....")
for snp_row in snp_list:
#    print("snp row: %s" % (snp_row))
    iterator += 1
    snp_row_list = snp_row.split(",")
    snp_of_interest = [snp_row_list[i] for i in indicies_of_interest]
    snp_of_interest = [int(i) for i in snp_of_interest]
    if sum(snp_of_interest) < 2 or sum(snp_of_interest) > 28:
        #only include if MAF > 0.05
        pass
    else:
        out_mat.append(snp_of_interest)
    if (iterator % 1000000) == 0:
        print("Finished iteration: %s" % (iterator))

print("Sites passing MAF > 0.05: %s" % (len(out_mat)))
#add column names and row names then output
out_df = pd.DataFrame(out_mat, columns = id_ext_list)

out_df.to_csv('srr_sub_snp_matrix.txt', index = True)





