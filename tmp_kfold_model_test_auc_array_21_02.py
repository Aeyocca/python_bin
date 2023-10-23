#!/bin/python
#01-22-2021
#tmp_kfold_model_test_auc_array_21_02.py
#Alan E. Yocca

#Dummy script just to get svc to run on diff features

import csv
import sys
import argparse
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split, StratifiedKFold
from sklearn.naive_bayes import GaussianNB
from sklearn.metrics import accuracy_score
from sklearn import preprocessing
from sklearn.inspection import permutation_importance
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score
import sklearn
from sklearn.model_selection import cross_val_score
from sklearn.metrics import roc_curve, auc
import statistics
#not sure if all the above are necessary but hey, we got em

parser = argparse.ArgumentParser()
parser.add_argument("--slurm_array_id", required=True, help = "slurm_task_array_id")

args = parser.parse_args()


def balance_df(ml_df = []):
	core_sub = ml_df[ml_df.Membership == 1]
	ncore = len(core_sub)
	disp_sub = ml_df[ml_df.Membership == 0]
	ndisp = len(disp_sub)
	
	minority_class = ""
	if ncore < ndisp:
		minority_class = 1
	else:
		minority_class = 0
	
	minority_subset = ml_df[ml_df.Membership == minority_class]
	nmin = len(minority_subset)
	majority_subset = ml_df[ml_df.Membership == abs(minority_class - 1)].sample(n=nmin, random_state=1)
	ml_df_balanced = minority_subset.append(majority_subset)
	return(ml_df_balanced)

def calc_cv(ml_df = "", cv = 10, model = "", feature_list = []):
	#extract features and scale
	target=np.ravel(ml_df['Membership'])
	features = ml_df[feature_list]
	features_scale = preprocessing.scale(features)
	
	n_splits = cv
	skf = StratifiedKFold(n_splits=n_splits)
	cv = skf.split(features_scale, target)
	
	roc_auc = []
	for train, test in cv:
		probas_ = model.fit(features_scale[train], target[train]).predict_proba(features_scale[test])
		# Compute ROC curve and area the curve
		fpr, tpr, thresholds = roc_curve(target[test], probas_[:, 1])
		roc_auc.append(auc(fpr, tpr))
	
	auc_stdev = [statistics.mean(roc_auc), statistics.stdev(roc_auc)]
	return auc_stdev

def calc_auc(train_ml_df = "", test_ml_df = "", model = "", feature_list = ""):
	#extract features and scale
	train_target=np.ravel(train_ml_df['Membership'])
	train_features = train_ml_df[feature_list]
	train_features_scale = preprocessing.scale(train_features)
	
	test_target=np.ravel(test_ml_df['Membership'])
	test_features = test_ml_df[feature_list]
	test_features_scale = preprocessing.scale(test_features)
	
	#train
	probas_ = model.fit(train_features, train_target).predict_proba(test_features_scale)
	# Compute ROC curve and area the curve
	fpr, tpr, thresholds = roc_curve(test_target, probas_[:, 1])
	roc_auc = auc(fpr, tpr)
	
	#compare
	auc_stdev = [roc_auc, 0]
	return auc_stdev

def calc_auc_split(train_ml_df = "", test_ml_df = "", model = "", feature_list = ""):
	#extract features and scale
	train_target=np.ravel(train_ml_df['Membership'])
	train_features = train_ml_df[feature_list]
	train_features_scale = preprocessing.scale(train_features)
	
	test_target=np.ravel(test_ml_df['Membership'])
	test_features = test_ml_df[feature_list]
	test_features_scale = preprocessing.scale(test_features)
	
	#ugh I need to train test split otherwise hella over training on balance v unbal on sames species
	data_train_test, data_test_test, target_train_test, target_test_test = train_test_split(test_features_scale, 
		test_target, test_size = 0.10, random_state = 10)
	
	data_train_train, data_test_train, target_train_train, target_test_train = train_test_split(train_features_scale, 
		train_target, test_size = 0.10, random_state = 10)
	
	probas_ = model.fit(data_train_train, target_train_train).predict_proba(data_test_test)
	# Compute ROC curve and area the curve
	fpr, tpr, thresholds = roc_curve(target_test_test, probas_[:, 1])
	roc_auc = auc(fpr, tpr)
	
	#compare
	auc_stdev = [roc_auc, 0]
	return auc_stdev

def loop_model_spec(df_dict = "", tmp_array = 0):
	#model_list = ["SVC","GNB","RFT"]
	model_list = ["SVC"]
	species_list = ["Osat_Bal","Osat","Bdis_Bal","Bdis"]
	#feature_list = [ "GC_Per", "Ortho_Ka_Ks", "Ortho_Ka", "Ortho_Ks", "Para_Ka_Ks", "Para_Ka",
	#			   "Para_Ks", "Length", "Exon_Count", "Exon_Length", "Intron_Length", 
	#			   "AA", "AC", "AG", "AT", "CA", "CC", "CG", "CT", "GA", "GC", "GG", "GT", 
	#			   "TA", "TC", "TG", "TT", "Dup_Type"]
	feature_list = [ "GC_Per", "Para_Ka_Ks", "Para_Ka",
				   "Para_Ks", "Length", "Exon_Count", "Exon_Length", "Intron_Length", 
				   "AA", "AC", "AG", "AT", "CA", "CC", "CG", "CT", "GA", "GC", "GG", "GT", 
				   "TA", "TC", "TG", "TT"]
	out_dict = dict()
	model = ""
	for model_tag in model_list:
		#initialize output
		try:
			tmp = out_dict[model_tag]
		except KeyError:
			out_dict[model_tag] = dict()
		if model_tag == "SVC":
			model = SVC(kernel='linear', C=1, probability = True)
		elif model_tag == "GNB":
			model = GaussianNB()
		else:
			model = RandomForestClassifier(n_estimators=100)
			
		func_loop = 0
		for train_species in species_list:
			try:
				tmp = out_dict[model_tag][train_species]
			except KeyError:
				out_dict[model_tag][train_species] = dict()
			for test_species in species_list:
				func_loop += 1
				#next if not loop
				if func_loop != tmp_array:
					continue
			
				print("Starting " + model_tag + " " + train_species + " " + test_species)
				acc_stdev = []
				if train_species.split("_")[0] == test_species.split("_")[0]:
					if train_species == test_species:
						acc_stdev = calc_cv(ml_df = df_dict[train_species], 
										cv = 10, model = model, 
										feature_list = feature_list)
					else:
						#train test split to avoid overtraining
						acc_stdev = calc_auc_split(train_ml_df = df_dict[train_species],
										 		   test_ml_df = df_dict[test_species],
										 		   model = model,
										 		   feature_list = feature_list)				
				else:
					#cross data, can't do cv, just list acc, stdev = 0
					acc_stdev = calc_auc(train_ml_df = df_dict[train_species],
										 test_ml_df = df_dict[test_species],
										 model = model,
										 feature_list = feature_list)
				out_dict[model_tag][train_species][test_species] = acc_stdev
	return out_dict
	
def output_table(out_dict = dict(), filename = ""):
	#incase this is broken somehow lets put to stdout
	print(out_dict)
	header = ["Model","Train","Test","AUC-ROC"]
	with open(filename, "w") as out_fh:
		for i in range(len(header) - 1):
			out_fh.write(header[i] + "\t")
		out_fh.write(header[-1] + "\n")
		for model in out_dict.keys():
			for train in out_dict[model].keys():
				for test in out_dict[model][train].keys():
					out_fh.write(model + "\t" + train + "\t" + test + "\t")
					out_fh.write(str(out_dict[model][train][test][0]) + " +/- " + 
									 str(out_dict[model][train][test][1]) + "\n")

#want to get a nice output table of train, test, model, accuracy with +/- stdev
if __name__ == "__main__":
	loop = int(args.slurm_array_id) - 1

	#load in osat and bdis data
	osat_ml_df = pd.read_csv("osat_meta_21_02.txt", delimiter = "\t")
	osat_ml_df = osat_ml_df.dropna()
	osat_ml_df['Membership'] = osat_ml_df['Membership'].map({'Core': 1, 'Dispensable': 0})
	
	bdis_ml_df = pd.read_csv("bdis_meta_21_02.txt", delimiter = "\t")
	bdis_ml_df = bdis_ml_df.dropna()
	bdis_ml_df['Membership'] = bdis_ml_df['Membership'].map({'Core': 1, 'Dispensable': 0})

	#balance
	osat_ml_df_bal = balance_df(ml_df = osat_ml_df)
	bdis_ml_df_bal = balance_df(ml_df = bdis_ml_df)
	
	df_dict = {"Osat_Bal" : osat_ml_df_bal,
			   "Osat" : osat_ml_df,
			   "Bdis_Bal" : bdis_ml_df_bal,
			   "Bdis" : bdis_ml_df}
	out_dict = loop_model_spec(df_dict = df_dict, tmp_array = int(args.slurm_array_id))
	
	#output
	output_table(out_dict = out_dict, filename = "tmp_osat_bdis_cv_auc_table_" + str(loop) + ".txt")


