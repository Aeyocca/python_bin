#!/bin/python
#01-22-2021
#bdis_model_test_21_01.py
#Alan E. Yocca


import csv
import sys
import argparse
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.naive_bayes import GaussianNB
from sklearn.metrics import accuracy_score
from sklearn import preprocessing
from sklearn.inspection import permutation_importance
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score
import sklearn
from sklearn.model_selection import cross_val_score

parser = argparse.ArgumentParser(prog='PROG')
parser.add_argument("--slurm_array_id", required=True, help = "slurm_task_array_id")
args = parser.parse_args()

#not sure if all the above are necessary but hey, we got em

def load_ml_df():
	#read in dataframe
	#bdis_meta_osat_cds.txt
	wkdir = "/mnt/research/edgerpat_lab/AlanY/00_collab/02_core_pred/"
	
	bdis_cds_ml_df = pd.read_csv(wkdir + "/04_brachy/bdis_meta_bdis_cds.txt", delimiter = "\t")
	
	#drop any row with NA value
	bdis_cds_ml_df = bdis_cds_ml_df.dropna()
	
	#convert membership to binary (helps with testing accuracies later)
	bdis_cds_ml_df['Membership'] = bdis_cds_ml_df['Membership'].map({'Core': 1, 'Dispensable': 0})
	
	#load in para_cds, against bdis... yea...
	#want to extract ka/ks from this, rename, cbind to above
	bdis_osat_cds_ml_df = pd.read_csv(wkdir + "/04_brachy/bdis_meta_osat_cds.txt", delimiter = "\t")
	#df.rename(columns={'oldName1': 'newName1', 'oldName2': 'newName2'})
	bdis_osat_cds_ml_df = bdis_osat_cds_ml_df.rename(columns={"Ka_Ks" : "Ortho_Ka_Ks", 
															  "Ka" : "Ortho_Ka",
															  "Ks" : "Ortho_Ks"})
	#add these three columns
	ortho_ka_ks_sub = bdis_osat_cds_ml_df[["Gene", "Ortho_Ka_Ks", "Ortho_Ka", "Ortho_Ks"]]
	bdis_ml_df = bdis_cds_ml_df.merge(ortho_ka_ks_sub)
	
	osat_cds_ml_df = pd.read_csv(wkdir + "/06_rice/osat_meta_osat_cds.txt", delimiter = "\t")
	
	#drop any row with NA value
	osat_cds_ml_df = osat_cds_ml_df.dropna()
	
	#convert membership to binary (helps with testing accuracies later)
	osat_cds_ml_df['Membership'] = osat_cds_ml_df['Membership'].map({'Core': 1, 'Dispensable': 0})
	
	#load in para_cds, against bdis... yea...
	#want to extract ka/ks from this, rename, cbind to above
	osat_bdis_cds_ml_df = pd.read_csv(wkdir + "/06_rice/osat_meta_bdis_cds.txt", delimiter = "\t")
	#df.rename(columns={'oldName1': 'newName1', 'oldName2': 'newName2'})
	osat_bdis_cds_ml_df = osat_bdis_cds_ml_df.rename(columns={"Ka_Ks" : "Ortho_Ka_Ks", 
															  "Ka" : "Ortho_Ka",
															  "Ks" : "Ortho_Ks"})
	#add these three columns
	ortho_ka_ks_sub = osat_bdis_cds_ml_df[["Gene", "Ortho_Ka_Ks", "Ortho_Ka", "Ortho_Ks"]]
	osat_ml_df = osat_cds_ml_df.merge(ortho_ka_ks_sub)
	
	df_dict = {"osat_df" : osat_ml_df,
			   "bdis_df" : bdis_ml_df}
	return(df_dict)

def replace_arb_value(ml_df = [], replace_value = 0, arb_value = 99):
	#replace all arbitrary values with something different
	#test model accuracy on three models
	#and report the feature importance of the six
	#ka/ks values
	
	#df.replace({'A': {0: 100, 4: 400}})
	ka_ks_list = ["Ka","Ks","Ka_Ks","Ortho_Ka","Ortho_Ks","Ortho_Ka_Ks"]
	print("Substituting %s for %s" % (arb_value, replace_value))
	for column in ka_ks_list:
		ml_df = ml_df.replace({column : {arb_value : replace_value}})
	return(ml_df)

def subset_scale_split(ml_df = []):	
	#subset features
	feature_names = [ "GC_Per", "Ka_Ks", "Ka", "Ks", "Ortho_Ka_Ks", "Ortho_Ka", "Ortho_Ks", "Length", "Exon_Count", "Exon_Length", "Intron_Length", "AA", "AC", "AG", "AT", "CA", "CC", "CG", "CT", "GA", "GC", "GG", "GT", "TA", "TC", "TG", "TT"]
	feat_df = ml_df[feature_names]
	
	#setup test/train	
	target=np.ravel(ml_df['Membership'])
	feat_scale = preprocessing.scale(feat_df)
	data_train, data_test, target_train, target_test = train_test_split(feat_scale, target, test_size = 0.10, random_state = 10)
	
	split_dict = {"data_train" : data_train,
				  "data_test" : data_test,
				  "target_train" : target_train,
				  "target_test" : target_test}
	return split_dict
	
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


def train_test_models(tag = "", data_train = [], 
						data_test = [], target_train = [], target_test = []):
	
	print("Training models: %s" % (tag))
	feature_names = [ "GC_Per", "Ka_Ks", "Ka", "Ks", "Ortho_Ka_Ks", "Ortho_Ka", "Ortho_Ks", "Length", "Exon_Count", "Exon_Length", "Intron_Length", "AA", "AC", "AG", "AT", "CA", "CC", "CG", "CT", "GA", "GC", "GG", "GT", "TA", "TC", "TG", "TT"]
	
	svc_model = SVC(kernel='linear', C=1, probability = True)
	trained_svc = svc_model.fit(data_train, target_train)
	
	gnb_model = GaussianNB()
	trained_gnb = gnb_model.fit(data_train, target_train)
	
	rft_model = RandomForestClassifier(n_estimators=100)
	trained_rft = rft_model.fit(data_train,target_train)
	
	pred_svc = trained_svc.predict(data_test)
	pred_gnb = trained_gnb.predict(data_test)
	pred_rft = trained_rft.predict(data_test)
	
	
	print("SVC bdis arb - accuracy : %s" % ( 
		accuracy_score(target_test, pred_svc, normalize = True)))
	print("GNB bdis arb - accuracy : %s" % ( 
		accuracy_score(target_test, pred_gnb, normalize = True)))
	print("RFT bdis arb - accuracy : %s" % ( 
		accuracy_score(target_test, pred_rft, normalize = True)))
	
	#get feature importances
	rft_os_fi = trained_rft.feature_importances_
	feat_imp_os_dict = {"Features" : feature_names,
						"Importance" : rft_os_fi}
	feat_imp_os_df = pd.DataFrame (feat_imp_os_dict, columns = ['Features','Importance'])
	
	#just print out the features we want
	ka_ks_list = ["Ka","Ks","Ka_Ks","Ortho_Ka","Ortho_Ks","Ortho_Ka_Ks"]
	print(feat_imp_os_df[feat_imp_os_df['Features'].isin(ka_ks_list)])



if __name__ == "__main__":
	arb_test_list = [0,0.25,0.5,0.75,1,1.5,2,3,5,10,20,50,99]
	
	df_dict = load_ml_df()
	
	#new idea, different thread for each arb value
	#how to keep track of outputs
	#think there is flag when running python we can just output to the slurm eo file
	
	arb_value = arb_test_list[(int(args.slurm_array_id) - 1)]

	osat_arb = replace_arb_value(ml_df = df_dict["osat_df"], replace_value = arb_value)
	bdis_arb = replace_arb_value(ml_df = df_dict["bdis_df"], replace_value = arb_value)

	osat_ml_df_bal = balance_df(ml_df = osat_arb)
	bdis_ml_df_bal = balance_df(ml_df = bdis_arb)
	
	osat_split_dict = subset_scale_split(ml_df = osat_arb)
	osat_bal_split_dict = subset_scale_split(ml_df = osat_ml_df_bal)
	bdis_split_dict = subset_scale_split(ml_df = bdis_arb)
	bdis_bal_split_dict = subset_scale_split(ml_df = bdis_ml_df_bal)
	
	#train osat test osat
	train_test_models(tag = "Train Osat Test Osat " + str(arb_value), 
					  data_train = osat_split_dict["data_train"], 
					  data_test = osat_split_dict["data_test"], 
					  target_train = osat_split_dict["target_train"], 
					  target_test = osat_split_dict["target_test"])
	
	#train osat_bal test osat_bal
	train_test_models(tag = "Train Osat_bal Test Osat_bal " + str(arb_value),
					  data_train = osat_bal_split_dict["data_train"], 
					  data_test = osat_bal_split_dict["data_test"], 
					  target_train = osat_bal_split_dict["target_train"], 
					  target_test = osat_bal_split_dict["target_test"])
	
	#train osat test bdis
	train_test_models(tag = "Train Osat Test Bdis " + str(arb_value),
					  data_train = osat_split_dict["data_train"], 
					  data_test = bdis_split_dict["data_test"], 
					  target_train = osat_split_dict["target_train"], 
					  target_test = bdis_split_dict["target_test"])
	
	#train osat_bal test bdis_bal
	train_test_models(tag = "Train Osat_bal Test Bdis_bal " + str(arb_value),
					  data_train = osat_bal_split_dict["data_train"], 
					  data_test = bdis_bal_split_dict["data_test"], 
					  target_train = osat_bal_split_dict["target_train"], 
					  target_test = bdis_bal_split_dict["target_test"])
	
	####################################################################
	#train bdis test bdis
	train_test_models(tag = "Train Bdis Test Bdis " + str(arb_value), 
					  data_train = bdis_split_dict["data_train"], 
					  data_test = bdis_split_dict["data_test"], 
					  target_train = bdis_split_dict["target_train"], 
					  target_test = bdis_split_dict["target_test"])
	
	#train bdis_bal test bdis_bal
	train_test_models(tag = "Train Bdis_bal Test Bdis_bal " + str(arb_value),
					  data_train = bdis_bal_split_dict["data_train"], 
					  data_test = bdis_bal_split_dict["data_test"], 
					  target_train = bdis_bal_split_dict["target_train"], 
					  target_test = bdis_bal_split_dict["target_test"])
	
	#train bdis test bdis
	train_test_models(tag = "Train Bdis Test Osat " + str(arb_value),
					  data_train = bdis_split_dict["data_train"], 
					  data_test = osat_split_dict["data_test"], 
					  target_train = bdis_split_dict["target_train"], 
					  target_test = osat_split_dict["target_test"])
	
	#train bdis_bal test osat_bal
	train_test_models(tag = "Train Bdis_bal Test Osat_bal " + str(arb_value),
					  data_train = bdis_bal_split_dict["data_train"], 
					  data_test = osat_bal_split_dict["data_test"], 
					  target_train = bdis_bal_split_dict["target_train"], 
					  target_test = osat_bal_split_dict["target_test"])

	
	#ugh I could easily split these and finish them all today....

"""	
	for i in arb_test_list:
		osat_arb = replace_arb_value(ml_df = df_dict["osat_ml_df"])
		bdis_arb = replace_arb_value(ml_df = df_dict["bdis_ml_df"])
	
		osat_ml_df_bal = balance_df(ml_df = osat_arb)
		bdis_ml_df_bal = balance_df(ml_df = bdis_arb)
		
		osat_split_dict = subset_scale_split(ml_df = osat_arb)
		osat_bal_split_dict = subset_scale_split(ml_df = osat_ml_df_bal)
		bdis_split_dict = subset_scale_split(ml_df = bdis_arb)
		bdis_bal_split_dict = subset_scale_split(ml_df = bdis_ml_df_bal)
		
		#train osat test osat
		
		#train osat_bal test osat_bal
		
		#train osat test bdis
		
		#train osat_bal test bdis_bal
		
		#ugh I could easily split these and finish them all today....
"""






