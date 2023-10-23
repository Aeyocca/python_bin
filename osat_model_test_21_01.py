#!/bin/python
#01-22-2021
#osat_model_test_21_01.py
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

#not sure if all the above are necessary but hey, we got em

#read in dataframe
#osat_meta_osat_cds.txt
osat_cds_ml_df = pd.read_csv("osat_meta_osat_cds.txt", delimiter = "\t")

#drop any row with NA value
osat_cds_ml_df = osat_cds_ml_df.dropna()

#convert membership to binary (helps with testing accuracies later)
osat_cds_ml_df['Membership'] = osat_cds_ml_df['Membership'].map({'Core': 1, 'Dispensable': 0})

#load in para_cds, against bdis... yea...
#want to extract ka/ks from this, rename, cbind to above
osat_bdis_cds_ml_df = pd.read_csv("osat_meta_bdis_cds.txt", delimiter = "\t")
#df.rename(columns={'oldName1': 'newName1', 'oldName2': 'newName2'})
osat_bdis_cds_ml_df = osat_bdis_cds_ml_df.rename(columns={"Ka_Ks" : "Ortho_Ka_Ks", 
														  "Ka" : "Ortho_Ka",
														  "Ks" : "Ortho_Ks"})
#add these three columns
ortho_ka_ks_sub = osat_bdis_cds_ml_df[["Gene", "Ortho_Ka_Ks", "Ortho_Ka", "Ortho_Ks"]]
osat_ml_df = osat_cds_ml_df.merge(ortho_ka_ks_sub)


#hmm need to vary arbitrary value,
#lets make a function to run these?
#whats our criteria? Feature importance and accuracy
#can have the function spit these out
#looping a few arbitrary values


#what else? we just want accuracy and feature importance right?
#Also want to see what happens if we weight GC percentage

def test_arb_value(ml_df = [], replace_value = 0, arb_value = 99):
	#replace all arbitrary values with something different
	#test model accuracy on three models
	#and report the feature importance of the six
	#ka/ks values
	
	#df.replace({'A': {0: 100, 4: 400}})
	ka_ks_list = ["Ka","Ks","Ka_Ks","Ortho_Ka","Ortho_Ks","Ortho_Ka_Ks"]
	print("Substituting %s for %s" % (arb_value, replace_value))
	for column in ka_ks_list:
		ml_df = ml_df.replace({column : {arb_value : replace_value}})
	
	#subset features
	feature_names = [ "GC_Per", "Ka_Ks", "Ka", "Ks", "Ortho_Ka_Ks", "Ortho_Ka", "Ortho_Ks", "Length", "Exon_Count", "Exon_Length", "Intron_Length", "AA", "AC", "AG", "AT", "CA", "CC", "CG", "CT", "GA", "GC", "GG", "GT", "TA", "TC", "TG", "TT"]
	osat_feat_df = ml_df[feature_names]
	
	#setup test/train	
	osat_target=np.ravel(ml_df['Membership'])
	osat_feat_scale = preprocessing.scale(osat_feat_df)
	data_train_osat, data_test_osat, target_train_osat, target_test_osat = train_test_split(osat_feat_scale, osat_target, test_size = 0.10, random_state = 10)
	
	print("Training models")
	svc_model = SVC(kernel='linear', C=1, probability = True)
	trained_osat_svc = svc_model.fit(data_train_osat, target_train_osat)
	
	gnb_model = GaussianNB()
	trained_osat_gnb = gnb_model.fit(data_train_osat, target_train_osat)
	
	rft_model = RandomForestClassifier(n_estimators=100)
	trained_osat_rft = rft_model.fit(data_train_osat,target_train_osat)
	
	pred_osat_svc = trained_osat_svc.predict(data_test_osat)
	pred_osat_gnb = trained_osat_gnb.predict(data_test_osat)
	pred_osat_rft = trained_osat_rft.predict(data_test_osat)
	
	
	print("SVC osat arb %s - accuracy : %s" % (replace_value, 
		accuracy_score(target_test_osat, pred_osat_svc, normalize = True)))
	print("RFT osat arb %s - accuracy : %s" % (replace_value, 
		accuracy_score(target_test_osat, pred_osat_gnb, normalize = True)))
	print("RFT osat arb %s - accuracy : %s" % (replace_value, 
		accuracy_score(target_test_osat, pred_osat_rft, normalize = True)))
	
	#get feature importances
	rft_os_fi = trained_osat_rft.feature_importances_
	feat_imp_os_dict = {"Features" : feature_names,
						"Importance" : rft_os_fi}
	feat_imp_os_df = pd.DataFrame (feat_imp_os_dict, columns = ['Features','Importance'])
	
	#just print out the features we want
	print(feat_imp_os_df[feat_imp_os_df['Features'].isin(ka_ks_list)])



arb_test_list = [0,0.25,0.5,0.75,1,1.5,2,3,5,10,20,50,99]

for i in arb_test_list:
	test_arb_value(ml_df = osat_ml_df, replace_value = i, arb_value = 99)



