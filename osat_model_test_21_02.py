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


######LOAD OSAT DATA
osat_ml_df = pd.read_csv("osat_meta_21_02.txt", delimiter = "\t")
osat_ml_df = osat_ml_df.dropna()
	
osat_ml_df['Membership'] = osat_ml_df['Membership'].map({'Core': 1, 'Dispensable': 0})

#"""
feature_list = [ "GC_Per", "Ortho_Ka_Ks", "Ortho_Ka", "Ortho_Ks", "Para_Ka_Ks", "Para_Ka",
				   "Para_Ks", "Length", "Exon_Count", "Exon_Length", "Intron_Length", 
				   "AA", "AC", "AG", "AT", "CA", "CC", "CG", "CT", "GA", "GC", "GG", "GT", 
				   "TA", "TC", "TG", "TT", "Dup_Type"]
#"""
#No dup type or ortho ka/ks
"""
feature_list = [ "GC_Per", "Para_Ka_Ks", "Para_Ka", "Para_Ks", "Length", "Exon_Count", 
				 "Exon_Length", "Intron_Length", "AA", "AC", "AG", "AT", "CA", "CC", "CG", 
				 "CT", "GA", "GC", "GG", "GT", "TA", "TC", "TG", "TT"]
"""

osat_features = osat_ml_df[feature_list]

osat_target=np.ravel(osat_ml_df['Membership'])

osat_features_scale = preprocessing.scale(osat_features)

data_train_osat, data_test_osat, target_train_osat, target_test_osat = train_test_split(osat_features_scale, osat_target, 
                                                       test_size = 0.10, random_state = 10)                                                       
                                                       
######LOAD BDIS DATA
bdis_ml_df = pd.read_csv("bdis_meta_21_02.txt", delimiter = "\t")
bdis_ml_df = bdis_ml_df.dropna()
	
bdis_ml_df['Membership'] = bdis_ml_df['Membership'].map({'Core': 1, 'Dispensable': 0})

bdis_features = bdis_ml_df[feature_list]

bdis_target=np.ravel(bdis_ml_df['Membership'])

bdis_features_scale = preprocessing.scale(bdis_features)

data_train_bdis, data_test_bdis, target_train_bdis, target_test_bdis = train_test_split(bdis_features_scale, bdis_target, 
                                                       test_size = 0.10, random_state = 10)

#######Train models
osat_svc_model = SVC(kernel='linear', C=1, probability = True)
bdis_svc_model = SVC(kernel='linear', C=1, probability = True)
trained_osat_svc = osat_svc_model.fit(data_train_osat, target_train_osat)
trained_bdis_svc = bdis_svc_model.fit(data_train_bdis, target_train_bdis)

osat_gnb_model = GaussianNB()
bdis_gnb_model = GaussianNB()
trained_osat_gnb = osat_gnb_model.fit(data_train_osat, target_train_osat)
trained_bdis_gnb = bdis_gnb_model.fit(data_train_bdis, target_train_bdis)

osat_rft_model = RandomForestClassifier(n_estimators=100)
bdis_rft_model = RandomForestClassifier(n_estimators=100)
trained_osat_rft = osat_rft_model.fit(data_train_osat,target_train_osat)
trained_bdis_rft = bdis_rft_model.fit(data_train_bdis,target_train_bdis)


#####Test models
osat_pred_osat_svc = trained_osat_svc.predict(data_test_osat)
bdis_pred_bdis_svc = trained_bdis_svc.predict(data_test_bdis)
bdis_pred_osat_svc = trained_osat_svc.predict(data_test_bdis)
osat_pred_bdis_svc = trained_bdis_svc.predict(data_test_osat)

osat_pred_osat_gnb = trained_osat_gnb.predict(data_test_osat)
bdis_pred_bdis_gnb = trained_bdis_gnb.predict(data_test_bdis)
bdis_pred_osat_gnb = trained_osat_gnb.predict(data_test_bdis)
osat_pred_bdis_gnb = trained_bdis_gnb.predict(data_test_osat)

osat_pred_osat_rft = trained_osat_rft.predict(data_test_osat)
bdis_pred_bdis_rft = trained_bdis_rft.predict(data_test_bdis)
bdis_pred_osat_rft = trained_osat_rft.predict(data_test_bdis)
osat_pred_bdis_rft = trained_bdis_rft.predict(data_test_osat)


print("SVC Train osat, Test osat - accuracy : ",accuracy_score(target_test_osat, osat_pred_osat_svc, normalize = True))
print("SVC Train bdis, Test bdis - accuracy : ",accuracy_score(target_test_bdis, bdis_pred_bdis_svc, normalize = True))
print("SVC Train osat, Test bdis - accuracy : ",accuracy_score(target_test_bdis, bdis_pred_osat_svc, normalize = True))
print("SVC Train bdis, Test osat - accuracy : ",accuracy_score(target_test_osat, osat_pred_bdis_svc, normalize = True))

print("GNB Train osat, Test osat - accuracy : ",accuracy_score(target_test_osat, osat_pred_osat_gnb, normalize = True))
print("GNB Train bdis, Test bdis - accuracy : ",accuracy_score(target_test_bdis, bdis_pred_bdis_gnb, normalize = True))
print("GNB Train osat, Test bdis - accuracy : ",accuracy_score(target_test_bdis, bdis_pred_osat_gnb, normalize = True))
print("GNB Train bdis, Test osat - accuracy : ",accuracy_score(target_test_osat, osat_pred_bdis_gnb, normalize = True))

print("RFT Train osat, Test osat - accuracy : ",accuracy_score(target_test_osat, osat_pred_osat_rft, normalize = True))
print("RFT Train bdis, Test bdis - accuracy : ",accuracy_score(target_test_bdis, bdis_pred_bdis_rft, normalize = True))
print("RFT Train osat, Test bdis - accuracy : ",accuracy_score(target_test_bdis, bdis_pred_osat_rft, normalize = True))
print("RFT Train bdis, Test osat - accuracy : ",accuracy_score(target_test_osat, osat_pred_bdis_rft, normalize = True))


###Balanced training, balanced testing
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

def subset_scale_split(ml_df = [], feature_list = []):	
	
	feat_df = ml_df[feature_list]
	
	#setup test/train	
	target=np.ravel(ml_df['Membership'])
	feat_scale = preprocessing.scale(feat_df)
	data_train, data_test, target_train, target_test = train_test_split(feat_scale, target, test_size = 0.10, random_state = 10)
	
	split_dict = {"data_train" : data_train,
				  "data_test" : data_test,
				  "target_train" : target_train,
				  "target_test" : target_test}
	return split_dict
	

osat_ml_df_bal = balance_df(ml_df = osat_ml_df)
bdis_ml_df_bal = balance_df(ml_df = bdis_ml_df)
osat_bal_split_dict = subset_scale_split(ml_df = osat_ml_df_bal, feature_list = feature_list)
bdis_bal_split_dict = subset_scale_split(ml_df = bdis_ml_df_bal, feature_list = feature_list)

osat_bal_svc_model = SVC(kernel='linear', C=1, probability = True)
bdis_bal_svc_model = SVC(kernel='linear', C=1, probability = True)
trained_osat_bal_svc = osat_bal_svc_model.fit(osat_bal_split_dict["data_train"], osat_bal_split_dict["target_train"])
trained_bdis_bal_svc = bdis_bal_svc_model.fit(bdis_bal_split_dict["data_train"], bdis_bal_split_dict["target_train"])

osat_bal_gnb_model = GaussianNB()
bdis_bal_gnb_model = GaussianNB()
trained_osat_bal_gnb = osat_bal_gnb_model.fit(osat_bal_split_dict["data_train"], osat_bal_split_dict["target_train"])
trained_bdis_bal_gnb = bdis_bal_gnb_model.fit(bdis_bal_split_dict["data_train"], bdis_bal_split_dict["target_train"])

osat_bal_rft_model = RandomForestClassifier(n_estimators=100)
bdis_bal_rft_model = RandomForestClassifier(n_estimators=100)
trained_osat_bal_rft = osat_bal_rft_model.fit(osat_bal_split_dict["data_train"], osat_bal_split_dict["target_train"])
trained_bdis_bal_rft = bdis_bal_rft_model.fit(bdis_bal_split_dict["data_train"], bdis_bal_split_dict["target_train"])

#####Test models
osat_bal_pred_osat_bal_svc = trained_osat_svc.predict(osat_bal_split_dict["data_test"])
bdis_bal_pred_bdis_bal_svc = trained_bdis_svc.predict(bdis_bal_split_dict["data_test"])
bdis_bal_pred_osat_bal_svc = trained_osat_svc.predict(bdis_bal_split_dict["data_test"])
osat_bal_pred_bdis_bal_svc = trained_bdis_svc.predict(osat_bal_split_dict["data_test"])

osat_bal_pred_osat_bal_gnb = trained_osat_gnb.predict(osat_bal_split_dict["data_test"])
bdis_bal_pred_bdis_bal_gnb = trained_bdis_gnb.predict(bdis_bal_split_dict["data_test"])
bdis_bal_pred_osat_bal_gnb = trained_osat_gnb.predict(bdis_bal_split_dict["data_test"])
osat_bal_pred_bdis_bal_gnb = trained_bdis_gnb.predict(osat_bal_split_dict["data_test"])

osat_bal_pred_osat_bal_rft = trained_osat_rft.predict(osat_bal_split_dict["data_test"])
bdis_bal_pred_bdis_bal_rft = trained_bdis_rft.predict(bdis_bal_split_dict["data_test"])
bdis_bal_pred_osat_bal_rft = trained_osat_rft.predict(bdis_bal_split_dict["data_test"])
osat_bal_pred_bdis_bal_rft = trained_bdis_rft.predict(osat_bal_split_dict["data_test"])


print("SVC Train osat_bal, Test osat_bal - accuracy : ",
	accuracy_score(osat_bal_split_dict["target_test"], osat_bal_pred_osat_bal_svc, normalize = True))
#SVC Train osat_bal, Test osat_bal - accuracy :  0.6487223168654174
print("SVC Train bdis_bal, Test bdis_bal - accuracy : ",
	accuracy_score(bdis_bal_split_dict["target_test"], bdis_bal_pred_bdis_bal_svc, normalize = True))
#SVC Train bdis_bal, Test bdis_bal - accuracy :  0.6670247046186896
print("SVC Train osat_bal, Test bdis_bal - accuracy : ",
	accuracy_score(bdis_bal_split_dict["target_test"], bdis_bal_pred_osat_bal_svc, normalize = True))
#SVC Train osat_bal, Test bdis_bal - accuracy :  0.673469387755102
print("SVC Train bdis_bal, Test osat_bal - accuracy : ",
	accuracy_score(osat_bal_split_dict["target_test"], osat_bal_pred_bdis_bal_svc, normalize = True))
#SVC Train bdis_bal, Test osat_bal - accuracy :  0.606473594548552

print("GNB Train osat_bal, Test osat_bal - accuracy : ",
	accuracy_score(osat_bal_split_dict["target_test"], osat_bal_pred_osat_bal_gnb, normalize = True))
#GNB Train osat_bal, Test osat_bal - accuracy :  0.6313458262350937
print("GNB Train bdis_bal, Test bdis_bal - accuracy : ",
	accuracy_score(bdis_bal_split_dict["target_test"], bdis_bal_pred_bdis_bal_gnb, normalize = True))
#GNB Train bdis_bal, Test bdis_bal - accuracy :  0.7352309344790547
print("GNB Train osat_bal, Test bdis_bal - accuracy : ",
	accuracy_score(bdis_bal_split_dict["target_test"], bdis_bal_pred_osat_bal_gnb, normalize = True))
#GNB Train osat_bal, Test bdis_bal - accuracy :  0.7030075187969925
print("GNB Train bdis_bal, Test osat_bal - accuracy : ",
	accuracy_score(osat_bal_split_dict["target_test"], osat_bal_pred_bdis_bal_gnb, normalize = True))
#GNB Train bdis_bal, Test osat_bal - accuracy :  0.6248722316865417

print("RFT Train osat_bal, Test osat_bal - accuracy : ",
	accuracy_score(osat_bal_split_dict["target_test"], osat_bal_pred_osat_bal_rft, normalize = True))
#RFT Train osat_bal, Test osat_bal - accuracy :  0.8255536626916524
print("RFT Train bdis_bal, Test bdis_bal - accuracy : ",
	accuracy_score(bdis_bal_split_dict["target_test"], bdis_bal_pred_bdis_bal_rft, normalize = True))
#RFT Train bdis_bal, Test bdis_bal - accuracy :  0.5988184747583244
print("RFT Train osat_bal, Test bdis_bal - accuracy : ",
	accuracy_score(bdis_bal_split_dict["target_test"], bdis_bal_pred_osat_bal_rft, normalize = True))
#RFT Train osat_bal, Test bdis_bal - accuracy :  0.6535982814178303
print("RFT Train bdis_bal, Test osat_bal - accuracy : ",
	accuracy_score(osat_bal_split_dict["target_test"], osat_bal_pred_bdis_bal_rft, normalize = True))
#RFT Train bdis_bal, Test osat_bal - accuracy :  0.5117546848381601


#Train balanced, test unbalanced. No sense in reciprocal thankfully
osat_bal_pred_osat_svc = trained_osat_svc.predict(data_test_osat)
bdis_bal_pred_bdis_svc = trained_bdis_svc.predict(data_test_bdis)
bdis_bal_pred_osat_svc = trained_osat_svc.predict(data_test_bdis)
osat_bal_pred_bdis_svc = trained_bdis_svc.predict(data_test_osat)

osat_bal_pred_osat_gnb = trained_osat_gnb.predict(data_test_osat)
bdis_bal_pred_bdis_gnb = trained_bdis_gnb.predict(data_test_bdis)
bdis_bal_pred_osat_gnb = trained_osat_gnb.predict(data_test_bdis)
osat_bal_pred_bdis_gnb = trained_bdis_gnb.predict(data_test_osat)

osat_bal_pred_osat_rft = trained_osat_rft.predict(data_test_osat)
bdis_bal_pred_bdis_rft = trained_bdis_rft.predict(data_test_bdis)
bdis_bal_pred_osat_rft = trained_osat_rft.predict(data_test_bdis)
osat_bal_pred_bdis_rft = trained_bdis_rft.predict(data_test_osat)

print("SVC Train osat_bal, Test osat - accuracy : ",
	accuracy_score(target_test_osat, osat_bal_pred_osat_svc, normalize = True))
print("SVC Train bdis_bal, Test bdis - accuracy : ",
	accuracy_score(target_test_bdis, bdis_bal_pred_bdis_svc, normalize = True))
print("SVC Train osat_bal, Test bdis - accuracy : ",
	accuracy_score(target_test_bdis, bdis_bal_pred_osat_svc, normalize = True))
print("SVC Train bdis_bal, Test osat - accuracy : ",
	accuracy_score(target_test_osat, osat_bal_pred_bdis_svc, normalize = True))

print("GNB Train osat_bal, Test osat - accuracy : ",
	accuracy_score(target_test_osat, osat_bal_pred_osat_gnb, normalize = True))
print("GNB Train bdis_bal, Test bdis - accuracy : ",
	accuracy_score(target_test_bdis, bdis_bal_pred_bdis_gnb, normalize = True))
print("GNB Train osat_bal, Test bdis - accuracy : ",
	accuracy_score(target_test_bdis, bdis_bal_pred_osat_gnb, normalize = True))
print("GNB Train bdis_bal, Test osat - accuracy : ",
	accuracy_score(target_test_osat, osat_bal_pred_bdis_gnb, normalize = True))

print("RFT Train osat_bal, Test osat - accuracy : ",
	accuracy_score(target_test_osat, osat_bal_pred_osat_rft, normalize = True))
print("RFT Train bdis_bal, Test bdis - accuracy : ",
	accuracy_score(target_test_bdis, bdis_bal_pred_bdis_rft, normalize = True))
print("RFT Train osat_bal, Test bdis - accuracy : ",
	accuracy_score(target_test_bdis, bdis_bal_pred_osat_rft, normalize = True))
print("RFT Train bdis_bal, Test osat - accuracy : ",
	accuracy_score(target_test_osat, osat_bal_pred_bdis_rft, normalize = True))



