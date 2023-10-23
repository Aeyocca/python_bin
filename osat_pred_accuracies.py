#!/bin/python
#09-01-20
#Alan E. Yocca
#osat_pred_accuracies.py

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
parser.add_argument('--input', required=True, help='input meta file')
parser.add_argument('--feature_list', required=False, 
	help='list of features from meta file to train on',
	default = "GC_Per,Ka_Ks,Ka,Ks,Length,Exon_Count,Exon_Length,Intron_Length")
parser.add_argument('--output', required=True, help='output tsv file')

args = parser.parse_args()

if __name__ == "__main__":

	ml_df = pd.read_csv(args.input, delimiter = "\t")

	#drop NA
	ml_df = ml_df.dropna()
	#Code binary
	ml_df['Membership'] = ml_df['Membership'].map({'Core': 1, 'Dispensable': 0})

	#how to code feature subsets??
	#right now, I want: All, certain Ka/ks subsets, dinucleotide + all??	
	#lets feed features on CMD line? we should just run this separately with different features
	feature_list = args.feature_list.split(",")

	#subset the features
	features = ml_df[ feature_list ]

	#np.ravel to get pd df column into an array
	target=np.ravel(ml_df['Membership'])

	#scale features for comparison between features
	features_scale = preprocessing.scale(features)

	#don't need train_test_split because running cross validation
	#data_train, data_test, target_train, target_test = train_test_split(features_scale, target, 
    #                                                   test_size = 0.10, random_state = 10)

	model_dict = {"SVC" : {"Accuracy" : {"Mean" : 0, "Std" : 0},
						   "AUC-ROC" : {"Mean" : 0, "Std" : 0}},
				  "GNB" : {"Accuracy" : {"Mean" : 0, "Std" : 0},
						   "AUC-ROC" : {"Mean" : 0, "Std" : 0}},
				  "RF" : {"Accuracy" : {"Mean" : 0, "Std" : 0},
						   "AUC-ROC" : {"Mean" : 0, "Std" : 0}}}
	
	svc_model = SVC(kernel='linear', C=1)
	scores = cross_val_score(svc_model, features_scale, target, cv=10)
	model_dict["SVC"]["Accuracy"]["Mean"] = scores.mean()
	model_dict["SVC"]["Accuracy"]["Std"] = scores.std() * 2
	
	gnb_model = GaussianNB()
	scores = cross_val_score(gnb_model, features_scale, target, cv=10)
	model_dict["GNB"]["Accuracy"]["Mean"] = scores.mean()
	model_dict["GNB"]["Accuracy"]["Std"] = scores.std() * 2
	
	rf_model=RandomForestClassifier(n_estimators=100)
	scores = cross_val_score(rf_model, features_scale, target, cv=10)
	model_dict["RF"]["Accuracy"]["Mean"] = scores.mean()
	model_dict["RF"]["Accuracy"]["Std"] = scores.std() * 2
	
	with open(args.output, "w") as output:
		output.write("Model\tAccuracy_Mean\tAccuracy_Std\n")
		for model in model_dict.keys():
			output.write(model)
			output.write("\t")
			output.write(str(model_dict[model]["Accuracy"]["Mean"]))
			output.write("\t")
			output.write(str(model_dict[model]["Accuracy"]["Std"]))
			output.write("\n")


"""
#copy and paste to cmd line:

#ml_df = pd.read_csv("osat_meta.txt", delimiter = "\t", names = ["Gene", "Membership", "Ka_Ks", "Ka", "Ks", "Length", "Exon_Count", "Exon_Length", "Intron_Length"])

ml_df = pd.read_csv("osat_meta_dinuc.txt", delimiter = "\t")
#names = ["Gene", "Membership", "Ka_Ks", "Ka", "Ks", "Length", "Exon_Count", "Exon_Length", "Intron_Length", "TPM", "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", "X"]

ml_df = ml_df.dropna()
	
ml_df['Membership'] = ml_df['Membership'].map({'Core': 1, 'Dispensable': 0})

#features = ml_df[[ "GC_Per", "Ka_Ks", "Ka", "Ks", "Length", "Exon_Count", "Exon_Length", "Intron_Length", "TPM", "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", "X"]]
#features = ml_df[[ "GC_Per", "Ka_Ks", "Ka", "Ks", "Length", "Exon_Count", "Exon_Length", "Intron_Length", "TPM"]]

features = ml_df[[ "GC_Per", "Ka_Ks", "Ka", "Ks", "Length", "Exon_Count", "Exon_Length", "Intron_Length", "TPM", "AA", "AC", "AG", "AT", "CA", "CC", "CG", "CT", "GA", "GC", "GG", "GT", "TA", "TC", "TG", "TT", "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", "X"]]


target=np.ravel(ml_df['Membership'])

features_scale = preprocessing.scale(features)

#data_train, data_test, target_train, target_test = train_test_split(features_scale, target, 
#                                                       test_size = 0.30, random_state = 10)

data_train, data_test, target_train, target_test = train_test_split(features_scale, target, 
                                                       test_size = 0.10, random_state = 10)


#data_train_scale = preprocessing.scale(data_train)
#data_test_scale = preprocessing.scale(data_test)


#svc_model = SVC(probability = True)
svc_model = SVC(kernel='linear', C=1, probability = True)

scores = cross_val_score(svc_model, features_scale, target, cv=10)

print("SVC= Accuracy: %0.2f (+/- %0.2f)" % (scores.mean(), scores.std() * 2))

#pred = svc_model.fit(data_train, target_train).predict(data_test)

#print("SVC- accuracy : ",accuracy_score(target_test, pred, normalize = True))
#SVC- accuracy :  0.7011472634944518

gnb = GaussianNB()

scores = cross_val_score(gnb, features_scale, target, cv=10)
print("gnb= Accuracy: %0.2f (+/- %0.2f)" % (scores.mean(), scores.std() * 2))

#pred = gnb.fit(data_train, target_train).predict(data_test)
#print("gnb - accuracy : ",accuracy_score(target_test, pred, normalize = True))
#gnb - accuracy :  0.6507429001316531

weights = [1 if x == 1 else 1.5 for x in target_train]

clf=RandomForestClassifier(n_estimators=100)

scores = cross_val_score(clf, features_scale, target, cv=10)
print("RF Accuracy: %0.2f (+/- %0.2f)" % (scores.mean(), scores.std() * 2))

#Train the model using the training sets y_pred=clf.predict(X_test)
#clf.fit(data_train,target_train)

#pred=clf.predict(data_test)

#print("rf - accuracy : ",accuracy_score(target_test, pred, normalize = True))
#rf - accuracy :  0.6658830167387625


predict golfer's fantasy score, then come up with a list of score/$ on sites
but also would need results of contest to see what good scores are


pred = svc_model.fit(data_train_scale, target_train).predict(data_test_scale)
print("SVC- accuracy : ",accuracy_score(target_test, pred, normalize = True))
SVC- accuracy :  0.6782502810041214

roc_auc_score(target_test,clf.fit(data_train,data_test).predict(data_test))
"""