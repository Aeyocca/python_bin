#!/bin/python
#Alan E. Yocca
#09-16-20
#osat_auc_roc_curve_bal_test_bdis_bal.py

import pandas as pd
import numpy as np
from scipy import interp
from sklearn.model_selection import train_test_split, StratifiedKFold
from sklearn import preprocessing
from sklearn.naive_bayes import GaussianNB
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_curve, auc
import matplotlib.pyplot as plt
from sklearn.metrics import accuracy_score
from statistics import mean
from statistics import stdev
import argparse


def auc_roc_curves(ml_df = "pd.data.frame()", features = "", tag = "", model = ""):
	print("Running tag %s" % (tag))
	#drop NA
	ml_df = ml_df.dropna()
	
	#binarize
	ml_df['Membership'] = ml_df['Membership'].map({'Core': 1, 'Dispensable': 0})
	
	#true labels
	target=np.ravel(ml_df['Membership'])
	
	#scale values (normalize)
	features_scale = preprocessing.scale(features)

def split_and_process(ml_df = "", majority_var = 1):
	#train balanced? test balanced??? hmmmmm
	#train balanced, test balanced for now
	
	ml_df = ml_df.dropna()
	ml_df['Membership'] = ml_df['Membership'].map({'Core': 1, 'Dispensable': 0})
	
	min_sub = ml_df[ml_df.Membership == abs(majority_var - 1)]
	n_min = len(min_sub)
	maj_subset = ml_df[ml_df.Membership == majority_var].sample(n=n_min, random_state=1)
	ml_df_balanced = maj_subset.append(min_sub)
	
	
	feature_list = [ "GC_Per", "Ka_Ks", "Ka", "Ks", "Length", "Exon_Count", 
					   "Exon_Length", "Intron_Length", "AA", "AC", "AG", 
					   "AT", "CA", "CC", "CG", "CT", "GA", "GC", "GG", "GT", "TA", 
					   "TC", "TG", "TT", "A", "C", "D", "E", "F", "G", "H", "I", "K", 
					   "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]
	features = ml_df_balanced[ feature_list ]
	target=np.ravel(ml_df_balanced['Membership'])
	features_scale = preprocessing.scale(features)
	
	return features_scale, target


def plot_auc_roc(target = "", probas_ = "", pred = "", tag = ""):
	tpr = 0.0
	fpr = np.linspace(0, 1, 100)
	all_tpr = []
	all_acc = []
	all_pre = []
	all_rec = []
	
	# Compute ROC curve and area the curve
	fpr, tpr, thresholds = roc_curve(target, probas_[:, 1])
	roc_auc = auc(fpr, tpr)
	plt.plot(fpr, tpr, lw=1, label='Area = %0.2f)' % (roc_auc))
	plt.plot([0, 1], [0, 1], '--', color=(0.6, 0.6, 0.6), label='Luck')
	plt.xlim([-0.05, 1.05])
	plt.ylim([-0.05, 1.05])
	plt.xlabel('False Positive Rate')
	plt.ylabel('True Positive Rate')
	plt.title('AUC_ROC %s' % (tag))
	plt.legend(loc="lower right")
	plt.savefig('auc_roc_%s.pdf' % (tag)) 
	
	#Also keep track of accuracy, precision, and recall
	pred = model.fit(features_scale[train], target[train]).predict(features_scale[test])
	print("Accuracy for %s: %0.2f" % (tag, accuracy_score(target[test], pred, normalize = True)))
	print("Recall for %s: %0.2f" % (tag, len([x for x in range(0,len(target[test])) 
										if target[test][x] == 1 and pred[x] == 1])/ sum(target[test])))
	print("Precision for %s: %0.2f" % (tag, len([x for x in range(0,len(target[test])) 
										if target[test][x] == 1 and pred[x] == 1])/ sum(pred)))


if __name__ == "__main__":
	#load data
	wkdir = "/mnt/research/edgerpat_lab/AlanY/00_collab/02_core_pred/"
	ml_df_osat = pd.read_csv(wkdir + "/06_rice/osat_meta_dinuc.txt", delimiter = "\t")
	ml_df_bdis = pd.read_csv(wkdir + "/04_brachy/bdis_meta_ortho.txt", delimiter = "\t")
	
	#split into features/classes, preprocess features
	osat_feat_norm, osat_class = split_and_process(ml_df = ml_df_osat, majority_var = 1)
	bdis_feat_norm, bdis_class = split_and_process(ml_df = ml_df_bdis, majority_var = 0)
	
	#train three models with all you got
	svc = SVC(kernel='linear', C=1, probability = True)
	gnb = GaussianNB()
	rft = RandomForestClassifier(n_estimators=100)

	svc_osat = svc.fit(osat_feat_norm, osat_class)
	gnb_osat = gnb.fit(osat_feat_norm, osat_class)
	rft_osat = rft.fit(osat_feat_norm, osat_class)
	
	#need acc, pre, rec, tpr, fpr
	print("Predicting probs")
	svc_osat_probas_bdis = svc_osat.predict_proba(bdis_feat_norm)
	svc_osat_pred_bdis = svc_osat.predict(bdis_feat_norm)
	gnb_osat_probas_bdis = gnb_osat.predict_proba(bdis_feat_norm)
	gnb_osat_pred_bdis = gnb_osat.predict(bdis_feat_norm)
	rft_osat_probas_bdis = rft_osat.predict_proba(bdis_feat_norm)
	rft_osat_pred_bdis = rft_osat.predict(bdis_feat_norm)
	
	print("Plotting")
	plot_auc_roc(target = bdis_class, probas_ = svc_osat_probas_bdis, 
					pred = svc_osat_pred_bdis, tag = "SVC_osat_train_bdis_test")
	plot_auc_roc(target = bdis_class, probas_ = gnb_osat_probas_bdis, 
					pred = gnb_osat_pred_bdis, tag = "GNB_osat_train_bdis_test")
	plot_auc_roc(target = bdis_class, probas_ = rft_osat_probas_bdis, 
					pred = rft_osat_pred_bdis, tag = "RFT_osat_train_bdis_test")
