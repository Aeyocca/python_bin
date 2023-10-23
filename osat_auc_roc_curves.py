#!/bin/python
#Alan E. Yocca
#09-16-20
#osat_auc_roc_curves.py

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


def auc_roc_curves(ml_df = "pd.data.frame()", features = "", tag = ""):
	print("Running tag %s" % (tag))
	#drop NA
	ml_df = ml_df.dropna()
	
	#binarize
	ml_df['Membership'] = ml_df['Membership'].map({'Core': 1, 'Dispensable': 0})
	
	#true labels
	target=np.ravel(ml_df['Membership'])
	
	#scale values (normalize)
	features_scale = preprocessing.scale(features)
	
	#cut up data
	n_splits = 10
	skf = StratifiedKFold(n_splits=n_splits)
	cv = skf.split(features_scale, target)
	
	svc_model = SVC(kernel='linear', C=1, probability = True)
	gnb_model = GaussianNB()
	rft_model = RandomForestClassifier(n_estimators=100)
	
	cv_auc_roc_curve(model = svc_model, features_scale = features_scale, 
					target = target, tag = tag + "_svc_model", cv = cv)
	cv_auc_roc_curve(model = gnb_model, features_scale = features_scale, 
					target = target, tag = tag + "_gnb_model", cv = cv)
	cv_auc_roc_curve(model = rft_model, features_scale = features_scale, 
					target = target, tag = tag + "_rft_model", cv = cv)

def cv_auc_roc_curve(model = "", features_scale = "", target = "", tag = "", cv = "cv"):
	print("Cross validating %s" % (tag))
	mean_tpr = 0.0
	mean_fpr = np.linspace(0, 1, 100)
	all_tpr = []
	all_acc = []
	all_pre = []
	all_rec = []
	
	fold = 0
	for train, test in cv:
		fold += 1
		probas_ = model.fit(features_scale[train], target[train]).predict_proba(features_scale[test])
		# Compute ROC curve and area the curve
		fpr, tpr, thresholds = roc_curve(target[test], probas_[:, 1])
		mean_tpr += interp(mean_fpr, fpr, tpr)
		mean_tpr[0] = 0.0
		roc_auc = auc(fpr, tpr)
		plt.plot(fpr, tpr, lw=1, label='ROC fold %d (area = %0.2f)' % (fold, roc_auc))
		#Also keep track of accuracy, precision, and recall
		pred = model.fit(features_scale[train], target[train]).predict(features_scale[test])
		all_acc.append(accuracy_score(target[test], pred, normalize = True))
		all_rec.append(len([x for x in range(0,len(target[test])) if target[test][x] == 1 and pred[x] == 1])
						/ sum(target[test]))
		all_pre.append(len([x for x in range(0,len(target[test])) if target[test][x] == 1 and pred[x] == 1])
						/ sum(pred))
	"""
cv = skf.split(features, target)
fold = 0
all_acc = []
all_pre = []
all_rec = []
for train, test in cv:
	fold += 1
	print(fold)
	#probas_ = model.fit(features_scale[train], target[train]).predict_proba(features_scale[test])
	# Compute ROC curve and area the curve
	#fpr, tpr, thresholds = roc_curve(target[test], probas_[:, 1])
	#mean_tpr += interp(mean_fpr, fpr, tpr)
	#mean_tpr[0] = 0.0
	#roc_auc = auc(fpr, tpr)
	#plt.plot(fpr, tpr, lw=1, label='ROC fold %d (area = %0.2f)' % (fold, roc_auc))
	#Also keep track of accuracy, precision, and recall
	pred = model.fit(features_scale[train], target[train]).predict(features_scale[test])
	all_acc.append(accuracy_score(target[test], pred, normalize = True))
	all_rec.append(len([x for x in range(0,len(target[test])) if target[test][x] == 1 and pred[x] == 1])
					/ sum(target[test]))
	all_pre.append(len([x for x in range(0,len(target[test])) if target[test][x] == 1 and pred[x] == 1])
					/ sum(pred))

print("Accuracy: %0.2f (+/- %0.2f)" % (mean(all_acc), stdev(all_acc) * 2))
print("Precision: %0.2f (+/- %0.2f)" % (mean(all_pre), stdev(all_pre) * 2))
print("Recall: %0.2f (+/- %0.2f)" % (mean(all_rec), stdev(all_rec) * 2))

	"""


	print("Accuracy: %0.2f (+/- %0.2f)" % (mean(all_acc), stdev(all_acc) * 2))
	print("Precision: %0.2f (+/- %0.2f)" % (mean(all_pre), stdev(all_pre) * 2))
	print("Recall: %0.2f (+/- %0.2f)" % (mean(all_rec), stdev(all_rec) * 2))
	
	plt.plot([0, 1], [0, 1], '--', color=(0.6, 0.6, 0.6), label='Luck')
	
	mean_tpr /= len(list(cv))
	mean_tpr[-1] = 1.0
	mean_auc = auc(mean_fpr, mean_tpr)
	plt.plot(mean_fpr, mean_tpr, 'k--',
	         label='Mean ROC (area = %0.2f)' % mean_auc, lw=2)
	
	plt.xlim([-0.05, 1.05])
	plt.ylim([-0.05, 1.05])
	plt.xlabel('False Positive Rate')
	plt.ylabel('True Positive Rate')
	plt.title('AUC_ROC %s' % (tag))
	plt.legend(loc="lower right")
	plt.savefig('auc_roc_%s.pdf' % (tag)) 

if __name__ == "__main__":
	
	#load data
	ml_df = pd.read_csv("osat_meta_dinuc.txt", delimiter = "\t")
	
	feature_list = [ "GC_Per", "Ka_Ks", "Ka", "Ks", "Length", "Exon_Count", 
					   "Exon_Length", "Intron_Length", "TPM", "AA", "AC", "AG", 
					   "AT", "CA", "CC", "CG", "CT", "GA", "GC", "GG", "GT", "TA", 
					   "TC", "TG", "TT", "A", "C", "D", "E", "F", "G", "H", "I", "K", 
					   "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", "X"]
	
	auc_roc_curves(ml_df = ml_df, features = ml_df[feature_list], tag = "all")
	auc_roc_curves(ml_df = ml_df, features = ml_df[feature_list[9:24]], tag = "dinuc")
	auc_roc_curves(ml_df = ml_df, features = ml_df[feature_list[25:45]], tag = "aa")
	auc_roc_curves(ml_df = ml_df, features = ml_df[feature_list[0:8]], tag = "gene_feat")
		

