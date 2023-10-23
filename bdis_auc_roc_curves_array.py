#!/bin/python
#Alan E. Yocca
#09-16-20
#bdis_auc_roc_curves_array.py

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

parser = argparse.ArgumentParser()
parser.add_argument("--slurm_array_id", required=True, help = "slurm_task_array_id")

args = parser.parse_args()


def auc_roc_curves(ml_df = "pd.data.frame()", features = "", tag = "", model = ""):
	print("Running tag %s" % (tag))
	
	#true labels
	target=np.ravel(ml_df['Membership'])
	
	#scale values (normalize)
	features_scale = preprocessing.scale(features)
	
	#cut up data
	n_splits = 10
	skf = StratifiedKFold(n_splits=n_splits)
	cv = skf.split(features_scale, target)
	
	model_obj = ""
	if model == "svc":
		model_obj = SVC(kernel='linear', C=1, probability = True)
	elif model == "gnb":
		model_obj = GaussianNB()
	elif model == "rft":
		model_obj = RandomForestClassifier(n_estimators=100)
	
	cv_auc_roc_curve(model = model_obj, features_scale = features_scale, 
					target = target, tag = tag + "_" + model, cv = cv)

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
		print(fold)
		probas_ = model.fit(features_scale[train], target[train]).predict_proba(features_scale[test])
		# Compute ROC curve and area the curve
		fpr, tpr, thresholds = roc_curve(target[test], probas_[:, 1])
		mean_tpr += np.interp(mean_fpr, fpr, tpr)
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
	
	print("Accuracy: %0.2f (+/- %0.2f)" % (mean(all_acc), stdev(all_acc) * 2))
	print("Precision: %0.2f (+/- %0.2f)" % (mean(all_pre), stdev(all_pre) * 2))
	print("Recall: %0.2f (+/- %0.2f)" % (mean(all_rec), stdev(all_rec) * 2))
	
	plt.plot([0, 1], [0, 1], '--', color=(0.6, 0.6, 0.6), label='Luck')
	
	mean_tpr /= fold
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

if __name__ == "__main__":
	
	#load data
	ml_df = pd.read_csv("bdis_meta_21_02.txt", delimiter = "\t")
	#drop na
	ml_df = ml_df.dropna()
	#binarize
	ml_df['Membership'] = ml_df['Membership'].map({'Core': 1, 'Dispensable': 0})

	feature_list = [ "Para_Ka_Ks", "Para_Ka", "Para_Ks", "Ortho_Ka_Ks", "Ortho_Ka", 
					 "Ortho_Ks", "Dup_Type", "GC_Per", "Length", "Exon_Count", 
					 "Exon_Length", "Intron_Length", "AA", "AC", "AG", "AT", 
					 "CA", "CC", "CG", "CT", "GA", "GC", "GG", "GT", "TA", "TC", 
					 "TG", "TT"]
	
	comb_list=[slice(0,28),slice(0,28),slice(0,28),slice(12,28),slice(12,28),slice(12,28)
				,slice(7,12),slice(7,12),slice(7,12)]
	tag_list=["all","all","all","dinuc","dinuc","dinuc","gene_feat","gene_feat","gene_feat"]
	model_list=["svc","gnb","rft","svc","gnb","rft","svc","gnb","rft","svc","gnb","rft"]
	loop = int(args.slurm_array_id) - 1
	
	auc_roc_curves(ml_df = ml_df, features = ml_df[feature_list[comb_list[loop]]], 
					tag = "bdis_ortho_dup_" + tag_list[loop], model = model_list[loop])
	#balance ml_df, then make curves also
	ml_df_balanced = balance_df(ml_df = ml_df)
	auc_roc_curves(ml_df = ml_df_balanced, features = ml_df_balanced[feature_list[comb_list[loop]]], 
					tag = "bdis_bal_ortho_dup_" + tag_list[loop], model = model_list[loop])

