#!/bin/python
#03-28-20
#Alan E. Yocca
#test_ml_core_bol.py

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

parser = argparse.ArgumentParser(prog='PROG')
parser.add_argument('--input', required=True, help='input meta file')
#parser.add_argument('--output', required=True, help='output tsv file')

args = parser.parse_args()



if __name__ == "__main__":

    ml_df = pd.read_csv(args.input, delimiter = "\t", names = ["Gene", "Membership", "Ka_Ks", "Ka", "Ks", "Length", "Exon_Count", "Exon_Length", "Intron_Length", "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", "GO_terms"])
	
    aa_array = ["A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y"]
	
    le = preprocessing.LabelEncoder()
    
    ml_df = ml_df.fillna("0")
    ml_df["GO_terms"] = le.fit_transform(ml_df["GO_terms"])
    
    #split go terms into an array in the column? hmm could test on first column
    
    #ml_df["AED"] = ml_df["AED"]
    
    
	#add header

	#now split into train / testing sets?
    features = ml_df[[ "Ka_Ks", "Ka", "Ks", ]]
	#features = ml_df[["Ka_Ks"]]
    
    #features = ml_df[aa_array]
    target = ml_df[["Membership"]]
	
    data_train, data_test, target_train, target_test = train_test_split(features, target, 
														test_size = 0.10, random_state = 10)
	
	#print(np.where(np.isnan(data_test)))
	#print(data_test.iloc[119])
    gnb = GaussianNB()
    model = gnb.fit(data_train, target_train)
    pred = gnb.fit(data_train, target_train).predict(data_test)

    print("Naive-Bayes accuracy : ",accuracy_score(target_test, pred, normalize = True))

"""
	r = permutation_importance(pred, data_test, target_test,
                            n_repeats=30,
                            random_state=0)
                            
    for i in r.importances_mean.argsort()[::-1]:
      if r.importances_mean[i] - 2 * r.importances_std[i] > 0:
        print(f"{list(features)[i]:<8}"
          f"{r.importances_mean[i]:.3f}"
          f" +/- {r.importances_std[i]:.3f}")                            
"""