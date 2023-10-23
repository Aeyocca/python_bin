#!/bin/python
#07-21-20
#Alan E. Yocca
#test_ml_core_osat.py

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

parser = argparse.ArgumentParser(prog='PROG')
parser.add_argument('--input', required=True, help='input meta file')
#parser.add_argument('--output', required=True, help='output tsv file')

args = parser.parse_args()



if __name__ == "__main__":

    ml_df = pd.read_csv("osat_meta.txt", delimiter = "\t", names = ["Gene", "Membership", "Length", "Exon_Count", "Exon_Length", "Intron_Length"])
	
    #aa_array = ["A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y", "X"]
	    
    ml_df['Membership'] = ml_df['Membership'].map({'Core': 1, 'Dispensable': 0})
    
	#now split into train / testing sets?
	features = ml_df[[ "Ka_Ks", "Ka", "Ks", "Length", "Exon_Count", "Exon_Length", "Intron_Length", "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", "X", "GO_terms" ]]
	target=np.ravel(ml_df['Membership'])
	
	features = ml_df[[ "Length", "Exon_Count", "Exon_Length", "Intron_Length" ]]
	target=np.ravel(ml_df['Membership'])

    data_train, data_test, target_train, target_test = train_test_split(features, target, 
                                                       test_size = 0.10, random_state = 10)

    svc_model = SVC()
    pred = svc_model.fit(data_train, target_train).predict(data_test)

    print("SVC- accuracy : ",accuracy_score(target_test, pred, normalize = True))

    gnb = GaussianNB()
    pred = gnb.fit(data_train, target_train).predict(data_test)

    print("gnb - accuracy : ",accuracy_score(target_test, pred, normalize = True))

