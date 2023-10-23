#!/bin/python
#Alan E. Yocca
#bbpan_snp_pca.py

import pandas as pd
import numpy as np
#import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
import argparse

parser = argparse.ArgumentParser(prog='PROG')
parser.add_argument('--input', required=True, help='input vcf file')
parser.add_argument('--output', required=True, help='output txt file')

args = parser.parse_args()

snp_map_maf = pd.read_csv(args.input, sep='\t')

#need to transpose this?
snp_mat_t = snp_map_maf.T

print("Transposed matrix")

x = snp_mat_t.iloc[1:,]

pca = PCA(n_components=28)
principalComponents = pca.fit_transform(x)

print("Finished calculating PCA")

principalDf = pd.DataFrame(data = principalComponents)

spec_np = np.array(list(snp_mat_t.index)[1:])

principalDf["Accession"] = spec_np

principalDf.to_csv(args.output, sep = "\t")



