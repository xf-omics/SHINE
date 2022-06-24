# train on insertions using 1 aa only
# validate using 2-3 aa

# import random
# from collections import Counter
# from tpdm import tqdm
# import itertools
import argparse
import pathlib
import torch
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
# import seaborn as sns
import os.path
from os import path
import esm

import scipy
from sklearn.model_selection import GridSearchCV, train_test_split
from sklearn.decomposition import PCA
from sklearn.neighbors import KNeighborsClassifier, KNeighborsRegressor
from sklearn.svm import SVC, SVR
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from sklearn.naive_bayes import GaussianNB
from sklearn.linear_model import LogisticRegression, SGDRegressor, LinearRegression
from sklearn.linear_model import ElasticNet
from sklearn.ensemble import GradientBoostingRegressor

# from sklearn.preprocessing import OneHotEncoder
from sklearn.metrics import mean_squared_error, r2_score

from sklearn.metrics import roc_curve, auc
from sklearn.metrics import roc_auc_score


def create_parser():
    parser = argparse.ArgumentParser(
        description="Evaluate representations for pathogenicity prediction"  # noqa
    )

    parser.add_argument(
        "input_csv",
        type=pathlib.Path,
        help="variant file on which to extract representations",
    )
 
    parser.add_argument(
        "test_csv",
        type=pathlib.Path,
        help="variant file on which to extract representations",
    )

    parser.add_argument(
        "result",
        type=pathlib.Path,
        help="result file on predictions",
    )

    parser.add_argument(
        "msa_dir",
        type=pathlib.Path,
        help="directory for extracted representations from ESM-MSA",
    )

    parser.add_argument(
        "esm1b_dir",
        type=pathlib.Path,
        help="directory for extracted representations from ESM-1b",
    )
    return parser



parser = create_parser()
args = parser.parse_args()

# load AAIndex, O and U are not available
# aaindex = pd.read_csv("/share/terra/Users/xf2193/resource/AAIndex/AAIndex_imputed.txt",sep='\t',header=0)
# aaindex = aaindex.to_dict('list')

INPUT_PATH = args.input_csv
esm1b_PATH = args.esm1b_dir # Path to directory of embeddings for fasta
msa_PATH = args.msa_dir
EMB_LAYER = 33
MSA_LAYER = 12


def dataset(INPUT_PATH,esm1b_PATH,msa_PATH):
# Load embeddings (Xs) and target effects (ys)
    df = pd.read_csv(INPUT_PATH,sep='\t',header=0,keep_default_na=False)
    ys = []
    Xs = []
    ids = []

# the initial token is not saved, so 0-based
    for idx, line in df.iterrows():
        a = line['AA_ref']
        if len(a)>1:
            print(line)
            a = a[0]
        pos = line['AA_pos'] - 1 # 1-based to 0-based
        pos_msa = line['AA_msa_pos'] # 0-based
# representation from ESM-1b
        transcript = line['Transcript']
        fn = f'{esm1b_PATH}/{transcript}.pt'
        if not path.exists(fn):
            print(fn)
            continue
        embs = torch.load(fn)
        esm1b = embs['representations'][EMB_LAYER] # 1280 features
        if esm1b.shape[0] <= pos:
            print(line)
            ids.append(idx)
            continue
        esm1b = esm1b[pos,]
# representation from ESM-MSA
        gene = line['Gene']
        fn = f'{msa_PATH}/{gene}.pt'
        if not path.exists(fn):
            print(fn)
            continue
        embs = torch.load(fn)
        msa = embs['representations'][MSA_LAYER] # 768 features
        if msa.shape[0] <= pos_msa:
            print(line)
            ids.append(idx)
            continue
        msa = msa[pos_msa,]
# AAIndex
#        if a in aaindex:
#            b = aaindex[a]
#            b = torch.tensor(b)
#        else:
#            print(line)
#            b = torch.zeros(544) # 544 features
# combine features
        embs = torch.cat((esm1b,msa))
#        embs = torch.cat((embs,b))
        Xs.append(embs)
        ys.append(line['Label'])
        ids.append(line['ID'])

    Xs = torch.stack(Xs, dim=0).numpy()
    ids = np.array(ids).reshape(-1,1)
    ys = np.array(ys).reshape(-1,1)
    ys = np.append(ids,ys,axis=1)
#    ys = np.unique(ys,axis=0)
    return Xs, ys


Xs_train, ys_train = dataset(INPUT_PATH,esm1b_PATH,msa_PATH)

INPUT_PATH = args.test_csv
Xs_test, ys_test = dataset(INPUT_PATH,esm1b_PATH,msa_PATH)

# PCA
n = Xs_train.shape[1]
ys_train = ys_train[:,1].astype(int)

num = 60
pca = PCA(num)
Xs_train_pca = pca.fit_transform(Xs_train)
# evaluation
Xs_test_pca = pca.transform(Xs_test)

#   ML model
# reg = GradientBoostingRegressor()
reg = ElasticNet(alpha=0.5,l1_ratio=0.1)
model = reg.fit(Xs_train_pca, ys_train)
preds = model.predict(Xs_test_pca)
# score = reg.score(Xs_train_pca, ys_train)
print(reg.sparse_coef_)
# print(score)


preds = np.append(ys_test,preds.reshape(-1,1),axis=1)
preds = pd.DataFrame(preds,columns = ['ID','Label','ESM-MSA'])
preds["Label"] = pd.to_numeric(preds["Label"])
preds["ESM-MSA"] = pd.to_numeric(preds["ESM-MSA"])
preds = preds.groupby('ID').max()

fpr, tpr, _ = roc_curve(preds["Label"], preds["ESM-MSA"])
roc_auc = auc(fpr, tpr)
print(roc_auc)
print(f'{scipy.stats.spearmanr(preds["Label"], preds["ESM-MSA"])}')
print('\n', '-' * 80, '\n')

# np.savetxt(args.result,ids,fmt='%s',delimiter="\t")
preds.to_csv(args.result,sep="\t")
