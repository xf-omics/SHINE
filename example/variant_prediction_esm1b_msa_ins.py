# train on insertions using 1 aa only

import argparse
import pathlib
import torch
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
# import seaborn as sns
import os.path
from os import path
# import esm
import pickle as pk
import scipy

from sklearn.decomposition import PCA
from sklearn.linear_model import ElasticNet
from sklearn.metrics import mean_squared_error, r2_score

from sklearn.metrics import roc_curve, auc
from sklearn.metrics import roc_auc_score


def create_parser():
    parser = argparse.ArgumentParser(
        description="Evaluate representations for pathogenicity prediction"  # noqa
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
        help="output directory for extracted representations",
    )

    parser.add_argument(
        "esm1b_dir",
        type=pathlib.Path,
        help="output directory for extracted representations",
    )
    return parser



parser = create_parser()
args = parser.parse_args()


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
        a = line['AA_alt']
        c = len(a)
        pos = line['AA_pos'] - 1 # 1-based to 0-based
        pos_msa = line['AA_msa_pos'] # 0-based
# representation from ESM-1b
        varid = line['VarID']
        fn = f'{esm1b_PATH}/{varid}.pt'
        if not path.exists(fn):
            print(fn)
            continue
        embs = torch.load(fn)
        esm1b = embs['representations'][EMB_LAYER] # 1280 features
        if esm1b.shape[0] <= pos:
            print(line)
            ids.append(idx)
            continue
# representation from ESM-MSA
        gene = line['Gene']
        fn = f'{msa_PATH}/{gene}.pt'
        if not path.exists(fn):
            print(fn)
            continue
        embs = torch.load(fn)
        msa = embs['representations'][MSA_LAYER]
        if msa.shape[0] <= pos_msa:
            print(line)
            ids.append(idx)
            continue
        msa = msa[pos_msa,]
# AAIndex
        for i in range(c):
            esm1b_pos = esm1b[pos,]
            embs = torch.cat((esm1b_pos,msa))
            pos = pos+1
            Xs.append(embs)
            ys.append(line['Label'])
            ids.append(line['ID'])

    Xs = torch.stack(Xs, dim=0).numpy()
    ids = np.array(ids).reshape(-1,1)
    ys = np.array(ys).reshape(-1,1)
    ys = np.append(ids,ys,axis=1)
    return Xs, ys


INPUT_PATH = args.test_csv
Xs_test, ys_test = dataset(INPUT_PATH,esm1b_PATH,msa_PATH)

# PCA
pca = pk.load(open("model/SHINE_pca_ins.pkl",'rb'))
Xs_test_pca = pca.transform(Xs_test)

# ML model
model = pk.load(open("model/SHINE_model_ins.pkl",'rb'))
preds = model.predict(Xs_test_pca)

# maximum
preds = np.append(ys_test,preds.reshape(-1,1),axis=1)
preds = pd.DataFrame(preds,columns = ['ID','Label','ESM-MSA'])
preds["Label"] = pd.to_numeric(preds["Label"])
preds["ESM-MSA"] = pd.to_numeric(preds["ESM-MSA"])
preds = preds.groupby('ID').max()

# evaluation
fpr, tpr, _ = roc_curve(preds["Label"], preds["ESM-MSA"])
roc_auc = auc(fpr, tpr)
print(roc_auc)
print(f'{scipy.stats.spearmanr(preds["Label"], preds["ESM-MSA"])}')
print('\n', '-' * 80, '\n')

# np.savetxt(args.result,ids,fmt='%s',delimiter="\t")
preds.to_csv(args.result,sep="\t")

