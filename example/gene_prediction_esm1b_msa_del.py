# train on insertions using 1 aa only
# validate using 2-3 aa

# import random
# from collections import Counter
# from tpdm import tqdm
# import itertools
import argparse
import pathlib
import torch
# import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
# import seaborn as sns
import os.path
from os import path
# import esm
import pickle as pk
import scipy
# from sklearn.model_selection import GridSearchCV, train_test_split
from sklearn.decomposition import PCA

from sklearn.linear_model import ElasticNet


def create_parser():
    parser = argparse.ArgumentParser(
        description="Evaluate representations for pathogenicity prediction"  # noqa
    )

    parser.add_argument(
        "test_csv",
        type=pathlib.Path,
        help="gene file on which to extract representations",
    )

    parser.add_argument(
        "result",
        type=pathlib.Path,
        help="directory for prediction results",
    )

    parser.add_argument(
        "seq_dir",
        type=pathlib.Path,
        help="directory for segmented protein sequences (inputs of ESM-MSA)",
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

seq_PATH = args.seq_dir
esm1b_PATH = args.esm1b_dir # Path to directory of embeddings for fasta
msa_PATH = args.msa_dir
EMB_LAYER = 33
MSA_LAYER = 12

def findnot(s,ch):
    st=""
    l=[]
    for i,ltr in enumerate(s):
        if ltr!=ch:
            l.append(i-1)
            st = st + ltr
    return l,st

def dataset(gene,transcript,esm1b_PATH,msa_PATH,seq_PATH):
# Load embeddings (embs)
# representation from ESM-1b
    fn = f'{esm1b_PATH}/{transcript}.pt'
    if path.exists(fn):
        embs = torch.load(fn)
        esm1b = embs['representations'][EMB_LAYER] # 1280 features
    else:
        i = 0
        while i<100:
            fn = f'{esm1b_PATH}/{transcript}_{i}.pt'
            if not path.exists(fn):
                break
            embs = torch.load(fn)
            embs = embs['representations'][EMB_LAYER] # 1280 features
            if i==0:
                esm1b = embs
            else:
                esm1b = torch.cat((esm1b,embs))
            i = i + 1
# representation from ESM-MSA
    i = 0
    while i<100:
        fn = f'{seq_PATH}/{gene}_{i}.a3m'
        if not path.exists(fn):
            break
        file1 = open(fn,'r')
        seq = file1.readline() # header
        seq = file1.readline() # sequence
        seq = seq.rstrip() # remove \n
        file1.close()
        a,seq = findnot(seq,'-')
        fn = f'{msa_PATH}/{gene}_{i}.pt'
        if not path.exists(fn):
            print(fn)
            break
        embs = torch.load(fn)
        embs = embs['representations'][MSA_LAYER] # 768 features
        embs = embs[a,:]
        if i==0:
            msa = embs
            protein = seq
        else:
            msa = torch.cat((msa,embs),0)
            protein = protein + seq
        i = i + 1
# combine features
    if 'msa' in locals():
        if esm1b.shape[0]==msa.shape[0]:
            embs = torch.cat((esm1b,msa),1)
            flag = 0
        else:
            flag = 1
    else:
        protein=""
        flag = 1
    return embs,protein,flag


INPUT_PATH = args.test_csv
df = pd.read_csv(INPUT_PATH,sep='\t',header=0,keep_default_na=False)
# PCA preload
pca = pk.load(open("model/SHINE_pca_del.pkl",'rb'))
# model preload
model = pk.load(open("model/SHINE_model_del.pkl",'rb'))

for idx, line in df.iterrows():
    gene = line['Gene.stable.ID']
    transcript = line['Transcript.stable.ID']
    protein = line['Protein.stable.ID']
    Xs_test,seq,flag = dataset(gene,transcript,esm1b_PATH,msa_PATH,seq_PATH)
    if flag==1:
        print(line)
        continue
    # PCA transformation
    with torch.no_grad():
        Xs_test_pca = pca.transform(Xs_test)
    # ML model
    preds = model.predict(Xs_test_pca)
    start = np.arange(1,len(seq)+1)
    seq = np.array([*seq]) # string to list of characters to numpy.array
    seq = np.append(start.reshape(-1,1),seq.reshape(-1,1),axis=1)
    preds = np.append(seq,preds.reshape(-1,1),axis=1)
    preds = pd.DataFrame(preds,columns = ['Position','Amino.Acid','SHINE.score'])
    fn = f'{args.result}/{protein}.txt'
    # np.savetxt(args.result,ids,fmt='%s',delimiter="\t")
    preds.to_csv(fn,sep="\t",index=False)

