import sys, os, glob

module_root = os.path.dirname(os.path.realpath('__file__'))
module_root = '\\'.join(module_root.split('\\')[0:-1])

import pandas as pd
import torch
from torch.utils.data import DataLoader

from data import PolyAromaticDataset, ReducedPolyAromaticDataset
from model import GraphConv

files = glob.glob(module_root+"\\_DATA\\*.csv")

data0  = []
smiles0 = []

for file in files:
    if "smiles" in file:
        smiles0.append(file)
    else:data0.append(file)

def read_concat_data(data1,smiles1,smiles_index = "id"):
    data1 = pd.concat(
        [pd.read_csv(file,index_col="No") for file in data1]
    )
    smiles1 = pd.concat([
        pd.read_csv(file,index_col=smiles_index) for file in smiles1
    ]).rename({smiles_index:"No"})

    data1 = pd.concat([data1,smiles1],axis=1)
    return data1

raw_data = read_concat_data(data0,smiles0).reset_index()

data = ReducedPolyAromaticDataset(raw_data)
dataloader = DataLoader(
    data,batch_size=2,shuffle=False)
