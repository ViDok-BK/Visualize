#%%
import os, sys, glob

module_root = os.path.dirname(os.path.realpath('__file__'))
module_root = '\\'.join(module_root.split('\\')[0:-1])
sys.path.append(module_root)

import numpy as np
import pandas as pd
import torch
from torch.utils.data import Dataset, DataLoader

from poly_rings.rings import PolyRingGraph
#%%
files = glob.glob(module_root + "\\_DATA\\*.csv")
for file in files:
    if "smiles" in file:
        try: 
            smiles = pd.concat(
                [smiles,pd.read_csv(file,index_col="id")])
        except: 
            smiles = pd.read_csv(file,index_col="id")
    else:
        try: 
            data = pd.concat(
                [data,pd.read_csv(file,index_col="No")])
        except: 
            data = pd.read_csv(file,index_col="No")

smiles = smiles.rename({"id":"No"},axis=1)
dataframe = pd.concat([data,smiles],axis=1).drop("DPO equation",axis=1)
del data, smiles

max_num_seg = 0

for smiles in dataframe.loc[:,"smiles"]:
    graph = PolyRingGraph(smiles)
    if len(graph.segments) > max_num_seg:
        max_num_seg = len(graph.segments)
#%%
def repr_vector(graph,size=max_num_seg):
    V = np.zeros((size,2))
    for i,seg in enumerate(graph.segments):
        ring_count = [len(cyc) for cyc in seg.cycles]
        V[i,0] = ring_count.count(6)
        V[i,1] = ring_count.count(5)
    return np.stack([V,V,V])

def adjacency_matrix(graph,size=max_num_seg):
    if isinstance(graph,str):
        graph = PolyRingGraph(graph)
    n = size
    A1 = np.zeros((n,n))
    A2 = np.zeros((n,n))
    A3 = np.zeros((n,n))
    _dict = {"120":A1,"0":A2,"60":A3}
    for i,segment in enumerate(graph.segments):
        for orien, segs in segment.adjacents.items():
            segs_index = list(map(
                lambda x: graph.segments.index(x),
                segs
            ))
            if not segs_index:
                continue
            segs_index = np.array(segs_index)
            A = _dict[orien]
            A[i,segs_index] = np.ones(len(segs_index))
    return np.stack([A1,A2,A3])

class PolyAromaticDataset(Dataset):
    def __init__(self,datasets):
        self.datasets = datasets

    def __len__(self):
        return len(self.datasets)

    def __getitem__(self, index):
        target = self.datasets.loc[:,"BG"].iloc[index]
        smiles = self.datasets.loc[:,"smiles"].iloc[index]
        graph = PolyRingGraph(smiles)
        A = torch.tensor(adjacency_matrix(graph))
        V = torch.tensor(repr_vector(graph))
        return A,V, target

data = PolyAromaticDataset(dataframe)
print(data[1])

dataloader = DataLoader(
    data,batch_size=2,shuffle=False)

# %%
"""
Test
"""
sample = dataframe.loc[:,"smiles"].iloc[150]
sample2 = dataframe.loc[:,"smiles"].iloc[180]

graph = PolyRingGraph(sample)
graph2 = PolyRingGraph(sample2)

H1 = torch.tensor(repr_vector(graph,5))
print("H1",H1.size()) #H1 torch.Size([3, 5, 2])
A1 = torch.tensor(adjacency_matrix(graph,5))
print("A1",A1.size()) #A1 torch.Size([3, 5, 5])

H2 = torch.tensor(repr_vector(graph2,5))
A2 = torch.tensor(adjacency_matrix(graph2,5))
# %%
H = torch.stack([H1,H2],dim = 0)
print("H",H.size()) #H torch.Size([2, 3, 5, 2])
A = torch.stack([A1,A2],dim = 0)
print("A",A.size()) #A torch.Size([2, 3, 5, 5])
# %%
result1 = torch.matmul(A1,H1)
print("single matmul",result1.size()) 
#single matmul torch.Size([3, 5, 2])
result = torch.matmul(A,H)
print("batch matmul",result.size())
#batch matmul torch.Size([2, 3, 5, 2])
# %%
data_iter = iter(dataloader)
sample_batch = next(data_iter)
# %%
A,V,y = sample_batch
# %%
torch.matmul(A,V)
# %%
#custom torch module
import torch.nn as nn
import torch.optim as optim
import torch.optim.lr_scheduler as lr_scheduler

X = np.concatenate(
    [
        np.arange(0,100).reshape(-1,1),
        np.arange(1,101).reshape(-1,1)
    ],axis=1)
W = np.array([[1.5],[0.75]])
Y = np.dot(X,W)

class Linear(nn.Module):
    def __init__(self):
        super().__init__()
        self.W = torch.nn.Parameter(
            torch.randn((2,1))
        )
    def forward(self,X):
        X = torch.matmul(X,self.W)
        return X

X = torch.tensor(X,dtype=torch.float32)
Y = torch.tensor(Y,dtype=torch.float32)
loss_f = nn.MSELoss()

model = Linear()
optimizer = optim.Adam(model.parameters())

optimizer2 = optim.Adagrad(
    model.parameters(),lr = 0.1
)

scheduler1 = lr_scheduler.ReduceLROnPlateau(optimizer2)
scheduler2 = lr_scheduler.CosineAnnealingLR(
    optimizer2,100)


for i in range(10000):
    Y_hat = model(X)
    loss = loss_f(Y_hat,Y)
    optimizer2.zero_grad()
    loss.backward()
    optimizer2.step()
    if i%100 == 0:
        print(loss)
    scheduler1.step(loss)
    scheduler2.step()

# %%
