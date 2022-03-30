import os, sys, glob

import numpy as np
import torch
from torch.utils.data import Dataset, DataLoader

module_root = os.path.dirname(os.path.realpath('__file__'))
module_root = '\\'.join(module_root.split('\\')[0:-1])
sys.path.append(module_root)

from poly_rings.rings import PolyRingGraph

def findMaxNumSegs(data):
    max_num_seg = 0
    for smiles in data.loc[:,"smiles"]:
        graph = PolyRingGraph(smiles)
        if len(graph.segments) > max_num_seg:
            max_num_seg = len(graph.segments)
    return max_num_seg

def repr_vector(graph,size):
    V = np.zeros((size,2))
    for i,seg in enumerate(graph.segments):
        ring_count = [len(cyc) for cyc in seg.cycles]
        V[i,0] = ring_count.count(6)
        V[i,1] = ring_count.count(5)
    return np.stack([V,V,V])

def adjacency_matrix(graph,size):
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
                lambda x: graph.segments.index(x),segs))
            if not segs_index:
                continue
            segs_index = np.array(segs_index)
            A = _dict[orien]
            A[i,segs_index] = np.ones(len(segs_index))
    return np.stack([A1,A2,A3])

def reduce_adj_matrix(graph,size):    
    if isinstance(graph,str):
        graph = PolyRingGraph(graph)
    n = size
    A1 = np.zeros((n,n))
    A2 = np.zeros((n,n))
    A3 = np.zeros((n,n))
    _dict = {"120":A1,"0":A2,"60":A3}
    longest = max(map(lambda x: len(x),graph.segments))

    for i,segment in enumerate(graph.segments):
        if len(segment) != longest:
            continue
        for orien, segs in segment.adjacents.items():
            segs_index = list(map(
                lambda x: graph.segments.index(x),segs))
            if not segs_index:
                continue
            segs_index = np.array(segs_index)
            A = _dict[orien]
            A[i,segs_index] = np.ones(len(segs_index))
    return np.stack([A1,A2,A3])

class PolyAromaticDataset(Dataset):
    def __init__(self,datasets):
        self.datasets = datasets
        self.max_num_seg = findMaxNumSegs(datasets)

    def __len__(self):
        return len(self.datasets)

    def __getitem__(self, index):
        target = self.datasets.loc[index,"BG"]
        smiles = self.datasets.loc[index,"smiles"]
        graph = PolyRingGraph(smiles)
        A = adjacency_matrix(graph,size = self.max_num_seg)
        A = torch.tensor(A)
        V = repr_vector(graph,size = self.max_num_seg)
        V = torch.tensor(V)
        return A,V, target


class ReducedPolyAromaticDataset(Dataset):
    def __init__(self,datasets):
        self.datasets = datasets
        self.max_num_seg = findMaxNumSegs(datasets)

    def __len__(self):
        return len(self.datasets)

    def __getitem__(self, index):
        target = self.datasets.loc[index,"BG"]
        smiles = self.datasets.loc[index,"smiles"]
        graph = PolyRingGraph(smiles)
        A = reduce_adj_matrix(graph,size = self.max_num_seg)
        A = torch.tensor(A)
        V = repr_vector(graph,size = self.max_num_seg)
        V = torch.tensor(V)
        return A,V, target

