import numpy as np
import pandas as pd
import torch
import torch.nn as nn

class GraphConv(nn.Module):
    """
    Layer for graph convolution.
    Input: 
        + X \in R^(N,input_size) is feature vector of sample
            where N is the number of nodes 
        + A \in R^(N,N) is adjacency matrix
    Output:
        + Y \in R^(N,output_size) 
    Formula:
        + Y = A.X.W
            where W \in R^(input_size,output_size)
    """
    def __init__(self,input_size,output_size):
        super().__init__()
        self.W = nn.Parameter(
            torch.randn(input_size,output_size,dtype= torch.float64))
    def forward(self,A,X):
        X = torch.matmul(A,X)
        X = torch.matmul(X,self.W)
        return X

"""
class BaseLineModel(nn.Module):
    def __init__(self):
        super().__init__()
        self.parmaeter = torch.rand
"""
