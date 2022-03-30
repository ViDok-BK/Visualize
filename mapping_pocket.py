import glob
import os

current_dir = os.getcwd()

def get_mapping():
    pocket = glob.glob(current_dir + '/data/Data_Vidok_Infer/*_pocket.pdb')[0]
    lines = open(pocket,'r').readlines()
    mapping_index = {}
    index = 0
    for line in lines:
        t=[i for i in line.split(' ') if i!='']
        if t[0]=='ATOM':
            mapping_index[index] = int(t[1])
            index+=1

    return mapping_index

    
