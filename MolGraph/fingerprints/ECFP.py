#%%
import os, sys
from zlib import crc32

from numpy import double

module_root = os.path.dirname(os.path.realpath('__file__'))
module_root = '\\'.join(module_root.split('\\')[0:-1])
sys.path.append(module_root)

from data_structures.essential_data import *
from data_structures.node import AtomNode
from data_structures.tree import Tree
from data_structures.molecular_graph import MolecularGraph
# %%
def sum_bonds(connectivity):
    numbers_bond = 0
    double_flag = False
    for bond_order in connectivity.values():
        if bond_order in range(0,10):
            numbers_bond += bond_order
        elif bond_order == 10:
            if double_flag:
                numbers_bond += 2
                double_flag = False
            else:
                numbers_bond += 1
                double_flag = True
        else:
            raise Exception("unknown bond type",bond_order,"in",connectivity)
    return numbers_bond

class Atom(AtomNode):
    def __init__(self, value):
        super().__init__(value)
        self.identifier = None
        self.all_valence = sorted(
            VALENCE[self.value.title()])
        
    def get_valence(self,total_bonds):
        for val in self.all_valence:
            if val >= total_bonds:
                break
        return val

    def hash(self,array):
        array = str(array)
        return crc32(bytes(array,"utf8"))

    def sort_key(self,atom):
        return self.connectivity[atom]*10**11 + atom.identifier

    def init_identifier(self):
        #get total number of attached H
        all_bonds = sum_bonds(self.connectivity)
        implicit_H = self.get_valence(all_bonds) - all_bonds
        total_H = max(implicit_H,self.explicit_H)
        
        #immediate adjacent non-hydrogen atoms
        self.adj_atoms = [atom for atom in 
            self.connectivity.keys() if atom.value != "H"]
        
        p1 = len(self.adj_atoms)
        #p2 = the valence minus the number of hydrogen
        p3 = self.atomic_number
        p4 = self.atomic_mass
        p5 = self.charge
        #number of attached hydrogen
        p6 = total_H
        assert p5 >= 0
        p7 = int(self.in_cycles)
        self.array = [p1,p3,p4,p5,p6,p7]
        self.identifier = self.hash(self.array)
    
    def new_identifier(self):
        self.adj_atoms = sorted(
            self.adj_atoms,key = self.sort_key)
        array = list(map(
            lambda x:x.identifier,self.adj_atoms
        ))
        return self.hash(array)

    def update_identifier(self,identifier):
        self.identifier = identifier

class ECFPGraph(MolecularGraph):
    def __init__(self, smiles, node_class = Atom):
        super().__init__(node_class=node_class)
        self.from_smiles(smiles)
        self.find_cycles()
        self.traverse(f = lambda x:x.init_identifier())
        self.identifiers = []

    def update_identifier(self):
        new_identifiers = {}
        self.traverse(f = lambda x: new_identifiers.update({x:x.new_identifier()}))
        self.traverse(f = lambda x: x.update_identifier(new_identifiers[x]))
        return new_identifiers 

    def returnECFP(self,radius):
        self.update_identifier(n = radius)
        self.traverse(f = lambda x:print(x.identifier))
# %%

sample = "c1ccc2c(cccc2O)c1"
ecfp_graph = ECFPGraph(smiles= sample,node_class=Atom)
print("initial identier")
ecfp_graph.traverse(f = lambda x: print(x,x.identifier))

ecfp = ecfp_graph.update_identifier()
print("Updated identifer")
for atom,identifier in ecfp.items():
    print(atom,identifier)

ecfp = ecfp_graph.update_identifier()
print("Updated identifer")
for atom,identifier in ecfp.items():
    print(atom,identifier)

ecfp = ecfp_graph.update_identifier()
print("Updated identifer")
for atom,identifier in ecfp.items():
    print(atom,identifier)

