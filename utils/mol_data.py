import numpy as np

from utils.utilities import sqrt_dist, dist
from utils.structures import Atom, ReceptorAtom, Residue, Protein, Molecule
from utils.general_data import MAXIMUM_BOND_LENGTH

import MolGraph.__init__
from MolGraph.data_structures.molecular_graph import MolecularGraph

from utils.PDBReader import ComplexReader
from utils.SDFReader import SDFReader
from utils.fingerprint import AdjacentFP

class LigandGraph(Molecule):
    def __init__(self, mol, keep_hydrogen = False):
        super().__init__()
        self.__dict__ = mol.__dict__
        if keep_hydrogen == False:
            atoms = []
            for atom in self.atoms:
                if atom.atomic_number == 1:
                    continue
                atoms.append(atom)
            self.atoms = atoms
        self.map_adj_atoms()
        if not self.atoms[0].connectivity:
            self.bond_assume()

    def map_adj_atoms(self):
        for atom_order,atom in enumerate(self.atoms):
            atom.adjacents = {}
            for atom2 in self.atoms[:atom_order]:
                d12 = sqrt_dist(atom.coord,atom2.coord)
                if d12 < MAXIMUM_BOND_LENGTH(atom.element,atom2.element)**2:
                    atom.adjacents.update({atom2:np.sqrt(d12)})
                    atom2.adjacents.update({atom:np.sqrt(d12)})

    def bond_assume(self):
        for atom in self.atoms:
            atom.connectivity = atom.adjacents

    def list_all_atoms(self):
        return self.atoms

    def update_connectivity(self,input_graph):
        fp = AdjacentFP(self)
        fp.maximize_unique_FP()

        # for x,y in fp.id2atom.items():
        #     print(x)
        #     for m in y:
        #         m.printAtom()
        #from pdbqt

        for k,v in fp.id2atom.items():
            if len(v) >= 3:
                print("3 or more!" + str(len(v)))

        input_fp = AdjacentFP(input_graph)
        input_fp.maximize_unique_FP()
        # for x,y in input_fp.id2atom.items():
        #     print(x)
        #     for m in y:
        #         m.printAtom()
        #from sdf

        a2a = Atom2Atom(fp,input_fp)
        #1 pdb, 2 sdf
        connect = []
        result = ''
        for atom in self.atoms:
            connectivity = a2a.atom1_to_2[atom].connectivity
            
            connectivity = {a2a.atom2_to_1[k]:v for k,v in connectivity.items()}
            atom.connectivity = connectivity
            # print('Atom: {}'.format(atom.seq_num))
            # print(atom.printAtom())
            ele,co,conn,atomic,seq,charge = atom.returnAtom()
            
            for u,v in conn.items():
                if (u.seq_num,atom.seq_num, v) not in connect:
                    connect.append((atom.seq_num, u.seq_num, v))
            result += "{:>10s}{:>10s}{:>10s} {:4s}0  0  0  0  0  0  0  0  0  0  0  0\n".format("{:.4f}".format(co[0]), "{:.4f}".format(co[1]),"{:.4f}".format(co[2]), str(ele))

        for con in connect:
            if con == connect[-1]:
                result += "{:>3s}{:>3s}{:>3s}  0  0  0  0".format(str(con[0]), str(con[1]), str(con[2]))
            else:
                result += "{:>3s}{:>3s}{:>3s}  0  0  0  0\n".format(str(con[0]), str(con[1]), str(con[2]))
        return result, a2a.mapping

        
            
        
def read_smi(smi_path):
    smiles = open(smi_path,"r").readline()
    mol_graph = MolecularGraph()
    mol_graph.from_smiles(smiles)
    return mol_graph

class Atom2Atom:
    def __init__(self,fp1,fp2):
        fp1_ids = sorted(list(fp1.atom2id.values()))
        fp2_ids = sorted(list(fp2.atom2id.values()))
        assert fp1_ids == fp2_ids
        # print(fp1_ids)
        fp1_atoms = list(fp1.atom2id.keys())
        # for i in fp1_atoms:
        #     i.printAtom()
        fp1_atoms = sorted(fp1_atoms, key = lambda i: fp1.atom2id[i])

        fp2_atoms = list(fp2.atom2id.keys())
        fp2_atoms = sorted(fp2_atoms, key = lambda i: fp2.atom2id[i])

        self.atom1_to_2 = dict(zip(fp1_atoms,fp2_atoms))
        self.atom2_to_1 = dict(zip(fp2_atoms,fp1_atoms))
        self.mapping = {}
        for k,v in self.atom1_to_2.items():
            _,_,_,_,seq_1,_ = k.returnAtom()
            _,_,_,_,seq_2,_ = v.returnAtom()
            # print("{}---{}".format(seq_1,seq_2+1))
            if seq_2+1 not in self.mapping:
                self.mapping[seq_2+1] = seq_1
        
        # for w,y in self.mapping.items():
        #     print("{}---{}".format(w,y))
        # print('-*10')
        # print(self.atom2_to_1)

# reader = ComplexReader()
# reader.parseComplex("_DATA/ZINC000001697411_4.pdbqt")

# ligand, receptor = reader.returnLigandReceptor()

# lig_graph = LigandGraph(ligand)
# fp = AdjacentFP(lig_graph)
# fp.maximize_unique_FP()
