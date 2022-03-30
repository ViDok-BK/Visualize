import numpy as np
from utils.general_data import PERIODIC_TABLE
from utils.structures import Atom, ReceptorAtom, Residue, Protein, Molecule

class SDFReader:
    def __init__(self):
        self.molecule = Molecule()
        self.number_of_bonds = 0

    def read(self,path):
        start = ""
        end = ""
        with open(path,"r") as  o:
            lines = o.readlines()
            self.n = 0
            for i, line in enumerate(lines):
                if i < 3:
                    start+=line
                    continue
                sparsed_line = line.split()
                if "END" in line or "CHG" in line or "M" in line:
                    end+=line
                    break
                elif i == 3:
                    start += line
                    n_atoms, n_bonds = tuple(sparsed_line[:2])
                    n_atoms = int(n_atoms)
                    n_bonds = int(n_bonds)
                #elif i > 3 and i < 3 + n_atoms + 1: 
                elif len(sparsed_line) > 10:
                    atom = self.parse_atom_line(sparsed_line,i - 4)
                    self.molecule.append(atom)
                #elif i > 3 + n_atoms and i < 3 + n_atoms + n_bonds + 1:
                else:
                    self.parse_bond_line(sparsed_line)
        return start, end

    def parse_atom_line(self,line,seq_num):
        """1.1378    0.0901   -0.8134 C   0  0  0  0  0  0  0  0  0  0  0  0"""
        atom = Atom()
        atom.coord = np.array([float(i) for i in line[:3]])
        atom.symbol = line[3]
        atom.element = line[3].strip().title()
        atom.atomic_number = PERIODIC_TABLE[atom.element][0]
        atom.seq_num = seq_num
        return atom

    def parse_bond_line(self,line):
        self.number_of_bonds += 1
        idx1, idx2, bond_type = tuple(line[:3])
        idx1 = int(idx1) - 1
        idx2 = int(idx2) - 1
        bond_type = int(bond_type)
        self.molecule[idx1].connectivity.update({self.molecule[idx2]:bond_type})
        self.molecule[idx2].connectivity.update({self.molecule[idx1]:bond_type})

"""
class PDBComplexReader:
    def __init__(self):
        self.receptor = None
        self.ligand = None
        
    def read(self,file_path):
        with open(file_path,"r") as o:
            for line in o.readlines():
                self.parse(line)

    def parse(self,line):
        #sample line with index
        #          1         2         3         4         5         6         7
        #0     6     2   67   12       0       8       6   0         0         0     6         
        #ATOM   1768  CA  ARG A 188     -14.327  13.087  76.166  1.00 46.10     0.159 C
        #ATOM    366 1HH1 ARG A  40     -20.806  20.348  79.130  1.00  0.00     0.174 HD
        #ATOM    367 2HH1 ARG A  40     -20.116  21.607  78.100  1.00  0.00     0.174 HD
        self.constructor = ReceptorConstructor()

        terms = line.split()
        if terms[0] == "ATOM" or terms[0] == "HETATOM":
            atom = Atom()

            atom.seq_num = terms[1]
            atom.element = terms[2]
            atom.coord = np.array([float[i] for i in terms[6:9]])
            atom.charge = float(terms[11])

            residue = terms[3]
            res_seq_num = terms[5]
            symbol2 = terms[12]

            constructor.append(atom,residue,re_seq_num,symbol2)

        elif "VINA" in line:
            self.swap_constructor()

    def swap_constructor(self):
        self.constructor = LigandConstructor()
            
class ResidueConstructor:
    def __init__(self,atom,residue,re_seq_num,symbol2):
        self.residue = Residue()
        self.residue.amino_acid = residue
        self.residue.seq_num = re_seq_num

    def append(self,atom,residue,re_seq_num,symbol2):
        if re_seq_num != self.residue.seq_num:
            return 0
        
class ReceptorConstructor:
    def __init__(self):
        self.receptor = Protein()
        
    def append(self,atom,residue,re_seq_num,symbol2):
        atom = ReceptorAtom(atom)
        

class LigandConstructor:
    def __init__(self):
        
    def append(self):
"""
