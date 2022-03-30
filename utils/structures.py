import numpy as np

class Atom:
    def __init__(self):
        """
        Class for atom
        """
        self.element = None
        self.atomic_number = None
        self.seq_num = None

        self.coord = None
        self.connectivity = {}
        self.bound_stereo = {}

        self.charge = 0

    def parse_symbol(self,symbol):
        self.element = symbol[:2].strip().title()

    def total_connectivity(self):
        return sum(self.connectivity.values())

    def printAtom(self):
        print(self.element, self.coord, self.connectivity, self.atomic_number, self.seq_num, self.charge)
    
    def returnAtom(self):
        return self.element, self.coord, self.connectivity, self.atomic_number, self.seq_num, self.charge
    def get_coords(self):
        return list(self.coord)

class ReceptorAtom(Atom):
    def __init__(self):
        super().__init__()
        #position or remoteness of atom 
        #alpha (A) C bond with acid and amine group 
        #connect to beta (B), gamma (G), delta (D), epsilon (E), zeta (Z), ???(H)
        #last digit distinguish atoms of the same remoteness
        self.position = None
        
        self.residue = None #put residue sequential number here

    def parse_symbol(self,symbol):
        """
        if len(symbol) == 1:
            self.element = symbol
        elif len(symbol) < 4:
            self.position = symbol[1:]
            self.element = symbol[0] 
        elif len(symbol) == 4:
            self.position = symbol[2:]
            self.element = symbol[1] 
        else:
            raise Exception(f"Method<ReceptorAtom.parse_symbol>: unpredicted case {symbol}")
        """
        self.element = symbol[1].strip().title()
        self.position = symbol[2:]

class Molecule:
    def __init__(self):
        self.atoms = []

    def append(self,atom):
        self.atoms.append(atom)

    def append_atom(self,atom,res_seq_num,residue_symbol):
        self.append(atom)    

    def __len__(self):
        return len(self.atoms)

    def __getitem__(self,idx):
        return self.atoms[idx]

    def printMolecule(self):
        print(self.atoms)

class Residue(Molecule):
    def __init__(self,amino_acid=None,seq_num=None):
        super().__init__()
        self.amino_acid = amino_acid
        self.seq_num = seq_num

class Protein:
    def __init__(self):
        self.residues = []

    def append(self,residue):
        self.residues.append(residue)

    def append_atom(self,atom,res_seq_num,residue_symbol):
        if not self.residues or res_seq_num != self[-1].seq_num:
            self.append(Residue(amino_acid = residue_symbol,seq_num=res_seq_num))
        elif residue_symbol != self[-1].amino_acid:
            raise Exception()
        self[-1].append(atom)
        atom.residue = self[-1]

    def __len__(self):
        return len(self.residues)

    def __getitem__(self,idx):
        return self.residues[idx]

    def printRes(self):
        print(self.residues)
