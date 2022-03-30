from copy import deepcopy
from zlib import crc32
import numpy as np

from utils.PDBReader import ComplexReader

class BaseFP:
    def __init__(self,mol):
        self.atoms = mol.list_all_atoms()
        self.reset_FP()
        self.fp = None

    def reset_FP(self):
        self.atom2id = {atom:self.init_identifier(atom) for atom in self.atoms}

class AdjacentFP(BaseFP):
    """
    Refer to 
        + atom2id (dict) to retrieve identifier from atom
        + id2atom (dict) to retrieve list of atom(s) from identifier
    """
    def __init__(self,mol):
        super().__init__(mol)

    def hash_f(self,id_list):
        """
        hash_f: list -> bytes_string
        """
        if not isinstance(id_list,list):
            assert(isinstance(id_list,int))
            return crc32(bytes([id_list]))
        assert(isinstance(id_list,list))
        id_list = sorted(id_list)
        id_list = np.array(id_list)
        return crc32(id_list.tobytes())

    def init_identifier(self,atom):
        """create the initial identifier for each node"""
        return self.hash_f(atom.atomic_number)

    def update_identifiers(self,atom):
        """
        update identifier of each node in one loop
        """
        id_list = [self.atom2id[atom]]
        id_list += [self.atom2id[a] for a in atom.adjacents.keys()]
        return self.hash_f(id_list)

    def update_FP(self):
        """
        update all identifiers of all node in one loop, effectively, execute a loop
        """
        fp = {} 
        for atom in self.atoms:
            identifier = self.update_identifiers(atom)
            fp.update({atom:identifier})
        self.atom2id = fp

    def len_unique_FP(self):
        """Compute number of unique identifiers"""
        return len(set(self.atom2id.values()))

    def compute_FP(self,cycles):
        """
        Compute the identifiers of all atoms in the molecule after a number of cycles
        """
        for i in range(cycles):
            self.update_FP()

    def maximize_unique_FP(self):
        """
        Iteratively update the identifiers of all atoms until the number 
        of unique identifiers does not increase.
        This method come with the final creation of id2atom 
        """
        self.reset_FP()
        no_unique_FP = self.len_unique_FP()
        while True:
            self.update_FP()
            if self.len_unique_FP() == no_unique_FP:
                break
            no_unique_FP = self.len_unique_FP()
        self.create_id2atom()

    def create_id2atom(self):
        """
        Create the id2atom dict
        """
        self.id2atom = {f:[] for f in self.atom2id.values()}
        for k,v in self.atom2id.items():
            self.id2atom[v].append(k)
