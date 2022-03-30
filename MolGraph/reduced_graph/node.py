import os, sys
from copy import deepcopy

from personal_custom_pkgs.adv_import import adv_import
adv_import(explicit_path = "C:\\Users\\hoang\\Dropbox\\coding\\Cheminformatic\\MolGraph")

from data_structures.node import AtomNode
from data_structures.essential_data import VALENCE, get_valence

class MarkedAtomNode(AtomNode):
    def __init__(self,atom_node):
        self.__dict__ = deepcopy(atom_node.__dict__)
        self.label = None

def HydrogenDonor(atom_node):
    if atom_node.atomic_number == 6:  
        return False
    flag = False
    if sum(atom_node.connectivity.values()) == get_valence(atom_node):
        return False
    return True

def HydrogenAcceptor(atom_node):
    if atom_node.atomic_number in [6,16,17,35,53]:
        return False
    return True

def AcceptorDonorCheck(atom_node):
    result = []
    if HydrogenAcceptor(atom_node):
        result += "A"
    if HydrogenDonor(atom_node):
        result += "D"
    if result:
        return result#" & ".join(result)
    else:
        return None

class RingNode:
    def __init__(self):
        self.label = None
        self.atoms = None
        self.connectivity = None

class LinkerNode:
    def __init__(self):
        self.label = None
        self.atoms = None
        self.connectivity = None

class FeatureNode:
    def __init__(self):
        self.label = None
        self.atoms = None
        self.connectivity = None
