from node import MarkedAtomNode

from personal_custom_pkgs.adv_import import adv_import
adv_import(explicit_path = "C:\\Users\\hoang\\Dropbox\\coding\\Cheminformatic\\MolGraph")

from data_structures.tree import DFS
from data_structures.cycles import CycleFinding, Cycle

from reduced_graph.node import MarkedAtomNode, AcceptorDonorCheck

class ReducedGraph:
    """
    Class for classic reduced graph.
    """
    def __init__(self,mol_graph):
        atom_nodes = mol_graph.list_all_atoms()
        self.atom_nodes = [MarkedAtomNode(node) for node in atom_nodes]
        for i,node in enumerate(self.atom_nodes):
            self.atom_nodes[i].label = AcceptorDonorCheck(node)

    def preliminary(self):
        """
        Delete terminal node
        """
        pass
    
    def mark_ring(self):
        cycle_finder = CycleFinding()
        cycles = cycle_finder.find_cycle(self.nodes_list[0])

    def mark_linker(self):
        pass

    def mark_feature(self):
        pass

    def collapse(self):
        pass

