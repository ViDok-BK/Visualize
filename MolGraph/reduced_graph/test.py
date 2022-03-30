from personal_custom_pkgs.adv_import import adv_import
adv_import(explicit_path = "C:\\Users\\hoang\\Dropbox\\coding\\Cheminformatic\\MolGraph")

from data_structures.molecular_graph import MolecularGraph
from reduced_graph.node import AcceptorDonorCheck
from reduced_graph.graph import ReducedGraph

sample = "CC1(C)C(C=C2C(CC3=CC=C(C=C3)C(OC)=O)CO)=C(C=C2O)C(C)(C)CC1" 

molgraph = MolecularGraph()
molgraph.from_smiles(sample)

#reduced_graph = ReducedGraph(molgraph)
molgraph.traverse(lambda x: print(AcceptorDonorCheck(x)))

regraph = ReducedGraph(molgraph)
