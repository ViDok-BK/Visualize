class RingNode:
    def __init__(self, label):
        self.label = label
        self.connectivity = {}

    def add_connection(self,node,fused_atoms):
        assert isinstance(node,RingNode)
        self.connectivity[node] = fused_atoms
        node.connectivity[self] = fused_atoms

class ReducedGraph:
    def __init__(self):
        self.generic_node_classes = None
        self.nodes = []

    def from_mol_graph(self,mol_graph):
        pass


class ReducedRingGraph(ReducedGraph):
    
