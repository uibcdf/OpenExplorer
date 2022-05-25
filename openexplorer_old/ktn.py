from networkx import DiGraph

class KTN():

    _topology = None
    _system = None

    network = DiGraph()

    def __init__(self, topology, system):

        self._topology = topology
        self._system = system

    def reset(self):

        self.network = DiGraph()

    def add_transition(self, origin=None, end=None):

        if origin not in self.network:
            self.network.add_node(origin)

        if end in self.network[origin]:
            self.network[origin][end]['weight']+=1
        else:
            self.network.add_edge(origin, end, weight=1)

