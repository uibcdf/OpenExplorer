from simtk.unit import Quantity
import simtk.unit as u
from networkx import Graph

class PES():

    minima = []
    saddle_points = []
    potential_energy_minima = []*u.kilojoules_per_mole
    potential_energy_saddle_points = []*u.kilojoules_per_mole
    n_minima = 0
    n_saddle_points = 0
    basins_transition_network = Graph()
    disconnectivity_transition_network = Graph()

    def __init__(self):


    def collect_minimum(self, explorer, similarity_criterion='least_rmsd', similarity_threshold=Quantity(0.01, u.angstroms)):

        new_minimum = True
        inherent_structure_index = None

        if similarity_criterion=='least_rmsd':
            similarity=explorer.distance.least_rmsd
        else:
            raise NotImplementedError

        for minimum_index in reversed(range(self.n_minima)):
            if similarity(self.minima[minimum_index]) < similarity_threshold:
                new_minimum=False
                inherent_structure_index=minimum_index
                break

        if new_minimum:
            self.minima.append(explorer.get_coordinates())
            self.potential_energy_minima.append(explorer.get_potential_energy())
            inherent_structure_index=self.n_minima
            self.basins_transition_network.add_node(inherent_structure_index)
            self.n_minima+=1

        return inherent_structure_index

    def add_transition_between_basins(self, from=None, to=None):

        if to in self.basins_transition_network[from]:
            self.basins_transition_network[from][to]['weight']+=1
        else:
            self.basins_transition_network.add_edge(from, to, weight=1)

