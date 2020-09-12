from simtk.unit import Quantity
import simtk.unit as u
from networkx import Graph
from numpy import infty

class PES():

    _topology = None
    _system = None

    minima = []
    saddle_points = []
    potential_energy_minima = []*u.kilojoules_per_mole
    potential_energy_saddle_points = []*u.kilojoules_per_mole
    n_minima = 0
    n_saddle_points = 0
    basins_network = Graph()
    disconnectivity_transition_network = Graph()
    global_minimum_index = None
    global_minimum_potential_energy = infty*u.kilojoules_per_mole

    def __init__(self, topology, system):

        self._topology = topology
        self._system = system

    def reset(self):

        self.minima = []
        self.saddle_points = []
        self.potential_energy_minima = []*u.kilojoules_per_mole
        self.potential_energy_saddle_points = []*u.kilojoules_per_mole
        self.n_minima = 0
        self.n_saddle_points = 0
        self.basins_network = Graph()
        self.disconnectivity_transition_network = Graph()
        self.global_minimum_index = None
        self.global_minimum_potential_energy = infty*u.kilojoules_per_mole

    def collect_minimum(self, explorer, coordinates=None, potential_energy=None, similarity_criterion='least_rmsd', similarity_threshold=Quantity(0.01, u.angstroms)):

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

            if coordinates is None:
                coordinates = explorer.get_coordinates()
            if potential_energy is None:
                potential_energy = explorer.get_potential_energy()

            self.minima.append(coordinates)
            self.potential_energy_minima.append(potential_energy)
            inherent_structure_index=self.n_minima
            self.basins_network.add_node(inherent_structure_index)
            self.n_minima+=1

            if self.global_minimum_potential_energy > potential_energy:
                self.global_minimum_index = inherent_structure_index
                self.global_minimum_potential_energy = potential_energy

        return inherent_structure_index

    def add_pair_of_neighbor_minima(self, minimum_i=None, minimum_j=None):

        self.basins_network.add_edge(minimum_i, minimum_j)

