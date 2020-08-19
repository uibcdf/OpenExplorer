from simtk.unit import Quantity
import simtk.unit as u

class PES():

    minima = []
    saddle_points = []
    potential_energy_minima = []*u.kilojoules_per_mole
    potential_energy_saddle_points = []*u.kilojoules_per_mole
    n_minima = 0
    n_saddle_points = 0

    def __init__(self):

        self.minima = []
        self.n_minima = 0

    def collect(self, explorer, similarity_criterion='least_rmsd', similarity_threshold=Quantity(0.25, u.angstroms)):

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
            self.n_minima+=1

        return inherent_structure_index

