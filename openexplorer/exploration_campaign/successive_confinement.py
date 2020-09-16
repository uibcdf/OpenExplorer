import numpy as np
import simtk.unit as u
from simtk.unit import Quantity
from openexplorer.tools.quench import L_BFGS
from openexplorer.tools.md import Langevin

class SuccessiveConfinement():

    explorer = None

    quench = None
    md = None

    n_confinements_per_basin = None
    md_time_before_quench = None
    md_steps_before_quench = None
    similarity_threshold = None

    pes = None
    ktn = None
    time_per_confinement = {}
    trajectory_inherent_structures = {}
    basins_to_be_explored = set()
    basins_explored = set()

    initialized = False

    def __init__(self, explorer, n_confinements_per_basin=250, md_time_before_quench=Quantity(1.0, u.picoseconds),
                 similarity_threshold=Quantity(0.01, u.angstroms), quench=L_BFGS, md=Langevin):

        from openexplorer import PES, KTN

        self.explorer = explorer

        self.quench = quench(self.explorer)
        self.md = md(self.explorer)

        self.n_confinements_per_basin = n_confinements_per_basin

        self.md_time_before_quench = md_time_before_quench
        md_timestep = self.md.get_parameters()['timestep']
        self.md_steps_before_quench = int(self.md_time_before_quench/md_timestep)
        self.similarity_threshold = similarity_threshold

        self.pes=PES(self.explorer.topology, self.explorer.context.getSystem())
        self.ktn=KTN(self.explorer.topology, self.explorer.context.getSystem())

        self.reset()

    def reset(self):

        self.pes.reset()
        self.ktn.reset()
        self.time_per_confinement = {}
        self.trajectory_inherent_structures = {}
        basins_to_be_explored = []

        self.initialized = False

    def run(self, n_basins=None, progress_bar=False, verbose=False):

        temperature = self.md.get_parameters()['temperature']
        n_basins_explored = 0

        # With n_basins in [0, None] there should be a convergence criterium
        # With n_confinements_per_basin in [0, None] there should also be a convergence criterium

        if n_basins is None:
            n_basins = np.inf

        if not self.initialized:
            self.quench()
            basin_index = self.pes.collect_minimum(self.explorer)
            self.basins_to_be_explored.add(basin_index)
            self.initialized = True
        else:
            if not len(self.basins_to_be_explored):
                self.quench()
                basin_index = self.pes.collect_minimum(self.explorer)
                self.basins_to_be_explored.add(basin_index)

        while len(self.basins_to_be_explored) and (n_basins_explored < n_basins):

            current_basin_index = self.basins_to_be_explored.pop()
            current_basin_coordinates = self.pes.minima[current_basin_index]

            if current_basin_index not in self.time_per_confinement:
                self.time_per_confinement[current_basin_index]=[]*u.nanoseconds
                self.trajectory_inherent_structures[current_basin_index]=[]

            pre_coordinates = current_basin_coordinates

            if progress_bar:

                from tqdm import tqdm
                iterator = tqdm(range(self.n_confinements_per_basin))

            else:

                iterator = range(self.n_confinements_per_basin)

            for _ in iterator:

                n_quenchs = 0
                is_in = True

                self.explorer.set_coordinates(pre_coordinates)
                self.explorer.set_velocities_to_temperature(temperature)

                while is_in:

                    self.md(self.md_steps_before_quench)

                    coordinates = self.explorer.get_coordinates()
                    velocities = self.explorer.get_velocities()

                    self.quench()
                    n_quenchs+=1

                    if self.explorer.distance.least_rmsd(current_basin_coordinates) <= self.similarity_threshold:

                        self.explorer.set_coordinates(coordinates)
                        self.explorer.set_velocities(velocities)
                        basin_index = current_basin_index
                        pre_coordinates = coordinates

                    else:

                        is_in = False

                        basin_index = self.pes.collect_minimum(self.explorer)
                        if (basin_index not in self.basins_explored) or (basin_index not in self.basins_to_be_explored):
                            self.basins_to_be_explored.add(basin_index)
                            self.pes.add_pair_of_neighbor_minima(current_basin_index, basin_index)

                        self.time_per_confinement[current_basin_index].append(n_quenchs*self.md_time_before_quench)

                    self.trajectory_inherent_structures[current_basin_index].append(basin_index)
                    self.ktn.add_transition(current_basin_index, basin_index)

            n_basins_explored+=1
            self.basins_explored.add(current_basin_index)

            if verbose:
                print('Basin {} explored. {} basins detected.'.format(current_basin_index, self.pes.n_minima))

