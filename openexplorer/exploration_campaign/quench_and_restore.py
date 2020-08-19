import numpy as np
import simtk.unit as u
from simtk.unit import Quantity
from openmmtools.constants import kB
from tqdm import tqdm

class QuenchAndRestore():

    explorer = None

    temperature = Quantity(500.0, u.kelvin)
    collision_rate = Quantity(1.0, u.picoseconds**-1)
    md_timestep = Quantity(2.0, u.femtoseconds)
    time_iteration = Quantity(50.0, u.picoseconds)
    md_steps = int(time_iteration/md_timestep)

    similarity_threshold = Quantity(0.25, u.angstroms)

    pes = None
    time = []*u.nanoseconds
    trajectory_inherent_structures = []

    def __init__(self, explorer):

        from openexplorer import PES

        self.explorer = explorer
        self.pes = PES()

        self.explorer.md.langevin.set_parameters(temperature=self.temperature, timestep=self.md_timestep,
                                                 collision_rate=self.collision_rate)

    def set_parameters(self, temperature = Quantity(500.0, u.kelvin), collision_rate = Quantity(1.0, u.picoseconds**-1),
                       md_timestep = Quantity(2.0, u.femtoseconds), time_iteration = Quantity(50.0, u.picoseconds),
                       similarity_threshold = Quantity(0.25, u.angstroms)):

        self.temperature = temperature
        self.collision_rate = collision_rate
        self.md_timestep = md_timestep
        self.time_iteration = time_iteration
        self.md_steps = int(time_iteration/md_timestep)
        self.similarity_threshold = similarity_threshold

        self.explorer.md.langevin.set_parameters(temperature=self.temperature, timestep=self.md_timestep,
                                                 collision_rate=self.collision_rate)

    def reset(self):

        from openexplorer import PES

        self.pes = PES()

    def run(self, n_iterations=1, time=None, verbose=False):

        if time is not None:

            n_iterations = int(time/self.time_iteration)

        if verbose:

            iterator = tqdm(range(n_iterations))

        else:

            iterator = range(n_iterations)


        for iteration_index in iterator:

            self.explorer.md.langevin(self.md_steps)
            coordinates = self.explorer.get_coordinates()
            velocities = self.explorer.get_velocities()
            self.explorer.quench.l_bfgs()
            inherent_structure_index = self.pes.collect(self.explorer, similarity_threshold=self.similarity_threshold)
            self.trajectory_inherent_structures.append(inherent_structure_index)
            self.time.append(self.explorer.md.langevin.get_time())
            self.explorer.set_coordinates(coordinates)
            self.explorer.set_velocities(velocities)

