import numpy as np
import simtk.unit as u
from simtk.unit import Quantity
from openmmtools.constants import kB
from tqdm import tqdm

class QuenchAndRestore():

    explorer = None

    md_time_before_quench = Quantity(1.0, u.picoseconds)
    similarity_threshold = Quantity(0.010, u.angstroms)

    temperature = Quantity(500.0, u.kelvin)
    collision_rate = Quantity(1.0, u.picoseconds**-1)
    md_timestep = Quantity(2.0, u.femtoseconds)

    md_steps_before_quench = int(md_time_before_quench/md_timestep)

    pes = None
    time = []*u.nanoseconds
    trajectory_inherent_structures = []

    initialized = False

    def __init__(self, explorer):

        from openexplorer import PES

        self.explorer = explorer

        self.set_parameters()
        self.reset()

    def set_parameters(self, md_time_before_quench = Quantity(1.0, u.picoseconds), temperature = Quantity(500.0, u.kelvin),
            collision_rate = Quantity(1.0, u.picoseconds**-1), md_timestep = Quantity(2.0, u.femtoseconds),
            similarity_threshold = Quantity(0.010, u.angstroms)):

        self.md_time_before_quench = md_time_before_quench
        self.similarity_threshold = similarity_threshold

        self.temperature = temperature
        self.collision_rate = collision_rate
        self.md_timestep = md_timestep

        self.md_steps_before_quench = int(self.md_time_before_quench/self.md_timestep)

        self.explorer.md.langevin.set_parameters(temperature=self.temperature, timestep=self.md_timestep,
                                                 collision_rate=self.collision_rate)

    def reset(self):

        from openexplorer import PES

        self.pes = PES()
        self.time = []*u.nanoseconds
        self.trajectory_inherent_structures = []

        self.initialized = False

    def run(self, n_quenchs=1, time=None, verbose=False):

        if not self.initialized:
            coordinates = self.explorer.get_coordinates()
            velocities = self.explorer.get_velocities()
            self.explorer.quench.l_bfgs()
            inherent_structure_index = self.pes.collect(self.explorer, similarity_threshold=self.similarity_threshold)
            self.trajectory_inherent_structures.append(inherent_structure_index)
            self.time.append(self.explorer.md.langevin.get_time())
            self.explorer.set_coordinates(coordinates)
            self.explorer.set_velocities(velocities)
            self.initialized = True

        else:
            inherent_structure_index = self.trajectory_inherent_sturctures[-1]

        if time is not None:
            n_iterations = int(time/self.md_time_before_quench)

        if verbose:
            iterator = tqdm(range(n_quenchs))
        else:
            iterator = range(n_quenchs)


        for iteration_index in iterator:

            self.explorer.md.langevin(self.md_steps_before_quench)
            coordinates = self.explorer.get_coordinates()
            velocities = self.explorer.get_velocities()
            self.explorer.quench.l_bfgs()
            new_inherent_structure_index = self.pes.collect(self.explorer, similarity_threshold=self.similarity_threshold)
            self.pes.add_transition_between_minima(from=inherent_structure_index, to=new_inherent_structure_index)
            self.trajectory_inherent_structures.append(new_inherent_structure_index)
            self.time.append(self.explorer.md.langevin.get_time())
            self.explorer.set_coordinates(coordinates)
            self.explorer.set_velocities(velocities)
            inherent_structure_index=new_inherent_structure_index

