import numpy as np
import simtk.unit as u
from simtk.unit import Quantity
from openmmtools.constants import kB
from .montecarlo import Acceptance_Metropolis_Hastings

class SuccesiveConfinement():

    explorer = None

    n_confinements_per_basin = 250
    md_time_period_before_quench = Quantity(1.0, u.picoseconds)
    similarity_threshold = Quantity(0.01, u.angstroms)

    temperature = Quantity(500.0, u.kelvin)
    collision_rate = Quantity(1.0, u.picoseconds**-1)
    md_timestep = Quantity(2.0, u.femtoseconds)

    pes = None
    time_per_confinement = []*u.nanoseconds
    trajectory_inherent_structures = []
    n_basins_explored = 0

    initialized = False

    def __init__(self, explorer):

        self.explorer = explorer

        self.set_parameters()
        self.reset()

    def set_parameters(self, n_confinements_per_basin = 250, md_time_period_before_quench = Quantity(1.0, u.picoseconds),
            temperature = Quantity(500.0, u.kelvin), collision_rate = Quantity(1.0, u.picoseconds**-1),
            md_timestep = Quantity(2.0, u.femtoseconds), similarity_threshold = Quantity(0.01, u.angstroms)):

        self.n_confinements_per_basin = n_confinements_per_basin
        self.md_time_period_before_quench = md_time_period_before_quench

        self.temperature = temperature
        self.collision_rate = collision_rate
        self.md_timestep = md_timestep

        self.explorer.md.langevin.set_parameters(temperature=self.temperature,
                timestep=self.md_timestep, collision_rate=self.collision_rate)

    def reset(self):

        self.pes = PES()
        self.time_per_confinement = []*u.nanoseconds
        self.trajectory_inherent_structures = []
        self.n_basins_explored = 0

        self.initialized = False

    def run(self, n_basins=None):

        if n_basins is None:
            n_basins = np.inf

        if not self.initialized:
            self.explorer.quench.l_bfgs()
            current_basin_index = self.pes.collect(self.explorer, similarity_threshold=self.similarity_threshold)
            self.initialized = True
        else:
            if self.n_basins_explored < self.pes.n_minima:
                current_basin_index = self.n_basins_explored
            else:
                self.explorer.quench.l_bfgs()
                current_basin_index = self.pes.collect(self.explorer, similarity_threshold=self.similarity_threshold)
            pass

        current_basin_coordinates = self.pes.minima[current_basin_index]
        self.explorer.set_coordinates(current_basin_coordinates)
        self.explorer.set_velocities_to_temperature(self.temperature)

        n_basins_detected = 0

        while (current_minimum_index < self.pes.n_minima) or (n_basins_detected < n_basins):



