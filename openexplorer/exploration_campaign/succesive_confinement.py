import numpy as np
import simtk.unit as u
from simtk.unit import Quantity
from openmmtools.constants import kB
from .montecarlo import Acceptance_Metropolis_Hastings
from openexplorer.tools.quench import L_BFGS
from openexplorer.tools.md import Langevin

class SuccesiveConfinement():

    explorer = None

    quench = None
    md = None

    n_confinements_per_basin = None
    md_time_period_before_quench = None
    md_steps_before_quench = None

    pes = None
    time_per_confinement = []*u.nanoseconds
    trajectory_inherent_structures = []

    initialized = False

    def __init__(self, explorer, n_confinements_per_basin=250, md_time_before_quench=Quantity(1.0, u.picoseconds), quench=L_BFGS, md=Langevin):

        self.explorer = explorer

        self.quench = quench(self.explorer)
        self.md = md(self.explorer)

        self.md_time_before_quench = md_time_before_quench
        md_timestep = self.md.get_parameters()['timestep']
        self.md_steps_before_quench = int(self.md_time_before_quench/md_timestep)

        self.pes=PES(self.explorer.topology, self.explorer.context.getSystem())

        self.reset()

    def reset(self):

        self.pes = PES()
        self.time_per_confinement = []*u.nanoseconds
        self.trajectory_inherent_structures = []

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



