import numpy as np
import simtk.unit as u
from simtk.unit import Quantity
from openmmtools.constants import kB
from openexplorer.tools.quench import L_BFGS
from openexplorer.tools.md import Langevin

class QuenchAndRestore():

    explorer = None

    quench = None
    md = None


    md_time_before_quench = None
    md_steps_before_quench = None

    pes = None
    time = []*u.nanoseconds
    trajectory_inherent_structures = []

    initialized = False

    def __init__(self, explorer, md_time_before_quench = Quantity(1.0, u.picoseconds), quench=L_BFGS, md=Langevin):

        from openexplorer import PES

        self.explorer = explorer

        self.quench = quench(self.explorer)
        self.md = md(self.explorer)

        self.md.set_parameters(temperature=Quantity(500.0, u.kelvin))

        self.md_time_before_quench = md_time_before_quench
        md_timestep = self.md.get_parameters()['timestep']
        self.md_steps_before_quench = int(self.md_time_before_quench/md_timestep)

        self.pes=PES(self.explorer.topology, self.explorer.context.getSystem())

        self.reset()


    def reset(self):

        self.pes.reset()
        self.trajectory_inherent_structures = []
        self.time = []*u.nanoseconds

        self.initialized = False

    def run(self, n_steps=1, time=None, progress_bar=False):

        if time is not None:
            n_steps = int(time/self.md_time_before_quench)

        if not self.initialized:
            coordinates = self.explorer.get_coordinates()
            velocities = self.explorer.get_velocities()
            self.quench()
            isi = self.pes.collect_minimum(self.explorer)
            self.trajectory_inherent_structures.append(isi)
            self.time.append(self.md.get_time())
            self.explorer.set_coordinates(coordinates)
            self.explorer.set_velocities(velocities)
            self.initialized = True
        else:
            isi = self.trajectory_inherent_sturctures[-1]

        if progress_bar:
            from tqdm import tqdm
            iterator = tqdm(range(n_steps))
        else:
            iterator = range(n_steps)


        for iteration_index in iterator:

            self.md(self.md_steps_before_quench)
            coordinates = self.explorer.get_coordinates()
            velocities = self.explorer.get_velocities()
            self.quench()
            new_isi = self.pes.collect_minimum(self.explorer)
            self.pes.add_transition_between_minima(origin=isi, end=new_isi)
            self.trajectory_inherent_structures.append(new_isi)
            self.time.append(self.md.get_time())
            self.explorer.set_coordinates(coordinates)
            self.explorer.set_velocities(velocities)
            isi=new_isi

