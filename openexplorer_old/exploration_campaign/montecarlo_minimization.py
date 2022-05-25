import numpy as np
import simtk.unit as u
from openmmtools.constants import kB
from openexplorer.tools.quench import L_BFGS
from openexplorer.tools.move import DihedralShifts
from openexplorer.tools.acceptance import MetropolisHastings

class MonteCarloMinimization():

    explorer = None

    quench = None
    move = None
    acceptance = None

    pes = None
    trajectory_inherent_structures = []
    global_minimum_inherent_structures = []
    global_minimum_potential_energies = []*u.kilojoules_per_mole

    _previous_coordinates = None
    _previous_potential_energy = None

    initialized = False

    def __init__(self, explorer, quench=L_BFGS, move=DihedralShifts, acceptance=MetropolisHastings):

        self.explorer = explorer

        self.quench = quench(self.explorer)
        self.move = move(self.explorer)
        self.acceptance = acceptance(self.explorer)

        from openexplorer import PES

        self.pes=PES(self.explorer.topology, self.explorer.context.getSystem())

        self.reset()

    def reset(self):

        self.acceptance.reset()
        self.pes.reset()
        self.trajectory_inherent_structures = []

    def run(self, n_steps=1, progress_bar=False):

        if progress_bar:

            from tqdm import tqdm
            iterator=tqdm(range(n_steps))

        else:

            iterator=range(n_steps)

        if not self.initialized:

            self.quench()
            self._previous_coordinates = self.explorer.get_coordinates()
            self._previous_potential_energy = self.explorer.get_potential_energy()
            isi = self.pes.collect_minimum(self.explorer, coordinates=self._previous_coordinates,
                                           potential_energy=self._previous_potential_energy)
            gmis = self.pes.global_minimum_index
            gmpe = self.pes.global_minimum_potential_energy
            self.trajectory_inherent_structures.append(isi)
            self.global_minimum_inherent_structures.append(gmis)
            self.global_minimum_potential_energies.append(gmpe)
            self.initialized = True

        isi = self.trajectory_inherent_structures[-1]
        gmis = self.global_minimum_inherent_structures[-1]
        gmpe = self.global_minimum_potential_energies[-1]

        for _ in iterator:

            self.move()
            self.quench()
            self.acceptance(previous_coordinates=self._previous_coordinates,
                            previous_potential_energy=self._previous_potential_energy)

            if self.acceptance.accepted:

                self._previous_coordinates = self.acceptance.coordinates
                self._previous_potential_energy = self.acceptance.potential_energy

                isi = self.pes.collect_minimum(self.explorer, coordinates= self._previous_coordinates, potential_energy= self._previous_potential_energy)
                gmis = self.pes.global_minimum_index
                gmpe = self.pes.global_minimum_potential_energy

            self.trajectory_inherent_structures.append(isi)
            self.global_minimum_inherent_structures.append(gmis)
            self.global_minimum_potential_energies.append(gmpe)

