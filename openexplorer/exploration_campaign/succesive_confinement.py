import numpy as np
import simtk.unit as unit
from openmmtools.constants import kB
from .montecarlo import Acceptance_Metropolis_Hastings

class BasinHopping():

    def __init__(self, explorer, movement, minimizer='L-BFGS', tolerance=1.0*unit.kilojoules_per_mole/unit.nanometers, max_iterations=0,
                 acceptance='Metropolis-Hastings', temperature=300.0*unit.kelvin, seed=None, storage=None):

        self.minimizer = minimizer
        self.minimizer_tolerance = tolerance
        self.minimizer_max_iterations = max_iterations
        self.explorer = explorer
        self.movement = movement

        self.explorer.quench(minimizer=self.minimizer, tolerance=self.minimizer_tolerance,
                             max_iterations=self.minimizer_max_iterations)

        self.coordinates = self.explorer.get_energy()
        self.energy = self.explorer.get_energy()

        self.temperature = temperature
        self.seed = seed

        self.acceptance = None
        if acceptance=='Metropolis-Hastings':
            self.acceptance = Acceptance_Metropolis_Hastings(temperature, seed)

        self.n_tries = 0
        self.n_acceptances = 0
        self.accepted = False

        self.storage = storage

    def reset(self, temperature=300.0*unit.kelvin, seed=None):

        self.n_tries = 0
        self.n_acceptances = 0
        self.accepted = False
        self.acceptance.reset(temperature)

    def run(self, n_steps=1):


        for _ in range(n_steps):

            self.n_tries += 1

            self.movement(self.explorer)
            self.explorer.quench(minimizer=self.minimizer, tolerance=self.minimizer_tolerance,
                                 max_iterations=self.minimizer_max_iterations)

            new_coordinates = self.explorer.get_coordinates()
            new_energy = self.explorer.get_energy()

            self.accepted = self.acceptance(self.energy, new_energy)

            if self.accepted:

                self.explorer.set_coordinates(new_coordinates)
                self.coordinates=new_coordinates
                self.energy=new_energy
                self.n_acceptances += 1

                if self.storage:
                    self.storage.dump(new_coordinates)

