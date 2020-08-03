import numpy as np
import simtk.unit as unit
from openmmtools.constants import kB

class MonteCarlo ():

    def __init__(self, explorer, movement, acceptance='Metropolis-Hastings',
                 temperature=300.0*unit.kelvin, seed=None, reporter=None):

        self.explorer = explorer
        self.movement = movement
        self.coordinates = explorer.get_coordinates()
        self.potential_energy = explorer.get_potential_energy()

        self.temperature = temperature
        self.seed = seed

        self.acceptance = None
        if acceptance=='Metropolis-Hastings':
            self.acceptance = Acceptance_Metropolis_Hastings(temperature, seed)

        self.step = 0
        self.n_tries = 0
        self.n_acceptances = 0
        self.accepted = False

        self.reporter = reporter

        if self.reporter:
            self.reporter.report(self)

    def reset(self, temperature=300.0*unit.kelvin, seed=None):

        self.n_tries = 0
        self.n_acceptances = 0
        self.accepted = False
        self.acceptance.reset(temperature)

    def run(self, n_steps=1):

        for _ in range(n_steps):

            self.n_tries += 1
            self.step += 1

            self.movement(self.explorer)

            new_coordinates = self.explorer.get_coordinates()
            new_potential_energy = self.explorer.get_potential_energy()

            self.accepted = self.acceptance(self.potential_energy, new_potential_energy)

            if self.accepted:

                self.n_acceptances += 1
                self.coordinates = new_coordinates
                self.potential_energy = new_potential_energy

            else:

                self.explorer.set_coordinates(self.coordinates)

            if self.reporter:
                self.reporter.report(self)


class Acceptance_Metropolis_Hastings():

    def __init__(self, temperature, seed=None):

        self.kT = kB*temperature
        self.random_generator = np.random.default_rng(seed)

    def decide(self, energy, new_energy):

        accept=False

        if new_energy < energy:
            accept=True
        else:
            accept = True
            weight = np.exp(-(new_energy - energy) / self.kT)
            w = min(1.0, weight)
            rand = self.random_generator.random()
            if rand > w:
                accept = False

        return accept

    def __call__(self, energy, new_energy):
        return self.decide(energy, new_energy)

    def reset(self, temperature, seed=None):

        self.kT = kB*temperature
        self.random_generator = np.random.default_rng(seed)

