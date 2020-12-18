from simtk.unit import Quantity
from simtk import unit as u
import numpy as np
from openmmtools.constants import kB

class MetropolisHastings():

    _explorer = None
    _initialized = False

    temperature = None
    _kT = None

    n_tries = 0
    n_accepted = 0
    accepted = False

    coordinates = None
    previous_coordinates = None
    potential_energy = None
    previous_potential_energy = None

    weight = None
    random = None

    _rnd_gen = None

    def __init__(self, explorer):

        self._explorer = explorer

    def _initialize(self):

        self.reset()
        self._initialized = True

    def set_parameters(self, temperature=Quantity(298.0, u.kelvin)):

        self.temperature = temperature.in_units_of(u.kelvin)
        self._kT = (kB*temperature).in_units_of(u.kilojoules_per_mole)

    def reset(self):

        self.n_tries = 0
        self.n_accepted = 0
        self.accepted = False
        self._rnd_gen = np.random.default_rng()

        self.coordinates = None
        self.previous_coordinates = None
        self.potential_energy = None
        self.previous_potential_energy = None

        self.weight = None
        self.random = None

    def run(self, previous_coordinates=None, previous_potential_energy=None, coordinates=None,
            potential_energy=None, update_explorer=True):

        if not self._initialized:

            self.set_parameters()
            self._initialize()

        self.accepted=False

        current_coordinates = self._explorer.get_coordinates()

        self.previous_coordinates = previous_coordinates

        if previous_potential_energy is None:

            if previous_coordinates is None:

                raise ValueError('The input argument "previous_coordinates" is necessary if no previous potential energy is provided')

            else:

                self._explorer.set_coordinates(self.previous_coordinates)
                self.previous_potential_energy = self._explorer.get_potential_energy()
                self._explorer.set_coordinates(current_coordinates)

        else:

            self.previous_potential_energy = previous_potential_energy

        if coordinates is None:

            self.coordinates = current_coordinates

        else:

            self.coordinates = coordinates


        if potential_energy is None:

            self.potential_energy = self._explorer.get_potential_energy()

        else:

            self.potential_energy = potential_energy

        if self.potential_energy <= self.previous_potential_energy:

            self.accepted = True

            self.weight = 1.0
            self.random = None

        else:

            self.weight = np.exp(-(self.potential_energy - self.previous_potential_energy) / self._kT)
            self.random = self._rnd_gen.random()

            if self.random <= self.weight:
                self.accepted = True
            else:
                self.accepted = False

        if update_explorer:

            if not self.accepted:

                self._explorer.set_coordinates(self.previous_coordinates)


        self.n_tries += 1

        if self.accepted:

            self.n_accepted += 1


    def __call__(self, *args, **kwargs):

        return self.run(*args, **kwargs)

