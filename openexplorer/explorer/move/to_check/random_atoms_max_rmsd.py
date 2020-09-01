import numpy as np
import simtk.unit as u
from simtk.unit import Quantity

class RandomAtomsMaxRMSD():

    _explorer = None
    _initialized = False

    _n_atoms = 0
    _n_dof = 0
    _stepsize = Quantity(value=0.1, unit=u.nanometer)

    def __init__(self, explorer):

        self._explorer = explorer

    def _initialize(self):

        self._n_atoms = self._explorer.n_atoms
        self._n_dof = 3*self._explorer.n_atoms
        self._initialized = True

    def set_parameters(self, stepsize = Quantity(value=1.0, unit=u.nanometer)):

        self._stepsize = stepsize.in_units_of(u.nanometers)
        if not self._initialized:
            self._initialize()

    def replicate_parameters(self, explorer):

        stepsize = explorer.move.random_atoms_max_rmsd._stepsize

        self.set_parameters(stepsize)

    def _coordinates_to_explorer(self, coordinates):

        return self._explorer.set_coordinates(coordinates)

    def _coordinates_from_explorer(self):

        return self._explorer.get_coordinates()

    def run(self):

        if not self._initialized:

            self._initialize()

        coordinates = self._coordinates_from_explorer()

        v = np.random.normal(0.0,1.0, self._n_dof)
        r = np.random.uniform(0.0, 1.0, 1)**(1./self._n_dof)
        norm = np.linalg.norm(v)
        uv = v/norm
        shift = self._stepsize*r*uv.reshape(self._n_atoms, 3)
        coordinates+= shift

        self._coordinates_to_explorer(coordinates)

    def __call__(self, *args, **kwargs):

        return self.run(*args, **kwargs)

