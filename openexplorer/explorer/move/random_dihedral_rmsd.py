import numpy as np
from simtk.unit import Quantity
import simtk.unit as u
from molsysmt import set_dihedral_angles
from .dihedral_tools import get_quartets_and_blocks, get_blocks

class RandomDihedralRMSD():

    _explorer = None
    _initialized = False

    _stepsize = Quantity(value=5.0, unit=u.degrees)
    _quartets = None
    _n_quartets = None
    _blocks = None
    _dihedral_angles = 'all'

    def __init__(self, explorer):

        self._explorer = explorer

    def _initialize(self):

        if self._quartets is None:
            self._quartets, self._blocks = get_quartets_and_blocks(self._explorer, self._dihedral_angles)
            self._n_quartets = self._quartets.shape[0]

        self._initialized = True

    def set_parameters(self, stepsize=Quantity(value=5.0, unit=u.degrees), dihedral_angles='all', quartets=None, blocks=None):

        self._stepsize = stepsize.in_units_of(u.degrees)

        if quartets is not None:
            self._quartets = quartets
            if blocks is not None:
                self._blocks = blocks
        else:
            self._quartets = None
            self._blocks = None
            if dihedral_angles != self._dihedral_angles and self._initialized==True:
                self._dihedral_angles=dihedral_angles
                self._initialize()

        if not self._initialized:
            self.initialize()

    def run(self):

        if not self._initialized:

            self._initialize()

        v = np.random.uniform(-1.0,1.0, self._n_quartets)
        norm = np.linalg.norm(v)
        uv = v/norm
        shifts = self._stepsize * uv
        set_dihedral_angles(self._explorer, quartets=self._quartets, angles_shifts=shifts,
                blocks=self._blocks, pbc=self._explorer.pbc)

    def __call__(self, *args, **kwargs):

        return self.run(*args, **kwargs)

