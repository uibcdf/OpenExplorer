import numpy as np
from simtk.unit import Quantity
import simtk.unit as u
from molsysmt import set_dihedral_angles
from .dihedral_tools import get_quartets_and_blocks, get_blocks

class RandomDihedralShifts():

    _explorer = None
    _initialized = False

    mode = 'all'
    stepsize = Quantity(value=5.0, unit=u.degrees)
    dihedral_angles = 'all'
    quartets = None
    n_quartets = None
    blocks = None

    def __init__(self, explorer):

        self._explorer = explorer

    def _initialize(self):

        if self._quartets is None:
            self._quartets, self._blocks = get_quartets_and_blocks(self._explorer, self._dihedral_angles)
            self._n_quartets = self._quartets.shape[0]

        self._initialized = True

    def set_parameters(self, mode='all', dihedral_angles='all', stepsize=Quantity(value=5.0, unit=u.degrees), quartets=None, blocks=None):

        self.stepsize = stepsize.in_units_of(u.degrees)
        self.mode = mode

        if quartets is not None:
            self.dihedral_angles = None
            self.quartets = quartets
            self.n_quartets = self._quartets.shape[0]
            if blocks is not None:
                self.blocks = blocks
            else:
                self.blocks = get_blocks(quartets)
        else:
            self._quartets = None
            self._blocks = None
            if dihedral_angles != self._dihedral_angles and self._initialized==True:
                self._dihedral_angles=dihedral_angles
                self._initialize()

        if not self._initialized:
            self.initialize()

    def replicate_parameters(self, explorer):

        stepsize = explorer.move.random_dihedral_max_rmsd._stepsize
        dihedral_angles = explorer.move.random_dihedral_max_rmsd._dihedral_angles
        quartets = explorer.move.random_dihedral_max_rmsd._quartets
        blocks = explorer.move.random_dihedral_max_rmsd._blocks

        self.set_parameters(stepsize, dihedral_angles, quartets, blocks)

    def run(self):

        if not self._initialized:

            self._initialize()

        shifts = np.random.choice([-1,1], size=self._n_quartets)*self._stepsize
        set_dihedral_angles(self._explorer, quartets=self._quartets, angles_shifts=shifts,
                blocks=self._blocks, pbc=self._explorer.pbc)

    def __call__(self, *args, **kwargs):

        return self.run(*args, **kwargs)

