import numpy as np
from simtk.unit import Quantity
import simtk.unit as u
from molsysmt import set_dihedral_angles
from molsysmt import covalent_dihedral_quartets, covalent_blocks
import numpy as np

class DihedralShifts():

    _explorer = None
    _initialized = False

    dihedral_angle = 'all'
    mode_angles = 'random' # 'all', 'random'
    n_random_angles = 1
    mode_steps = 'random' # 'fixed', 'random', 'random_sign'
    stepsize = Quantity(value=180.0, unit=u.degrees)
    quartets = None
    n_quartets = None
    blocks = None

    _rnd_gen_angles = None
    _rnd_gen_steps = None

    quartets_moved = None
    shifts_moved = None

    def __init__(self, explorer):

        self._explorer = explorer

    def _initialize(self):

        if self.quartets is None:
            self.quartets = covalent_dihedral_quartets(self._explorer, self.dihedral_angle)
            self.blocks = []
            for quartet in self.quartets:
                tmp_blocks = covalent_blocks(self._explorer, remove_bonds=[quartet[1], quartet[2]])
                self.blocks.append(tmp_blocks)
            self.blocks = np.array(self.blocks)
            self.n_quartets = self.quartets.shape[0]

        self._rnd_gen_angles = np.random.default_rng()
        self._rnd_gen_steps = np.random.default_rng()
        self._initialized = True

    def set_parameters(self, dihedral_angle='all', quartets=None, blocks=None, mode_angles='random', n_random_angles=1,
                       stepsize=Quantity(value=180.0, unit=u.degrees), mode_steps='random'):

        self.dihedral_angle = dihedral_angle
        self.mode_angles = mode_angles
        self.n_random_angles = n_random_angles
        self.stepsize = stepsize.in_units_of(u.degrees)
        self.mode_steps = mode_steps

        if quartets is not None:
            self.quartets = quartets
            self.n_quartets = self.quartets.shape[0]
            if blocks is not None:
                self.blocks = []
                for quartet in self.quartets:
                    tmp_blocks = covalent_blocks(item, remove_bonds=[quartet[1], quartet[2]])
                    self.blocks.append(tmp_blocks)
                self.blocks = np.array(self.blocks)
            else:
                self.blocks = blocks

        self._initialize()

    def replicate_parameters(self, explorer):

        self.dihedral_angle = explorer.move.dihedral_shifts.dihedral_angle
        self.mode_angles = explorer.move.dihedral_shifts.mode_angles
        self.n_random_angles = explorer.move.dihedral_shifts.n_random_angles
        self.stepsize = explorer.move.dihedral_shifts.stepsize
        self.mode_steps = explorer.move.dihedral_shifts.mode_steps
        self.quartets = explorer.move.dihedral_shifts.quartets
        self.n_quartets = explorer.move.dihedral_shifts.n_quartets
        self.blocks = explorer.move.dihedral_shifts.blocks
        if explorer.move.dihedral_shifts._rnd_gen_angles is not None:
            self._rnd_gen_angles = np.random.default_rng()
        if explorer.move.dihedral_shifts._rnd_gen_steps is not None:
            self._rnd_gen_steps = np.random.default_rng()
        self._initialized = explorer.move._initialized

    def run(self):

        if not self._initialized:

            self._initialize()

        if self.mode_angles == 'all':

            if self.mode_steps == 'fixed':

                self.shifts_moved = np.repeat(self.stepsize, self.n_quartets)

            elif self.mode_steps == 'random_sign':

                self.shifts_moved = self._rnd_gen_steps.choice([-1, 1], size=self.n_quartets)*self.stepsize

            elif self.mode_steps == 'random':

                self.shifts_moved = self._rnd_gen_steps.uniform([-1.0, 1.0], size=self.n_quartets)*self.stepsize

            set_dihedral_angles(self._explorer, quartets=self.quartets, angles_shifts=self.shifts_moved, blocks=self.blocks,
                                pbc=self._explorer.pbc)

        elif self.mode_angles == 'random':

            self.quartets_moved = self._rnd_gen_angles.choice(self.n_quartets, self.n_random_angles, replace=False, shuffle=False)
            self.quartets_moved.sort()

            if self.mode_steps == 'fixed':

                self.shifts_moved = np.repeat(self.stepsize, self.n_random_angles)

            elif self.mode_steps == 'random_sign':

                self.shifts_moved = self._rnd_gen_steps.choice([-1, 1], size=self.n_random_angles)*self.stepsize

            elif self.mode_steps == 'random':

                self.shifts_moved = self._rnd_gen_steps.uniform([-1.0, 1.0], size=self.n_random_angles)*self.stepsize

            set_dihedral_angles(self._explorer, quartets=self.quartets[self.quartets_moved], angles_shifts=self.shifts_moved,
                                blocks=self.blocks[self.quartets_moved], pbc=self._explorer.pbc)

    def __call__(self, *args, **kwargs):

        return self.run(*args, **kwargs)

