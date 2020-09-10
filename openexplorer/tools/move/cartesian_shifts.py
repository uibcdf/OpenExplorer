import numpy as np
from simtk.unit import Quantity
import simtk.unit as u
import numpy as np

class CartesianShifts():

    _explorer = None
    _initialized = False

    atom_indices = 'all'
    n_atoms = None
    mode_atoms = 'random' # 'all', 'random'
    n_random_atoms = 1
    mode_steps = 'random' # 'random', 'random_orientation', 'rmsd', 'random_rmsd'
    stepsize = Quantity(value=0.25, unit=u.nanometers)

    _rnd_gen_atoms = None
    _rnd_gen_steps = None

    atoms_moved = None
    shifts_moved = None

    def __init__(self, explorer):

        self._explorer = explorer

    def _initialize(self):

        if self.atom_indices is 'all':
            self.atom_indices = np.arange(self._explorer.n_atoms)
            self.n_atoms = self._explorer.n_atoms
        else:
            self.n_atoms = len(self.atom_indices)

        self._rnd_gen_atoms = np.random.default_rng()
        self._rnd_gen_steps = np.random.default_rng()

        self._initialized = True

    def set_parameters(self, atom_indices='all', selection=None, mode_atoms='random', n_random_atoms=1,
                       stepsize=Quantity(value=0.25, unit=u.nanometers), mode_steps='random', syntaxis='MolSysMT'):

        if selection is not None:
            from molsysmt import select
            atom_indices =  select(self._explorer, selection=selection, syntaxis=syntaxis)

        self.atom_indices = atom_indices
        self.mode_atoms = mode_atoms
        self.n_random_atoms = n_random_atoms
        self.stepsize = stepsize.in_units_of(u.nanometers)
        self.mode_steps = mode_steps

        self._initialize()

    def replicate_parameters(self, explorer):

        self.atom_indices = explorer.move.dihedral_shifts.atom_indices
        self.mode_atoms = explorer.move.dihedral_shifts.mode_atoms
        self.n_random_atoms = explorer.move.dihedral_shifts.n_random_atoms
        self.stepsize = explorer.move.dihedral_shifts.stepsize
        self.mode_steps = explorer.move.dihedral_shifts.mode_steps
        if explorer.move.cartesian_shifts._rnd_gen_atoms is not None:
            self._rnd_gen_atoms = np.random.default_rng()
        if explorer.move.cartesian_shifts._rnd_gen_steps is not None:
            self._rnd_gen_steps = np.random.default_rng()
        self._initialized = explorer.move.cartesian_shifts._initialized

    def run(self):

        from molsysmt import translate

        if not self._initialized:

            self._initialize()


        if self.mode_atoms == 'all':

            self.atoms_moved = self.atom_indices
            self.shifts_moved = np.zeros([self._explorer.n_atoms, 3], dtype=float)*u.nanometers

        elif self.mode_atoms == 'random':

            self.atoms_moved = self._rnd_gen_atoms.choice(self.atom_indices, self.n_random_atoms, replace=False, shuffle=False)
            self.atoms_moved.sort()
            self.shifts_moved = np.zeros([self.n_random_atoms, 3], dtype=float)*u.nanometers

        n_atoms_moved = self.atoms_moved.shape[0]

        if self.mode_steps == 'random':

            for ii in range(n_atoms_moved):

                v = self._rnd_gen_steps.normal(0.0, 1.0, 3)
                r = self._rnd_gen_steps.uniform(0.0, 1.0, 1)**(1./3)
                norm = np.linalg.norm(v)
                uv = v/norm
                self.shifts_moved[ii,:] = self.stepsize*r*uv

        elif self.mode_steps == 'random_orientation':

            for ii in range(n_atoms_moved):

                v = self._rnd_gen_steps.normal(0.0,1.0,3)
                norm = np.linalg.norm(v)
                uv = v/norm
                self.shifts_moved[ii,:] = self.stepsize*uv

        elif self.mode_steps == 'rmsd':

            v = self._rnd_gen_steps.normal(0.0,1.0, n_atoms_moved*3)
            norm = np.linalg.norm(v)
            uv = v/norm
            self.shifts_moved[self.atoms_moved] = self.stepsize*uv.reshape(n_atoms_moved, 3)

        elif self.mode_steps == 'random_rmsd':

            v = self._rnd_gen_steps.normal(0.0, 1.0, n_atoms_moved*3)
            r = self._rnd_gen_steps.uniform(0.0, 1.0, 1)**(1./n_atoms_moved*3)
            norm = np.linalg.norm(v)
            uv = v/norm
            self.shifts_moved = self.stepsize*r*uv.reshape(self.atoms_moved.shape[0], 3)

        translate(self._explorer, translation=self.shifts_moved, selection=self.atoms_moved)

    def __call__(self, *args, **kwargs):

        return self.run(*args, **kwargs)

