from IPython.display import Markdown, display
import numpy as np
import simtk.unit as unit
from pandas import Series
import molsysmt as msm

class CartesianShifts():

    _explorer = None
    _initialized = False

    atom_indices = 'all'
    n_atoms = None
    mode_atoms = 'random' # 'all', 'random'
    n_random_atoms = 1
    mode_steps = 'random' # 'random', 'random_orientation', 'rmsd', 'random_rmsd'
    step_size = 0.25*unit.nanometers

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
                       step_size= 0.25*unit.nanometers, mode_steps='random', syntaxis='MolSysMT'):

        if selection is not None:
            atom_indices =  msm.select(self._explorer, selection=selection, syntaxis=syntaxis)

        self.atom_indices = atom_indices
        self.mode_atoms = mode_atoms
        if self.mode_atoms=='all':
            self.n_random_atoms = 0
        else:
            self.n_random_atoms = n_random_atoms
        self.step_size = step_size.in_units_of(unit.nanometers)
        self.mode_steps = mode_steps

        self._initialize()

    def show_parameters(self, verbose=True):

        tmp_dict = {
                'mode_atoms': self.mode_atoms,
                'mode_steps': self.mode_steps,
                'atom_indices': self.atom_indices,
                'n_random_atoms': self.n_random_atoms,
                'step_size': self.step_size
                }

        if verbose:
            display(Markdown(Series(tmp_dict, name='parameters').to_markdown()))
        else:
            return tmp_dict

    def replicate_parameters(self, explorer):

        self.atom_indices = explorer.move.dihedral_shifts.atom_indices
        self.mode_atoms = explorer.move.dihedral_shifts.mode_atoms
        self.n_random_atoms = explorer.move.dihedral_shifts.n_random_atoms
        self.step_size = explorer.move.dihedral_shifts.step_size
        self.mode_steps = explorer.move.dihedral_shifts.mode_steps
        if explorer.move.cartesian_shifts._rnd_gen_atoms is not None:
            self._rnd_gen_atoms = np.random.default_rng()
        if explorer.move.cartesian_shifts._rnd_gen_steps is not None:
            self._rnd_gen_steps = np.random.default_rng()
        self._initialized = explorer.move.cartesian_shifts._initialized

    def run(self):

        if not self._initialized:

            self._initialize()


        if self.mode_atoms == 'all':

            self.atoms_moved = self.atom_indices
            self.shifts_moved = np.zeros([self.atoms_moved.shape[0], 3], dtype=float)*unit.nanometers

        elif self.mode_atoms == 'random':

            self.atoms_moved = self._rnd_gen_atoms.choice(self.atom_indices, self.n_random_atoms, replace=False, shuffle=False)
            self.atoms_moved.sort()
            self.shifts_moved = np.zeros([self.n_random_atoms, 3], dtype=float)*unit.nanometers

        n_atoms_moved = self.atoms_moved.shape[0]

        if self.mode_steps == 'random':

            for ii in range(n_atoms_moved):

                v = self._rnd_gen_steps.normal(0.0, 1.0, 3)
                r = self._rnd_gen_steps.uniform(0.0, 1.0, 1)**(1./3)
                norm = np.linalg.norm(v)
                uv = v/norm
                self.shifts_moved[ii,:] = self.step_size*r*uv

        elif self.mode_steps == 'random_orientation':

            for ii in range(n_atoms_moved):

                v = self._rnd_gen_steps.normal(0.0,1.0,3)
                norm = np.linalg.norm(v)
                uv = v/norm
                self.shifts_moved[ii,:] = self.step_size*uv

        elif self.mode_steps == 'rmsd':

            v = self._rnd_gen_steps.normal(0.0,1.0, n_atoms_moved*3)
            norm = np.linalg.norm(v)
            uv = v/norm
            self.shifts_moved = self.step_size*uv.reshape(n_atoms_moved, 3)

        elif self.mode_steps == 'random_rmsd':

            v = self._rnd_gen_steps.normal(0.0, 1.0, n_atoms_moved*3)
            r = self._rnd_gen_steps.uniform(0.0, 1.0, 1)**(1./n_atoms_moved*3)
            norm = np.linalg.norm(v)
            uv = v/norm
            self.shifts_moved = self.step_size*r*uv.reshape(self.atoms_moved.shape[0], 3)

        msm.translate(self._explorer, translation=self.shifts_moved, selection=self.atoms_moved)

    def __call__(self, *args, **kwargs):

        return self.run(*args, **kwargs)

