import numpy as np
from simtk.unit import Quantity
import simtk.unit as u
import numpy as np

class DihedralShifts():

    _explorer = None
    _initialized = False

    mode_angles = 'random' # 'all', 'random'
    n_random_angles = 1
    mode_steps = 'random' # 'random', 'random_orientation', 'rmsd', 'random_rmsd'
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

        self.set_parameters()

    def set_parameters(self, dihedral_angle='all', selection='all', quartets=None, blocks=None, mode_angles='random', n_random_angles=1,
                       stepsize=Quantity(value=180.0, unit=u.degrees), mode_steps='random', syntaxis='MolSysMT'):

        from molsysmt import covalent_blocks, covalent_dihedral_quartets

        self.mode_angles = mode_angles
        self.n_random_angles = n_random_angles
        self.stepsize = stepsize.in_units_of(u.degrees)
        self.mode_steps = mode_steps

        if quartets is not None:
            self.quartets = quartets
            if blocks is None:
                self.blocks = []
                for quartet in self.quartets:
                    tmp_blocks = covalent_blocks(self._explorer, remove_bonds=[quartet[1], quartet[2]])
                    self.blocks.append(tmp_blocks)
                self.blocks = np.array(self.blocks)
            else:
                self.blocks = blocks
        else:

            if selection is 'all':
                sel='all'
                sel_not_pro='group_name!="PRO"'
            else:
                sel=selection
                sel_not_pro=selection+'and group_name!="PRO"'

            if dihedral_angle == 'all':

                phi_q, phi_blk = covalent_dihedral_quartets(self._explorer, dihedral_angle='phi',
                                                            selection=sel_not_pro,
                                                            with_blocks=True, syntaxis=syntaxis)
                psi_q, psi_blk = covalent_dihedral_quartets(self._explorer, dihedral_angle='psi',
                                                            selection=sel,
                                                            with_blocks=True, syntaxis=syntaxis)
                chi1_q, chi1_blk = covalent_dihedral_quartets(self._explorer, dihedral_angle='chi1',
                                                              selection=sel_not_pro,
                                                              with_blocks=True, syntaxis=syntaxis)
                chi2_q, chi2_blk = covalent_dihedral_quartets(self._explorer, dihedral_angle='chi2',
                                                              selection=sel_not_pro,
                                                              with_blocks=True, syntaxis=syntaxis)
                chi3_q, chi3_blk = covalent_dihedral_quartets(self._explorer, dihedral_angle='chi3',
                                                              selection=sel,
                                                              with_blocks=True, syntaxis=syntaxis)
                chi4_q, chi4_blk = covalent_dihedral_quartets(self._explorer, dihedral_angle='chi4',
                                                              selection=sel,
                                                              with_blocks=True, syntaxis=syntaxis)
                chi5_q, chi5_blk = covalent_dihedral_quartets(self._explorer, dihedral_angle='chi5',
                                                              selection=sel,
                                                              with_blocks=True, syntaxis=syntaxis)

                aux_q = []
                aux_blk = []
                if phi_q.shape[0]>0: aux_q.append(phi_q), aux_blk.append(phi_blk)
                if psi_q.shape[0]>0: aux_q.append(psi_q), aux_blk.append(psi_blk)
                if chi1_q.shape[0]>0: aux_q.append(chi1_q), aux_blk.append(chi1_blk)
                if chi2_q.shape[0]>0: aux_q.append(chi2_q), aux_blk.append(chi2_blk)
                if chi3_q.shape[0]>0: aux_q.append(chi3_q), aux_blk.append(chi3_blk)
                if chi4_q.shape[0]>0: aux_q.append(chi4_q), aux_blk.append(chi4_blk)
                if chi5_q.shape[0]>0: aux_q.append(chi5_q), aux_blk.append(chi5_blk)

                self.quartets = np.vstack(aux_q)
                self.blocks = np.vstack(aux_blk)

                del(aux_q, aux_blk)
                del(phi_q, phi_blk, psi_q, psi_blk)
                del(chi1_q, chi1_blk, chi2_q, chi2_blk, chi3_q, chi3_blk, chi4_q, chi4_blk, chi5_q, chi5_blk)

            elif dihedral_angle == 'backbone':

                phi_q, phi_blk = covalent_dihedral_quartets(self._explorer, dihedral_angle='phi',
                                                            selection=sel_not_pro,
                                                            with_blocks=True, syntaxis=syntaxis)
                psi_q, psi_blk = covalent_dihedral_quartets(self._explorer, dihedral_angle='psi',
                                                            selection=sel,
                                                            with_blocks=True, syntaxis=syntaxis)

                aux_q = []
                aux_blk = []
                if phi_q.shape[0]>0: aux_q.append(phi_q), aux_blk.append(phi_blk)
                if psi_q.shape[0]>0: aux_q.append(psi_q), aux_blk.append(psi_blk)

                self.quartets = np.vstack(aux_q)
                self.blocks = np.vstack(aux_blk)

                del(aux_q, aux_blk)
                del(phi_q, phi_blk, psi_q, psi_blk)

            elif dihedral_angle == 'sidechains':

                chi1_q, chi1_blk = covalent_dihedral_quartets(self._explorer, dihedral_angle='chi1',
                                                              selection=sel_not_pro,
                                                              with_blocks=True, syntaxis=syntaxis)
                chi2_q, chi2_blk = covalent_dihedral_quartets(self._explorer, dihedral_angle='chi2',
                                                              selection=sel_not_pro,
                                                              with_blocks=True, syntaxis=syntaxis)
                chi3_q, chi3_blk = covalent_dihedral_quartets(self._explorer, dihedral_angle='chi3',
                                                              selection=sel,
                                                              with_blocks=True, syntaxis=syntaxis)
                chi4_q, chi4_blk = covalent_dihedral_quartets(self._explorer, dihedral_angle='chi4',
                                                              selection=sel,
                                                              with_blocks=True, syntaxis=syntaxis)
                chi5_q, chi5_blk = covalent_dihedral_quartets(self._explorer, dihedral_angle='chi5',
                                                              selection=sel,
                                                              with_blocks=True, syntaxis=syntaxis)

                aux_q = []
                aux_blk = []
                if chi1_q.shape[0]>0: aux_q.append(chi1_q), aux_blk.append(chi1_blk)
                if chi2_q.shape[0]>0: aux_q.append(chi2_q), aux_blk.append(chi2_blk)
                if chi3_q.shape[0]>0: aux_q.append(chi3_q), aux_blk.append(chi3_blk)
                if chi4_q.shape[0]>0: aux_q.append(chi4_q), aux_blk.append(chi4_blk)
                if chi5_q.shape[0]>0: aux_q.append(chi5_q), aux_blk.append(chi5_blk)

                self.quartets = np.vstack(aux_q)
                self.blocks = np.vstack(aux_blk)

                del(aux_q, aux_blk)
                del(chi1_q, chi1_blk, chi2_q, chi2_blk, chi3_q, chi3_blk, chi4_q, chi4_blk, chi5_q, chi5_blk)


            else:

                self.quartets, self.blocks = covalent_dihedral_quartets(self._explorer,
                                                                        dihedral_angle=dihedral_angle,
                                                                        selection=selection, syntaxis=syntaxis)

        self.n_quartets = self.quartets.shape[0]

        self._rnd_gen_angles = np.random.default_rng()
        self._rnd_gen_steps = np.random.default_rng()
        self._initialized = True

    def replicate_parameters(self, explorer):

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

        from molsysmt import set_dihedral_angles

        if not self._initialized:

            self._initialize()

        if self.mode_angles == 'all':

            if self.mode_steps == 'random_orientation':

                self.shifts_moved = self._rnd_gen_steps.choice([-1, 1], size=self.n_quartets)*self.stepsize

            elif self.mode_steps == 'random':

                self.shifts_moved = self._rnd_gen_steps.uniform([-1.0, 1.0], size=self.n_quartets)*self.stepsize

            elif self.mode_steps == 'rmsd':

                v = self._rnd_gen_steps.uniform(-1.0,1.0, self.n_quartets)
                norm = np.linalg.norm(v)
                uv = v/norm
                self.shifts_moved = self.stepsize * uv

            elif self.mode_steps == 'random_rmsd':

                v = self._rnd_gen_steps.uniform(-1.0,1.0, self.n_quartets)
                norm = np.linalg.norm(v)
                uv = v/norm
                r = np.random.uniform(0.0, 1.0, 1)
                self.shifts_moved = self.stepsize * r * uv


            set_dihedral_angles(self._explorer, quartets=self.quartets, angles_shifts=self.shifts_moved, blocks=self.blocks,
                                pbc=self._explorer.pbc)

        elif self.mode_angles == 'random':

            self.quartets_moved = self._rnd_gen_angles.choice(self.n_quartets, self.n_random_angles, replace=False, shuffle=False)
            self.quartets_moved.sort()

            if self.mode_steps == 'random_orientation':

                self.shifts_moved = self._rnd_gen_steps.choice([-1, 1], size=self.n_random_angles)*self.stepsize

            elif self.mode_steps == 'random':

                self.shifts_moved = self._rnd_gen_steps.uniform([-1.0, 1.0], size=self.n_random_angles)*self.stepsize

            elif self.mode_steps == 'rmsd':

                v = self._rnd_gen_steps.uniform(-1.0,1.0, self.n_random_angles)
                norm = np.linalg.norm(v)
                uv = v/norm
                self.shifts_moved = self.stepsize * r * uv

            elif self.mode_steps == 'random_rmsd':

                v = self._rnd_gen_steps.uniform(-1.0,1.0, self.n_random_angles)
                norm = np.linalg.norm(v)
                uv = v/norm
                r = np.random.uniform(0.0, 1.0, 1)
                self.shifts_moved = self.stepsize * r * uv

            set_dihedral_angles(self._explorer, quartets=self.quartets[self.quartets_moved], angles_shifts=self.shifts_moved,
                                blocks=self.blocks[self.quartets_moved], pbc=self._explorer.pbc)


    def __call__(self, *args, **kwargs):

        return self.run(*args, **kwargs)

