import numpy as np
from simtk import unit
from .movement import Movement
import molsysmt as msm

def _quartets_and_blocks(explorer, dihedral_angles):

    quartets = None
    blocks = None

    if dihedral_angles=='all':

        phi_q, phi_blk = msm.covalent_dihedral_quartets(explorer.topology, dihedral_angle='phi', selection='group_name!="PRO"', with_blocks=True)
        psi_q, psi_blk = msm.covalent_dihedral_quartets(explorer.topology, dihedral_angle='psi', with_blocks=True)
        chi1_q, chi1_blk = msm.covalent_dihedral_quartets(explorer.topology, dihedral_angle='chi1', selection='group_name!="PRO"', with_blocks=True)
        chi2_q, chi2_blk = msm.covalent_dihedral_quartets(explorer.topology, dihedral_angle='chi2', selection='group_name!="PRO"', with_blocks=True)
        chi3_q, chi3_blk = msm.covalent_dihedral_quartets(explorer.topology, dihedral_angle='chi3', with_blocks=True)
        chi4_q, chi4_blk = msm.covalent_dihedral_quartets(explorer.topology, dihedral_angle='chi4', with_blocks=True)
        chi5_q, chi5_blk = msm.covalent_dihedral_quartets(explorer.topology, dihedral_angle='chi5', with_blocks=True)


        aux_q = []
        aux_blk = []
        if phi_q.shape[0]>0: aux_q.append(phi_q), aux_blk.append(phi_blk)
        if psi_q.shape[0]>0: aux_q.append(psi_q), aux_blk.append(psi_blk)
        if chi1_q.shape[0]>0: aux_q.append(chi1_q), aux_blk.append(chi1_blk)
        if chi2_q.shape[0]>0: aux_q.append(chi2_q), aux_blk.append(chi2_blk)
        if chi3_q.shape[0]>0: aux_q.append(chi3_q), aux_blk.append(chi3_blk)
        if chi4_q.shape[0]>0: aux_q.append(chi4_q), aux_blk.append(chi4_blk)
        if chi5_q.shape[0]>0: aux_q.append(chi5_q), aux_blk.append(chi5_blk)

        quartets = np.vstack(aux_q)
        blocks = np.vstack(aux_blk)

        del(aux_q, aux_blk)
        del(phi_q, phi_blk, psi_q, psi_blk)
        del(chi1_q, chi1_blk, chi2_q, chi2_blk, chi3_q, chi3_blk, chi4_q, chi4_blk, chi5_q, chi5_blk)

    elif dihedral_angles=='backbone':

        phi_q, phi_blk = msm.covalent_dihedral_quartets(explorer.topology, dihedral_angle='phi', selection='group_name!="PRO"', with_blocks=True)
        psi_q, psi_blk = msm.covalent_dihedral_quartets(explorer.topology, dihedral_angle='psi', with_blocks=True)

        aux_q = []
        aux_blk = []
        if phi_q.shape[0]>0: aux_q.append(phi_q), aux_blk.append(phi_blk)
        if psi_q.shape[0]>0: aux_q.append(psi_q), aux_blk.append(psi_blk)

        quartets = np.vstack(aux_q)
        blocks = np.vstack(aux_blk)

        del(aux_q, aux_blk)
        del(phi_q, phi_blk, psi_q, psi_blk)
 
    elif dihedral_angles=='side_chains':

        chi1_q, chi1_blk = msm.covalent_dihedral_quartets(explorer.topology, dihedral_angle='chi1', selection='group_name!="PRO"', with_blocks=True)
        chi2_q, chi2_blk = msm.covalent_dihedral_quartets(explorer.topology, dihedral_angle='chi2', selection='group_name!="PRO"', with_blocks=True)
        chi3_q, chi3_blk = msm.covalent_dihedral_quartets(explorer.topology, dihedral_angle='chi3', with_blocks=True)
        chi4_q, chi4_blk = msm.covalent_dihedral_quartets(explorer.topology, dihedral_angle='chi4', with_blocks=True)
        chi5_q, chi5_blk = msm.covalent_dihedral_quartets(explorer.topology, dihedral_angle='chi5', with_blocks=True)

        aux_q = []
        aux_blk = []
        if chi1_q.shape[0]>0: aux_q.append(chi1_q), aux_blk.append(chi1_blk)
        if chi2_q.shape[0]>0: aux_q.append(chi2_q), aux_blk.append(chi2_blk)
        if chi3_q.shape[0]>0: aux_q.append(chi3_q), aux_blk.append(chi3_blk)
        if chi4_q.shape[0]>0: aux_q.append(chi4_q), aux_blk.append(chi4_blk)
        if chi5_q.shape[0]>0: aux_q.append(chi5_q), aux_blk.append(chi5_blk)

        quartets = np.vstack(aux_q)
        blocks = np.vstack(aux_blk)

        del(aux_q, aux_blk)
        del(chi1_q, chi1_blk, chi2_q, chi2_blk, chi3_q, chi3_blk, chi4_q, chi4_blk, chi5_q, chi5_blk)

    return quartets, blocks

def _blocks(explorer, quartets):

    n_atoms = explorer.n_atoms
    n_quartets = quartets.shape[0]

    blocks = _np.zeros([n_quartets, n_atoms], dtype=int)

    for quartet_index in range(n_quartets):

        quartet = quartets[quartet_index]
        blocks[quartet_index,:] = msm.covalent_blocks(item, remove_bonds=[quartet[1], quartet[2]], output_form='array')

    return blocks

class RandomDihedralShifts(Movement):

    def __init__(self, stepsize=5.0*unit.degrees, dihedral_angles='all', quartets=None, blocks=None):

        super().__init__(stepsize=stepsize)

        self.quartets = None
        self.n_quartets = None
        self.blocks = None
        self._dihedral_angles = dihedral_angles
        if quartets is not None:
            self.quartets = quartets
            self.n_quartets = self.quartets.shape[0]
        if blocks is not None:
            self.blocks = blocks

    def move(self, explorer):

        if self.quartets is None:
            self.quartets, self.blocks = _quartets_and_blocks(explorer, self._dihedral_angles)
            self.n_quartets = self.quartets.shape[0]
        if self.blocks is None:
            self.blocks = _blocks(explorer, quartets=quartets)

        shifts = np.random.choice([-1,1], size=self.n_quartets)*self.stepsize
        msm.set_dihedral_angles(explorer, quartets=self.quartets, angles_shifts=shifts, blocks=self.blocks, pbc=explorer.pbc)


class RandomDihedralMaxShifts(Movement):

    def __init__(self, stepsize=5.0*unit.degrees, dihedral_angles='all', quartets=None, blocks=None):

        super().__init__(stepsize=stepsize)

        self.quartets = None
        self.n_quartets = None
        self.blocks = None
        self._dihedral_angles = dihedral_angles
        if quartets is not None:
            self.quartets = quartets
            self.n_quartets = self.quartets.shape[0]
        if blocks is not None:
            self.blocks = blocks

    def move(self, explorer):

        if self.quartets is None:
            self.quartets, self.blocks = _quartets_and_blocks(explorer, self._dihedral_angles)
            self.n_quartets = self.quartets.shape[0]
        if self.blocks is None:
            self.blocks = _blocks(explorer, quartets=quartets)

        shifts = np.random.uniform(-1.0, 1.0, size=self.n_quartets)*self.stepsize
        msm.set_dihedral_angles(explorer, quartets=self.quartets, angles_shifts=shifts, blocks=self.blocks, pbc=explorer.pbc)


class RandomDihedralRMSD(Movement):

    def __init__(self, stepsize=5.0*unit.degrees, dihedral_angles='all', quartets=None, blocks=None):

        super().__init__(stepsize=stepsize)

        self.quartets = None
        self.n_quartets = None
        self.blocks = None
        self._dihedral_angles = dihedral_angles
        if quartets is not None:
            self.quartets = quartets
            self.n_quartets = self.quartets.shape[0]
        if blocks is not None:
            self.blocks = blocks

    def move(self, explorer):

        if self.quartets is None:
            self.quartets, self.blocks = _quartets_and_blocks(explorer, self._dihedral_angles)
            self.n_quartets = self.quartets.shape[0]
        if self.blocks is None:
            self.blocks = _blocks(explorer, quartets=quartets)

        u = np.random.uniform(-1.0,1.0, self.n_quartets)
        norm = np.linalg.norm(u)
        u = u/norm
        shifts = self.stepsize * u
        msm.set_dihedral_angles(explorer, quartets=self.quartets, angles_shifts=shifts, blocks=self.blocks, pbc=explorer.pbc)


class RandomDihedralMaxRMSD(Movement):

    def __init__(self, stepsize=5.0*unit.degrees, dihedral_angles='all', quartets=None, blocks=None):

        super().__init__(stepsize=stepsize)

        self.quartets = None
        self.n_quartets = None
        self.blocks = None
        self._dihedral_angles = dihedral_angles
        if quartets is not None:
            self.quartets = quartets
            self.n_quartets = self.quartets.shape[0]
        if blocks is not None:
            self.blocks = blocks

    def move(self, explorer):

        if self.quartets is None:
            self.quartets, self.blocks = _quartets_and_blocks(explorer, self._dihedral_angles)
            self.n_quartets = self.quartets.shape[0]
        if self.blocks is None:
            self.blocks = _blocks(explorer, quartets=quartets)

        u = np.random.uniform(-1.0,1.0, self.n_quartets)
        norm = np.linalg.norm(u)
        u = u/norm
        r = np.random.uniform(0.0, 1.0, 1)
        shifts = self.stepsize * r * u
        msm.set_dihedral_angles(explorer, quartets=self.quartets, angles_shifts=shifts, blocks=self.blocks, pbc=explorer.pbc)

