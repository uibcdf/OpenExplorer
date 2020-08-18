import numpy as np
import molsysmt as msm
import simtk.unit as unit
from simtk.openmm import app
from simtk.openmm import Context, Platform
from simtk.openmm import VerletIntegrator

class Explorer():

    topology = None
    context = None
    pbc = False
    n_atoms = 0
    md = None
    quench = None
    move = None

    def __init__(self, topology=None, system=None, pbc=False, platform='CUDA'):

        from .md import MD
        from .quench import Quench
        from .move import Move
        from .distance import Distance

        if topology is None:
            raise ValueError('topology is needed')

        if system is None:
            raise ValueError('system is needed')

        step_size  = 1.0*unit.femtoseconds
        integrator = VerletIntegrator(step_size)
        integrator.setConstraintTolerance(0.00001)

        if platform=='CUDA':
            platform = Platform.getPlatformByName('CUDA')
            properties = {'CudaPrecision': 'mixed'}
        elif platform =='CPU':
            platform = Platform.getPlatformByName('CPU')
            properties = {}

        self.topology = topology
        self.context = Context(system, integrator, platform, properties)
        self.n_atoms = msm.get(self.context, target='system', n_atoms=True)

        self.pbc = pbc

        if self.pbc:
            raise NotImplementedError

        self.md = MD(self)
        self.quench = Quench(self)
        self.move = Move(self)
        self.distance = Distance(self)

    def _copy(self):

        topology = self.topology
        coordinates = self.get_coordinates()
        system = self.context.getSystem()
        platform = self.context.getPlatform().getName()
        pbc = self.pbc

        tmp_explorer = Explorer(topology, system, pbc, platform)
        tmp_explorer.set_coordinates(coordinates)

        for ii,jj in vars(tmp_explorer.md).items():
            if not ii.startswith('_'):
                if jj._initialized:
                    jj.replicate_parameters(self)

        for ii,jj in vars(tmp_explorer.quench).items():
            if not ii.startswith('_'):
                if jj._initialized:
                    jj.replicate_parameters(self)

        for ii,jj in vars(tmp_explorer.move).items():
            if not ii.startswith('_'):
                if jj._initialized:
                    jj.replicate_parameters(self)

        return tmp_explorer

    def replicate(self, times=1):

        from copy import deepcopy

        if times==1:

            output = self._copy()

        else:

            output = [self._copy() for ii in range(times)]

        return output

    def set_coordinates(self, coordinates):

        self.context.setPositions(coordinates)

    def get_coordinates(self):

        return self.context.getState(getPositions=True).getPositions(asNumpy=True)

    def set_velocities(self, velocities):

        self.context.setVelocities(velocities)

    def set_velocities_to_temperature(self, temperature):

        self.context.setVelocitiesToTemperature(temperature)

    def get_velocities(self):

        return self.context.getState(getVelocities=True).getVelocities(asNumpy=True)

    def get_potential_energy(self):

        energy = self.context.getState(getEnergy=True).getPotentialEnergy()
        return energy

    def get_potential_energy_gradient(self):

        gradient = -self.context.getState(getForces=True).getForces(asNumpy=True)
        gradient = gradient.ravel()*gradient.unit
        return gradient

    def get_potential_energy_hessian(self, mass_weighted=False, symmetric=True):

        """OpenMM single frame hessian evaluation
        Since OpenMM doesnot provide a Hessian evaluation method, we used finite difference on forces

        from: https://leeping.github.io/forcebalance/doc/html/api/openmmio_8py_source.html

        Returns
        -------
        hessian: np.array with shape 3N x 3N, N = number of "real" atoms
            The result hessian matrix.
            The row indices are fx0, fy0, fz0, fx1, fy1, ...
            The column indices are x0, y0, z0, x1, y1, ..
            The unit is kilojoule / (nanometer^2 * mole * dalton) => 10^24 s^-2
        """

        n_dof = self.n_atoms*3
        pos = self.get_coordinates()
        hessian = np.empty((n_dof, n_dof), dtype=float)*unit.kilojoules_per_mole / (unit.nanometers**2)
        # finite difference step size
        diff = 0.0001*unit.nanometer
        coef = 1.0 /  (2.0*diff) # 1/2h

        for i in range(self.n_atoms):
            # loop over the x, y, z coordinates
            for j in range(3):
                # plus perturbation
                pos[i][j] += diff
                self.set_coordinates(pos)
                grad_plus = self.get_potential_energy_gradient()
                # minus perturbation
                pos[i][j] -= 2*diff
                self.set_coordinates(pos)
                grad_minus = self.get_potential_energy_gradient()
                # set the perturbation back to zero
                pos[i][j] += diff
                # fill one row of the hessian matrix
                hessian[i*3+j] = (grad_plus - grad_minus) * coef

        if mass_weighted:
            mass = np.array([self.context.getSystem().getParticleMass(k).value_in_unit(unit.dalton) for k in range(self.n_atoms)])*unit.dalton
            mass_weight = 1.0/np.sqrt(mass) * (mass.unit**-0.5)
            mass_weight = np.repeat(mass_weight,3)*mass_weight.unit
            hessian = np.multiply(hessian, mass_weight)*hessian.unit*mass_weight.unit
            hessian = np.multiply(hessian, mass_weight[:,np.newaxis])*hessian.unit*mass_weight.unit

        # make hessian symmetric by averaging upper right and lower left
        if symmetric:
            hessian += hessian.T*hessian.unit
            hessian *= 0.5

        # recover the original position
        self.set_coordinates(pos)
        return hessian

#     def normal_modes(self, shot=0, optimize=True):
#        """OpenMM Normal Mode Analysis
#        from: https://leeping.github.io/forcebalance/doc/html/api/openmmio_8py_source.html
#        Since OpenMM doesnot provide a Hessian evaluation method, we used finite difference on forces
#        
#        Parameters
#        ----------
#        shot: int
#            The frame number in the trajectory of this target
#        optimize: bool, default True
#            Optimize the geometry before evaluating the normal modes
#        
#        Returns
#        -------
#        freqs: np.array with shape (3N - 6) x 1, N = number of "real" atoms
#            Harmonic frequencies, sorted from smallest to largest, with the 6 smallest removed, in unit cm^-1
#        normal_modes: np.array with shape (3N - 6) x (3N), N = number of "real" atoms
#            The normal modes corresponding to each of the frequencies, scaled by mass^-1/2.
#        """
#        if self.precision == 'single':
#            warn_once("Single-precision OpenMM engine used for normal mode analysis - recommend that you use mix or double precision.")
#        if not optimize:
#            warn_once("Asking for normal modes without geometry optimization?")
#        # step 0: check number of real atoms
#        noa = len(self.realAtomIdxs)
#        if noa < 2:
#            error('normal mode analysis not suitable for system with one or less atoms')
#        # step 1: build a full hessian matrix
#        hessian_matrix = self.build_mass_weighted_hessian(shot=shot, optimize=optimize)
#        # step 2: diagonalize the hessian matrix
#        eigvals, eigvecs = np.linalg.eigh(hessian_matrix)
#        # step 3: convert eigenvalues to frequencies
#        coef = 0.5 / np.pi * 33.3564095 # 10^12 Hz => cm-1
#        negatives = (eigvals >= 0).astype(int) * 2 - 1 # record the negative ones
#        freqs = np.sqrt(np.abs(eigvals)) * coef * negatives
#        # step 4: convert eigenvectors to normal modes
#        # re-arange to row index and shape
#        normal_modes = eigvecs.T.reshape(noa*3, noa, 3)
#        # step 5: Remove mass weighting from eigenvectors
#        massList = np.array(self.AtomLists['Mass'])[self.realAtomIdxs] # unit in dalton
#        for i in range(normal_modes.shape[0]):
#            mode = normal_modes[i]
#            mode /= np.sqrt(massList[:,np.newaxis])
#            mode /= np.linalg.norm(mode)
#        # step 5: remove the 6 freqs with smallest abs value and corresponding normal modes
#        n_remove = 5 if len(self.realAtomIdxs) == 2 else 6
#        larger_freq_idxs = np.sort(np.argpartition(np.abs(freqs), n_remove)[n_remove:])
#        freqs = freqs[larger_freq_idxs]
#        normal_modes = normal_modes[larger_freq_idxs]
#        return freqs, normal_modes

