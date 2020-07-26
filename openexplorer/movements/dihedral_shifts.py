import numpy as np
from simtk import unit
from .movement import Movement
from molsysmt import dihedral_angles

class RandomDihedralShifts(Movement):

    def __init__(self, stepsize=5.0*unit.degrees, covalent_chains=None):

        super().__init__(stepsize=stepsize)

        self.covalent_chains = covalent_chains

    def move(self, explorer):

        angles = dihedral_angles(explorer, quartets=self.covalent_chains)
        shifts = np.random.choice([-1,1],size=n_angles)
        coordinates = explorer.get_coordinates()



        for ii in range(explorer.n_atoms):

            u = np.random.normal(0.0,1.0,3)
            norm = np.linalg.norm(u)
            u = u/norm
            shift = self.stepsize*u
            coordinates[ii,:]+=shift

        explorer.set_coordinates(coordinates)


class RandomDihedralMaxShifts(Movement):

    def __init__(self, stepsize=5.0*unit.degrees, covalent_chains=None):

        super().__init__(stepsize=stepsize)

        self.covalent_chains = covalent_chains

    def move(self, explorer):

        pass

