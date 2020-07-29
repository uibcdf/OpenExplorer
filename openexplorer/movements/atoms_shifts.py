import numpy as np
from simtk import unit
from .movement import Movement

class RandomAtomsShifts(Movement):

    def __init__(self, stepsize=0.1*unit.nanometers):

        super().__init__(stepsize=stepsize)


    def move(self, explorer):

        coordinates = explorer.get_coordinates()

        for ii in range(explorer.n_atoms):

            u = np.random.normal(0.0,1.0,3)
            norm = np.linalg.norm(u)
            u = u/norm
            shift = self.stepsize*u
            coordinates[ii,:]+=shift

        explorer.set_coordinates(coordinates)


class RandomAtomsMaxShifts(Movement):

    def __init__(self, stepsize=0.1*unit.nanometers):

        super().__init__(stepsize=stepsize)


    def move(self, explorer):

        coordinates = explorer.get_coordinates()

        for ii in range(explorer.n_atoms):

            u = np.random.normal(0.0, 1.0, 3)
            r = np.random.uniform(0.0, 1.0, 1)**(1./3)
            norm = np.linalg.norm(u)
            u = u/norm
            shift = self.stepsize*r*u
            coordinates[ii,:]+=shift

        explorer.set_coordinates(coordinates)


class RandomAtomsRMSD(Movement):

    def __init__(self, stepsize=0.1*unit.nanometers):

        super().__init__(stepsize=stepsize)


    def move(self, explorer):

        coordinates = explorer.get_coordinates()

        u = np.random.normal(0.0,1.0, 3*explorer.n_atoms)
        norm = np.linalg.norm(u)
        u = u/norm
        shift = self.stepsize * u.reshape(explorer.n_atoms, 3)
        coordinates+= shift

        explorer.set_coordinates(coordinates)


class RandomAtomsMaxRMSD(Movement):

    def __init__(self, stepsize=0.1*unit.nanometers):

        super().__init__(stepsize=stepsize)


    def move(self, explorer):

        coordinates = explorer.get_coordinates()

        n_dof = 3*explorer.n_atoms
        u = np.random.normal(0.0,1.0, n_dof)
        r = np.random.uniform(0.0, 1.0, 1)**(1./n_dof)
        norm = np.linalg.norm(u)
        u = u/norm
        shift = self.stepsize *r*u.reshape(explorer.n_atoms, 3)
        coordinates+= shift

        explorer.set_coordinates(coordinates)


