import numpy as np
import simtk.unit as u
from simtk.unit import Quantity
from openmmtools.integrators import GradientDescentMinimizationIntegrator
from simtk.openmm import Context

class Gradient_descent():

    _explorer = None
    _initialized = False
    _context = None
    _integrator = None

    _initial_step_size = Quantity(0.01, u.angstroms)
    _tolerance = Quantity(1.0, u.kilojoules_per_mole)

    def __init__(self, explorer):

        self._explorer=explorer

    def _initialize(self):

        system = self._explorer.context.getSystem()
        platform = self._explorer.context.getPlatform()
        properties = {}
        if platform.getName()=='CUDA':
            properties['CudaPrecision'] = 'mixed'

        self._integrator = GradientDescentMinimizationIntegrator(initial_step_size=self._initial_step_size)
        self._context = Context(system, self._integrator, platform, properties)
        self._initialized = True

    def set_parameters(self, tolerance=Quantity(1.0, u.kilojoules_per_mole),
            initial_step_size=Quantity(0.01, u.angstroms)):

        if not self._initialized:
            self._initialize()

        self._tolerance = tolerance.in_units_of(u.kilojoules_per_mole)
        self._initial_step_size = initial_step_size.in_units_of(u.nanometers)
        self._integrator.setGlobalVariableByName('step_size', self._initial_step_size._value)

    def _set_coordinates(self, coordinates):

        self._context.setPositions(coordinates)

    def _get_coordinates(self):

        return self._context.getState(getPositions=True).getPositions(asNumpy=True)

    def _coordinates_to_explorer(self):

        self._explorer.set_coordinates(self._get_coordinates())

    def _coordinates_from_explorer(self):

        self._set_coordinates(self._explorer.get_coordinates())

    def run(self, steps=0):

        if not self._initialized:

            self._initialize()

        self._coordinates_from_explorer()

        try:
            if steps == 0:
                delta=np.infty
                while delta > self._tolerance._value:
                    self._integrator.step(50)
                    delta=self._integrator.getGlobalVariableByName('delta_energy')
            else:
                integrator.step(steps)

            self._coordinates_to_explorer()
        except Exception as e:
            if str(e) == 'Particle coordinate is nan':
                print('NaN encountered in gradient descent minimizer; falling back to L-BFGS after resetting positions')
                self._explorer.quench.l_bfgs()
            else:
                raise e

        self._initialize()

    def __call__(self, *args, **kwargs):

        return self.run(*args, **kwargs)

