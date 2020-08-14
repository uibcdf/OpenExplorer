import numpy as np
import simtk.unit as unit
from openmmtools.integrators import GradientDescentMinimizationIntegrator

class Gradient_descent():

    _explorer = None
    _initialized = False
    _context = None
    _integrator = None

    _initial_step_size = None
    _tolerance = None

    def __init__(self, explorer):

        self._explorer=explorer

    def _initialize(self):

        system = self._explorer.context.getSystem()
        platform = self._explorer.context.getPlatform()
        if platform.getName()=='CUDA':
            properties['CudaPrecision'] = 'mixed'

        self._integrator = GradientDescentMinimizationIntegrator()
        self._context = Context(system, self._integrator, platform, properties)
        self.set_parameters()
        self._initialized = True

    def set_parameters(self, tolerance=1.0*unit.kilojoules_per_mole, initial_step_size=0.01 * unit.angstroms):

        self._tolerance = tolerance.value_in_unit(unit.kilojoules_per_mole)
        self._initial_step_size = initial_step_size

        self._integrator.setGlobalVariableByName('step_size', initial_step_size / unit.nanometers)

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
                while delta > self._tolerance:
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

    def __call__(self, *args, **kwargs):

        return self.run(*args, **kwargs)

