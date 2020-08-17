import simtk.unit as u
from simtk.unit import Quantity
from openmmtools.integrators import FIREMinimizationIntegrator
from simtk.openmm import Context

class FIRE():

    _explorer = None
    _initialized = False
    _context = None
    _integrator = None

    _timestep = Quantity(0.001, u.picoseconds)
    _tolerance = None
    _alpha = 0.1
    _dt_max = Quantity(0.01, u.picoseconds)
    _f_inc = 1.1
    _f_dec = 0.5
    _f_alpha = 0.99
    _N_min = 5

    def __init__(self, explorer):

        self._explorer=explorer

    def _initialize(self):

        system = self._explorer.context.getSystem()
        platform = self._explorer.context.getPlatform()
        properties = {}
        if platform.getName()=='CUDA':
            properties['CudaPrecision'] = 'mixed'

        self._integrator = FIREMinimizationIntegrator(timestep=self._timestep,
                tolerance=self._tolerance, alpha=self._alpha,
                dt_max=self._dt_max, f_inc=self._f_inc, f_dec=self._f_dec, f_alpha=self._f_alpha,
                N_min=self._N_min)

        self._context = Context(system, self._integrator, platform, properties)

        self._initialized = True

    def set_parameters(self, timestep=Quantity(1.0, u.femtoseconds), tolerance=None, alpha=0.1,
            dt_max=Quantity(10.0, u.femtoseconds), f_inc=1.1, f_dec=0.5, f_alpha=0.99, N_min=5):

        self._timestep = timestep.in_units_of(u.picoseconds)
        self._tolerance = tolerance
        self._alpha = alpha
        self._dt_max = dt_max.in_units_of(u.picoseconds)
        self._f_inc = f_inc
        self._f_dec = f_dec
        self._f_alpha = f_alpha
        self._N_min = N_min

        self._initialize()

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
                while self._integrator.getGlobalVariableByName('converged') < 1:
                        self._integrator.step(50)
            else:
                self._integrator.step(steps)
            self._coordinates_to_explorer()
        except Exception as e:
            if str(e) == 'Particle coordinate is nan':
                print('NaN encountered in FIRE minimizer; falling back to L-BFGS after resetting positions')
                self._explorer.quench.l_bfgs()
            else:
                raise e

        self._initialize()

    def __call__(self, *args, **kwargs):

        return self.run(*args, **kwargs)

