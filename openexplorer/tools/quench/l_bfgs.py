from simtk.openmm import LocalEnergyMinimizer
import simtk.unit as u
from simtk.unit import Quantity

class L_BFGS():

    _explorer = None
    _initialized = False

    _tolerance = Quantity(1.0, u.kilojoules_per_mole/u.nanometers)
    _max_iter = 0

    def __init__(self, explorer):

        self._explorer=explorer

    def _initialize(self):

        self._initialized = True

    def set_parameters(self, tolerance=Quantity(1.0, u.kilojoules_per_mole/u.nanometers), max_iter=0):

        self._tolerance = tolerance.in_units_of(u.kilojoules_per_mole/u.nanometers)
        self._max_iter = max_iter

        self._initialize()

    def replicate_parameters(self, explorer):

        tolerance = explorer.quench.l_bfgs._tolerance
        max_iter = explorer.quench.l_bfgs._max_iter

        self.set_parameters(tolerance, max_iter)

    def run(self):

        if not self._initialized:

            self._initialize()

        LocalEnergyMinimizer.minimize(self._explorer.context, self._tolerance, self._max_iter)

    def __call__(self, *args, **kwargs):

        return self.run(*args, **kwargs)

