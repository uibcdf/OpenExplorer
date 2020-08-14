from simtk.openmm import LocalEnergyMinimizer

class L_BFGS():

    _explorer = None
    _initialized = True

    _tolerance = None
    _max_iter = None

    def __init__(self, explorer):

        self._explorer=explorer

    def _initialize(self):

        self.set_parameters()
        self._initialized = True

    def set_parameters(self, tolerance=1.0*unit.kilojoules_per_mole/unit.nanometers, max_iter=0):

        self._tolerance = tolerance
        self._max_iter = max_iter

    def run(self):

        if not self._initialized:

            self._initialize()

        LocalEnergyMinimizer.minimize(self._explorer.context, self._tolerance, self._max_iter)

    def __call__(self, *args, **kwargs):

        return self.run(*args, **kwargs)

