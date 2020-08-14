from .l_bfgs import L_BFGS
from .fire import FIRE
from .gradient_descent import Gradient_descent

class Quench():

    _explorer = None

    l_bfgs = None
    fire = None
    gradient_descent = None

    def __init__(self, explorer):

        self._explorer = explorer

        self.l_bfgs = L_BFGS(explorer)
        self.fire = FIRE(explorer)
        self.gradient_descent = Gradient_descent(explorer)


