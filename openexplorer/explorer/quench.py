class Quench():

    _explorer = None

    l_bfgs = None
    fire = None
    gradient_descent = None

    def __init__(self, explorer):

        from openexplorer.tools import quench

        self._explorer = explorer

        self.l_bfgs = quench.L_BFGS(explorer)
        self.fire = quench.FIRE(explorer)
        self.gradient_descent = quench.GradientDescent(explorer)

