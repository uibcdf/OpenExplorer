from .langevin import Langevin

class MD():

    _explorer = None

    langevin = None
    verlet = None

    def __init__(self, explorer):

        self._explorer = explorer

        self.langevin = Langevin(explorer)


