
class MD():

    _explorer = None

    langevin = None

    def __init__(self, explorer):

        from openexplorer.tools import md

        self._explorer = explorer

        self.langevin = md.Langevin(explorer)


