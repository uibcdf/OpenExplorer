
class Acceptance():

    _explorer = None

    metropolis_hastings = None

    def __init__(self, explorer):

        from openexplorer.tools import acceptance

        self._explorer = explorer

        self.metropolis_hastings = acceptance.MetropolisHastings(explorer)

