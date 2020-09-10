from .metropolis_hastings import MetropolisHastings

class Acceptance():

    _explorer = None

    metropolis_hastings = None

    def __init__(self, explorer):

        self._explorer = explorer

        self.metropolis_hastings = MetropolisHastings(explorer)

