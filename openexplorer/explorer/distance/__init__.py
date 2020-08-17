from .least_rmsd import Least_RMSD
from .rmsd import RMSD
from .angular_rmsd import Angular_RMSD

class Distance():

    _explorer = None

    least_rmsd = None
    rmsd = None
    angular_rmsd = None

    def __init__(self, explorer):

        self._explorer = explorer

        self.least_rmsd = Least_RMSD(explorer)
        self.rmsd = RMSD(explorer)
        self.angular_rmsd = Angular_RMSD(explorer)

