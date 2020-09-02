from .dihedral_shifts import DihedralShifts
from .cartesian_shifts import CartesianShifts

class Move():

    _explorer = None

    dihedral_shifts = None
    cartesian_shifts = None

    def __init__(self, explorer):

        self._explorer = explorer

        self.dihedral_shifts = DihedralShifts(explorer)
        self.cartesian_shifts = CartesianShifts(explorer)

