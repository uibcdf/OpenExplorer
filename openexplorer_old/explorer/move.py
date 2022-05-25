
class Move():

    _explorer = None

    dihedral_shifts = None
    cartesian_shifts = None

    def __init__(self, explorer):

        from openexplorer.tools import move

        self._explorer = explorer

        self.dihedral_shifts = move.DihedralShifts(explorer)
        self.cartesian_shifts = move.CartesianShifts(explorer)

