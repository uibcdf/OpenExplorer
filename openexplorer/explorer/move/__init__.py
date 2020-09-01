from .dihedral_shifts import DihedralShifts

class Move():

    _explorer = None

    dihedral_shifts = None

    def __init__(self, explorer):

        self._explorer = explorer

        self.dihedral_shifts = DihedralShifts(explorer)

