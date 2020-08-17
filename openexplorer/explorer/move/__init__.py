from .random_atoms_shifts import RandomAtomsShifts
from .random_atoms_max_shifts import RandomAtomsMaxShifts
from .random_atoms_rmsd import RandomAtomsRMSD
from .random_atoms_max_rmsd import RandomAtomsMaxRMSD
from .random_dihedral_shifts import RandomDihedralShifts
from .random_dihedral_max_shifts import RandomDihedralMaxShifts
from .random_dihedral_rmsd import RandomDihedralRMSD
from .random_dihedral_max_rmsd import RandomDihedralMaxRMSD

class Move():

    _explorer = None

    random_atoms_shifts = None
    random_atoms_max_shifts = None
    random_atoms_rsmd = None
    random_atoms_max_rsmd = None

    random_dihedral_shifts = None
    random_dihedral_max_shifts = None
    random_dihedral_rmsd = None
    random_dihedral_max_rmsd = None

    def __init__(self, explorer):

        self._explorer = explorer

        self.random_atoms_shifts = RandomAtomsShifts(explorer)
        self.random_atoms_max_shifts = RandomAtomsMaxShifts(explorer)
        self.random_atoms_rsmd = RandomAtomsRMSD(explorer)
        self.random_atoms_max_rsmd = RandomAtomsMaxRMSD(explorer)

        self.random_dihedral_shifts = RandomDihedralShifts(explorer)
        self.random_dihedral_max_shifts = RandomDihedralMaxShifts(explorer)
        self.random_dihedral_rmsd = RandomDihedralRMSD(explorer)
        self.random_dihedral_max_rmsd = RandomDihedralMaxRMSD(explorer)

