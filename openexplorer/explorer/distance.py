class Distance():

    _explorer = None

    least_rmsd = None
    rmsd = None
    angular_rmsd = None

    def __init__(self, explorer):

        from openexplorer.tools import distance

        self._explorer = explorer

        self.least_rmsd = distance.Least_RMSD(explorer)
        self.rmsd = distance.RMSD(explorer)
        self.angular_rmsd = distance.Angular_RMSD(explorer)

