class RMSD():

    _explorer = None
    _initialized = False

    _selection = None
    _syntaxis = None
    _atom_indices = None

    def __init__(self, explorer):

        self._explorer=explorer

    def set_parameters(self, selection='atom_type!="H"', syntaxis='MolSysMT'):

        from molsysmt import select

        self._selection = selection
        self._syntaxis = syntaxis
        self._atom_indices = select(self._explorer, selection=selection, syntaxis=syntaxis)

    def run(self, reference, selection='atom_type!="H"', syntaxis='MolSysMT'):

        from molsysmt import rmsd

        if selection!=self._selection:
            self.set_parameters(selection=selection, syntaxis=syntaxis)

        return rmsd(self._explorer, reference_item=reference, selection=self._atom_indices)[0]

    def __call__(self, *args, **kwargs):

        return self.run(*args, **kwargs)

