import simtk.unit as unit
from molsysmt.native.molsys import MolSys
from molsysmt.native.trajectory import Trajectory
import molsysmt as msm
import time

class OpenExplorerReporter():

    def __init__(self, reportInterval, selection='all', syntaxis='MolSysMT',
                 step=True, coordinates=True, boxVectors=True, potentialEnergy=False):

        self._initialized = False
        self._selection = selection
        self._syntaxis = syntaxis

        self._reportInterval = reportInterval
        self._step = step
        self._coordinates = coordinates
        self._boxVectors = boxVectors
        self._potentialEnergy = potentialEnergy

        self.topology = None
        self.step = None
        self.coordinates = None
        self.box = None
        self.potential_energy = None

        self._needsPositions = self._coordinates
        self._needsVelocities = False
        self._needsForces = False
        self._needsEnergy = self._potentialEnergy

    #def describeNextReport(self, simulation):
    #    steps = self._reportInterval - simulation.currentStep%self._reportInterval
    #    return (steps, self._needsPositions, self._needsVelocities, self._needsForces, self._needsEnergy)

    def _initialize(self, exploration_campaign):

        topology = exploration_campaign.explorer.topology

        if self._selection is not 'all':
            self._atom_indices = msm.select(topology, selection=self._selection, syntaxis=self._syntaxis)
            self.topology = msm.convert(topology, to_form='molsysmt.Topology', selection=self._atom_indices)
        else:
            self._atom_indices = 'all'
            self.topology = msm.convert(topology, to_form='molsysmt.Topology')

        if self._step:
            self.step = []

        if self._coordinates:
            self.coordinates = []*unit.nanometers

        if self._boxVectors:
            self.box = []*unit.nanometers

        if self._potentialEnergy:
            self.potential_energy = []*unit.kilojoules_per_mole

        self._initialized=True

    def report(self, exploration_campaign):

        if not self._initialized:
            self._initialize(exploration_campaign)

        if self._step:
            value = exploration_campaign.step
            self.step.append(value)

        if self._coordinates:
            value=exploration_campaign.explorer.get_coordinates()
            if self._atom_indices is 'all':
                self.coordinates.append(value)
            else:
                value = value[self._atom_indices,:]
                self.coordinates.append(value)

        if self._boxVectors:
            value=exploration_campaign.explorer.getState().getPeriodicBoxVectors(asNumpy=True)
            self.box.append(value)

        if self._potentialEnergy:
            value=exploration_campaign.explorer.get_potential_energy()
            self.potential_energy.append(value)

