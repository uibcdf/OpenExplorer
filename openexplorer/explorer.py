import numpy as np

class Explorer():

    system = None

    openmm_system = None
    openmm_topology = None
    openmm_context = None

    def __init__(self, system=None):

        self.system = system

