class Movement():

    def __init__(self, stepsize):

        self.stepsize = stepsize

    def scale(self, factor):

        self.stepsize *= factor

    def move(self, coordinates):

        pass

    def __call__(self, *args, **kwargs):

        return self.move(*args, **kwargs)

    #def updateStep(self, accepted, **kwargs):
    #    """feedback from basin hopping if last step was accepted

    #    Parameters
    #    ----------

    #    accepted : boolean
    #        True if the last step was accepted, otherwise false
    #    """
    #    pass

