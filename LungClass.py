from math import sqrt

class EasyLung(object):
    """
    Description of the class:
    This class aims to implement a very basic model of the lung based on Weibel architecture
    and Lambert model for mechanical properties of the airway.
    In future a choice of the model might be implemented in the constructor.
    """
    totGenerations = 27

    def __init__(self, numberOfGen):
        self.numberOfGen = int(numberOfGen) # We may be interested in creating only the first N generations
        self.tracheaDiameter = 1.8 # cm
        self.h = 0.5**(1.0/3) # h can be 0.82 in humans if taken as constant or 0.5**(1.0/3) ~ 0.79
        #self.h = 0.82
        self.r0 = 0.0
        self.r1 = 0.0
        self.r2 = 0.0
        self.c1 = 0.0
        self.c2 = 0.0
        self.createModel()

    def createModel(self):
        R = [self.r0]
        C = [self.c0]
        for gen in range(self.numberOfGen):
            R.append
    
    def getDiameterFromGen(self, numGen):
        return self.d0 * self.h**numGen

    def getResistanceFromGen(self, numGen):
        """Poiseuille Resistance, Weibel Model"""
        d = self.getDiameterFromGen(numGen)
        return (pi*self.P*((d/2.0)**4))/(8.0*eta*self.l)

    def nothing(self, nothing):
        pass
