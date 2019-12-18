from math import sqrt, pi
from BronchoClass import EasyBroncho
from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt

class EasyLung(EasyBroncho):
    """
    Description of the class:
    This class should only expose the methods for solving the differential equations of the model of the EasyBroncho class.
    """
    totGenerations = 27

    def __init__(self, numberOfGen):
        if numberOfGen < 0 or numberOfGen > EasyLung.totGenerations:
            raise ValueError
        super().__init__(generationNumber = numberOfGen, paramsFromJson = True)
    
    def setModelParams(self, t, Vg, Vc0):
        self.t = t
        self.Vg = Vg
        self.Vc0 = Vc0
        return True
    
    def model(self, Vc, t):
        self.dVcdt = (1/(self.getTau()))*(self.Vg - Vc)
        return self.dVcdt

    def solveModel(self):
        sol = odeint(self.model, self.Vc0, self.t)
        return sol

    def getDiameterFromGen(self, numGen):
        """Description"""
        return self.d0 * self.h**numGen

    def getResistanceFromGen(self, numGen):
        """Poiseuille Resistance, Weibel Model"""
        d = self.getDiameterFromGen(numGen)
        return (pi*self.P*((d/2.0)**4))/(8.0*eta*self.l)


    def getTau(self):
        return self.resistance*self.compliance

if __name__ == "__main__":
    lung = EasyLung(2)
    lung.compliance = 0.8
    t0 = 0
    t_int = 1
    tf = 100
    t = np.linspace(t0, t_int, tf)
    Vg0 = 2
    f = 1
    Vg = Vg0*np.sin(2*pi*f*t) - Vg0
    Vc0 = np.zeros((len(t), ))
    lung.setModelParams(t, Vg, Vc0)
    y = lung.solveModel()

    # plot results
    plt.plot(t,y,'r-',label='Output (y(t))')
    plt.ylabel('values')
    plt.xlabel('time')
    plt.legend(loc='best')
    plt.show()