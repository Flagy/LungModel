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
    
    def model(self, y, t):
        iter = np.nditer(self.Vg)
        Vg = next(iter)
        dydt = -(1/self.Tau)*(y - Vg)
        return dydt

    def solveModel(self):
        self.getTau()
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
        self.Tau = self.resistance*self.compliance
        return self.Tau

if __name__ == "__main__":
    lung = EasyLung(2)
    lung.compliance = 0.8
    t0 = 0
    tf = 20000
    t = np.linspace(t0, tf)
    Vg0 = 200
    f = 100
    Vg = Vg0*np.sin(2*pi*f*t)
    Vc0 = 5
    lung.setModelParams(t, Vg, Vc0)
    y = lung.solveModel()

    # plot results
    plt.plot(t,y,'r-',label='Output (y(t))')
    plt.ylabel('values')
    plt.xlabel('time')
    plt.legend(loc='best')
    plt.show()