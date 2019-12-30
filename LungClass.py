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
        self.bronchi = []
        for gen in range(numberOfGen + 1):
            self.bronchi.append(EasyBroncho(generationNumber = numberOfGen, paramsFromJson = True))
    
    def setModelParams(self, t, Vg, Vc0):
        self.t = t
        self.Vg = Vg
        self.Ig = self.Vg / self.bronchi[0].resistance
        self.Vc0 = [Vc0, Vc0]
        return True
    
    def model(self, y, t):
        # Creating the iterators
        iterVg = np.nditer(self.Vg)
        iterIg = np.nditer(self.Ig)

        # Using the iterators. This is for how it works the function odeint of scipy.integrate
        Vg = next(iterVg)
        Ig = next(iterIg)

        # Definition of the model
        dVc1dt = 1/self.bronchi[1].getTau()*(Vg - self.bronchi[0].resistance*Ig) - 1/self.bronchi[1].getTau()*y
        dVc2dt = 1/self.bronchi[2].getTau()*(Vg - self.bronchi[0].resistance*Ig) - 1/self.bronchi[2].getTau()*y
        dVcdt = [dVc1dt[0], dVc2dt[0]]
        return dVcdt

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

if __name__ == "__main__":
    lung = EasyLung(2)
    lung.bronchi[1].compliance = 0.8
    lung.bronchi[2].compliance = 0.7
    t0 = 0
    tf = 20000
    t = np.linspace(t0, tf)
    Vg0 = 200
    f = 100
    Vg = Vg0*np.sin(2*pi*f*t)
    Vc0 = 5
    lung.setModelParams(t, Vg, Vc0)
    dVcdt = lung.solveModel()

    # plot results
    
    plt.plot(t,dVcdt[:,0],'r-',label='dVc1')
    plt.plot(t,dVcdt[:,1],'b-',label='dVc2')
    plt.ylabel('dVc/dt',)
    plt.xlabel('time [s]')

    plt.legend(loc='best')
    plt.grid(linestyle='dashed', linewidth=0.5)
    plt.show()