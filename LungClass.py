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
    
    def setModelParams(self, t, Ic0):
        self.t = t
        self.i = 0
        self.Ic0 = [Ic0, Ic0]
        return True
    
    def model_dIdt(self, y0, t, f):
        # Definition of the model
        self.i+=1
        print(self.i)
        Ig = 5*np.sin(2*pi*10*t)
        dIgdt = 5*np.cos(2*pi*10*t)
        C1 = self.bronchi[1].compliance
        C2 = self.bronchi[2].compliance
        R1 = self.bronchi[1].resistance
        R2 = self.bronchi[2].resistance
        dI1dt = ((Ig-y0[0])/C2 - (Ig-y0[1])/C1 + dIgdt*R2)/(R1+R2)
        dI2dt = ((Ig-y0[1])/C1 - (Ig-y0[0])/C2 + dIgdt*R1)/(R1+R2)
        dIdt = [dI1dt, dI2dt]
        #print(y0[0], y0[1], t)
        return dIdt

    def solveModel(self):
        # Creating the iterators
        f = 10
        sol = odeint(self.model_dIdt, self.Ic0, self.t, args=(f, ))
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
    t = np.linspace(0, 1000, 1000)
    lung.setModelParams(t, 5)
    dIdt = lung.solveModel()

    # plot results
    
    plt.plot(t,dIdt[:,0],'r-',label=r'$\frac{dVc_1}{dt}$')
    plt.plot(t,dIdt[:,1],'b-',label=r'$\frac{dVc_2}{dt}$')
    plt.ylabel("dVc/dt")
    plt.xlabel('time [s]')

    plt.legend(loc='best')
    plt.grid(linestyle='dashed', linewidth=0.5)
    plt.show()