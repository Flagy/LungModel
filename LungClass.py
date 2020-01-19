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
    eta = 1.81e-5 #Pa*s https://www.engineersedge.com/physics/viscosity_of_air_dynamic_and_kinematic_14483.htm
    def __init__(self, numberOfGen):
        if numberOfGen < 0 or numberOfGen > EasyLung.totGenerations:
            raise ValueError
        self.bronchi = []
        for gen in range(numberOfGen + 1):
            self.bronchi.append(EasyBroncho(generationNumber = numberOfGen, paramsFromJson = True))
    
    def setModelParams(self, t, f, initConds):
        self.f = f # 15 respiri al minuto --> 1 respiro ogni 4 secondi --> 0.25 Hz
        self.t = t
        self.Ic0 = initConds # Initial conditions in the 2 currents
        return True
    
    def model_dIdt(self, y0, t):
        # Definition of the model
        Ig = 5*np.sin(2*pi*self.f*t)
        dIgdt = 2*pi*self.f*5*np.cos(2*pi*self.f*t)
        C1 = self.bronchi[1].compliance
        C2 = self.bronchi[2].compliance
        R1 = self.bronchi[1].resistance
        R2 = self.bronchi[2].resistance
        R1 = 100000
        R2 = 900000
        # print(R2*dIgdt)
        # print(R1*dIgdt)
        dI1dt = (R2*dIgdt + Ig/C2 - y0[0]*(1/C1 + 1/C2))/(R1+R2)
        dI2dt = (R1*dIgdt + Ig/C1 - y0[1]*(1/C1 + 1/C2))/(R1+R2)
        print(dI1dt)
        print(dI2dt)
        dIdt = [dI1dt, dI2dt]
        return dIdt

    def solveModel(self):
        # Creating the iterators
        sol = odeint(self.model_dIdt, self.Ic0, self.t)
        return sol

    def getDiameterFromGen(self, numGen):
        """Description"""
        return self.diameter * self.h**numGen

    def getResistanceFromGen(self, numGen):
        """Poiseuille Resistance, Weibel Model"""
        d = self.getDiameterFromGen(numGen)
        return (pi*self.P*((d/2.0)**4))/(8.0*EasyBroncho.eta*self.l)

if __name__ == "__main__":
    lung = EasyLung(2)
    lung.bronchi[1].compliance = 0.05e-5
    lung.bronchi[2].compliance = 0.05e-5
    f = 0.25 # 15 respiri al minuto --> 1 respiro ogni 4 secondi --> 0.25 Hz
    startTime = 0
    stopTime = 30
    incrementTime = 0.01 # Be sure this is okay with the frequency
    t = np.arange(startTime, stopTime, incrementTime)
    lung.setModelParams(t, f, initConds=(0, 0))
    dIdt = lung.solveModel()
    # print(lung.getResistanceFromGen(2))
    # plot results
    
    plt.plot(t, 5*np.sin(2*pi*0.25*t), 'k', label=r'$I_g(t)$')
    plt.plot(t,dIdt[:,0],'r-',label=r'$\frac{dI_1}{dt}$')
    plt.plot(t,dIdt[:,1],'b-',label=r'$\frac{dI_2}{dt}$')
    plt.ylabel(r"$\frac{dI_c}{dt}$", rotation=0)
    plt.xlabel('time [s]')

    plt.legend(loc='best')
    plt.grid(linestyle='dashed', linewidth=0.5)
    plt.show()