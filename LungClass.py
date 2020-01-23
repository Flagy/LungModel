from math import sqrt, pi
from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt

class EasyLung(object):

    def __init__(self, R, C):
        (self.R1, self.R2) = R
        (self.C1, self.C2) = C
    
    def setModelParams(self, t, inputSignal, initConds):
        self.t = t
<<<<<<< HEAD
        self.Ic0 = initConds # Initial conditions in the 2 currents
        return True
    
    def model_dIdt(self, y0, t):
        # Definition of the model
        Ig = 5*np.sin(2*pi*self.f*t)
        dIgdt = 5*2*pi*self.f*np.cos(2*pi*self.f*t)
        C1 = self.bronchi[1].compliance
        C2 = self.bronchi[2].compliance
        R1 = self.bronchi[1].resistance
        R2 = self.bronchi[2].resistance
        R1 = 1e-5
        R2 = 2e-5
        dI1dt = ((Ig-y0[0])/C2 - (Ig-y0[1])/C1 + dIgdt*R2)/(R1+R2)
        dI2dt = ((Ig-y0[1])/C1 - (Ig-y0[0])/C2 + dIgdt*R1)/(R1+R2)
        dIdt = [dI1dt, dI2dt]
        return dIdt

    def solveModel(self):
        # Creating the iterators
        sol = odeint(self.model_dIdt, self.Ic0, self.t)
        return sol
=======
        self.inputSignal = inputSignal
        self.derivInputSignal = np.gradient(self.inputSignal, self.t)
        self.initConds = initConds # Initial conditions of the charges
>>>>>>> f52d8480c7c666253289730bf8c6068b9cd66cf5

    def model(self, z, t, s, ds):
        Req = self.R1 + self.R2
        Ceq = 1/self.C1 + 1/self.C2
        ds1dt = -Ceq/Req*z[0] + (self.R2/Req)*ds + s/(Req*self.C2)
        ds2dt = -Ceq/Req*z[1] + (self.R1/Req)*ds + s/(Req*self.C1)
        dsdt = [ds1dt, ds2dt]
        return dsdt

<<<<<<< HEAD
    def getResistanceFromGen(self, numGen):
        """Poiseuille Resistance, Weibel Model"""
        d = self.getDiameterFromGen(numGen)
        return (pi*self.P*((d/2.0)**4))/(8.0*EasyBroncho.eta*self.l)

if __name__ == "__main__":
    lung = EasyLung(2)
    lung.bronchi[1].compliance = 0.01e-5
    lung.bronchi[2].compliance = 0.01e-5
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
=======
    def SolveModel(self, model):
        z0 = self.initConds
        sol1 = np.zeros_like(self.t)
        sol2 = np.zeros_like(self.t)
        sol1[0] = z0[0]
        sol2[0] = z0[1]
        for i in range(1, len(self.t)):
            # time span necessary for derivative inside odeint
            tspan = [self.t[i-1], self.t[i]]
            sol = odeint(model, z0, tspan, args=(self.inputSignal[i], self.derivInputSignal[i]))
            sol1[i] = sol[1][0]
            sol2[i] = sol[1][1]
            z0 = sol[1]
        return (sol1, sol2)
>>>>>>> f52d8480c7c666253289730bf8c6068b9cd66cf5

    