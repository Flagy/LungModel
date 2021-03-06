from math import sqrt, pi
from scipy.integrate import odeint
from scipy.special import softmax
import numpy as np
import matplotlib.pyplot as plt
import cmath

class EasyLung(object):

    def __init__(self, R, C):
        (self.R1, self.R2) = R
        (self.C1, self.C2) = C
    
    def setModelParams(self, t, inputSignal, initConds):
        self.t = t
        self.inputSignal = inputSignal
        self.derivInputSignal = np.gradient(self.inputSignal, self.t)
        self.initConds = initConds # Initial conditions of the charges

    def linearModel(self, z, t, s, ds):
        Req = self.R1 + self.R2
        Ceq = 1/self.C1 + 1/self.C2
        ds1dt = -Ceq/Req*z[0] + (self.R2/Req)*ds + s/(Req*self.C2)
        ds2dt = -Ceq/Req*z[1] + (self.R1/Req)*ds + s/(Req*self.C1)
        dsdt = [ds1dt, ds2dt]
        return dsdt

    def notLinearResistanceModel(self, z, t, s, ds):
        R1 = self.R1 + np.tanh(self.R1*s**2)
        R2 = self.R2*np.cos(self.R2*s**2)
        Req = R1 + R2
        Ceq = 1/self.C1 + 1/self.C2
        ds1dt = -Ceq/Req*z[0] + (R2/Req)*ds + s/(Req*self.C2)
        ds2dt = -Ceq/Req*z[1] + (R1/Req)*ds + s/(Req*self.C1)
        dsdt = [ds1dt, ds2dt]
        return dsdt

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

    def LaplaceSolution(self, f, externalForce):
        w = 2*pi*f
        s = complex(0, 1)*w
        externalForce = np.fft.fft(externalForce)
        Z0 = 843 # (Pa*s)/m^3 it's the resistance of the trachea
        Z1 = 1/(s*self.C1) + self.R1
        Z2 = 1/(s*self.C2) + self.R2
        Zeq = 1/(Z0 + 1/(1/Z1 + 1/Z2))
        i0 = externalForce * Zeq
        i1 = i0*Z2/(Z1+Z2)
        i2 = i0*Z1/(Z1+Z2)
        return(i0, i1, i2)