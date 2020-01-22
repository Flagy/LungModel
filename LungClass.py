from math import sqrt, pi
# from BronchoClass import EasyBroncho
from scipy.integrate import odeint
from scipy import signal
import numpy as np
import matplotlib.pyplot as plt

class EasyLung(object):

    def __init__(self, R, C):
        (self.R1, self.R2) = R
        (self.C1, self.C2) = C
    
    def setModelParams(self, t, inputSignal, initConds):
        self.t = t
        self.inputSignal = inputSignal
        self.derivInputSignal = np.gradient(self.inputSignal, self.t)
        self.initConds = initConds # Initial conditions of the charges

    def model(self, z, t, s, ds):
        Req = self.R1 + self.R2
        Ceq = 1/self.C1 + 1/self.C2
        ds1dt = -Ceq/Req*z[0] + (self.R2/Req)*ds + s/(Req*self.C2)
        ds2dt = -Ceq/Req*z[1] + (self.R1/Req)*ds + s/(Req*self.C1)
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

if __name__ == "__main__":

    R = (5e-5, 2e-5)
    C = (1e-5, 1000e-5)
    lung = EasyLung(R, C)
    f = 0.25 # 15 respiri al minuto --> 1 respiro ogni 4 secondi --> 0.25 Hz
    t = np.arange(0, 30, 0.01) # start, stop, step: be sure step is ok with frequency
    
    ########## Volume partitions ##########
    Qg = 5*np.sin(2*pi*f*t)
    lung.setModelParams(t, Qg, initConds=(0, 0))
    Q = lung.SolveModel(lung.model)
    plt.subplot(211)
    plt.plot(t, Qg, 'k', label=r'$Q_0(t)$')
    plt.plot(t,Q[0],'r-',label=r'$Q_1(t)$')
    plt.plot(t,Q[1],'b-',label=r'$Q_2(t)$')
    plt.plot(t,Q[0] + Q[1],'g--',label=r'$Q_1 + Q_2$')
    plt.ylabel(r"$Q(t)$", rotation=0)
    plt.xlabel('time [s]')
    plt.legend(loc='best')
    plt.grid(linestyle='dashed', linewidth=0.5)
    plt.title(r"Volume partition")
    
    ########## Flow partitions ###########
    Ig = 10*pi*f*np.cos(2*pi*f*t)
    derivInput = np.gradient(Ig, t)
    initConds = tuple(Ig[0] * np.array((R[1]/(R[0] + R[1]), R[0]/(R[0] + R[1]))))
    lung.setModelParams(t, Ig, initConds)
    I = lung.SolveModel(lung.model)

    plt.subplot(212)
    plt.plot(t, Ig, 'k', label=r'$I_0(t)$')
    plt.plot(t,I[0],'r-',label=r'$I_1(t)$')
    plt.plot(t,I[1],'b-',label=r'$I_2(t)$')
    plt.plot(t,I[0] + I[1],'g--',label=r'$I_1 + I_2$')
    plt.ylabel(r"$I(t)$", rotation=0)
    plt.xlabel('time [s]')
    plt.legend(loc='best')
    plt.grid(linestyle='dashed', linewidth=0.5)
    plt.title(r"Flow partition")
    plt.show()
