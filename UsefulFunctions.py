import numpy as np
import matplotlib.pyplot as plt
from math import pi

def solveAndPlot(f, t, lung, func):
    ########## Volume partitions ##########
    Qg = 5*np.sin(2*pi*f*t)
    lung.setModelParams(t, Qg, initConds=(0, 0))
    Q = lung.SolveModel(func)

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
    Ig = np.gradient(Qg, t)
    derivInput = np.gradient(Ig, t)
    initConds = tuple(Ig[0] * np.array((lung.R2/(lung.R1 + lung.R2), lung.R1/(lung.R1 + lung.R2))))
    lung.setModelParams(t, Ig, initConds)
    I = lung.SolveModel(func)

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

    plt.suptitle(r"$R_1: %g, R_2: %g, C_1: %g, C_2: %g$" % (lung.R1, lung.R2, lung.C1, lung.C2), fontsize=12)
    return plt