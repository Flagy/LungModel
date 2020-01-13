import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from math import pi

def model(y0, t, f):
    # Definition of the model
    Ig = 5*np.sin(2*pi*10*t)
    dIgdt = 5*np.cos(2*pi*10*t)
    C1 = 2
    C2 = 3
    R1 = 10
    R2 = 12
    dI1dt = ((Ig-y0[0])/C2 - (Ig-y0[1])/C1 + dIgdt*R2)/(R1+R2)
    dI2dt = ((Ig-y0[1])/C1 - (Ig-y0[0])/C2 + dIgdt*R1)/(R1+R2)
    dIdt = [dI1dt, dI2dt]
    #print(y0[0], y0[1], t)
    return dIdt

def solveModel():
    # Creating the iterators
    f = 10
    Ic0 = [5, 5]
    t = np.arange(0, 10, 0.01)
    sol = odeint(model, Ic0, t, args=(f, ))
    plotSolution(t, sol)
    
def plotSolution(t, sol):
    plt.plot(t, sol)
    plt.xlabel("time")
    plt.ylabel("di/dt")
    plt.title("Result")
    plt.show()

if __name__ == "__main__":
    solveModel()