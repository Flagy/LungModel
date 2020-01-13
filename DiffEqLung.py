import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from math import pi

def model(y0, t, f, A, B, C):
    ft = np.sin(2*pi*f*t)
    gt = np.cos(2*pi*f*t)
    it = 2*t
    di0dt = A*ft + B*gt #+ C*it
    return di0dt

def modelDefinition():
    y0 = 10
    t = np.arange(0, 10, 0.1)
    f = 0.25
    A = 2.3
    B = 3.1
    C = -2.1
    sol = odeint(model, y0, t, args=(f, A, B, C))
    plotSolution(t, sol)
    
def plotSolution(t, sol):
    plt.plot(t, sol)
    plt.xlabel("time")
    plt.ylabel("di/dt")
    plt.title("Result")
    plt.show()

if __name__ == "__main__":
    modelDefinition()


