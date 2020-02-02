from scipy import signal
from LungClass import *
from UsefulFunctions import solveAndPlot

f = 0.25 # 15 breaths per minute --> 1 breath every 4 seconds --> 0.25 Hz
w = 2*pi*f
t = np.arange(0, 30, 0.01) # start, stop, step: be sure step is ok with frequency
extForce = 5*np.exp(1j*w*t)
R = (843, 1.211e5, 3e5)
C = (1e-5, 2.235e-5)
lung = EasyLung(R, C)

solveAndPlot(f, t, extForce, lung, lung.linearModelFromQ).show()
solveAndPlot(f, t, extForce, lung, lung.notLinearModelFromQ).show()


