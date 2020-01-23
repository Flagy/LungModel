from scipy import signal
from LungClass import *
from UsefulFunctions import solveAndPlot

R = (1.211e5, 3e5)
C = (1e-5, 2.235e-5)
lung = EasyLung(R, C)
f = 0.25 # 15 breaths per minute --> 1 breath every 4 seconds --> 0.25 Hz
t = np.arange(0, 30, 0.01) # start, stop, step: be sure step is ok with frequency
solveAndPlot(f, t, lung, lung.linearModel).show()
solveAndPlot(f, t, lung, lung.notLinearResistanceModel).show()




