from scipy import signal
from LungClass import *
from UsefulFunctions import solveAndPlot

R = (1.211e5, 3e5)
C = (1e-5, 22.35)
lung = EasyLung(R, C)
f = 0.25 # 15 breaths per minute --> 1 breath every 4 seconds --> 0.25 Hz
t = np.arange(0, 30, 0.01) # start, stop, step: be sure step is ok with frequency
forzante = 5*np.sin(2*pi*f*t + pi/3) # Voltage
i0, i1, i2 = lung.LaplaceSolution(f, forzante)
plt.subplot(211)
plt.plot(t, forzante)
plt.grid(linestyle='dashed', linewidth=0.5)
plt.subplot(212)
plt.plot(t, np.fft.ifft(i0))
plt.plot(t, np.fft.ifft(i1))
plt.plot(t, np.fft.ifft(i2))
plt.grid(linestyle='dashed', linewidth=0.5)
plt.show()
# solveAndPlot(f, t, lung, lung.linearModel).show()
# solveAndPlot(f, t, lung, lung.notLinearResistanceModel).show()




