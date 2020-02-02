from scipy import signal
from LungClass import *
from UsefulFunctions import solveAndPlot

f = 0.25 # 15 breaths per minute --> 1 breath every 4 seconds --> 0.25 Hz
t = np.arange(0, 30, 0.01) # start, stop, step: be sure step is ok with frequency
R = (843, 4.5e5, 4.5e5)
C = (1e-5, 1e-5)
lung = EasyLung(R, C)
# forzante = 5*(np.cos(2*pi*f*t) + 1j*np.sin(2*pi*f*t)) # Voltage
forzante = 5*np.exp(1j*2*pi*f*t) # complex sin # Voltage
i0, i1, i2 = lung.LaplaceSolution(f, forzante)
plt.subplot(211)
plt.plot(t, forzante)
plt.xlabel("time (s)")
plt.ylabel(r"V$_g$(t)")
plt.grid(linestyle='dashed', linewidth=0.5)
plt.subplot(212)
plt.plot(t, np.fft.ifft(i0)*1e6)
plt.plot(t, np.fft.ifft(i1)*1e6)
plt.plot(t, np.fft.ifft(i2)*1e6)
plt.plot(t, (np.fft.ifft(i1) + np.fft.ifft(i2))*1e6, '--')
plt.xlabel("time (s)")
plt.ylabel(r"i(t) ($\mu$A)")
plt.grid(linestyle='dashed', linewidth=0.5)
plt.show()
lung.BodeV0(f, forzante).show()
# solveAndPlot(f, t, lung, lung.linearModel).show()
# solveAndPlot(f, t, lung, lung.notLinearResistanceModel).show()




