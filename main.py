from scipy import signal
from LungClass import *
from UsefulFunctions import solveAndPlot

f = 0.25 # 15 breaths per minute --> 1 breath every 4 seconds --> 0.25 Hz
t = np.arange(0, 60, 0.01) # start, stop, step: be sure step is ok with frequency
R = (843, 700e3, 500e3)
C = (100e-6, 1e-6)
lung = EasyLung(R, C)
# forzante = 5*(np.cos(2*pi*f*t) + 1j*np.sin(2*pi*f*t)) # Voltage
phi = pi/2
forzante = 250*np.exp(-1j*2*pi*f*t + 1j*phi) # complex sin # Voltage
# forzante = 500*np.sin(2*pi*f*t + phi)
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




