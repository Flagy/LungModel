import numpy as np
from scipy.integrate import odeint
from scipy import signal
from math import pi
from matplotlib import pyplot as plt
"""RC model, solve i(t) in the frequency domain, knowing v(t).
        Complex domain"""
f = 5
w = 2*pi*f
t = np.linspace(0, 1, 1000)
phi = pi/3
Vcomplex = 5*(np.cos(w*t)+1j*np.sin(w*t))
# Vreal = 5*1j*np.sin(w*t)
R = 100
C = 1/w

Vf = np.fft.fft(Vcomplex)
Z = R + 1/(1j*w*C)
If = Vf/Z
It = np.fft.ifft(If)
plt.plot(t, Vcomplex, t, It)
plt.show()

"""Bode"""
H = 1/(1+1j*w*R*C)
sys = signal.lti([1], [R*C, 1])
w, mag, phase = sys.bode()

plt.figure()
plt.semilogx(w, mag)    # Bode magnitude plot
plt.figure()
plt.semilogx(w, phase)  # Bode phase plot
plt.show()