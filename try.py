import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from math import pi

# function that returns dz/dt
def model(z,t):
    R1 = 10
    R2 = 100
    C1 = 0.005
    C2 = 0.004
    f = 10
    A0 = 10
    Ig = A0*np.sin(2*pi*f*t)
    dIgdt = A0*np.cos(2*pi*f*t)
    # dI1dt = ((Ig-z0[2])/C2 - (Ig-z0[3])/C1 + dIgdt*R2)/(R1+R2)
    # dI2dt = ((Ig-z0[3])/C1 - (Ig-z0[2])/C2 + dIgdt*R1)/(R1+R2)
    dzdt = [Ig,dIgdt]
    return dzdt

# initial condition
z0 = [0,0]

# time points
t = np.linspace(0,1,500)

# solve ODE
z = odeint(model,z0,t)

# plot results
plt.plot(t,z[:,0],'b-',label=r'$Ig = 10sin(2\pi ft)$')
plt.plot(t,z[:,1],'r--',label=r'$\frac{dI_g}{dt} = 10cos(2\pi ft)$')
# plt.plot(t,z[:,2],'g-',label=r'$\frac{dI_1}{dt}$')
# plt.plot(t,z[:,3],'y--',label=r'$\frac{dI_2}{dt}$')
plt.ylabel('response')
plt.xlabel('time')
plt.legend(loc='best')
plt.show()