from scipy import signal
from LungClass import *

R = (1.211e5, 3e5)
C = (1e-5, 2.235e-5)
lung = EasyLung(R, C)
f = 0.25 # 15 respiri al minuto --> 1 respiro ogni 4 secondi --> 0.25 Hz
f1 = 0.6
t = np.arange(0, 30, 0.01) # start, stop, step: be sure step is ok with frequency
    
########## Volume partitions ##########
Qg = 5*np.sin(2*pi*f*t)
lung.setModelParams(t, Qg, initConds=(0, 0))
Q = lung.SolveModel(lung.notLinearResistanceModel)

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
initConds = tuple(Ig[0] * np.array((R[1]/(R[0] + R[1]), R[0]/(R[0] + R[1]))))
lung.setModelParams(t, Ig, initConds)
I = lung.SolveModel(lung.notLinearResistanceModel)

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

plt.suptitle(r"$R_1: %g, R_2: %g, C_1: %g, C_2: %g$" % (R[0], R[1], C[0], C[1]), fontsize=12)
plt.show()

