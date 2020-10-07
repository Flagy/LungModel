#!/usr/bin/env python
# coding: utf-8

# # Linear two-compartment models

# ## 7.1 Passive expiration
# 
# One source of experimental evidence for the presence of more than one compartment in the lung comes from *passive expiration*. This is the *homogeneous* equation that describes how the model behaves during a relaxed expiration when pressure at the mouth is zero.
# \begin{equation}
# \frac{dV(t)}{dt} = - V(t)\frac{1}{RC}
# \end{equation}

# In[1]:


import sympy as sym

sym.init_printing()
t, R, C, C1, V0 = sym.symbols('t R C C1 V0')
V = sym.Function('V')(t)
dVdt = V.diff(t)
expr = sym.Eq(dVdt, -V/(R*C))
sol = sym.dsolve(expr)
sol = sol.subs({C1: V0})
sol


# In[2]:


import numpy as np
from matplotlib import pyplot as plt
from sympy import lambdify

lam = lambdify([t, R, C, V0], sol.rhs, modules=['numpy'])
np_t = np.linspace(0, 2, 100)
np_R = 3 #cmH2O/L/s
np_C = 0.1 #L
np_V0 = 2.1
np_sol = lam(np_t, np_R, np_C, np_V0)

fig, ax = plt.subplots(1, 2, figsize=(15,5))
ax[0].plot(np_t, np_sol)
ax[0].set(xlabel='Time (s)', ylabel='Volume (L)')
ax[0].grid()
ax[1].plot(np_t, np.gradient(np_sol, np_t))
ax[1].set(xlabel='Time (s)', ylabel='Flow (L/s)')
ax[1].grid()
tau = np_R*np_C
# plt.suptitle(r'Solution mono-compartment linear model R = %d $\frac{cmH_2O}{L/s}$, C = %.2f, $\tau$ = %.2f') % np_R, np_C, tau
plt.suptitle(f'Solution mono-compartment linear model R = {np_R} C = {np_C}, tau = {tau:.1f}', fontsize=15)
plt.show()


# Based on initial condition we have $V(t) = V_0 = e^{-t/\tau}$ where $\tau = RC$. The process of fitting this equation to data is not as easy as for the single compartment model because the multiple linear regression technique can't be used due to the non linearity of the time dependent $V(t)$ to the time constant $\tau$.\
# In an experiment a dog was anesthetized, paralyzed and tracheostomized. Its lungs were filled with air through a syringe and the deflation behavior of the lung was studied. The monocompartment model has been trained on these data and the results were poor. This indicates that **the mono-compartment model has not sufficient degrees of freedom** to account to all the aspects of the data.\
# So what about 2 exponential curves? $$V(t) = V_1e^{-t/\tau_1} + V_2e^{-t/\tau_2}$$

# ## The parallel model
# ![ParallelModel.PNG](attachment:ParallelModel.PNG)
# The referral equation for this model is a *second-order linear differential equation*, knowing that $E = 1/C$
# $$(R_1 + R_2)\dot{P}(t) + (E_1 + E_2)P(t) = \Big[R_1R_2 +R_c(R_1 + R_2)\Big]\ddot{V}(t) + \Big[(R_2 + R_c)E_1 + (R_1 + R_c)E_2\Big]\dot{V}(t) + E_1E_2V(t)$$
# 
# Here, the *homogeneous* equation is:
# $$ \Big[R_1R_2 + (R_1 +R_2)R_c\Big]\ddot{V}(t) + \Big[E_1R_2 + E_2R_1 + (E_1 + E_2)R_c\Big]\dot{V}(t) + E_1E_2V(t) = 0$$

# In[3]:


R1, R2, Rc, E1, E2, C1, C2, A1, A2 = sym.symbols('R1 R2 Rc E1 E2 C1 C2 A1 A2')
# Relative values used for solving numerically the equations: PARAMETERS!
n_Rc = 0.25
n_R1 = n_R2 = 5.5
n_E1 = n_E2 = 50
n_P0 = 0.75
paralSol = sym.dsolve(V.diff(t, t)*(R1*R2 + (R1+R2)*Rc) + V.diff(t)*(E1*R2 + E2*R1 + (E1+E2)*Rc) + E1*E2*V)
paralSol = paralSol.subs({C1: A1, C2: A2})
paralSol


# In[4]:


5.5*5.5/11


# In[5]:


exp1 = paralSol.subs({A2: 0})
exp1


# In[6]:


exp2 = paralSol.subs({A1: 0})
exp2


# The initial volume $V(0)$ is the sum of the compartments' volumes, each determined by the inflation pressure **$P_0$** applied just prior to expiration. We have that:
# 
# $$V(0)=P_0\Big(\frac{1}{E_1}+\frac{1}{E_2}\Big) = A_1 + A_2$$
# 
# Also the initial flow is the sum of the 2 compartment's flows, each again determined by $P_0$ and the downstram resistances.
# 
# $$\frac{dV_0}{dt}=P_0\Big(\frac{1}{R_1+R_c}+\frac{1}{R_2+R_c}\Big)=-\frac{A_1}{\tau_1}-\frac{A_2}{\tau_2}$$
# 
# From here we obtain the values of $A_1$ and $A_2$ in function of $P_0, \tau_1, \tau_2, R_1, R_2, E_1, E_2$

# In[7]:


# Getting the time constants, we need them to calculate A1 and A2
tau1 = -sym.ln(exp1.rhs/A1)/t
n_tau1 = tau1.subs({t: 1, Rc:n_Rc, R1:n_R1, E1:n_E1, R2:n_R2, E2:n_E2}).evalf()
print(f'Tau 1 = {n_tau1:.3f}')
tau2 = -sym.ln(exp2.rhs/A2)/t
n_tau2 = tau2.subs({t: 1, Rc:n_Rc, R1:n_R1, E1:n_E1, R2:n_R2, E2:n_E2}).evalf()
print(f'Tau 2 = {n_tau2:.3f}')


# In[8]:


# Solving the initial conditions of A1 and A2
A1, A2, Tau1, Tau2, P0 = sym.symbols('A1 A2 tau1 tau2 P0')
A1 = sym.Eq(A1, Tau1/(Tau1 + Tau2)*(P0*(1/E1 + 1/E2) + Tau2*P0*(1/(R1+Rc)+1/(R2+Rc))))
n_A1 = A1.subs({Tau1: n_tau1, Tau2: n_tau2, P0: n_P0, Rc:n_Rc, R1:n_R1, E1:n_E1, R2:n_R2, E2:n_E2}).evalf()
print(f'From the following expression we get A1 = {n_A1.rhs:.3f}')
A1


# In[9]:


A2 = sym.Eq(A2, Tau2/(Tau1 + Tau2)*(P0*(1/E1 + 1/E2) + Tau1*P0*(1/(R1+Rc)+1/(R2+Rc))))
n_A2 = A2.subs({Tau1: n_tau1, Tau2: n_tau2, P0: n_P0, Rc:n_Rc, R1:n_R1, E1:n_E1, R2:n_R2, E2:n_E2}).evalf()
print(f'From the following expression we get A2 = {n_A2.rhs:.3f}')
A2


# Finally we can solve and visualize the 2 exponential curves caracterized by the 2 time constants $\tau_1$ and $\tau_2$.

# In[10]:


A1 = sym.symbols('A1')
lamExp1 = sym.lambdify([t, Rc, A1, R1, E1, R2, E2], exp1.rhs, 'numpy')
solExp1 = lamExp1(np_t, n_Rc, n_A1.rhs, n_R1, n_E1, n_R2, n_E2)
A2 = sym.symbols('A2')
lamExp2 = sym.lambdify([t, Rc, A2, R1, E1, R2, E2], exp2.rhs, 'numpy')
solExp2 = lamExp1(np_t, n_Rc, n_A2.rhs, n_R1, n_E1, n_R2, n_E2)
plt.plot(np_t, solExp1, '--', label=r'$\tau_1$')
plt.plot(np_t, solExp2, '--', label=r'$\tau_2$')
plt.plot(np_t, solExp1 + solExp2, label=r'$\tau_1 + \tau_2$')
plt.grid()
plt.xlabel("Time (s)")
plt.ylabel("Volume (L)")
plt.title(rf'Bicompartment linear model, homogenous solution, $P_0$={n_P0}cmH$_2$O')
plt.legend()
plt.show()


# ## Frequency Analysis
# The second order differential equation becomes in the Laplace domain:
# 
# $$sP(s)(R_1+R_2)+P(s)(E_1+E_2)=s^2V(s)\Big[R_1R_2 +R_c(R_1 + R_2)\Big] + sV(s)\Big[(R_2 + R_c)E_1 + (R_1 + R_c)E_2\Big] + E_1E_2V(s)$$
# 
# meaning that the transfre function $H(s)$ is:
# 
# $$H(s) = \frac{V(s)}{P(s)} = \frac{s(R_1+R_2)+(E_1+E_2)}{s^2\Big[R_1R_2 +R_c(R_1 + R_2)\Big] + s\Big[(R_2 + R_c)E_1 + (R_1 + R_c)E_2\Big] + E_1E_2}$$

# In[11]:


from scipy import signal
s = signal.lti([n_R1+n_R2, n_E1+n_E2], [n_R1*n_R2 + n_Rc*(n_R1+n_R2), (n_R2+n_Rc)*n_E1 + (n_R1+n_Rc)*n_E2, n_E1*n_E2])
w, mag, phase = s.bode()

fig, ax = plt.subplots(2, 1, figsize=(18,9))
ax[0].semilogx(w, mag)
ax[0].set(ylabel=r"$|H(\omega)|$")
ax[0].grid(linestyle="--")
plt.setp(ax[0], xticklabels=[])
ax[1].semilogx(w, phase)
ax[1].set(xlabel=r"$\omega$  $\frac{rad}{s}$", ylabel=r"$\phi(\omega)$  $rad$")
ax[1].grid(linestyle="--")
fig.suptitle(r"Transfer Function $H(\omega)$", fontsize=16)
pass


# ## Impulse response and step response
# Let's see how the system respond to an impulse stimulus.

# In[12]:


t, y = s.impulse()
fig, ax = plt.subplots(1, 2, figsize=(18,9))
ax[0].plot(t, y)
ax[0].set(xlabel=r"time [s]", ylabel="Amplitude", title="Impulse Response")
ax[0].grid(linestyle="--")

t, y = s.step()
ax[1].plot(t, y)
ax[1].set(xlabel=r"time [s]", ylabel=r"Amplitude", title="Step Response")
ax[1].grid(linestyle="--")


# ## Simulation of tidal breathing

# In[33]:


from scipy import signal
from scipy.signal.ltisys import TransferFunctionContinuous
import numpy as np
import matplotlib.pyplot as plt

def p(t, p0, maxP):
    t = t % 4
    if t < 0.6:
        pp = -maxP*np.sin(np.pi*t/1.2)
    elif t < 0.8:
        pp = -maxP
    elif t < 1:
        pp = -maxP*np.sin(np.pi*(t - 0.6)/0.4)
    elif t < 4:
        pp = 0
    else:
        raise ValueError("Something wrong in t")
    return pp + p0

tspan = np.linspace(0, 20, 1000)
PEEP = 0 # 5 cmH2O
parameters = {'Rc':0.25, 'R1':5.5, 'R2':5.5, 'E1':20, 'E2':20}
topPressure = 10

class ParallelModel(TransferFunctionContinuous):
    def __init__(self, **params):
        self.setParams(**params)
    
    def setParams(self, **params):
        self.n_Rc = params['Rc'] # Normally 0.25
        self.n_R1 = params['R1'] # Normally = 5.5
        self.n_R2 = params['R2'] # Normally = 5.5
        self.n_E1 = params['E1'] # Normally 20
        self.n_E2 = params['E2'] # Normally 20
        super().__init__([self.n_R1 + self.n_R2, self.n_E1 + self.n_E2], 
                         [self.n_R1 * self.n_R2 + self.n_Rc * (self.n_R1 + self.n_R2),
                         (self.n_R2 + self.n_Rc) * self.n_E1 + (self.n_R1 + self.n_Rc)*self.n_E2,
                          self.n_E1 * self.n_E2])
s = ParallelModel(**parameters)
w, mag, phase = s.bode()

fig, ax = plt.subplots(2, 1, figsize=(18,9))
ax[0].semilogx(w, mag)
ax[0].set(ylabel=r"$|H(\omega)|$")
ax[0].grid(linestyle="--")
plt.setp(ax[0], xticklabels=[])
ax[1].semilogx(w, phase)
ax[1].set(xlabel=r"$\omega$  $\frac{rad}{s}$", ylabel=r"$\phi(\omega)$  $rad$")
ax[1].grid(linestyle="--")
fig.suptitle(r"Transfer Function $H(\omega)$", fontsize=16)

breathingPressure = np.vectorize(p)(tspan, PEEP, topPressure)
tout, tidalVolume, yout = s.output(breathingPressure, tspan)
plt.figure(figsize=(16, 8))
plt.plot(tspan, breathingPressure, 'r', linewidth=1, label="input (pressure [cmH2O])")
plt.plot(tout, tidalVolume, 'k', linewidth=1.5, label="output (volume [L])")
plt.plot(tout, np.gradient(tidalVolume), 'b', linestyle='--', label="Flow obtained as gradient of volume")
plt.xlabel("Time [s]")
plt.title("System response to tidal breathing bressure")
plt.legend(loc="lower right")
plt.grid(linestyle='--')
plt.show()
print(f'Min tidal Volume: {min(tidalVolume):.4f}')
print(f'Max tidal Volume: {max(tidalVolume):.4f}')
print(f'Min tidal Flow: {min(np.gradient(tidalVolume)):.4f}')
print(f'Max tidal Flow: {max(np.gradient(tidalVolume)):.4f}')


# ## Visualizing Pressure-Volume and Flow-Volume responses

# In[32]:


fig, ax = plt.subplots(1, 2, figsize=(18,9))
ax[0].plot(breathingPressure, tidalVolume)
ax[0].set(xlabel="Pressure [cmH2O]", ylabel="Volume [L]", title="Volume over Pressure, tidal breathing")
ax[0].grid(linestyle='--')
ax[1].plot(breathingPressure, np.gradient(tidalVolume))
ax[1].set(xlabel="Pressure [cmH2O]", ylabel="Flow [L/s]", title="Flow over Pressure, tidal breathing")
ax[1].grid(linestyle='--')
plt.show()


# ## Getting data
# Changing all the parameters and save the response in a dataframe.

# In[60]:


# Let's try to create clusters of data related to each known group. For istance:
# FIBROSIS: Ceq = (50 +/- 10) ml/cmH2O AND Req = (3 +/- 0.5) cm/H2O
# ASTHMA: Ceq = (100 +/- 10) ml/cmH2O AND Req = (5 +/- 0.5) cm/H2O
# HEALTHY: Ceq = (100 +/- 10) ml/cmH2O AND Req = (3 +/- 0.5) cm/H2O
# Each class will have a total of 1000 examples.
# The noise is white gaussian.
# Req = Rc + R1*R2/(R1+R2)
# Ceq = 1/E1 + 1/E2 = C1 + C2
# Eeq = E1*E2/(E1+E2)
### FIBROSIS ###
def generateData(**kwargs):
    np.random.seed(0)
    e1 = kwargs['E1'] + np.random.normal(loc=0, scale=kwargs['stdE'], size=1000)
    e2 = kwargs['E2'] + np.random.normal(loc=0, scale=kwargs['stdE'], size=1000)
    rc = 0.25
    r1 = kwargs['R1'] + np.random.normal(loc=0, scale=kwargs['stdR'], size=1000)
    r2 = kwargs['R2'] + np.random.normal(loc=0, scale=kwargs['stdR'], size=1000)
    return({'E1':e1, 'E2':e2, 'Eeq':e1*e2/(e1+e2), 'Rc':rc, 'R1':r1, 'R2':r2, 'Req': rc+r1*r2/(r1+r2), 'label':kwargs['label']})

Fibrosis = generateData(E1=40, E2=40, stdE=5, R1=5.5, R2=5.5, stdR=0.5, label='Fibrosis')
Asthma = generateData(E1=20, E2=20, stdE=5, R1=9.5, R2=9.5, stdR=0.5, label='Asthma')
Healthy = generateData(E1=20, E2=20, stdE=5, R1=5.5, R2=5.5, stdR=0.5, label='Healthy')

# Create pandas table
import pandas as pd
dfF = pd.DataFrame(Fibrosis)
dfA = pd.DataFrame(Asthma)
dfH = pd.DataFrame(Healthy)
df = pd.concat([dfF,dfA,dfH],ignore_index=True).drop_duplicates().reset_index(drop=True)
df[df.label=='Fibrosis'].head()
df[df.label=='Asthma'].head()
df[df.label=='Healthy'].head()


# In[57]:


# df = pd.DataFrame(newdict)
df[df.label=='Healthy'].head()


# In[16]:


# Now it's time to decide the set of parameters to test.
import pandas as pd
from itertools import permutations
import numpy as np

Rc = np.linspace(0, 0.75, 4)
R1 = np.linspace(1, 10, 10)
R2 = np.linspace(1, 10, 10)
E1 = np.linspace(0, 30, 10)
E2 = np.linspace(0, 30, 10)

# Let's keep the same tidal input and create our dataset:

def dataCreation(Rc, R1, R2, E1, E2):
    # Creation of a dataframe from a list of dicts.
    matrix = []
    for rc in Rc:
        for r1 in R1:
            for r2 in R2:
                for e1 in E1:
                    for e2 in E2:
                        myDict = {'Rc':rc, 'R1':r1, 'R2':r2, 'E1':e1, 'E2':e2}
                        matrix.append(myDict)
    dataFrame = pd.DataFrame(matrix)
    return dataFrame

df = dataCreation(Rc, R1, R2, E1, E2)


# In[61]:


def p(t, p0):
    t = t % 4
    if t < 0.6:
        pp = -10*np.sin(np.pi*t/1.2)
    elif t < 0.8:
        pp = -10
    elif t < 1:
        pp = -10*np.sin(np.pi*(t - 0.6)/0.4)
    elif t < 4:
        pp = 0
    else:
        raise ValueError("Something wrong in t")
    return pp + p0

tspan = np.linspace(0, 20, 1000)
PEEP = 0
breathingPressure = np.vectorize(p)(tspan, PEEP)
output = []
for index, row in df.iterrows():
    s.setParams(**dict(row))
    output.append(s.output(breathingPressure, tspan)[1])
df['output_Volume'] = output

import pickle
outfile = open("Table.pkl",'wb')
pickle.dump(df, outfile)
outfile.close()


# In[52]:


import pickle
outfile = open("Table.pkl",'wb')
pickle.dump(df, outfile)
outfile.close()


# ## Simulation of Forced Expiration
# Let's see an example of real forced expiration:

# In[19]:


'''from scipy.io import loadmat
matfile = loadmat("ForcedExpirationHealthy.mat")
signal = matfile['data'][0]
totLen = matfile['dataend'][0][0]
t = np.linspace(0, int(totLen)/matfile['samplerate'][0][0], int(totLen))

plt.figure(figsize=(16, 8))
plt.plot(t, signal)
plt.grid(linestyle='--')
plt.xlabel('Time [s]')
plt.ylabel('Flow or volume?')
plt.title("Healthy Forced Expiration Maneuvre")
plt.show()'''


# In[20]:


'''from scipy.io import loadmat
matfile = loadmat("../Python Scripts/plotMFVLhealthyyoung.mat")
signal = matfile['data'][0]
totLen = matfile['dataend'][0][0]
t = np.linspace(0, int(totLen)/matfile['samplerate'][0][0], int(totLen))

plt.figure(figsize=(16, 8))
plt.plot(t[50000:1000000], signal[50000:1000000])
plt.grid(linestyle='--')
plt.xlabel('Time [s]')
plt.ylabel('Flow or volume?')
plt.title("Healthy Forced Expiration Maneuvre")
plt.show()'''


# We can actually see that from this subject we can even extract tidal breathing from the first 10 seconds of signal. We can then use this signal for train our model using multiple linear regression to discover the parameters of the model.
