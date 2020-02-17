import numpy as np
from matplotlib import pyplot as plt

R0 = 843
Rst = 500e3
Rc_clean = np.linspace(-20000, 20000, 100)
r = np.random.RandomState(seed = 1) 
Rc = Rc_clean + r.uniform(-10000,10000,100)

alpha = (1 + Rc)*(R0 + Rst)/(R0 + Rc + Rst)
classes = np.array([-1]*40 + [0]*20 + [1]*40)

plt.scatter(Rc_clean[0:40] , alpha[0:40], color='red', marker = 'o', label = 'opposite')
plt.scatter(Rc_clean[60:-1] , alpha[60:-1], color='blue', marker = 'o', label = 'that one')
plt.scatter(Rc_clean[40:60] , alpha[40:60], color='green', marker = 'x', label = 'healthy')
plt.xlabel(r"$\Omega$")
plt.ylabel(r"$\alpha$")
plt.title(r"$\alpha(R_c)$")
plt.grid()
plt.legend(loc = 'upper left')
plt.show()

ppn = Perceptron(eta = 0.1, random_state=r.seed)
ppn.fit(alpha, classes)

plt.plot(range(1, len(ppn.errors_) + 1), ppn.errors_, marker='o')
plt.xlabel('Epochs')
plt.ylabel('Number of updates')

# plt.savefig('images/02_07.png', dpi=300)
plt.show()