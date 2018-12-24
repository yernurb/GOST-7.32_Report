import matplotlib
import matplotlib.pyplot as plt
import numpy as np 


plt.style.use('ggplot')

data = np.loadtxt("energy.dat")
t = data[:, [0]].flatten()
x1 = data[:, [1]].flatten()
x2 = data[:, [2]].flatten()
x3 = data[:, [3]].flatten()
x4 = data[:, [4]].flatten()
x5 = data[:, [5]].flatten()

fig, ax = plt.subplots()
ax.plot(t, x1, t, x2, t, x3, t, x4, t, x5)
ax.set(xlabel='time', ylabel='energy', title='Energy evolution')
plt.savefig('energy_evo.png')
plt.show()
