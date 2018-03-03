# Plots a simple logistic map of the evolution of a population

import matplotlib.pyplot as plt

x = 0.75
mu = 4.*(1-2E-14)
xi = []
n = 200
N = range(0,n+1)

for n in N:
    x = mu*x*(1-x)
    xi.append(x)
    
plt.plot(N,xi)
plt.ylim([0,1])
plt.xlim([0,n])
plt.show()