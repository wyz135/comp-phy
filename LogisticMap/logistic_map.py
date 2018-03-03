# Plots a simple logistic map of the evolution of a population

import matplotlib.pyplot as plt

# Initial population and mu value
x = 0.75
mu = 4.*(1-2E-14)

n = 200
N = range(0,n+1)

# Obtain the population value over n iterations
xi = []
for n in N:
    x = mu*x*(1-x)
    xi.append(x)
 
# Plotting the evolution of population   
plt.plot(N,xi)
plt.ylim([0,1])
plt.ylabel("Population x")
plt.xlim([0,n])
plt.xlabel("Time, n")
plt.show()
