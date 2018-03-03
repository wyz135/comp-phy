import numpy as np
import matplotlib.pyplot as plt

err = 1E-4
muplot = []; xstarplot = []

def f(x,mu):
    return mu*x*(1-x)

n = 200
N = range(0,n)

for mu in np.arange(1.,4.0,0.001):
    x = 0.75 # Initial x value
    xi = []
    
    for i in N: # perform 200 iterations to wash out transient behavior
        x = f(x,mu)

    for i in N: # obtain steady state values
        x = f(x,mu)
        xi.append(x)

    xstar = [xi[0]]
    

    for x in xi: # Obtain attractor values by isolating unique numbers in sequence
        for xs in xstar:
            if abs(xs-x)<err:
                i = True
                break
            else:
                i = False
    
        if i == False:
            xstar.append(x)
        else:
            pass
    
    for xs in xstar:
        muplot.append(mu)
        xstarplot.append(xs)
            
   
plt.plot(muplot,xstarplot,'b,')
plt.ylim([0,1])
plt.xlim([0,4])
plt.show()