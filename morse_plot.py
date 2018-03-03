from pylab import *
from math import *

x = arange(0.5,2.5,0.01)
e = 2.7182818
req = 0.74166
V0 = 4.747
a = 0.256 # Optimal value of a
y = V0*(((1-e**((req-x)/a))**2)-1)

# Calculated energy levels
e = [-4.447,-3.960,-3.476,-3.023,-2.601,-2.211,-1.853,-1.526, -1.231,
    -0.968, -0.736, -0.536, -0.368, -0.126, -0.053]

# Plots all the energy levels.
n = 0
while n < 15:
    xmin = req - a*log(1+sqrt(1+e[n]/V0))
    xmax = req - a*log(1-sqrt(1+e[n]/V0))
    xn = arange(xmin,xmax,0.001)
    yn = e[n]*xn/xn
    plot(xn,yn,'-r')
    n = n + 1

xlabel('Molecular distance (A)')
ylabel('Energy (eV)')
text(0.95,-4, 'Energy levels', color='red')
text(0.55,5, 'Potential', color='green')
legend(loc='upper left')
plot(x, y, '-g')
grid(True)
show()