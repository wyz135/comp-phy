from pylab import *
from math import *

x = arange(0.4,1.2,0.01)
req = 0.74166
V0 = 4.747
a = 0.535 # Optimal value of a
y = 4*V0*((x/a)-1)*((x/a)-2)

# Calculated energy levels
e = [-4.485,-3.960,-3.436,-2.911,-2.387,-1.862,-1.338,-0.813,-0.288,0.236,
    0.761,1.285,1.810,2.334,2.859]

# Plots all energy levels
n = 0
while n < 15:
    xmin = (a/2)*(3-sqrt(1+e[n]/V0))
    xmax = (a/2)*(3+sqrt(1+e[n]/V0))
    xn = arange(xmin,xmax,0.001)
    yn = e[n]*xn/xn
    plot(xn,yn,'-r')
    n = n + 1

xlabel('Molecular distance (A)')
ylabel('Energy (eV)')
text(0.95,-4, 'Energy levels', color='red')
text(1.15,5, 'Potential', color='green')
legend(loc='upper left')
plot(x, y, '-g')
grid(True)
show()