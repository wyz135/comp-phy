from cmath import *
from pylab import *

def y(t):
        return cos(20*pi*t)
        #if t < 1.0 or t > 3.0:
        #    return 0.
        #else:
        #    return 1.

T = 4.
h = 0.2
tk = arange(0,T+h,h)
yk = []
for t in tk:
    yk.append(y(t))
print yk

N = int(T/h)
an = zeros((N+1),complex)

# Evaluate fourier coefficients
for n in range(0,N+1):
    for k in range(1,N+1):
        an[n] += exp(-2j*pi*n*k/N)*yk[k]
    an[n] /= sqrt(2*pi)

def y_reconstruct(t):
    y = 0
    for n in range(1,N+1):
        if n <= N/2:
            y += an[n]*exp(2j*pi*n*t/T)
        else:
            y += an[n]*exp(2j*pi*(n-N)*t/T)
    return y.real*sqrt(2*pi)/N


#plot(w,abs(an),'bo')
t = arange(0,T+h/10,h/10)
plot(t,y_reconstruct(t),'b')
plot(tk,yk,'ro')
show()

print an