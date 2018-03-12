from pylab import *
from cmath import *
from random import *

def y(t): # Signal generator
        return 4*cos(pi*t)+ 2*cos(3*pi*t) + 3*cos(5*pi*t) + 5*(0.5-random())
        #if t < 1.0 or t > 3.0:
        #    return -1. + 0.5*random()
        #else:
        #    return 1. + 0.5*random()

T   = 4.0 # Sampling period
h   = 0.01 # Sampling interval
N = int(T/h)
#sig = [0.1, 1.2, 0.9, 1.1, -0.1]
sig = []
for t in arange(0,T+h,h):
    sig.append(y(t))

#print sig 
plot(arange(0,T+h,h),sig,'g')

def auto_corr(sig): # Evaluates the autocorrelation for a real signal
    
    N = len(sig) # Number of sample points
    Ay = []
    
    for t in arange(0,T+h,h):
        n = int(t/h) # Units of shifting
        A = 0.
        for i in range(1,N):
            A += sig[i]*sig[(i+n)%N]
        A *= h
        Ay.append(A)
    
    return Ay

def DFT(sig):     # Performs a discrete fourier transform
    an = zeros((N+1),complex)
    for n in range(0,N+1):
        for k in range(1,N+1):
            an[n] += exp(-2j*pi*n*k/N)*sig[k]
        an[n] /= sqrt(2*pi)
    return an

Aw = DFT(auto_corr(sig)) # Power spectrum
print Aw
Sw = []

for a in DFT(auto_corr(sig)):
    Sw.append(sqrt(abs(a))/sqrt(2*pi))


def y_reconstruct(an,t): # Reconstructs signal
    y = 0.
    for n in range(1,N):
        if n <= N/2:
            y += an[n]*exp(2j*pi*n*t/T)
        else:
            y += an[n]*exp(2j*pi*(n-N)*t/T)
    return y.real*sqrt(2*pi)/N

yrecon = []
for t in arange(0,T+h,h):
    yrecon.append(y_reconstruct(Aw,t)/(2*pi)) # Why need to divide an extra 2pi?
    
plot(arange(0,T+h,h),yrecon,'r')
show()