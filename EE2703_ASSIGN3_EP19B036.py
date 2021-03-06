"""
        EE2703 Assignment 3 by Aniket Kukreti EP19B036
"""
from matplotlib.pyplot import *
import scipy as sp
import scipy.special as sps
import numpy as np
import sys
import random

#Loading the .dat file
try:
    t=np.loadtxt('fitting.dat', usecols=0)
    data=np.loadtxt('fitting.dat', usecols=np.arange(1,10))
except:
    print('Error while loading file. Please check your .dat file and try again.')
    sys.exit()

#Defining the function g
def g(t,A,B):
    return(A*sps.jv(2,t)+B*t)

sig=np.logspace(-1,-3,9)

#Plotting the noisy data along with true function
figure(1, (9.6,7.2))
for i in range(9):
    plot(t,data.T[i], label='$\u03C3_{{{ss}}}={v}$'.format(ss=(i+1),v=round(sig[i],3)))
plot(t,g(t,1.05,-0.105), color='black', label='$true\ value$',)
xlabel(r'$t$',size=15)
ylabel(r'$f(t)+noise$',size=15)
title(r'Plot of true function along with noisy data',size=18)
legend()
grid()
show()

#Plotting true func with error bars
col1=data.T[0]
figure(2, (9.6,7.2))
errorbar(t[::5], col1[::5], sig[0], fmt='ro', label='$error\ bar$')
plot(t,g(t,1.05,-0.105), color='black', label='$f(t)$')
xlabel(r'$t$',size=15)
title(r'Plot of true function with error bars for '+'$\u03C3=0.100$',size=18)
legend()
grid()
show()

#Creating matrix M
j1=[sps.jv(2,x) for x in t]
M=np.c_[j1,t]

#Checking if M*p=g for random values of A and B in a certain range
r1=random.uniform(0,2)
r2=-random.uniform(0,0.2)
p=np.array([r1,r2])
Mp=np.dot(M,p)
if g(t,r1,r2).all() == Mp.all():
    print('Equal.')
else:
    print('Not equal.')

#Generating the Mean Squared Error matrix
As=np.linspace(0,2,21)
Bs=np.linspace(-0.2,0,21)
eps=np.zeros((21,21))       
for i in range(len(As)):
    for j in range(len(Bs)):
        for k in range(len(t)):
            eps[i][j]+=(col1[k]-g(t[k],As[i],Bs[j]))**2
eps=eps/101

#Creating the contour plot of MS Error vs (A,B)
figure(3, (9.6,7.2))
cp=contour(As,Bs,eps, levels=np.linspace(0,0.4,17))
clabel(cp, inline=True)
plot(1.05,-0.105, marker='o', markersize=10, color='r')
annotate('Exact Location (1.05,-0.105)',xy=(1.05,-0.105))
xlabel(r'$A$',size=15)
ylabel(r'$B$',size=15)
title(r'Contour plot of '+'$\u03B5_{ij}$',size=20)
show()

#Solving for least square approx. and calculating error in A and B
U=[sp.linalg.lstsq(M,data.T[i])[0] for i in range(9)]
Aerr=[abs(U[i][0]-1.05) for i in range(9)]
Berr=[abs(U[i][1]+0.105) for i in range(9)]

#Plotting error in A and B versus noise
figure(4, (9.6,7.2))
plot(sig,Aerr, linestyle='--',  marker='o', color='r', label='$Aerr$')
plot(sig,Berr, linestyle='--',  marker='o', color='g', label='$Berr$')
xlabel('$Noise\ Standard\ Deviation\ (\u03C3_{n})$',size=15)
ylabel('$MS\ Error$',size=15)
title(r'Error in A and B versus noise',size=18)
legend()
grid()
show()

#Creating the above plot but with log scales
figure(5, (9.6,7.2))
loglog(sig,Aerr, linestyle='--',  marker='o', color='r', label='$Aerr$')
loglog(sig,Berr, linestyle='--',  marker='o', color='g', label='$Berr$')
xlabel('$Noise\ Standard\ Deviation\ (\u03C3_{n})\ (log)$',size=16)
ylabel('$MS\ Error\ (log)$',size=15)
title(r'Error in A and B versus noise in log scales',size=18)
legend()
grid()
show()


