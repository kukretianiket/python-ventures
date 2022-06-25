"""
            EE2703 Assignment 6.5 by Aniket Kukreti EP19B036
"""
from pylab import*
import scipy.signal as sp

#Question 1
def F(freq, decay):     #Laplace transform of f(t)
    num=poly1d([1,decay])
    den=polyadd(polymul(num, num), freq**2)
    return num, den

F1=F(1.5, 0.5)
X1=sp.lti(F1[0], polymul(F1[1], poly1d([1, 0, 2.25])))   #X(s)=H(s)F(s)
t1=linspace(0, 50, 1001)
t, x=sp.impulse(X1, None, t1)     #Impulse response of system

figure(1, (9.6,7.2))  
plot(t, x)
xlabel(r'$t$',size=16)
ylabel(r'$x$',size=16)
title(r'Time response of spring with $decay=0.5$',size=20)
grid()
show()

#Question 2
F2=F(1.5, 0.05)
X2=sp.lti(F2[0], polymul(F2[1], poly1d([1, 0, 2.25])))   #X(s)=H(s)F(s)
t, x=sp.impulse(X2, None, t1)     #Impulse response of system

figure(2, (9.6,7.2))  
plot(t, x)
xlabel(r'$t$',size=16)
ylabel(r'$x$',size=16)
title(r'Time response of spring with $decay=0.05$',size=20)
grid()
show()

#Question 3
def f(freq, t):
    return cos(freq*t)*exp(-0.05*t)

fRange=linspace(1.4, 1.6, 5)    
for i, fq in enumerate(fRange, start=1):
    H=sp.lti(1, poly1d([1, 0, 2.25]))
    t, y, _= sp.lsim(H, f(fq,t1), t1)
    
    figure(i+2, (9.6,7.2))
    plot(t, y, 'g')
    xlabel(r'$t$',size=16)
    ylabel(r'$x$',size=16)
    title(r'Time response of spring with frequency=${}$'.format(fq), size=20)
    grid()
    show()

#Question 4
X4=sp.lti([1,0,2],[1,0,3,0])
Y4=sp.lti(2,[1,0,3,0])
t4=linspace(0,20,1001)
tx, x=sp.impulse(X4, None, t4)
ty, y=sp.impulse(Y4, None, t4)

figure(8, (9.6,7.2))
plot(tx, x, 'r', label='$x(t)$')
plot(ty, y, 'b', label='$y(t)$')
xlabel(r'$t$',size=16)
ylabel(r'$x$',size=16)
legend()
title(r'Time evolution of coupled system', size=20)
grid()
show()

#Question 5
R=100   #Resistance of R
L=1e-6  #Inductance of L
C=1e-6  #Capacitance of C
H=sp.lti(1,[L*C,R*C,1])
w, S, phi=H.bode()

figure(9, (9.6,7.2))
semilogx(w, S)
xlabel('$\omega$', size=16)
ylabel('Magnitude   (dB)', size=16)
title('Magnitude response of transfer function', size=20)
grid()
show()

figure(10, (9.6,7.2))
semilogx(w, phi)
xlabel('$\omega$', size=16)
ylabel('Phase   (deg)', size=16)
title('Phase response of transfer function', size=20)
grid()
show()

#Question 6
def v_i(t):
    return cos(1e3*t)-cos(1e6*t)
# Fast response
t61=linspace(0, 3e-5, 301)
t, v_o, _=sp.lsim(H, v_i(t61), t61)

figure(11, (9.6,7.2))
plot(t, v_o)
xlabel(r'$t$',size=16)
ylabel(r'$v_{o}(t)$',size=16)
title(r'Output signal of two port network (Fast response)', size=20)
grid()
show()

# Slow response
t62=linspace(0, 0.01, 100001)
t, v_o, _=sp.lsim(H, v_i(t62), t62)

figure(12, (9.6,7.2))
plot(t, v_o)
xlabel(r'$t$',size=16)
ylabel(r'$v_{o}(t)$',size=16)
title(r'Output signal of two port network (Slow response)', size=20)
grid()
show()