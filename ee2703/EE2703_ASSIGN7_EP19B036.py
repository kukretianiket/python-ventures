"""
            EE2703 Assignment 7 by Aniket Kukreti EP19B036
"""
from sympy import *
import pylab as p
import scipy.signal as sp
init_session 

def lowpass(R1,R2,C1,C2,G,Vi):  #Low pass transfer function
    s=symbols('s')
    A=Matrix([[0,0,1,-1/G],
              [-1/(1+s*R2*C2),1,0,0],
              [0,-G,G,1],
              [-1/R1-1/R2-s*C1,1/R2,0,s*C1]])
    b=Matrix([0,0,0,-Vi/R1])
    V=A.inv()*b
    return (A,b,V)

def highpass(R1,R2,C1,C2,G,Vi):  #High pass transfer function
    s=symbols('s')
    A=Matrix([[0,0,G,-1],
              [(s*R2*C2)/(1+s*R2*C2),-1,0,0],
              [0,-G,G,1],
              [s*C1+s*C2+1/R1,-s*C2,0,-1/R1]])
    b=Matrix([0,0,0,s*C1*Vi])
    V=A.inv()*b
    return (A,b,V)

def sympyToTF(V):
    """
    Converts from sympy representation of the TF to scipy.signal usable form
    """
    s=symbols('s')
    num, den=fraction(V.simplify())
    n=[float(i) for i in Poly(num,s).all_coeffs()]
    d=[float(i) for i in Poly(den,s).all_coeffs()]
    return sp.lti(p.poly1d(n), p.poly1d(d))

#Unit step response of the low pass filter (Question 1)
s=symbols('s')
t1=p.linspace(0, 0.01, 100001)
A1, b1, V1=lowpass(10000, 10000, 1e-9, 1e-9, 1.586, 1/s)
X1=sympyToTF(V1[3])
t, v=sp.impulse(X1, None, t1)

p.figure(1, (9.6,7.2))
p.plot(t, v)
p.xlabel(r'$t$',size=16)
p.ylabel(r'$V_{o}(t)$',size=16)
p.title(r'Unit step response of the low-pass filter', size=20)
p.grid()
p.show()

#Response of low-pass filter with input vi(t) (Question 2)
vi=p.sin(2000*p.pi*t1)+p.cos(2e6*p.pi*t1)
A2, b2, V2=lowpass(10000, 10000, 1e-9, 1e-9, 1.586, 1)
H2=sympyToTF(V2[3])
t, y, _=sp.lsim(H2, vi, t1)

p.figure(2, (9.6,7.2))
p.plot(t, y)
p.xlabel(r'$t$',size=16)
p.ylabel(r'$V_{o}(t)$',size=16)
p.title(r'Response of the high-pass filter for $v_{i}(t)$', size=20)
p.grid()
p.show()

#Transfer function of the given high-pass filter (Question 3)
w=p.logspace(0,8,801)
ss=1j*w
s=symbols('s')
A3, b3, V3=highpass(10000, 10000, 1e-9, 1e-9, 1.586, 1)
V_o=V3[3]
hf=lambdify(s, V_o, 'numpy')
v3=hf(ss)
H3=sympyToTF(V_o)   #Transfer function of high-pass filter

p.figure(3, (9.6,7.2))
p.loglog(w, abs(v3), lw=2)
p.xlabel(r'$\omega$', size=16)
p.ylabel(r'Magnitude    $(dB)$', size=16)
p.title(r'Bode magnitude plot of highpass filter', size=18)
p.grid()
p.show()

#Response of damped sinusoid (Question 4)
t4=p.linspace(0, 1e-3, 100001)
vi_l=p.exp(-100*t4)*p.cos(2e3*p.pi*t4)    #low frequency input
t, yl, _=sp.lsim(H3, vi_l, t4)

p.figure(4, (9.6,7.2))
p.plot(t, yl)
p.xlabel(r'$t$',size=16)
p.ylabel(r'$V_{o}(t)$',size=16)
p.title(r'Response of the circuit for $v_{i}(t)=e^{-100t}\cos(2000\pi t)$', size=18)
p.grid()
p.show()

vi_h=p.exp(-100*t1)*p.cos(2e6*p.pi*t1)    #high frequency input
t, yh, _=sp.lsim(H3, vi_h, t1)

p.figure(5, (9.6,7.2))
p.plot(t, yh)
p.xlabel(r'$t$',size=16)
p.ylabel(r'$V_{o}(t)$',size=16)
p.title(r'Response of the circuit for $v_{i}(t)=e^{-100t}\cos(2\times10^6\pi t)$', size=18)
p.grid()
p.show()

#Unit step response of the circuit (Question 5)
A5, b5, V5=highpass(10000, 10000, 1e-9, 1e-9, 1.586, 1/s)
Y5=sympyToTF(V5[3])
t, v5=sp.impulse(Y5, None, t1)

p.figure(6, (9.6,7.2))
p.plot(t, v5)
p.xlabel(r'$t$',size=16)
p.ylabel(r'$V_{o}(t)$',size=16)
p.title(r'Unit step response of the high-pass filter', size=20)
p.grid()
p.show()



