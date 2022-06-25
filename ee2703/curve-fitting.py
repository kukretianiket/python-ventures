"""
        EE2703 Assignment 4 by Aniket Kukreti EP19B036
"""
from matplotlib.pyplot import *
import scipy as sp
from scipy import integrate
import numpy as np
import math as m

#Creating function f1(x)=exp(x)
def f1(x):
    return np.exp(x)

#Creating function f2(x)=cos(cos(x))
def f2(x):
    return np.cos(np.cos(x))

#Plotting Figure 1
x1=np.linspace(-2*m.pi,0,400, endpoint=False)
x2=np.linspace(0,2*m.pi,400, endpoint=False)
x3=np.linspace(2*m.pi,4*m.pi,400, endpoint=False)
x4=np.linspace(-2*m.pi,4*m.pi,1200, endpoint=False)
figure(1, (12,9))
semilogy(x1,f1(x2), color='red')
semilogy(x2,f1(x2), color='red')
semilogy(x3,f1(x2), color='red')
xlabel(r'$x$', size=15)
ylabel(r'$f_1(x)$    $(\log)$', size=15)
title(r'Plot of $f_1(x)$ in $semilog$ scale', size=18)
grid()
show()

#Plotting Figure 2
figure(2, (12,9))
plot(x4,f2(x4), color='green')
xlabel('$x$', size=15)
ylabel('$f_2(x)$', size=15)
title('Plot of $f_2(x)$', size=18)
grid()
show()

#Creating u(x,k) and v(x,k) for f1 and f2
def u1(x,k):
    return f1(x)*np.cos(k*x)
def v1(x,k):
    return f1(x)*np.sin(k*x)
def u2(x,k):
    return f2(x)*np.cos(k*x)
def v2(x,k):
    return f2(x)*np.sin(k*x)

#Creating vector containing coefficients for f1
clist1=[integrate.quad(u1,0,2*m.pi, args=(0))[0]]
for k in range(1,26):
    clist1.append(integrate.quad(u1,0,2*m.pi, args=(k))[0])
    clist1.append(integrate.quad(v1,0,2*m.pi, args=(k))[0])    
co1=np.c_[clist1]
co1/=m.pi
co1[0]/=2

#Creating vector containing coefficients for f2
clist2=[integrate.quad(u2,0,2*m.pi, args=(0))[0]]
for k in range(1,26):
    clist2.append(integrate.quad(u2,0,2*m.pi, args=(k))[0])
    clist2.append(integrate.quad(v2,0,2*m.pi, args=(k))[0])    
co2=np.c_[clist2]
co2/=m.pi
co2[0]/=2

#Plotting the coefficients for f1(x), (i) in semilog
n=np.arange(1,26)
figure(3, (12,9))
semilogy(0,abs(co1)[0], 'o', color='red')
semilogy(n,abs(co1)[1::2], 'o', color='red', label='$a_n$')
semilogy(n,abs(co1)[2::2], 'o', color='blue', label='$b_n$')
xlabel('$n$', size=15)
ylabel('Coefficients    $(\log)$', size=15)
title(r'Plot of coefficients of $e^x$ in $semilog$ scale', size=18)
legend()
grid()
show()
#(ii) in log-log
figure(4, (12,9))
loglog(n,abs(co1)[1::2], 'o', color='red', label='$a_n$')
loglog(n,abs(co1)[2::2], 'o', color='blue', label='$b_n$')
xlabel('$n$    $(\log)$', size=15)
ylabel('Coefficients    $(\log)$', size=15)
title(r'Plot of coefficients of $e^x$ in $\log\log$ scale', size=18)
legend()
grid()
show()

#Plotting the coefficients for f2(x), (i) in semilog
figure(5, (12,9))
semilogy(0,abs(co2)[0], 'o', color='red')
semilogy(n,abs(co2)[1::2], 'o', color='red', label='$a_n$')
semilogy(n,abs(co2)[2::2], 'o', color='blue', label='$b_n$')
xlabel('$n$', size=15)
ylabel('Coefficients    $(\log)$', size=15)
title(r'Plot of coefficients of $\cos(\cos(x))$ in $semilog$ scale', size=18)
legend()
grid()
show()
#(ii) in log-log
figure(6, (12,9))
loglog(n,abs(co2)[1::2], 'o', color='red', label='$a_n$')
loglog(n,abs(co2)[2::2], 'o', color='blue', label='$b_n$')
xlabel('$n$    $(\log)$', size=15)
ylabel('Coefficients    $(\log)$', size=15)
title(r'Plot of coefficients of $\cos(\cos(x))$ in $\log\log$ scale', size=18)
legend()
grid()
show()

#Generating A & b to find the least square estimate
xall=np.linspace(0,2*m.pi,400, endpoint=False)
b1=f1(xall)
A1=np.zeros([400,51])
A1[:,0]=1
for k in range(1,26):
    A1[:,2*k-1]=np.cos(k*xall)
    A1[:,2*k]=np.sin(k*xall)
c1=np.c_[sp.linalg.lstsq(A1,b1)[0]]

b2=f2(xall)
A2=np.zeros([400,51])
A2[:,0]=1
for k in range(1,26):
    A2[:,2*k-1]=np.cos(k*xall)
    A2[:,2*k]=np.sin(k*xall)
c2=np.c_[sp.linalg.lstsq(A2,b2)[0]]

#Plotting the estimated coefficients for f1(x), (i) in semilog
figure(7, (12,9))
semilogy(0,abs(c1)[0], 'o', color='green')
semilogy(n,abs(c1)[1::2], 'o', color='green', label='$a_n$')
semilogy(n,abs(c1)[2::2], 'o', color='orange', label='$b_n$')
xlabel('$n$', size=15)
ylabel('Coefficients    $(\log)$', size=15)
title(r'Plot of least square approx. coefficients of $e^x$ in $semilog$ scale', size=18)
legend()
grid()
show()
#(ii) in log-log
figure(8, (12,9))
loglog(n,abs(c1)[1::2], 'o', color='green', label='$a_n$')
loglog(n,abs(c1)[2::2], 'o', color='orange', label='$b_n$')
xlabel('$n$    $(\log)$', size=15)
ylabel('Coefficients    $(\log)$', size=15)
title(r'Plot of least square approx. coefficients of $e^x$ in $\log\log$ scale', size=18)
legend()
grid()
show()

#Plotting the estimated coefficients for f2(x), (i) in semilog
figure(9, (12,9))
semilogy(0,abs(c2)[0], 'o', color='green')
semilogy(n,abs(c2)[1::2], 'o', color='green', label='$a_n$')
semilogy(n,abs(c2)[2::2], 'o', color='orange', label='$b_n$')
xlabel('$n$', size=15)
ylabel('Coefficients    $(\log)$', size=15)
title(r'Plot of least squares approx. coefficients of $\cos(\cos(x))$ in $semilog$ scale', size=18)
legend()
grid()
show()
#(ii) in log-log
figure(10, (12,9))
loglog(n,abs(c2)[1::2], 'o', color='green', label='$a_n$')
loglog(n,abs(c2)[2::2], 'o', color='orange', label='$b_n$')
xlabel('$n$    $(\log)$', size=15)
ylabel('Coefficients    $(\log)$', size=15)
title(r'Plot of least squares approx. coefficients of $\cos(\cos(x))$ in $\log\log$ scale', size=18)
legend()
grid()
show()

#Calculating the error between least square and quad estimate
cerr1=abs(co1-c1) #For f1
cerr2=abs(co2-c2) #For f2
print(np.amax(cerr1)) #Maximum deviation for f1
print(np.amax(cerr2)) #Maximum deviation for f2

#Calculating the function f from the fourier series expansion
#Using least square
f1_ls=np.dot(A1,c1)
f2_ls=np.dot(A2,c2)
#Using quad
f1_fs=np.dot(A1,co1)
f2_fs=np.dot(A2,co2)

#Plotting all estimates of f1
figure(11, (12,9))
semilogy(xall,f1(xall), color='red', label='$f_1(x)$')
semilogy(xall,f1_fs, color='blue', label='$f_1(x):\ fourier\ series$')
semilogy(xall,f1_ls, 'o', ms=3, color='green', label='$f_1(x):\ least\ square$')
xlabel(r'$x$', size=15)
ylabel(r'$f_1(x)$    $(\log)$', size=15)
title(r'Plot of estimates of $e^x$ in $semilog$ scale', size=18)
legend()
grid()
show()

#Plotting all estimates of f2
figure(12, (12,9))
plot(xall,f2(xall), color='red', label='$f_2(x)$')
plot(xall,f2_fs, color='blue', label='$f_2(x):\ fourier\ series$')
plot(xall,f2_ls, 'o', ms=3, color='green', label='$f_2(x):\ least\ square$')
xlabel(r'$x$', size=15)
ylabel(r'$f_2(x)$', size=15)
title(r'Plot of estimates of $\cos(\cos(x))$', size=18)
legend()
grid()
show()