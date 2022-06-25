"""
            EE2703 Assignment 6 by Aniket Kukreti EP19B036
"""
from pylab import*
import sys

#Default parameter values
n=100       # spatial grid size.
M=5         # number of electrons injected per turn.
nk=500      # number of turns to simulate.
u0=5.0        # threshold velocity.
p=0.25      # probability that ionization will occur
Msig=2.0      #standard deviation of the injection distribution
model2=0    #initialising the variable that decides which model is used

#Reading command line arguments
if len(sys.argv) != 1:
    try:
        n=int(sys.argv[1])
        M=int(sys.argv[2])
        nk=int(sys.argv[3])
        u0=float(sys.argv[4])
        p=float(sys.argv[5])
        Msig=float(sys.argv[6])
        model2=int(sys.argv[7])
    except IndexError:
        print("Please input all the parameters in the following order: n, M, nk, u0, p, Msig, model2")
        sys.exit()
    except ValueError:
        print("Please input correct datatypes. (int, int, int, float, float, float, int)")
        sys.exit()

#Defining the vectors and lists
xx=zeros(n*M)
u=zeros(n*M)
dx=zeros(n*M)
I=[]
X=[]
V=[]

#Simulating loops
if model2 == 0:
    for i in range(nk):
        #Updating the vectors
        ii=where(xx>0)      #tuple containing indices of non zero entries
        dx[ii]=u[ii]+0.5
        xx[ii]=xx[ii]+dx[ii]
        u[ii]=u[ii]+1
        #Removing e-s that cross the anode
        aa=where(xx>n)      #tuple containing indices of e- that crossed the anode
        xx[aa]=0
        u[aa]=0
        dx[aa]=0
        #Checking which electrons collide to give emission
        kk=where(u>=u0)     #tuple containing indices of energetic e-
        ll=where(rand(len(kk[0]))<=p)   #tuple containing indices of kk that denote e- that collided
        kl=kk[0][ll]        #list containing indices of e- that collided
        rho=rand(len(kl))   #distribution of excitation positions
        xx[kl]=xx[kl]-dx[kl]*rho    #positions at which emissions take place
        u[kl]=0             #setting e- velocity to zero
        I.extend(xx[kl].tolist())   #storing emission positions
        #Injecting e-
        m=int(randn()*Msig+M)       #number of e- injected
        es=where(xx==0)
        Ni=min(m,len(es[0]))
        xx[es[0][:Ni]]=1
        #Storing positions and velocities
        ii=where(xx>0)
        X.extend(xx[ii].tolist())
        V.extend(u[ii].tolist())
elif model2 == 1:
    for i in range(nk):
        #Injecting e-
        m=int(randn()*Msig+M)
        es=where(xx==0)
        Ni=min(m,len(es[0]))
        xx[es[0][:Ni]]=1
        #Updating the vectors
        ii=where(xx>0) #tuple containing indices of non zero entries
        dx[ii]=u[ii]+0.5
        xx[ii]=xx[ii]+dx[ii]
        u[ii]=u[ii]+1
        #Removing e-s that cross the anode
        aa=where(xx>n) #tuple containing indices of e- that crossed the anode
        xx[aa]=0
        u[aa]=0
        dx[aa]=0
        #Checking which electrons collide to give emission
        kk=where(u>=u0)                  #tuple containing indices of energetic e-
        ll=where(rand(len(kk[0]))<=p)    #tuple containing indices of kk that denote e- that collided
        kl=kk[0][ll]                     #list containing indices of e- that collided
        dt=rand(len(kl))                 #distribution of times at which collision took place
        xx[kl]=(xx[kl]-dx[kl]+
                (u[kl]-1)*dt+0.5*(dt)**2) #positions at which emissions take place
        I.extend(xx[kl].tolist())        #storing the emission positions
        xx[kl]=xx[kl]+0.5*(1-dt)**2      #position of e- after ionization
        u[kl]=1-dt                       #velocity of e- after ionization
        #Storing positions and velocities
        ii=where(xx>0)
        X.extend(xx[ii].tolist())
        V.extend(u[ii].tolist()) 
else:
    print('To use model 1, please input 0 as 7th argument.\nTo use model 2, please input 1 as 7th argument.')
    sys.exit()
    
#Plotting the histograms
#Electron density distribution
figure(1, (9.6,7.2))    
hist(X, bins=linspace(0, n, n+1), ec='black')
title('Electron density   ($u_0=${} and $p=${})'.format(u0, p))
xlabel('$x$')
ylabel('Number of electrons')
xticks(arange(0, n+1, 10)) 
show()
#Light intensity distribution    
figure(2, (9.6,7.2))    
popCount, bins, rects=hist(I, bins=linspace(0, n, n+1), fill=False) 
xticks(arange(0, n+1, 10))
xlabel(r'$x$')
ylabel(r'$I$')
title(r'Population plot of Intensity $I$')
show()
#Electron phase space
figure(3, (9.6,7.2))    
scatter(X, V, c='green',s=12, marker='x')
xlabel(r'Position of electrons')
ylabel(r'Velocity of electrons')
xticks(arange(0, n+1, 10))
title(r'Phase space of the electrons')

#Printing intensity table
xpos=0.5*(bins[0:-1]+bins[1:])
print('Intensity Data:')
print('xpos\tcount')
data=c_[xpos, popCount]
for row in data[:]:
    print(' {} \t {} \n'.format(row[0],row[1]))
    
    

    