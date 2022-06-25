"""
        EE2703 Assignment by Aniket Kukreti EP19B036
"""
from pylab import*
import mpl_toolkits.mplot3d.axes3d as p3
import sys
from scipy.linalg import*

#Initialising parameters
Nx=25
Ny=25
radius=8
Niter=1500

#Reading command line arguments
if(len(sys.argv)!=1):
    try:
        Nx=int(sys.argv[1])
        Ny=int(sys.argv[2])
        radius=float(sys.argv[3])
        Niter=int(sys.argv[4])
        if 2*radius >= Nx or 2*radius >=Ny:
            print('Diameter cannot exceed the grid length. Please input a valid radius.')
            sys.exit()
    except IndexError:
        print("Please input parameters in the following order: Nx, Ny, radius, Niter")
        sys.exit()
    except ValueError:
        print("Please input correct datatypes.")
        sys.exit()

#Constructing potential array
phi=zeros((Ny,Nx))
x=linspace(-Nx/2+0.5, Nx/2-0.5, Nx)
y=linspace(-Ny/2+0.5, Ny/2-0.5, Ny)
X,Y=meshgrid(x,y)
ii=where(X*X+Y*Y <= radius*radius)
phi[ii]=1

#Plotting the potential
figure(1, (9.6,7.2))
cp1=contourf(X, Y, phi, cmap='plasma', alpha=0.35)
colorbar(cp1)
scatter(X, Y, c='red', s=18*phi)
gca().invert_yaxis()
xlabel(r'$x$', size=16)
ylabel(r'$y$', size=16)
title(r'Plot of potential', size=18)
show()

#Using algorithm to calculate solutions and estimating error b/w iterations.
errors=zeros(Niter)
for k in range(Niter):
    #Creating a copy of the old potential
    oldphi=phi.copy()
    #Updating the potential
    phi[1:-1,1:-1]=0.25*(phi[1:-1,0:-2]+phi[1:-1,2:]+phi[0:-2,1:-1]+phi[2:,1:-1])
    #Asserting boundary conditions
    phi[:,0]=phi[:,1]     #left
    phi[:,-1]=phi[:,-2]   #right
    phi[0,:]=phi[1,:]     #top 
    phi[-1,:]=0           #bottom
    phi[ii]=1             #centre
    #Allocating the errors vector
    errors[k]=(abs(phi-oldphi)).max()

#Plotting error vs no. of iterations in loglog scale
iterRange=linspace(1, Niter, Niter)  
figure(2, (9.6,7.2))
loglog(iterRange[::50], errors[::50], 'ro')
xlabel(r'$No.\ of\ iterations$    ($\log$)', size=16)
ylabel(r'$Error$    ($\log$)', size=16)
title(r'$\log\log$ plot of Error and Niter', size=18)
show()

#Fitting the entire vector using least square method
P1=c_[ones(Niter), iterRange]    # P*c=log(errors)
c1=lstsq(P1, log(errors))[0]
A1=exp(c1[0])
B1=c1[1]
print('Fit1 (A,B): ',A1,', ',B1,'\n')
errorsFit1=exp(dot(P1,c1))

#Fitting the portion beyond 500th iteration using least square method
P2=c_[ones(Niter-500), iterRange[500:]]
c2=lstsq(P2, log(errors[500:]))[0]
A2=exp(c2[0])
B2=c2[1]
print('Fit2 (A,B): ',A2,', ',B2)
errorsFit2=exp(dot(P2,c2))

#Plotting error vs Niter for the different fits
figure(3, (9.6,7.2))
semilogy(iterRange[::50], errors[::50], 'ro', label='errors')
semilogy(iterRange, errorsFit1, label='fit1')
semilogy(iterRange[500:],  errorsFit2, label='fit2')
xlabel(r'$No.\ of\ iterations\ (Niter)$', size=16)
ylabel(r'$Error$    ($\log$)', size=16)
title(r'Error vs. No. of iterations', size=18)
legend()
show()

#Creating a 3D surface plot of the potential
ax=p3.Axes3D(figure(4, (9.6,7.2)))
title('The 3-D surface plot of the potential', size=20)
surf = ax.plot_surface(X, Y, phi, rstride=1, cstride=1, cmap=cm.jet,
                       linewidth=0, antialiased=False)
ax.set_ylim(Ny/2,-Ny/2)
colorbar(surf,shrink=0.3,aspect=6)
ax.set_xlabel(r'$x$',size=18)
ax.set_ylabel(r'$y$',size=18)
show()

#Generating the current arrays
Jx=zeros((Ny,Nx))
Jy=zeros((Ny,Nx))
Jx[:,1:-1]=0.5*(phi[:,0:-2]-phi[:,2:])
Jy[1:-1,:]=0.5*(phi[0:-2,:]-phi[2:,:])
Jx[:,0]=0
Jx[:,-1]=0
Jy[0,:]=0

elec=np.zeros((Ny,Nx))      #Electrode location
elec[ii]=1

#Contour plot of potential along with vector plot of currents
figure(5, (9.6,7.2))
cp=contourf(X, Y, phi)
colorbar()
quiver(x, y, Jx, -Jy, scale=4)
clabel(cp, inline=True)
gca().invert_yaxis()
scatter(X, Y, c='red', s=16*elec, label='electrode')
xlabel(r'$x$', size=16)
ylabel(r'$y$', size=16)
title(r'Contour of potential along with current', size=18)
legend()
show()

#Finding numerical solution of temperature from the Poisson's equation
T=zeros((Ny,Nx))

for k in range(Niter):
    T[1:-1,1:-1]=0.25*(T[1:-1,0:-2]+T[1:-1,2:]+T[0:-2,1:-1]+T[2:,1:-1]
                       +square(Jx[1:-1,1:-1])+square(Jy[1:-1,1:-1]))
    #Boundary conditions
    T[ii]=300         #electrode
    T[-1,:]=300
    T[0,:]=T[1,:]     #top
    T[:,0]=T[:,1]     #left
    T[:,-1]=T[:,-2]   #right
    
#Creating a 3D surface plot of temperature
ax1=p3.Axes3D(figure(6, (9.6,7.2)))
title('The 3-D surface plot of the temperature $T$', size=20)
surf = ax1.plot_surface(Y, X, T, rstride=1, cstride=1, cmap=cm.jet,
                       linewidth=0, antialiased=False)
colorbar(surf,shrink=0.3,aspect=6)
ax1.set_xlabel(r'$y$',size=18)
ax1.set_ylabel(r'$x$',size=18)
show()

#Contour plot of temperature
figure(7, (9.6,7.2))
contourf(X, Y, T, cmap='coolwarm')
colorbar()
scatter(X, Y, c='red', s=16*elec, label='electrode')
gca().invert_yaxis()
xlabel(r'$x$', size=16)
ylabel(r'$y$', size=16)
title(r'Contour of temperature', size=18)
legend()
show()