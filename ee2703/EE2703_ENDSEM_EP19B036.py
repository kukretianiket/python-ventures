"""
            EE2703 Endsem by Aniket Kukreti EP19B036
"""
from pylab import *
import scipy.linalg as sp

#Initialising parameters
a=10        #Radius of loop antenna
sect=100    #Number of sections

#Creating the 3 by 3 by 1000 mesh
x=linspace(-1, 1, 3)
y=linspace(-1, 1, 3)
z=linspace(1, 1000, 1000)
xx, yy, zz=meshgrid(x, y, z)

#Breaking the loop into sections
phi=linspace(0, 2*pi, sect+1)   #Points on the loop
phi_l=0.5*(phi[:-1]+phi[1:])    #phi value of the current elements
#Array of all position vectors (x, y, z) of current elements (below)
r_=c_[a*cos(phi_l), a*sin(phi_l), np.zeros(sect)]   #(100, 3)
#Array of all dl vectors (below)
dl_=(a*2*pi/sect)*c_[-sin(phi_l), cos(phi_l)]   #(100, 2)
#Array of all position vectors (x, y, z) of points in the mesh (below) 
r=c_[xx.reshape(-1), yy.reshape(-1), zz.reshape(-1)]    #(9000, 3)
#Array containing current directions 
i=c_[-cos(phi_l)*sin(phi_l), (cos(phi_l))**2]   #(100, 2)

#Plotting the current elements of the ring
figure(1, (9.6, 9.6))
quiver(r_[:,0], r_[:,1], i[:,0], i[:,1], label='Current flow')
scatter(r_[:,0], r_[:,1], c='r', label='Current element', s=20)
xlabel('$x$', size=16)
ylabel('$y$', size=16)
title('Plot of the current elements', size=20)
legend()
grid()
show()

#Defining the calc function
def calc(l):
    """
    Calculates distance between all points in the mesh and current element 'l'.
    
    Returns a (3, 3, 1000) array.
    """
    R_=r[:]-r_[l]   #Difference between every vector r and chosen r'_l
    R=sqrt((R_[:, 0])**2+(R_[:, 1])**2+(R_[:, 2])**2)   #Computation of the magnitude
    return R.reshape((3, 3, 1000))

#Finding the magnetic potential through summation
Ax=np.zeros_like(calc(0), dtype=complex)    #Initialising x component of A
Ay=np.zeros_like(calc(0), dtype=complex)    #Initialising y component of A
#Using for loop to sum over all 100 elements
for l in range(sect):
    Ax+=cos(phi_l[l])*exp(-0.1*1j*calc(l))*dl_[l][0]/calc(l)
    Ay+=cos(phi_l[l])*exp(-0.1*1j*calc(l))*dl_[l][1]/calc(l)

#Using curl approximation to obtain B
B=0.25*(Ay[1,2,:] - Ay[1,0,:] - Ax[2,1,:] + Ax[0,1,:])      #Vectorised sum

#Estimating the parameters c & b
M=c_[ones(1000), log(z)]    #M*v=log(B)
v=sp.lstsq(M, log(abs(B)))[0]   #Least square approximation of v
c=exp(v[0])
b=v[1]

#Plotting magnetic field and its fit
figure(2, (14.4, 9.6))
loglog(z, abs(B), label='$B_z(z)$', linewidth=0.8)
loglog(z, c*z**b, label='$cz^b$')
xlabel('$z$    $(\log)$', size=16)
ylabel('$B_z$    $(\log)$', size=16)
title('$\log\log$ plot of $B_z$ vs $z$', size=19)
annotate('$c={}$\n$b={:.6f}$'.format(c, b), (2, 1e-15))
legend()
grid()
show()

#Static case
staticAx=np.zeros_like(calc(0))
staticAy=np.zeros_like(calc(0))
for l in range(sect):
    staticAx+=cos(phi_l[l])*dl_[l][0]/calc(l)   #Complex exponential is not present
    staticAy+=cos(phi_l[l])*dl_[l][1]/calc(l)   #Complex exponential is not present
  
staticB=0.25*(staticAy[1,2,:] - staticAy[1,0,:] - staticAx[2,1,:] + staticAx[0,1,:])

v2=sp.lstsq(M, log(abs(staticB)))[0]
c2=exp(v2[0])
b2=v2[1]

#Plotting magnetic field along with its fit for the static case
figure(3, (14.4, 9.6))
loglog(z, abs(staticB), label='$B_z(z)$', linewidth=0.8)
loglog(z, (c2)*z**(b2), label='$cz^b$')
xlabel('$z$    $(\log)$', size=16)
ylabel('$B_z$    $(\log)$', size=16)
title('$\log\log$ plot of $B_z$ vs $z$ (static case)', size=19)
annotate('$c={}$\n$b={:.6f}$'.format(c2, b2), (2, 7e-16))
legend()
grid()
show()

