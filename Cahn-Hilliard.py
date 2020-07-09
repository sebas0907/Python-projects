import numpy as N
import matplotlib.pyplot as plt
import random
T = 100 #final time
dx = 1.0
dy = dx
dt = 0.01
n=64
nx = ny = n
kappa = 0.5 #grad-coefficient
c0 = 0.4 #average cocentration
M = 4.0 #mobility
A = N.zeros((n,n))
L = N.zeros((n,n))
H = N.zeros((n,n)) #clear variables
f = N.zeros((n,n))
F = N.zeros((n,n))

A = c0 + 0.02*(0.5-N.random.rand(n,n)) #random noise
#boundary conditions
A[0,:] = A[n-1,:]
A[:,0] = A[:,n-1] 
f= A**3-A#2*A*(1-A)**2-2*A**2*(1-A) #variation of chemical energy
g= A**2-1 #second derivative of chemical energy
g0=c0**2-1 #second derivative evaluated at the average 
#fig = plt.figure()
#mfig = [1,100,200,300]
def laplace(Z,d_x,d_y):
    for i in range(n):
        for j in range(n):
            jp=j+1
            jm=j-1
            ip=i+1
            im=i-1
            if ip > nx-1:  
                ip = ip - nx
            if im < 0:
                im = im + nx
            if jp > ny-1:
                jp = jp - ny
            if jm < 0:
                jm = jm + ny
            uxx = (Z[ip,j] - 2*Z[i,j] + Z[im,j]) / dx**2
            uyy = (Z[i,jp] - 2*Z[i,j] + Z[i,jm]) / dy**2
            L[i,j] =  (uxx + uyy)
    return L
print('wait please...')
def del4(c,d_x,d_y):
    for i in range(n):
        for j in range(n):
            jp=j+1
            jm=j-1
            ip=i+1
            im=i-1
            jpp=j+2
            jmm=j-2
            ipp=i+2
            imm=i-2
            if ip > n-1:  
                ip = ip - n
            if im < 0:
                im = im + n
            if jp > n-1:
                jp = jp - n
            if jm < 0:
                jm = jm + n
            if ipp > n-1: 
                ipp = ipp - n
            if imm < 0:
                imm = imm + n
            if jpp > n-1:
                jpp = jpp - n
            if jmm < 0:
                jmm = jmm + n
            uxxxx = (c[imm,j]-4*c[im,j]+6*c[i,j]-4*c[ipp,j]+c[ipp,j])/dx**4 
            uyyyy = (c[i,jmm]-4*c[i,jm]+6*c[i,j]-4*c[i,jp]+c[i,jpp])/dy**4
            uxxyy = (c[ip,jp]-2*c[i,jp]+c[im,jp]-2*c[ip,j]+4*c[i,j]-2*c[im,j]+c[ip,jm]-2*c[i,jm]+c[im,jm])/(dy**2*dx**2)
            H[i,j] = uxxxx +uyyyy +2*uxxyy
    return H
num = 0
#Animation
for m in range(T):
    plt.cla()
    F = laplace(A,dx,dy)
    D = del4(A,dx,dy)
    An = A + M*dt*(g0*F - kappa*D) #Cahn-Hilliard equation
    A = An
##    if m % (T/4) == 0: #print four figures for equal time intervals
##        num += 1
##        ax = fig.add_subplot(220 + num) #2x2 grid, subplot: 0 + num
##        im = ax.imshow(An,interpolation='lanczos',cmap='bwr')
##fig.subplots_adjust(right=0.8)
##cbar_ax = fig.add_axes([0.9, 0.15, 0.03, 0.7])
##fig.colorbar(im,cax=cbar_ax)
##plt.show()
    plt.title('Spinoidal decomposition')
    im = plt.imshow(An,interpolation='lanczos',cmap='bwr')
    plt.pause(0.00001)
cbar = plt.colorbar()
cbar.set_label('Concentration')
plt.show()


