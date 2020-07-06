from math import sin, exp, pi, cos, sqrt, tanh
import matplotlib.pyplot as plt
import numpy as N
import random
L = 1.0
u = 0.5
h = 1.0
dt = 0.01 #timestep size
M = 1.0 #mobility
dx = dy = h #grid spacing size
A = 1.0 
kappa = 0.5 #gradient coefficient
nx = 64#(int)(L/dx) #number of grid points 
ny = nx
q = 20000#(int)(u/dt)
c0 = 0.4 #Average concentration value
U = N.zeros((nx,ny)) #concentration at next step
f0 = N.zeros((nx,ny))#variation of the chemical energy
c = N.zeros((nx,ny)) #initial concentration

c = c0 + 0.02*(0.5-N.random.rand(nx,ny)) #initial condition 
f0 = 2*c*(1-c)**2-2*c**2*(1-c)#c**3 - c #Chemical energy double-well potential
mfig = [1,20,50,100] #output figures at these timesteps
fig = plt.figure()
num = 0
for m in range(1,101):
    for i in range(nx):
        for j in range(ny):
            jp=j+1
            jm=j-1
            ip=i+1
            im=i-1
            jpp=j+2
            jmm=j-2
            ipp=i+2
            imm=i-2
            if ip > nx-1:  # periodic boundary condition
                ip = ip - nx
            if im < 0:
                im = im + nx
            if jp > ny-1:
                jp = jp - ny
            if jm < 0:
                jm = jm + ny
            if ipp > nx-1: 
                ipp = ipp - nx
            if imm < 0:
                imm = imm + nx
            if jpp > ny-1:
                jpp = jpp - ny
            if jmm < 0:
                jmm = jmm + ny
                Uxxxx = (c[imm,j]-4*c[im,j]+6*c[i,j]-4*c[ip,j]+c[ipp,j])/dx**4 
                Uyyyy = (c[i,jmm]-4*c[i,jm]+6*c[i,j]-4*c[i,jp]+c[i,jpp])/dy**4
                Uxxyy = (c[ip,jp]-2*c[i,jp]+c[im,jp]-2*c[ip,j]+4*c[i,j]-2*c[im,j]+c[ip,jm]-2*c[i,jm]+c[im,jm])/(dy**2*dx**2) #centered finite differences
                fxx = (f0[ip,j] - 2*f0[i,j] + f0[im,j]) / dx**2
                fyy = (f0[i,jp] - 2*f0[i,j] + f0[i,jm]) / dy**2   
                H = Uxxxx+Uyyyy+2*Uxxyy #second Laplacian of the concentration
                F = fxx + fyy #Laplacian of the derivative of the Chemical energy
                U = c + M*dt*(F - kappa*H) #Cahn-Hilliard equation
    c = U
    if m in mfig:
        num += 1
        ax = fig.add_subplot(220 + num) #2x2 grid, subplot: 0 + num
        im = ax.imshow(U,interpolation='lanczos',cmap='bwr')
fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.9, 0.15, 0.03, 0.7])
fig.colorbar(im,cax=cbar_ax)
plt.show()
