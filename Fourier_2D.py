from math import exp, sqrt
import numpy as N
import matplotlib.pyplot as plt
import random
import mpld3
#2D Heat diffusion equation 
t = 100 #final time 
dx = 1.0 #grid spacing size
dy = dx 
dt = 0.01 #time interval size
n=11 
nx = ny = n #number of grid points
T_tb = 10 #temperature at the ends
T_lr = 0
M = 4.0 #diffusivity
T = N.zeros((n,n))
L = N.zeros((n,n))

#boundary conditions:
T[0,:] = T[n-1,:] = T_tb #up and down
T[:,0] = T[:,n-1] = T_lr #left and right
#Definition of Laplacian:
#fig = plt.figure()
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
#Solution is displayed as animation
for m in range(t):
    #plt.cla()
    F = laplace(T,dx,dy)
    Tn = T + M*dt*F #Heat equation 
    T = Tn
plt.title('2D Heat Propagation')
plt.xlabel('X (m)')
plt.ylabel('Y (m)')
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
im = plt.imshow(Tn,interpolation='lanczos',cmap='hot')    
    #plt.pause(0.00001)
cbar = plt.colorbar(im)
cbar.set_label('T (K)')
#plt.show()
html_str = mpld3.fig_to_html(fig, no_extras=True, template_type='simple')
#Html_file= open("index.html","w")
#Html_file.write(html_str)
#Html_file.close()
print(html_str)


