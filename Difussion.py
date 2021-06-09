from math import sin, exp, pi, cos, sqrt
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np

#This is an explicit implementation of the Crank-Nicolson method for the heat equation in the case of a one-dimensional rod.
L = 1.0 #Lenght of the rod
u = 0.01 #final time
alpha = 1.0 #Diffusion coefficient
dt = 0.0001 #time interval
dx = 0.01 #grid interspacing size
K = alpha*dt/dx**2 
n = int(L/dx) #number of grid points
q = int(u/dt) #number of time intervals
x = np.linspace(0, L, n+1) #definition of axes
t = np.linspace(0, u, q+1) 
A = np.zeros((n+1,n+1)) #initial variables
B = np.zeros((n+1,n+1))
b = np.zeros(n+1)
T = np.zeros(n+1) #temperature at previous step
Tn = np.zeros(n+1) #Temperature at new step

def f(x):
    return np.sin(np.pi*x)

#initial condition:
T[:] = f(x)

for i in range(1,n):
    A[i,i-1] = -K/2
    A[i,i+1] = -K/2   #A and B are both tri-diagonal matrices and the
    A[i,i] = 1 + K    #initial conditions must be forced here for consistency
    A[0,0] = A[n,n] = 1
for i in range(1,n):
    B[i,i-1] = K/2
    B[i,i+1] = K/2
    B[i,i] = 1 - K
    B[0,0] = B[n,n] = 1

fig, ax = plt.subplots()
line, = ax.plot(x,f(x))
ax.set_title('1D heat propagation')
ax.set_xlabel('Distance')
ax.set_ylabel('Temperature')
ax.grid(True)

def animate(m):
    for m in range(0,q): #temporal evolution of the temperature for each T(x)
        b[:] = T
        b[0] = b[n] = 0
        Tn[:] = np.linalg.solve(A,np.dot(B,b)) #Solve equation A*T[n+1] = B*T[n] for T[n+1]
        T[:] = Tn
        line.set_ydata(T)
    return line

ani = animation.FuncAnimation(fig, animate, interval=100)
plt.show()
