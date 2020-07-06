from math import sin, exp, pi, cos, sqrt
import matplotlib.pyplot as plt
import numpy as N
#This is my very own implementation of the Crank-Nicolson method for the 1D heat equation
L = 1.0
u = 0.2
alpha = 4.0
dt = 0.001
dx = 0.01
K = alpha*dt/dx**2
n = (int)(L/dx)
q = (int)(u/dt)
x = N.linspace(0, L, n+1) #for some reason a stupid list can be plotted instead of an array.
t = N.linspace(0, u, q+1) 
A = N.zeros((n+1,n+1))
B = N.zeros((n+1,n+1))
b = N.zeros(n+1)
T = N.zeros(n+1)
Tn = N.zeros(n+1)
def f(x):
    return sin(pi*x)
for i in range(1,n):
        T[i] = f(x[i]) #initial condition
for i in range(1,n):
    A[i,i-1] = -K/2
    A[i,i+1] = -K/2   #A and B are both tri-diagonal matrices but the
    A[i,i] = 1 + K    #initial conditions must be forced for consistency
    A[0,0] = A[n,n] = 1
for i in range(1,n):
    B[i,i-1] = K/2
    B[i,i+1] = K/2
    B[i,i] = 1 - K
    B[0,0] = B[n,n] = 1
plt.figure(1)
for i in range(0,q): #temporal evolution of the temperature for each T(x)
    plt.cla()
    for j in range(1,n):
        b[j] = T[j]
    b[0] = b[n] = 0
    Tn[:] = N.linalg.solve(A,N.dot(B,b)) #Solve equation A*x[n+1] = B*x[n] for x[n+1]
    T[:] = Tn
    plt.xlabel('Distance')
    plt.ylabel('Temperature')
    plt.title('1D heat propagation')
    plt.axis([0,L,-0.1,1.0])
    plt.grid()
    plt.plot(x,abs(T))
    plt.pause(0.0001)
plt.show()


T = N.array(T)
print(T)
