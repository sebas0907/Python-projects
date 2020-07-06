from math import cos, sin, pi
import matplotlib.pyplot as plt
import numpy as N
import pylab
#definition of the differential equation dy/dx = F(x,y)
def F(x , y):
    return y*cos(x)
print("Numer of steps n: ")
n = int(input()) #number of steps. Warning! this is the time complexity
x = 0.0 #initial conditions
print("initial condition y0: ")
y = float(input())
print("Step size h: ")
h = float(input()) #step size
xn = [] #initial empty columns to store the obtained values after evaluation
yn = []
#implementation of the Runge-Kutta method
if h == 0:
    h = 0.01
for i in range(0,n):
        k1 = F(x, y)*h
        k2 = F(x + h/2, y + k1/2)*h
        k3 = F(x + h/2, y + k2/2)*h
        k4 = F(x + h, y + k3)*h
        y = y + (k1 + 2*k2 + 2*k3 + k4)/6
        x = x + h
        xn.append(x)
        yn.append(y)
        #print(x, y, sep = '\t\t')
#plt.ylim((0.0, 1.1))
plt.xlim((0.0, 2*pi))
plt.plot(xn,yn)
plt.xlabel('x')
plt.ylabel('y(x)')
plt.grid()
plt.show()

        
        

    
