import numpy as N
import matplotlib.pyplot as plt
import random
from math import pi
#define function:
def G(r):
    x = r[0]
    y = r[1]
    return (1-x)**2+100*(y-x**2)**2  #Rosenbrock function
 
#define gradient:
def gradG(z,h,d_x):
    return (G(N.add(z,h))-G(N.subtract(z,h)))/(2*d_x)

r0 = [-1.5,2.5]#initial guess
hx = [0.01,0.0]#step in the x-direction
hy = [0.0,0.01]#step in the y-direction
dx = 0.01#step size in x
dy = 0.01#step size in y
n = 30000#number of steps

print('initial guess: ')
print('x0 = ',r0[0])
print('y0 = ',r0[1])
plt.plot(r0[0],r0[1],'og',label= 'Guess')

#iteration using the steepest descent method:
for i in range(n):
    DxG,DyG = [gradG(r0,hx,dx),gradG(r0,hy,dy)]
    gradient = [DxG,DyG]
    delG = N.sqrt(N.dot(gradient,gradient))
    alpha = 0.0002
    beta = alpha/delG
    rn = r0 - N.dot(beta,gradient)
    r0 = rn
    alpha = beta
print("after {} iterations the series converged to: ".format(str(n)))
if rn[0] != 0:
    print('minimum located at: x = ',rn[0])
if rn[1] != 0:
    print('minimum located at: y = ',rn[1])
print("value at minimum: G(r0) = ",G(rn))

X = N.linspace(-2.0,2.0,100)
Y = N.linspace(-1.0,3.0,100)
x1,y1 = N.meshgrid(X,Y)
Z = (1-x1)**2+100*(y1-x1**2)**2
plt.contour(X,Y,Z,N.logspace(-0.5,3.5,20,base=10),cmap='gray')
plt.xlabel('x')
plt.ylabel('y')
plt.title('Gradient descent method for the Rosenbrock function. \n The theoretical minimum is f(1,1)=0')
plt.plot(rn[0],rn[1],'vr',label = 'Result')
plt.legend()
plt.show()
 



