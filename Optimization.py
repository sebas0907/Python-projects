import numpy as N
import matplotlib.pyplot as plt

#define function:
def G(r):
    x = r[0]
    y = r[1]
    return -x**2/2+x**4/4 #double-well potential
            #(1-x)**2+100*(y-x**2)**2  #Rosebrock function
            #(x-1)**2*N.exp(-y**2) + y*(y+2)*N.exp(-2*x**2)
        
#define gradient:
def gradG(z,h,d_x):
    return (G(N.add(z,h))-G(N.subtract(z,h)))/(2*d_x)

r0 = [-0.1,0.0] #initial guess
hx = [0.01,0.0] #step in the x-direction
hy = [0.0,0.01] #step in the y-direction
dx = 0.01 #step size in x
dy = 0.01 #step size in y
n = 10000 #number of steps

#the programm will only return values close to a non-zero guess , i.e. either x or y are empty
#for a one-dimensional function or both are non-empty for a two-dimensional function. An initial zero-vector will thus throw an error (work on progress.)
print('initial guess: ')
if r0[0] != 0:
    print('x0 = ',r0[0])
if r0[1] != 0:
    print('y0 = ',r0[1])
    
#iteration using the steepest descent method with adaptative step size:
for i in range(n):
    DxG,DyG = [gradG(r0,hx,dx),gradG(r0,hy,dy)]
    gradient = [DxG,DyG]
    delG = N.sqrt(N.dot(gradient,gradient))
    alpha = 0.001
    beta = alpha/delG
    rn = r0 - N.dot(beta,gradient)
    r0 = rn
if rn[0] != 0:
    print('minimum located at: x = ',rn[0])
if rn[1] != 0:
    print('minimum located at: y = ',rn[1])
print("value at minimum: G = ",G(rn))


 



