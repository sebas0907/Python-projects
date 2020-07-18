import numpy as N
import scipy
import matplotlib.pyplot as plt
import random

#Lest squares method to fit an exponential function of the form c*exp(-d*x) + k to some randomly generated data:

n = 50 #number of data points
c = 2.5
d = 1.3 #parameters can be given by the user
k = 0.0

#fake data:
Xi = N.linspace(0,4,n)
def g(Z,c1,d1,k1):
    return c1*N.exp(-d1*Z)+k1
Y = g(Xi,c,d,k)

Yi = Y + 0.01*N.random.normal(Xi) #add some noise

#Least squares method:
up = n-12 #lower and upper limits for the fit
lo = 0
#linearization:
xi = Xi
yi = N.log(Yi-k)
#averages:

xa = N.mean(xi[lo:up]) 
ya = N.mean(yi[lo:up])

#slope:
m = N.sum(yi[lo:up]*(xi[lo:up]-xa))/N.sum((xi[lo:up]*(xi[lo:up]-xa)))
#intercept:
b = ya - m*xa
#fit function:
def f(Z,m,b):
    return m*Z + b

#errors:
ssx = N.sum((xi[lo:up]-xa)*(xi[lo:up]-xa))
ssy = N.sum((yi[lo:up]-ya)*(yi[lo:up]-ya))
err = N.sum((yi[lo:up]-f(xi[lo:up],m,b))*(yi[lo:up]-f(xi[lo:up],m,b))) 
#standard Deviation: (there are n - 2 degrees of freedom because of the two extracted parameters)
sigma = N.sqrt(err/(len(range(lo,up))-2))
#error in slope:
me = sigma/N.sqrt(ssx)
#error in intercept:
be = sigma*N.sqrt(1/n + xa**2/ssx)
#coefficient of determination:
R2 = (ssy - err)/ssy
#reduced xi-square:
chi = err/(len(range(lo,up))-2)

#exponential fit:
def h(Z,c1,d1,k1):
    return N.exp(c1+d1*Z) + k1

print('Set of parameters obtained: ')
print('c = ',N.exp(b),'\u00B1',N.exp(b)*be)
print('d = ',N.abs(m),'\u00B1',me)
print('sigma =',sigma)
print('Coefficient of determination =',R2) #R\u00B2
print('Reduced chi-square = ',chi) #'\u03C7

#plot scattered data values and fit function:
plt.plot(Xi,Yi,'+',label='Data')
plt.plot(Xi,h(Xi,b,m,k),label='Fit')
#plt.plot(xi,yi,'+',label='Data')
#plt.plot(xi,f(xi,m,b),label='Fit')
plt.grid()
plt.legend()
plt.title('Least squares fit')
plt.xlabel('x')
plt.ylabel('y')
plt.show()
  

