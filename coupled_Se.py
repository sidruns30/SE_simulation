import numpy as np
from numpy import append
import time
import matplotlib.pyplot as plt

def RK4(x_0, f_0, g, h):
    k1 = h*g(x_0,f_0)
    k2 = h*g(x_0 + h/2, f_0 + k1/2)
    k3 = h*g(x_0 + h/2, f_0 + k2/2)
    k4 = h*g(x_0 + h, f_0 + k3)

    x_1 = x_0 + h
    f_1 = f_0 + (1/6)*(k1 + 2*(k2 + k3) + k4) 
    return x_1, f_1

def integrate(x_0, x_f, y_0, g, h):  #f and g are fucntions, h is step size
    (x,y) = *map(np.asarray, ([x_0], [y_0])),

    while(x[-1]<x_f):
        (a, b) = RK4(x[-1], y[-1], g, h)
        x = append(x, a)
        y = append(y, b)
    return x,y
    
def inner_prod(x_0, x_f, f, g, h): 
    
    y_ = lambda x, _: np.conj(f(x))*g(x)
    y_0 = y_(0, None) #y_ depends only of first para

    (x, y) = integrate(x_0, x_f, y_0, y_, h)    
    return y[-1]-y[0]
    
    
#Set boundaries   
x_0 = -np.pi
x_f = np.pi
y_0 = 0
h_bar = 1

#Create basis functions
e = [lambda x: np.exp(1j*k*x) for k in range(5)]

#Initialize Hamiltonian matrix
H_mat = np.empty((len(e),len(e)),dtype=complex)

#Set kinetic and potential energy functions, find hamiltonian
T = lambda k: -h_bar*1j*k**2 #Is like the derivative on exp basis
U = lambda x: 0
H = lambda x: T(x) + U(x)

ket = [lambda x: H(x)*e[i](x) for i in range(len(e))]  #Function on the right side
bra = e

for i in range(len(H_mat[:,0])):
    for j in range(len(H_mat[0,:])):
        H_mat[i,j] = inner_prod(x_0, x_f, bra[i], ket[j], h=0.01)

print(H_mat)


