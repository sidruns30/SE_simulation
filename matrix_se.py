import numpy as np
from numpy import append
import time
import matplotlib.pyplot as plt
import numpy.linalg as lalg
from matplotlib import cm


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
    
    y_ = lambda x, _: (f(x))*g(x)
    y_0 = y_(0, None) #y_ depends only of first para

    (x, y) = integrate(x_0, x_f, y_0, y_, h)    
    return 2*(y[-1]-y[0])/(x_f-x_0)
    
    
#Set boundaries   
x_0 = -np.pi
x_f = np.pi
l = x_f - x_0
y_0 = 0
h_bar = 1
n_basis = 12  #Number of basis vectors
h = 0.005
x = np.linspace(x_0,x_f,l/h)


#Create basis functions (sin(kx) if k odd or cos(kx) if k even)
def generate_basis(n):
    k = n*np.pi/l
    def f(x):
        if (n%2!=0):
            return np.sin((k+1)*x)
        else:
            return np.cos((k+1)*x)
    return f

e = [generate_basis(k) for k in range(n_basis)] 
bra = e

#Make the H matrix for kinetic energy
t1 = [(i+1)**2 for i in range(int(n_basis/2))]
t2 = [(i+1)**2 for i in range(n_basis - int(n_basis/2))]
temp = [None]*(n_basis)
temp[::2] = t2
temp[1::2] = t1

H_mat_K = np.diag(temp)

#Set potential
U = lambda x: x**10
ket_U = []

#Get kets functions for <i\U\j>
for i in range(n_basis):
    ket_U.append(lambda x: U(x)*e[i](x))  

#Get H matrix for potential energy term
H_mat_U = np.zeros((len(e),len(e))) 

start = time.time()
  
for i in range(n_basis):
    for j in range(n_basis):
        H_mat_U[j,i] = inner_prod(x_0, x_f, bra[j], ket_U[i], h)

end = time.time()
print("Elapsed calculation time: " + str(end-start) + ' s')

H = H_mat_K + H_mat_U

print('The Hamiltonian Matrix is: ')
print(np.round(H,1))


#Get eigenvalues and eigenvectors
E, v = lalg.eig(H) #Side note: e has complex values with small values for j part
ind = np.argsort(E)
E = E[ind]
v = v[ind]

print('Allowed energies are: ')
print(E)

fig = plt.figure(figsize=(20,20))
plt.plot(x, U(x), '--', label='$U(x)$')
#Plot wavefunctions
for j in range(3):
    psi = [bra[i](x)*v[j][i] for i in range(n_basis)]
    psi = np.sum(psi, axis=0)
    
    plt.title('Ground State Wavefunction', fontsize=30)
    plt.plot(x, psi, label='$\psi_{'+str(j)+'} (x)$'+ '; ' + '$ E = $' + str(round(E[j],2)))
    plt.grid()
    plt.legend(fontsize=25)