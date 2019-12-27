from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np 

#Defining global variables h bar and m
h_b = 1  
m = 1

class system():
    def __init__(self, x, V, psi_0, psi_1, E):
        '''
        x: array like- eveny sampled position array
        V: array like- potential array
        psi_0, psi_1: float- initial conditons for the wave function psi(0), psi(dx)
        E: float- energy of the system
        '''
        self.x, self.V = map(np.asarray, (x, V))
        self.N = self.x.shape[0]
        assert (self.V.shape[0] == self.N), 'V and x must have the same shape'
        self.E, self.psi_0, self.psi_1 = map(np.float, (E, psi_0, psi_1))
        self.psi = np.array([self.psi_0, self.psi_1])
        self.h = self.x[1] - self.x[0]
        
<<<<<<< HEAD
    def func(self, V):    
        return 2*m/(h_b**2)*(V-np.ones(self.N)*self.E)
=======
    def f(self, V):    
        return -2*m/(h_b**2)*(V-np.ones(self.N)*self.E)
>>>>>>> 2d792a8c779709eaacfef5b240db329870d5e941
    
    def numerov(self, k): #return k+1 step 
        y_1 = self.psi[-1]
        y_0 = self.psi[-2]
        h = self.h
        f = self.f(self.V)
        y_2 = (2*(1 - (5/12)*h**2*f[k])*y_1 - (1 + h**2/12*f[k-1])*y_0)/(1 + h**2/12*f[k+1])
        return y_2
    
    def integrate(self):  
        k = 0
        while(k<self.N-2):
            self.psi = np.append(self.psi, [self.numerov(k)], axis=-1)
            k+=1
        self.psi = self.psi/np.trapz(self.psi**2, dx=self.h)
        return self.psi


def main():
    x = np.linspace(-np.pi, np.pi, 1000)
    y = np.linspace(-np.pi, np.pi, 1000)
    
    def f_1(x):
        v = np.exp(-np.absolute(x))
        return v
    
    def f_2(y):
        v = np.exp(-np.absolute(y))
        return v
        
    V_x = f_1(x)
    V_y = f_2(y)
    psi_0 = 0
    psi_1 = 0.0001
    E = -2
    
    s_x = system(x, V_x, psi_0, psi_1, E/2)
    psi_x = s_x.integrate()
    s_y = system(y, V_y, psi_0, psi_1, E/2)
    psi_y = s_y.integrate()
    
    psi = np.array([k*psi_y for k in psi_x])
    V = np.array([k + V_y for k in V_x])
    xx, yy = np.meshgrid(x,y)
    
    fig = plt.figure()
    ax1 = ax1 = fig.add_subplot(121,projection='3d')
    surf = ax1.plot_surface(xx, yy, psi, cmap=cm.inferno,
                       linewidth=0, antialiased=False)
    plt.title('Wavefunction')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.legend()

   
    ax2 = ax1 = fig.add_subplot(122,projection='3d')
    surf = ax2.plot_surface(xx, yy, V, cmap=cm.viridis,
                       linewidth=0, antialiased=False)
    plt.title('Potential')
    plt.xlabel('X')
    plt.ylabel('Y')   
    
    plt.show()
    plt.savefig('plots/'+ 'exponential potenials.png')

    return None

if __name__ =="__main__":
	main()
