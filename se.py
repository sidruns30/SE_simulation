import numpy as np
import matplotlib.pyplot as plt

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
        assert (self.x.shape== self.V.shape), 'V and x must have the same shape'
        self.N = self.x.shape[0]
        self.E, self.psi_0, self.psi_1 = map(np.float, (E, psi_0, psi_1))
        self.psi = np.array([self.psi_0, self.psi_1])
        self.h = self.x[1] - self.x[0]
        
    def func(self):    
        return -2*m/(h_b**2)*(self.V-np.ones(self.N)*self.E)
    
    def numerov(self,k): #return k+1 step 
        y_1 = self.psi[-1]
        y_0 = self.psi[-2]
        h = self.h
        f = self.func()
        y_2 = (2*(1 - (5/12)*h**2*f[k])*y_1 - (1 + h**2/12*f[k-1])*y_0)/(1 + h**2/12*f[k+1])
        return y_2
    
    def integrate(self):  #does the actual integration
        k = 0
        while(k<self.N-2):
            self.psi = np.append(self.psi, [self.numerov(k)], axis=-1)
            k+=1
        self.psi = self.psi/max(self.psi)
        return self.psi

def main():
    x = np.linspace(-np.pi, np.pi, 1000)
    f = lambda x: np.sin(x)
    V = f(x)
    psi_0 = 0
    psi_1 = 0.01
    E = 5
    h_osc = system(x, V, psi_0, psi_1, E)
    p = h_osc.integrate()
    plt.plot(x, V, 'k')
    plt.plot(x,p, 'r')
    return None

if __name__ =="__main__":
	main()
