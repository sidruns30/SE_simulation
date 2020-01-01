# SE_simulation
This is a side project where I try to simulate Schrodinger's Equation for various potentials. The end goal is to make an animation of a surface representing the wavefunction subject to a time dependent potential.  

## se.py (v1)
### Seperable, time independent potentials
The file *se.py* is the generates a wavefunction for a seperable 2D potential. Enter V_X and V_y as the potentials in main() along with E (energy), psi_0 and psi_1 (boundary conditon of psi at 0, dx- don't worry about convergence as long as the values are reasonable!) and the x-y arrays. 
I use a [Numerov integrator](https://en.wikipedia.org/wiki/Numerov%27s_method) to integrate the wavefunction. Future versions of *se.py* may have an option of using different integrators to compute the wavefunction. 

## time_dep_se.py (v1)
### Seperable, time independent potential but the wavefunction showing the unitary evolution of the wavefunction
[animation taken from](https://pythonmatplotlibtips.blogspot.com/2018/11/animation-3d-surface-plot-funcanimation-matplotlib.html)

## matrix_se.py (v1)
*matrix_se.py* uses the matrix representation of the Schrodinger equation to solve for the wavefunction. It involves calculating the Hamiltonian matrix, where each element of the matrix H(i,j) is given by <i|H|j>, |i>, <j| being the kets and the bras respecitvely. I use the basis functions for the infinite square well: sin(kx), cos(kx), k=1,2,...n to calculate the Hamiltonian matrix. 

In the calculation, the Hamiltonain matrix is split into H_mat_K and H_mat_U, the matrices for the kinetic and potential energies respectively. H_mat_K is set to be a diagonal matrix with entries (1,1,4,4,9,9,...n^2) depending on what n is set to. 

The core of the file involves calulating H_mat_U where the inner product <i|U|j> is carried out. The program uses an RK4 integrator to perform this and the values are appended to H_mat_U. 

The full hamiltonian matrix (H) is then just H_mat_K + H_mat_U. The eigenvalues and eigenvectors of H are computed using numpy's [linalg module](https://docs.scipy.org/doc/numpy/reference/generated/numpy.linalg.eig.html) and the corresponding wavefunctions are then plotted. 

The code is stil under development and for now only solves one dimensional potentials.  
