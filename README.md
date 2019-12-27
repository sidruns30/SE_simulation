# SE_simulation
This is a side project where I try to simulate Schrodinger's Equation for various potentials. The end goal is to make an animation of a surface representing the wavefunction subject to a time dependent potential.  

## se.py (v1)
### Seperable, time independent potentials
The file *se.py* is the generates a wavefunction for a seperable 2D potential. Enter V_X and V_y as the potentials in main() along with E (energy), psi_0 and psi_1 (boundary conditon of psi at 0, dx- don't worry about convergence as long as the values are reasonable!) and the x-y arrays. 
I use a [Numerov integrator](https://en.wikipedia.org/wiki/Numerov%27s_method) to integrate the wavefunction. Future versions of *se.py* may have an option of using different integrators to compute the wavefunction. 

