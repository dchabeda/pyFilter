# pyFilter

tl;dr Submit a filter calculation with >> python main.py. The code takes inputs: conf.par, input.par, and potentials pot*.par.

This code calculates electronic states of an atomistic, pseudopotential-based, quasiparticle Hamiltonian (J. Chem. Phys. 157, 020901 (2022); Eq. 2). The calculation has 5 steps:

1. INITIALIZE: read in parameters and set macros, calculate the value of the potential on a basis of real space grid points, set up energy targets for filter diagonalization.
2. GET ENERGY RANGE: compute max(H) - min(H), the energy range of the Hamiltonian (algorithmic parameter for efficient filtering).
3. FILTER: apply the filter diagonalization procedure to obtain states within a Hilbert subspace having eigenvalues close to the target energies.
4. ORTHOGONALIZE: orthogonalize the filtered states using an SVD procedure to obtain a basis for subspace diagonalization.
5. DIAGONALIZE: diagonalize the Hamiltonian matrix constructed in the basis of orthogonal states, print the wavefunctions and energy eigenvalues.

The code is separated into multiple files by functionality, and the dependencies of each file are separated into separate file descriptors (fd_main.py, fd_init.py, etc). 
