
# Description

hubbard_complex is used to create 16-site square lattice Hamiltonian and to solve for the eigenvalues and eigenvectors by exact diagonalization.


# Compilation 
On TACC frontera, load intel compiler as well as PETSc and SLEPc modules:

- module  load  intel/18.0.5
- module  load  petsc/3.11-complexi64
- module  load  slepc/3.11-complexi64
- make
