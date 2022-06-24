#hubbard_complex is used to create Hamiltonian and to solve for the eigenvalues and eigenvectors.

##On TACC frontera, load intel compiler as well as PETSc and SLEPc modules:

module  load  intel/18.0.5
module  load  petsc/3.11-complexi64
module  load  slepc/3.11-complexi64

##Then type

make
