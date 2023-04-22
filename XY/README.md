# Classical XY model in any lattices #
---

## Hamiltonian ##
H = -J\sum_{ij} S_i J_{ij} S_j
J>0

## Representation ##
expand \exp[K\cos(\theta)] with the modiﬁed Bessel function

## Algorithms ##
worm algorithm


## Some symbol names for the preprocessor
- IFORT: when using ifort, this must be turn on
- GFORT: when using gfortran, this must be turn on
- MPI: mpi version
- THERMALIZATION: to print the evolution of energy to decide the thermalization time
- ACF: use this symbol to print observable data for calculating autocorrelation function
- RESTART: after serval blocks, restart the simulation from a random configurations
