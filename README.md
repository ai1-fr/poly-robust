# poly-robust
Matlab files for the computational experiment for the paper "On robust recovery of signals from indirect observations"

Uses CVX and Mosek optimization solver

Parameters of the simulation may be modified in the principal script is p2025sim.m
Choices of the identity nuisance matrix N (cf. the paper) and of the second conrast H'=I_m cannot be changed
Other parameters: problem dimensions m and p, noise std, confidence parameter eps, parameter L of the signal class, and nuisance sparsity s may be modified at will.

Solving hundreds of optimization problems takes time. Be patient...

by A. Juditsky and A. Nemirovski
