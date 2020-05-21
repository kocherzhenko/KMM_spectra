// Beyond KMM Eigensolver - (c) Aleksey Kocherzhenko, 2011

// Transformation of Hamiltonian to 3-diagonal form by a sequence of (SIZE-2) Householder reflections
// Adapted from "Numerical recepies in C" tred2.c
void hhtridiag (double *pmatr, double *pdiag, double *poffdiag, int UNITS);

// Diagonalization using QL factorization with implicit shifts
// Adapted from "Numerical recepies in C" tqli.c
int qlfact(double* pdiag, double* poffdiag, double *pmatr, int UNITS);

