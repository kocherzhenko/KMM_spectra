// Matrix Equation Solver for Finite System - (c) Aleksey Kocherzhenko, 2012

// Function to perform LU decomposition of a rowwise permutation of a[1..n][1..n]
// Adapted from "Numerical Recipies in C"
int ludcmp(double **a, int n, int *indx, double *d);

// Function to solve a set of linear equations
// Adapted from "Numerical Recipies in C"
int lubksb(double **a, int n, int *indx, double b[]);

