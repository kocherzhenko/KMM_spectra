// Matrix Equation Solver for Finite System - (c) Aleksey Kocherzhenko, 2012

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stddef.h>
#include "Initialize.h"
#include "Common.h"
#include "LinEqSolver.h"
#include "Setup.h"

#define PI 3.14159265358

FILE *fp;

// Main program starts here
int main()
{
  complex_t *pMatrA, *pVecX, *pVecB;
  double **pMatrIn, **pLUmatr;
  double **pMatrRhsRe, **pMatrRhsIm;
  double *pRhs, *pSolutions;
  double *pSolRe, *pSolIm;
  int n, *pIndx;
  double d, *pd = &d;
  double determinant, trace;
  int i, j, k; 
  double *pkplus, *pkminus;
  int UNITS;

  printf("Starting finiteEqs...\n");
  
  fp = fopen("Input.in", "r");
  fscanf(fp, " %d ", &UNITS);
  fclose(fp);
  
  n = UNITS*2;

  pMatrIn = dmatrix(1,n, 1,n);
  pLUmatr = dmatrix(1,n, 1,n);
  pMatrRhsRe = dmatrix(1,UNITS, 1,UNITS);
  pMatrRhsIm = dmatrix(1,UNITS, 1,UNITS);
  pRhs = dvector(1,n);
  pSolutions = dvector(1,n);
  pSolRe = dvector(1,UNITS); pSolIm = dvector(1,UNITS);
  pIndx = ivector(1,n);

  pkplus = dvector(1, UNITS);
  pkminus = dvector(1, UNITS);

  printf("Reading AMatrix.in and bVector.in...\n");

// Read in Matrix A
  fp = fopen("AMatrix.in", "r");
  for (i=1; i<=n; i+=2) { for (j=2; j<=n; j+=2) { 
    fscanf(fp, "%lf  %lf", &pMatrIn[i][j-1], &pMatrIn[i][j]);
    pMatrIn[i][j] = -pMatrIn[i][j];
    pMatrIn[i+1][j-1] = -pMatrIn[i][j]; pMatrIn[i+1][j] = pMatrIn[i][j-1];
  } }
  fclose(fp);

// Read in ALL Vectors B
  fp = fopen("bVector.in", "r");
  for (i=1; i<=UNITS; i++) { for (j=1; j<=UNITS; j++) { 
    fscanf(fp, "%lf  %lf", &pMatrRhsRe[i][j], &pMatrRhsIm[i][j]); 
  } }
  fclose(fp);

// Set up alpha and beta values
  i = AlphaBetaVals(pkplus, pkminus, UNITS);
//  for (i=1; i<=UNITS; i++) { printf(" %10.4f  %10.4f\n", pkplus[i], pkminus[i]); }

   printf("Factorizing matrix A = L*U... \n");

// Save input matrix pMatrIn, do all transformations for pLUmatr instead!
  for (i=1; i<=n; i++) { for (j=1; j<=n; j++) { pLUmatr[i][j] = pMatrIn[i][j]; } }

// Do LU factorization of pLUmatr
  i = ludcmp(pLUmatr, n, pIndx, pd);


// Check that we have 2M independent equations
  determinant = 1.0; for (i=1; i<=n; i++) { determinant *= pLUmatr[i][i]; }
  printf("Determinant = %10.5f \n", determinant);
  if (fabs(determinant) <= 1.e-10) { printf("Error: determinant is zero!\n"); return 1; }

// Solve set of eqs. sum_(alpha1) (A_{alpha1,beta1} * F2(alpha1,alpha2)) = B(beta1, alpha2)
// with alpha2 - fixed parameter, for all values of alpha2!
  printf("Solving set of linear equations... \n");
  fp = fopen("FinKsVsAlpha.out", "w");
  fclose(fp);

  for (k=1; k<=UNITS; k++)
  {

// Read in Vector B
    for (i=1; i<=UNITS; i++) 
    { 
      pRhs[2*i-1] = pMatrRhsRe[k][i];
      pRhs[2*i] = pMatrRhsIm[k][i]; 
    }

// Print vector pRhs
//    printf("Rhs: \n"); for (i=1; i<=n; i++) { printf("%20.10f  ", pRhs[i]); if (i % 2 == 0) {printf("\n");} } printf("\n"); 

// Save input vector pRhs, do all transformations for pSolutions instead!
    for (i=1; i<=n; i++) { pSolutions[i] = pRhs[i]; } 

// Solve system of linear  
    i = lubksb(pLUmatr, n, pIndx, pSolutions);

// Write solutions of the (complex) system of equations to file "Solutions.out"
// Four columns: alpha1, alpha2, Re(K_M), Im(K_M)
    fp = fopen("FinKsVsAlpha.out", "a");
    for (i=1; i<=UNITS; i++)
    {
      pSolRe[i] = pSolutions[2*i-1];
      pSolIm[i] = pSolutions[2*i];
//      printf("%10.5f   %10.5f \n", pSolRe[i], pSolIm[i]);
      fprintf(fp, "%20.10f   %20.10f  %20.10f  %20.10f \n", pkplus[i], pkplus[k], pSolRe[i], pSolIm[i]);
    }
    fprintf(fp, "\n");
    fclose(fp);
  
  }
  
  printf("FinKsVsAlpha.out successfully written...\n");

// End of main calculation here
// *****************************************************************************************************************************


// Compare to infinite size limit
  double **pInfKsRe, **pInfKsIm;
  double norm = 0;
  double temp1, temp2;

  printf("Comparing solutions to infinite volume limit... \n");

  pInfKsRe = dmatrix(1,UNITS, 1,UNITS);
  pInfKsIm = dmatrix(1,UNITS, 1,UNITS);
  
// Read in infinite size Ks
  fp = fopen("InfKsVsAlpha.out", "r");
  for (i=1; i<=UNITS; i++) { for (j=1; j<=UNITS; j++) { 
    fscanf(fp, "%lf %lf %lf  %lf", &temp1, &temp2, &pInfKsRe[i][j], &pInfKsIm[i][j]); 
//    printf("%f %f %f %f\n", temp1, temp2, pInfKsRe[i][j], pInfKsIm[i][j]);
  } }
  fclose(fp);

//  printf("\n");

// Read in UNITS size Ks
  fp = fopen("FinKsVsAlpha.out", "r");
  for (i=1; i<=UNITS; i++) { for (j=1; j<=UNITS; j++) { 
    fscanf(fp, "%lf %lf %lf %lf", &temp1, &temp2, &pMatrRhsRe[i][j], &pMatrRhsIm[i][j]); 
//    printf("%f %f %f %f\n", temp1, temp2, pMatrRhsRe[i][j], pMatrRhsIm[i][j]);
  } }
  fclose(fp);

// Compare infinite size and UNITS-size Ks using the metric: M^(-2) sum_(alpha1,alpha2) |K_inf(alpha1,alpha2) - K_M(alpha1,alpha2)|^2
  norm = 0.0;
  for (i=1; i<=UNITS; i++) { for (j=1; j<=UNITS; j++) { 
    pInfKsRe[i][j] -= pMatrRhsRe[i][j]; pInfKsIm[i][j] -= pMatrRhsIm[i][j];
    pInfKsRe[i][j] *= pInfKsRe[i][j]; pInfKsIm[i][j] *= pInfKsIm[i][j];
    pInfKsRe[i][j] += pInfKsIm[i][j]; 
//    printf(" %20.10f %20.10f \n", pInfKsRe[i][j], pInfKsIm[i][j]);
    norm += pInfKsRe[i][j]; 
  } }
  norm /= (UNITS*UNITS);
  printf("M^(-2) sum_{alpha1,alpha2} |K_inf(alpha1,alpha2)-K_M(alpha1,alpha2)|^2 = %20.10f \n", norm);
  
  free_dmatrix(pInfKsRe, 1,UNITS, 1,UNITS); 
  free_dmatrix(pInfKsIm, 1,UNITS, 1,UNITS); 


  free_dvector(pkplus, 1,UNITS);
  free_dvector(pkminus, 1,UNITS);
  free_dmatrix(pMatrIn, 1,n, 1,n);
  free_dmatrix(pLUmatr, 1,n, 1,n);
  free_dmatrix(pMatrRhsRe, 1,UNITS, 1,UNITS);
  free_dmatrix(pMatrRhsIm, 1,UNITS, 1,UNITS);
  free_dvector(pRhs, 1,n);
  free_dvector(pSolutions, 1,n);
  free_dvector(pSolRe, 1, UNITS);
  free_dvector(pSolIm, 1, UNITS);
  free_ivector(pIndx, 1,n);

  printf("finiteEqs successfully finished! \n\n");

  return 0;
}

