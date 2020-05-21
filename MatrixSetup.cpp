// Setup of the Input Matrix for Linear Equation Solver - (c) Aleksey Kocherzhenko, 2012

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Initialize.h"
#include "Common.h"
#include "Setup.h"

#define DEL 1.e-8

FILE *fp;


// Calculate elements of A matrix (l.h.s. of Eq. (14B))
complex_t getLHS(double alpha, double beta, int UNITS, double COUP)
{
  complex_t lhs, temp1, temp2, temp3;
  
  temp1.re = 0.0; temp1.im = alpha - beta;
  temp1 = cexp(temp1); temp1 = add(temp1, -1.0); temp1 = inv(temp1);
  temp2 = getTH(alpha, COUP); temp3 = getTH(beta, COUP); temp3 = inv(temp3); temp2 = prod(temp2, temp3); temp2 = add(temp2, 1.0);
  lhs = prod(temp1, temp2); lhs = prod(lhs, 1.0/UNITS);
  
  return lhs;
}

// Calculate elements of b vector(s) (r.h.s. of Eq. (14B))
complex_t getRHS(double alpha, double beta, double COUP)
{
  complex_t rhs, temp1, temp2, temp3;

  temp1.re = 0.0; temp1.im = alpha + beta;
  temp1 = cexp(temp1); temp3 = temp1; temp1 = add(temp1, -1.0); temp1 = inv(temp1); temp1 = prod(temp1, temp3);
  temp2 = getTH(alpha, COUP); temp3 = getTH(beta, COUP); temp2 = prod(temp2, temp3); temp2 = inv(temp2); temp2 = add(temp2, -1.0);
  rhs = prod(temp1, temp2);

  return rhs;
}

// Infinite volume params, for comparison
// 1/THETA(+)
complex_t getTHpi(complex_t z, double a, double A) 
{
  complex_t THpi, temp;

// for |2*b/epsilon| < 1
  if (fabs(A) > 1) { THpi = add(z, a); THpi = inv(THpi); THpi = prod(THpi, a); THpi = csqrt(THpi); THpi = inv(THpi); }
// for |2*b/epsilon| > 1
  if (fabs(A) < 1) { THpi = prod(z, a); THpi = add(THpi, 1.0); THpi = csqrt(THpi); THpi = inv(THpi); } 

  return THpi; 
}

// 1/THETA(-)
complex_t getTHmi(complex_t z, double ainv, double A)
{ 
  complex_t THmi, temp;

// for |2*b/epsilon| < 1
  if (fabs(A) > 1) { THmi = add(z, ainv); temp = inv(z); THmi = prod(THmi, temp); THmi = csqrt(THmi); THmi = inv(THmi); }
// for |2*b/epsilon| > 1
  if (fabs(A) < 1) { THmi = add(z, 1.0/ainv); THmi = inv(THmi); THmi = prod(z, THmi); THmi = csqrt(THmi); THmi = inv(THmi); } 

  return THmi;}

// K(z1,z2)
complex_t getK(complex_t z1, complex_t z2, double A)
{ complex_t temp1, temp2; complex_t K1; double a = A; double ainv = 1/A; 
  
//  printf("z1.re = %10.8f; z1.im = %10.8f; z2.re = %10.8f; z2.im = %10.8f; %10.8f\n", z1.re, z1.im, z2.re, z2.im, DEL);
  if ((z1.re == z2.re) && (z1.im == -z2.im)) { z2.re = z2.re*cos(DEL) - z2.im*sin(DEL); z2.im = z2.im*cos(DEL) - z2.re*sin(DEL); }

  temp1 = getTHpi(z1, a, A); temp2 = getTHmi(z2, ainv, A); K1 = prod(temp1, temp2);
  temp1 = getTHmi(z1, ainv, A); temp2 = getTHpi(z2, a, A); temp1 = prod(temp1, temp2); temp1 = neg(temp1); K1 = add(K1, temp1); 
  temp1 = prod(z1, z2); temp2 = add(temp1, -1.0); temp2 = inv(temp2); temp1 = prod(temp1, temp2); K1 = prod(temp1, K1); 

  return K1;
}

// Start of the main program
int main ()
{
  int i, j;
  int n;
  double *pkplus, *pkminus;
  double **pMatrOut;
  double **pInfKsRe, **pInfKsIm;
  
  double alpha;
  double temp;
  complex_t Theta;
  complex_t lhs, rhs;
  complex_t z1, z2;
  complex_t k1, k2;
  double ainv;
  int UNITS;
  double COUP, A;

  printf("Starting matrixSetup...\n");

  fp = fopen("Input.in", "r");
  fscanf(fp, " %d ", &UNITS);
  fscanf(fp, " %lf ", &COUP);
  fclose(fp);
 
  A = EPSILON/2.0/COUP;
  ainv = 1.0/A;
  
  n = UNITS*2;

//  pMatrOut = dmatrix(1,n, 1,n);
  pInfKsRe = dmatrix(1,UNITS, 1,UNITS);
  pInfKsIm = dmatrix(1,UNITS, 1,UNITS);

// Set k values

  pkplus = dvector(1,UNITS);
  pkminus = dvector(1,UNITS);

  i = AlphaBetaVals(pkplus, pkminus, UNITS);
  
// Set up A matrix and write to file
  printf("Setting up A matrix... \n");
  fp = fopen("AMatrix.in", "w");
  for (i=1; i<=UNITS; i++)
  {
    for (j=1; j<=UNITS; j++)
    {
      lhs = getLHS(pkplus[j], pkminus[i], UNITS, COUP);
      fprintf(fp, "%20.10f   %20.10f \n", lhs.re, lhs.im);
    }
    fprintf(fp, "\n");
  }
  fclose(fp);
  printf("AMatrix.in successfully written...\n");

// Set up B vectors (for all alpha2 values!) and write to file
  printf("Setting up b vector... \n");
  fp = fopen("bVector.in", "w");
  for (i=1; i<=UNITS; i++)
  {
    for (j=1; j<=UNITS; j++)
    {
      rhs = getRHS(pkplus[i], pkminus[j], COUP);
      fprintf(fp, "%20.10f   %20.10f \n", rhs.re, rhs.im);
    }
    fprintf(fp, "\n");
  }
  fclose(fp);
  printf("bVector.in successfully written...\n");
 
// Calculate infinite size limit Ks
// Four-column output file: alpha1, alpha2, Re(K_inf), Im(K_inf)
  printf("Calculating infinite volume Ks... \n");
  fp = fopen("InfKsVsAlpha.out", "w");  
  for (i=1; i<=UNITS; i++)
  {
    for (j=1; j<=UNITS; j++)
    {
      k1.re = k1.im = k2.re = k2.im = 0.0;
      z1.re = 0.0; z1.im = pkplus[i];
      z2.re = 0.0; z2.im = pkplus[j];
      z1 = cexp(z1); z2 = cexp(z2); 
      k1 = getK(z1, z2, A); k2 = getK(z2, z1, A);
      fprintf(fp, " %20.10f %20.10f %20.10f %20.10f \n", pkplus[j], pkplus[i], k2.re, k2.im);
    }
    fprintf(fp,"\n");
  }
  fclose(fp);
  printf("InfKsVsAlpha.out successfully written...\n");

// Calculate infinite size limit Js (same as K1!)
// Three-column output file: alpha, Re(J_inf), Im(J_inf)
  if (fabs(A) > 1)
  {
    printf("Calculating infinite volume Js (K1s)... \n");
    fp = fopen("InfJsVsAlpha.out", "w");
    for (i=1; i<=UNITS; i++)
    {
      z1.re = 0.0; z1.im = pkplus[i]; z1 = cexp(z1); 
      k1 = getTHmi(z1, ainv, A);
      fprintf(fp, " %20.10f %20.10f %20.10f \n", pkplus[i], k1.re, k1.im);
    }
    fclose(fp);
    printf("InfJsVsAlpha.out successfully written...\n");
  }

  free_dvector(pkplus, 1,UNITS); free_dvector(pkminus, 1,UNITS);

  printf("matrixSetup successfully finished!\n\n");

  return 0;
}

