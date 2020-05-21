// Vacuum overlap calculation - (c) Aleksey Kocherzhenko, 2013

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stddef.h>
#include "Initialize.h"
#include "Common.h"
#include "EigSolverU.h"
#include "Setup.h"

FILE *fp;

// Main program starts here
int main()
{
  double **pKsFinRe, **pKsFinIm;
  double **pExpRe, **pExpIm;
  double **pExpThRe, **pExpThIm;
  double *pkplus, *pkminus;
  double norm;
  
  double *pA, *pDiag, *pOffDiag;
  
  int n;
  double temp1, temp2;
  complex_t ctemp1, ctemp2;
  int i, j, k; 
  int UNITS;
  double COUP;

  printf("Starting vacCalc...\n");
  
  fp = fopen("Input.in", "r");
  fscanf(fp, " %d ", &UNITS);
  fscanf(fp, " %lf ", &COUP);
  fclose(fp);
  
  n = UNITS*2;


//################################################################################
// X(alpha1, alpha2) setup
//################################################################################

// Set up matrices/vectors
  pKsFinRe = dmatrix(1,UNITS, 1,UNITS); pKsFinIm = dmatrix(1,UNITS, 1,UNITS);
  pExpRe = dmatrix(1,UNITS, 1,UNITS); pExpIm = dmatrix(1,UNITS, 1,UNITS);
  pExpThRe = dmatrix(1,UNITS, 1,UNITS); pExpThIm = dmatrix(1,UNITS, 1,UNITS);
  pkplus = dvector(1, UNITS); pkminus = dvector(1, UNITS);

  printf("Reading FinKsVsAlpha.out...\n");
// Read in K(alpha1,alpha2) Re and Im parts (UNITS size)
  fp = fopen("FinKsVsAlpha.out", "r");
  for (i=1; i<=UNITS; i++) { for (j=1; j<=UNITS; j++) { fscanf(fp, "%lf  %lf  %lf  %lf", &temp1, &temp2, &pKsFinRe[i][j], &pKsFinIm[i][j]); } }
  fclose(fp);

// Set up alpha and beta values
  i = AlphaBetaVals(pkplus, pkminus, UNITS);

// WEAK COUPLING REGIME
  if (fabs(2*COUP) < EPSILON)
  {
  // Set up exp[-i*(alpha1+alpha2)/2]
    for (i=1; i<=UNITS; i++) { for (j=1; j<=UNITS; j++) { pExpRe[i][j] = cos((pkplus[i]+pkplus[j])/2.0); pExpIm[i][j] = -sin((pkplus[i]+pkplus[j])/2.0); } }

  // Set up exp[-i*theta(alpha1)]*exp[-i*theta(alpha2)]
    for (i=1; i<=UNITS; i++)
    {
      temp1 = EPSILON*EPSILON+4*COUP*COUP+4*EPSILON*COUP*cos(pkplus[i]);
      temp1 = sqrt(temp1);
      ctemp1.re = (EPSILON+2*COUP*cos(pkplus[i]))/temp1;
      ctemp1.im = -2*COUP*sin(pkplus[i])/temp1;
    
      for (j=1; j<=UNITS; j++)
      {
        temp2 = EPSILON*EPSILON+4*COUP*COUP+4*EPSILON*COUP*cos(pkplus[j]);
        temp2 = sqrt(temp2);
        ctemp2.re = (EPSILON+2*COUP*cos(pkplus[j]))/temp2;
        ctemp2.im = -2*COUP*sin(pkplus[j])/temp2;
     
        ctemp2 = prod(ctemp1, ctemp2);
        ctemp2 = csqrt(ctemp2);
         
        pExpThRe[i][j] = ctemp2.re; 
        pExpThIm[i][j] = ctemp2.im; 
      }
    } 

  // Calculate X(alpha1, alpha2)/F(0)
    printf("Setting up X(alpha1, alpha2) matrix...\n");
    for (i=1; i<=UNITS; i++)
    {
      for (j=1; j<=UNITS; j++) 
      {
        ctemp1.re = pExpRe[i][j]; ctemp1.im = pExpIm[i][j];
        ctemp2.re = pExpThRe[i][j]; ctemp2.im = pExpThIm[i][j];
        ctemp1 = prod(ctemp1, ctemp2);
 
        ctemp2.re = pKsFinRe[i][j]; ctemp2.im = pKsFinIm[i][j];
        ctemp1 = prod(ctemp1, ctemp2);

        ctemp2.re = 0.0; ctemp2.im = 1.0;
        ctemp1 = prod(ctemp1, ctemp2);
      
        ctemp1.re /= UNITS; ctemp1.im /= UNITS;
        pKsFinRe[i][j] = ctemp1.re; pKsFinIm[i][j] = ctemp1.im;
      }
    }
  }
  
// STRONG COUPLING REGIME
  else
  {
  // Set up exp[-i*theta(alpha1)]*exp[-i*theta(alpha2)]
    for (i=1; i<=UNITS; i++)
    {
      temp1 = EPSILON*EPSILON+4*COUP*COUP+4*EPSILON*COUP*cos(pkplus[i]);
      temp1 = sqrt(temp1);
      ctemp1.re = (EPSILON+2*COUP*cos(pkplus[i]))/temp1;
      ctemp1.im = -2*COUP*sin(pkplus[i])/temp1;
    
      for (j=1; j<=UNITS; j++)
      {
        temp2 = EPSILON*EPSILON+4*COUP*COUP+4*EPSILON*COUP*cos(pkplus[j]);
        temp2 = sqrt(temp2);
        ctemp2.re = (EPSILON+2*COUP*cos(pkplus[j]))/temp2;
        ctemp2.im = -2*COUP*sin(pkplus[j])/temp2;
     
        ctemp2 = prod(ctemp1, ctemp2);
        ctemp2 = csqrt(ctemp2);
         
        pExpThRe[i][j] = ctemp2.re; 
        pExpThIm[i][j] = ctemp2.im; 
      }
    } 

  // Calculate X(alpha1, alpha2)/F(0)
    printf("Setting up X(alpha1, alpha2) matrix...\n");
    for (i=1; i<=UNITS; i++)
    {
      for (j=1; j<=UNITS; j++) 
      {
        ctemp1.re = pExpThRe[i][j]; ctemp1.im = pExpThIm[i][j];
        ctemp2.re = pKsFinRe[i][j]; ctemp2.im = pKsFinIm[i][j];
        ctemp1 = prod(ctemp1, ctemp2);

        ctemp2.re = 0.0; ctemp2.im = 1.0;
        ctemp1 = prod(ctemp1, ctemp2);
      
        ctemp1.re /= UNITS; ctemp1.im /= UNITS;
        pKsFinRe[i][j] = ctemp1.re; pKsFinIm[i][j] = ctemp1.im;
      }
    }    
  }
    
//################################################################################
// X(alpha1, alpha2) eigenvalue/eigenvector calculation
//################################################################################

// Set up matrices for eigensolver 
  pA = dvector(0,4*UNITS*UNITS-1);
  pDiag = dvector(0, 2*UNITS-1);
  pOffDiag = dvector(0, 2*UNITS-1);
 
// Input matrix [ 0  -X ]
//              [ X   0 ] setup 
  for (i=0; i<UNITS; i++)
  {
    for (j=0; j<UNITS; j++) { pA[i*2*UNITS+j] = 0.0; } 
    for (j=UNITS; j<2*UNITS; j++) { pA[i*2*UNITS+j] = -pKsFinRe[i+1][j-UNITS+1]; }
  }
  for (i=UNITS; i<2*UNITS; i++)
  {
    for (j=0; j<UNITS; j++) { pA[i*2*UNITS+j] = pKsFinRe[i-UNITS+1][j+1]; }
    for (j=UNITS; j<2*UNITS; j++) { pA[i*2*UNITS+j] = 0.0; }
  }

// Diagonal/offdiagonal element vectors initially set to zero
  for (i=0; i<2*UNITS; i++) { pDiag[i] = 0.0; pOffDiag[i] = 0.0; }

// Eigenvalue/eigenvector calculation
  printf("Calculating the eigenvalues of X(alpha1, alpha2)\n");
  hhtridiag (pA, pDiag, pOffDiag, UNITS);
  qlfact(pDiag, pOffDiag, pA, UNITS);

// Write eigenvalues to "Xeigvals.out"
  printf("Writing the imaginary part of X(alpha1, alpha2) eigenvalues to Xeigvals.out\n");
  fp = fopen("Xeigvals.out", "w");
  for (i=0; i<2*UNITS; i+=2) { fprintf(fp, " %10.6f \n", pDiag[i]); } fprintf(fp, "\n");
  fclose(fp);
  
// Calculate |<PHI{-}|PHI{+}>|^2
  if (fabs(2*COUP) < EPSILON) { printf("Calculating <PHI{-}|PHI{+}>...\n"); }
  else { printf("Calculating <PHI{-}|sigma(1,x)|PHI{+}>...\n"); }
  norm = 1;
  for (i=0; i<2*UNITS; i++) { norm *= (1+pDiag[i]*pDiag[i]); }
  norm = sqrt(norm); norm = sqrt(norm);
  norm = 1/norm;

  if (fabs(2*COUP) < EPSILON) { printf("|<PHI{-}|PHI{+}>|^2 = %10.5f \n", norm); }
  else { printf("|<PHI{-}|sigma(1,x)|PHI{+}>|^2 = %10.5f \n", norm); }
  
// Write |<PHI{-}|PHI{+}>|^2 to "Overlap.out"
  if (fabs(2*COUP) < EPSILON) { printf("Writing |<PHI{-}|PHI{+}>|^2 to Overlap.out\n"); }
  else { printf("Writing |<PHI{-}|sigma(1,x)|PHI{+}>|^2 to Overlap.out\n"); }
  fp = fopen("Overlap.out", "w");
  fprintf(fp, "%lf", norm); 
  fclose(fp);
    
// Free memory
  free_dmatrix(pKsFinRe, 1,UNITS, 1,UNITS); free_dmatrix(pKsFinIm, 1,UNITS, 1,UNITS); 
  free_dmatrix(pExpRe, 1,UNITS, 1,UNITS); free_dmatrix(pExpIm, 1,UNITS, 1,UNITS);
  free_dmatrix(pExpThRe, 1,UNITS, 1,UNITS); free_dmatrix(pExpThIm, 1,UNITS, 1,UNITS);
  free_dvector(pkplus, 1,UNITS); free_dvector(pkminus, 1,UNITS);

  free_dvector(pA, 0,2*UNITS-1); free_dvector(pDiag, 0,2*UNITS-1); free_dvector(pOffDiag, 0,2*UNITS-1);

  printf("vacCalc successfully finished! \n\n");

  return 0;
}

