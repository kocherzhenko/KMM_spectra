// Parameter Calculation for the Weak Coupling Regime (EPSILON > |2*COUP|) - (c) Aleksey Kocherzhenko, 2012

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stddef.h>
#include "Initialize.h"
#include "Common.h"
#include "Setup.h"

FILE *fp;

// Get Js (K1s) for the finite-size (UNITS) system
complex_t getJ(double **pKsRe, double **pKsIm, double *pkplus, int j, int UNITS, double COUP)
{
  int i;
  double k1;
  complex_t temp1, temp2, theJ;
  
  theJ.re = 0.0; theJ.im = 0.0;

  for (i=1; i<=UNITS; i++)
  {
    k1 = pkplus[i];
    temp1 = getTH(k1, COUP);
    temp2.re = pKsRe[j][i]; temp2.im = pKsIm[j][i];
    temp1 = prod(temp1,temp2);
    theJ = add(theJ, temp1);
  }

  k1 = UNITS; theJ.re /= k1; theJ.im /= k1;
  k1 = pkplus[j]; temp1 = getTH(k1, COUP); temp1 = inv(temp1);
  theJ = add(theJ, temp1);
  
  return theJ;
}


// Get K(alpha1,alpha2,alpha3)
complex_t getK3(int i, int j, int k, double *pJre, double *pJim, double **pKre, double **pKim)
{ 
  complex_t temp1, temp2;
  complex_t K3V;
  
  K3V.re = K3V.im = 0.0;
  temp1.re = pJre[i]; temp1.im = pJim[i];
  temp2.re = pKre[j][k]; temp2.im = pKim[j][k];
  temp1 = prod(temp1, temp2);
  K3V = add(K3V, temp1);
  temp1.re = pJre[j]; temp1.im = pJim[j];
  temp2.re = pKre[i][k]; temp2.im = pKim[i][k];
  temp1 = prod(temp1, temp2);
  temp1 = neg(temp1);
  K3V = add(K3V, temp1);
  temp1.re = pJre[k]; temp1.im = pJim[k];
  temp2.re = pKre[i][j]; temp2.im = pKim[i][j];
  temp1 = prod(temp1, temp2);
  K3V = add(K3V, temp1);
  
  return K3V;
}



// Main program starts here
int main()
{
  double **pKsFinRe, **pKsFinIm;
  double **pKsInfRe, **pKsInfIm;
  double *pkplus, *pkminus;
  double *pJsFinRe, *pJsFinIm;
  double *pJsInfRe, *pJsInfIm;
  complex_t theJ;
  complex_t K123;
  double norm;

  int n;
  double temp1, temp2;
  int i, j, k; 
  int UNITS;
  double COUP, A;

  printf("Starting weakParams...\n");

  fp = fopen("Input.in", "r");
  fscanf(fp, " %d ", &UNITS);
  fscanf(fp, " %lf ", &COUP);
  fclose(fp);

  A = EPSILON/2.0/COUP;
  n = UNITS*2;

  if (fabs(A) > 1)
  {
// Set up matrices/vectors
    pKsFinRe = dmatrix(1,UNITS, 1,UNITS); pKsFinIm = dmatrix(1,UNITS, 1,UNITS);
    pKsInfRe = dmatrix(1,UNITS, 1,UNITS); pKsInfIm = dmatrix(1,UNITS, 1,UNITS);
    pJsFinRe = dvector(1, UNITS); pJsFinIm = dvector(1, UNITS);
    pJsInfRe = dvector(1, UNITS); pJsInfIm = dvector(1, UNITS);
    pkplus = dvector(1, UNITS); pkminus = dvector(1, UNITS);

    printf("Reading FinKsVsAlpha.out and InfKsVsAlpha.out...\n");
// Read in K(alpha1,alpha2) Re and Im parts (UNITS size)
    fp = fopen("FinKsVsAlpha.out", "r");
    for (i=1; i<=UNITS; i++) { for (j=1; j<=UNITS; j++) { fscanf(fp, "%lf  %lf  %lf  %lf", &temp1, &temp2, &pKsFinRe[i][j], &pKsFinIm[i][j]); } }
    fclose(fp);

// Read in K(alpha1,alpha2) Re and Im parts (infinite)
    fp = fopen("InfKsVsAlpha.out", "r");
    for (i=1; i<=UNITS; i++) { for (j=1; j<=UNITS; j++) { fscanf(fp, "%lf  %lf  %lf  %lf", &temp1, &temp2, &pKsInfRe[i][j], &pKsInfIm[i][j]); } }
    fclose(fp);

// Set up alpha and beta values
    i = AlphaBetaVals(pkplus, pkminus, UNITS);
 
// Calculate J(e^{i*alpha})
// 3-column file: alpha, Re(J), Im(J)
    printf("Calculating Js (K1s) for finite size system...\n");
    fp = fopen("FinJsVsAlpha.out", "w");
    for (i=1; i<=UNITS; i++) 
    { 
      theJ = getJ(pKsFinRe, pKsFinIm, pkplus, i, UNITS, COUP); 
      pJsFinRe[i] = theJ.re; pJsFinIm[i] = theJ.im; 
      fprintf(fp, " %20.10f  %20.10f  %20.10f \n", pkplus[i], pJsFinRe[i], pJsFinIm[i]); 
    }
    fclose(fp);
    printf("FinJsVsAlpha.out successfully written...\n");


// Compare J_UNITS(alpha) to infinite volume limit  
    printf("Comparing Js (K1s) to infinite volume limit...\n");

// Read in J_inf(alpha) Re and Im parts (infinite)
    printf("Reading InfJsVsAlpha.out...\n");
    fp = fopen("InfJsVsAlpha.out", "r");
    for (i=1; i<=UNITS; i++) { fscanf(fp, "%lf  %lf  %lf", &temp1, &pJsInfRe[i], &pJsInfIm[i]); }
    fclose(fp);

// Calculate metric M^(-1)*sum_alpha(|J_inf(alpha)-J_fin(alpha)|^2)
    norm = 0.0; 
    for (i=1; i<=UNITS; i++)
    {
      temp1 = pJsInfRe[i]-pJsFinRe[i];
      temp2 = pJsInfIm[i]-pJsFinIm[i];
      temp1 = temp1*temp1 + temp2*temp2;
      norm += temp1;
    }
    norm /= UNITS; 
    printf("M^(-1) sum_{alpha} |J_inf(alpha)-J_M(alpha)|^2 = %20.10f\n", norm);
    
    printf("Calculating K123 for finite size system...\n");
// Calculate K123(k1, k2, k3) for finite size (UNITS) system
    fp = fopen("FinK123s.out", "w");
    for (i=1; i<=UNITS; i++) { for (j=1; j<=UNITS; j++) { for (k=1; k<=UNITS; k++) {
      K123 = getK3(i, j, k, pJsFinRe, pJsFinIm, pKsFinRe, pKsFinIm);
      fprintf(fp, " %12.8f  %12.8f  %12.8f  %12.8f  %12.8f \n", pkplus[i], pkplus[j], pkplus[k], K123.re, K123.im);
    } fprintf(fp, "\n"); } fprintf(fp, "\n"); }
    fclose(fp);
    printf("FinK123s.out successfully written...\n");

    printf("Calculating K123 for infinite size system...\n");
// Calculate K123(k1, k2, k3) for infinite size system
    fp = fopen("InfK123s.out", "w");
    for (i=1; i<=UNITS; i++) { for (j=1; j<=UNITS; j++) { for (k=1; k<=UNITS; k++) {
      K123 = getK3(i, j, k, pJsInfRe, pJsInfIm, pKsInfRe, pKsInfIm);
      fprintf(fp, " %12.8f  %12.8f  %12.8f  %12.8f  %12.8f \n", pkplus[i], pkplus[j], pkplus[k], K123.re, K123.im);
    } fprintf(fp, "\n"); } fprintf(fp, "\n"); }
    fclose(fp);
    printf("InfK123s.out successfully written...\n");
    
    
// Free memory
    free_dmatrix(pKsFinRe, 1,UNITS, 1,UNITS); free_dmatrix(pKsFinIm, 1,UNITS, 1,UNITS); 
    free_dmatrix(pKsInfRe, 1,UNITS, 1,UNITS); free_dmatrix(pKsInfIm, 1,UNITS, 1,UNITS); 
    free_dvector(pJsFinRe, 1,UNITS); free_dvector(pJsFinIm, 1,UNITS);
    free_dvector(pJsInfRe, 1,UNITS); free_dvector(pJsInfIm, 1,UNITS);
    free_dvector(pkplus, 1,UNITS); free_dvector(pkminus, 1,UNITS);
  }
  else { printf("Weak coupling regime parameter calculation not needed for strong coupling case! \n"); }

  printf("weakParams successfully finished! \n\n");

  return 0;
}

