// Calculation of 4-excitation matrix elements for the Strong Coupling Regime (EPSILON < |2*COUP|) - (c) Aleksey Kocherzhenko, 2014

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stddef.h>
#include "Initialize.h"
#include "Common.h"
#include "LinEqSolver.h"
#include "Setup.h"

FILE *fp, *fp1, *fp2;

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


// Get vacuum energy
double getVac(double *pk, int UNITS, double COUP)
{ unsigned int i; double Evac = 0.0; double temp1, temp2; double k;
  for (i=1; i<=UNITS; i++) { k = pk[i]; temp1 = EPSILON + 2*COUP*cos(k); temp2 = sqrt( temp1*temp1 + 4*COUP*COUP*sin(k)*sin(k) ); Evac -= (temp2 - temp1)/2.0; }
  return Evac;
}


// Main program starts here
int main()
{
  complex_t **pKAAs;
  complex_t **pLHS, *pRHS;
  double *pkplus, *pkminus;
  double **pMatrIn, **pLUmatr;
  complex_t ctemp;
  complex_t ctemp1, ctemp2, ctemp3;
  double temp, temp1, temp2, temp3, temp4;
  int i, j, k, l;
  int ia2, ia3, ia4;

  int n, *pIndx; 
  double d, *pd = &d;
  double determinant, trace;
  double **pMatrRhsRe, **pMatrRhsIm;
  double *pRhs, *pSolutions;
  double *pSolRe, *pSolIm;

  double VacOverlap;
  double EvacP, EvacM;
  double Energy, Momentum;
  
  long totalnumber = 0;
  double sum = 0.0;
  int UNITS;
  double COUP, A;
  int FOURTHORDER;

  printf("Starting higherStrong...\n");

  fp = fopen("Input.in", "r");
  fscanf(fp, " %d ", &UNITS);
  fscanf(fp, " %lf ", &COUP);
  fscanf(fp, " %d", &FOURTHORDER);
  fclose(fp);

  A = EPSILON/2.0/COUP;

  if ((fabs(A) < 1) && (FOURTHORDER == 1))
  {
    pKAAs = (complex_t **) malloc(UNITS*sizeof(complex_t *));
    for (i=0; i<UNITS; i++) { pKAAs[i] = (complex_t *) malloc(UNITS*sizeof(complex_t)); }
    
    pLHS = (complex_t **) malloc(UNITS*sizeof(complex_t *));
    for (i=0; i<UNITS; i++) { pLHS[i] = (complex_t *) malloc(UNITS*sizeof(complex_t)); }
    pRHS = (complex_t *) malloc(UNITS*sizeof(complex_t));
   
    pkplus = dvector(1,UNITS);
    pkminus = dvector(1,UNITS);
    i = AlphaBetaVals(pkplus, pkminus, UNITS);
    
    n = 2*UNITS;
    pMatrIn = dmatrix(1,n, 1,n);
    pLUmatr = dmatrix(1,n, 1,n);
    pIndx = ivector(1,n);

    pMatrRhsRe = dmatrix(1,UNITS, 1,UNITS);
    pMatrRhsIm = dmatrix(1,UNITS, 1,UNITS);  
    pRhs = dvector(1,n);
    pSolutions = dvector(1,n);
    pSolRe = dvector(1,UNITS); pSolIm = dvector(1,UNITS);

// Read in vacuum overlap from file
    fp = fopen("Overlap.out", "r");
    fscanf(fp, " %lf ", &VacOverlap);
    fclose(fp);

// Get vacuum energies
    EvacP = getVac(pkminus, UNITS, COUP); EvacM = getVac(pkplus, UNITS, COUP);
    printf(" Minus vacuum = %25.20f \n Plus vacuum  = %25.20f \n", EvacM, EvacP);

// Read in two-excitation matrix elements 
    printf("Reading 2-excitation matrix elements from FinKsVsAlpha.out ...\n");
    fp = fopen("FinKsVsAlpha.out", "r");
    for (i=0; i<UNITS; i++) { for (j=0; j<UNITS; j++) { 
      fscanf(fp, "%*f %*f %lf %lf", &(pKAAs[i][j].re), &(pKAAs[i][j].im));
//      printf(" %20.10f  %20.10f\n", pKAAs[i][j].re, pKAAs[i][j].im);
    } }
    fclose(fp);

// Read in LHS matrix as real matrix of size 2*UNITS x 2*UNITS
    fp = fopen("AMatrix.in", "r");
    for (i=1; i<=n; i+=2) { for (j=2; j<=n; j+=2) { 
      fscanf(fp, "%lf  %lf", &pMatrIn[i][j-1], &pMatrIn[i][j]);
      pMatrIn[i][j] = -pMatrIn[i][j];
      pMatrIn[i+1][j-1] = -pMatrIn[i][j]; pMatrIn[i+1][j] = pMatrIn[i][j-1];
    } }
    fclose(fp);

// Perform LU-factorization for the LHS matrix (expanded as real 2*UNITS x 2*UNITS matrix)
    printf("Factorizing matrix A = L*U... \n");

// Save input matrix pMatrIn, do all transformations for pLUmatr instead!
    for (i=1; i<=n; i++) { for (j=1; j<=n; j++) { pLUmatr[i][j] = pMatrIn[i][j]; } }

// Do LU factorization of pLUmatr
    i = ludcmp(pLUmatr, n, pIndx, pd);

// Check that we have 2M independent equations
    determinant = 1.0; for (i=1; i<=n; i++) { determinant *= pLUmatr[i][i]; }
    printf("Determinant = %10.5f \n", determinant);
    if (fabs(determinant) <= 1.e-10) { printf("Error: determinant is zero!\n"); return 1; }

// Calculate RHS of the equations for four-excitation matrix elements
    printf("Calculating the right-hand-side of the equations for 4-excitation matrix elements ...\n");
    fp = fopen("HigherRHScolumns.in", "w");
    for (ia4=1; ia4<=UNITS; ia4++) { for (ia3=1; ia3<=UNITS; ia3++) { for (ia2=1; ia2<=UNITS; ia2++) {
      for (j=1; j<=UNITS; j++)
      {
        ctemp1 = getRHS(pkplus[ia2], pkminus[j], COUP);
        ctemp1 = prod(ctemp1, pKAAs[ia3-1][ia4-1]);
        ctemp2 = getRHS(pkplus[ia3], pkminus[j], COUP);
        ctemp2 = prod(ctemp2, -1.0);
        ctemp2 = prod(ctemp2, pKAAs[ia2-1][ia4-1]);
        ctemp3 = getRHS(pkplus[ia4], pkminus[j], COUP);
        ctemp3 = prod(ctemp3, pKAAs[ia2-1][ia3-1]);
        ctemp1 = add(ctemp1, ctemp2);
        ctemp1 = add(ctemp1, ctemp3);

        fprintf(fp, " %20.10f  %20.10f  \n", ctemp1.re, ctemp1.im);
        totalnumber++;
        if ((totalnumber) % (UNITS*UNITS*UNITS*UNITS/100) == 0) { printf("%6ld percent completed\n", totalnumber/(UNITS*UNITS*UNITS*UNITS/100)); }
      }
      fprintf(fp, "\n");
    } } }
    fclose(fp);
    printf("File HigherRHScolumns.in successfully written!\n");

// Read in RHS matrix for fixed alpha4, alpha3, varied alpha2, alpha1
// Solve system of linear equations to get K(alpha1, alpha2, alpha3, alpha4)
    printf("Solving linear equations...\n");
    totalnumber = 0; 
    fp = fopen("HigherRHScolumns.in", "r");
    fp1 = fopen("4FinIntens.out", "w");
    fp2 = fopen("FinK1234s.out", "w");
    for (ia4=1; ia4<=UNITS; ia4++) { for (ia3=1; ia3<=UNITS; ia3++)
    {
      for (i=1; i<=UNITS; i++) { for (j=1; j<=UNITS; j++) { 
        fscanf(fp, "%lf  %lf", &pMatrRhsRe[i][j], &pMatrRhsIm[i][j]); 
      } }

      for (k=1; k<=UNITS; k++)
      {
      // Read in Vector B
        for (i=1; i<=UNITS; i++) 
        { 
          pRhs[2*i-1] = pMatrRhsRe[k][i];
          pRhs[2*i] = pMatrRhsIm[k][i]; 
        }
                
      // Save input vector pRhs, do all transformations for pSolutions instead!
        for (i=1; i<=n; i++) { pSolutions[i] = pRhs[i]; } 

      // Solve system of linear  
        i = lubksb(pLUmatr, n, pIndx, pSolutions);

      // Print the calculated matrix elements squared to file "FinK1234.out"
        for (i=1; i<=UNITS; i++)
        {
          pSolRe[i] = pSolutions[2*i-1];
          pSolIm[i] = pSolutions[2*i];
          if ((i<k) && (k<ia3) && (ia3<ia4))
          {
            fprintf(fp2, " %20.10f  %20.10f  %20.10f  %20.10f  %20.10f  %20.10f \n", pkplus[i], pkplus[k], pkplus[ia3], pkplus[ia4], pSolRe[i], pSolIm[i]);

            temp = pSolRe[i]*pSolRe[i] + pSolIm[i]*pSolIm[i];
            temp /= (UNITS*UNITS*UNITS*UNITS)*VacOverlap;

            Energy = getEn(pkplus[i], COUP) + getEn(pkplus[k], COUP) + getEn(pkplus[ia3], COUP) + getEn(pkplus[ia4], COUP) + EvacM - EvacP;
            Momentum = pkplus[i] + pkplus[k] + pkplus[ia3] + pkplus[ia4];
          
            fprintf(fp1, " %+25.20lf   %+25.20lf  %+25.20lf \n", Energy, Momentum, temp);
            sum += temp;
          }
          totalnumber++;
          if ((totalnumber) % (UNITS*UNITS*UNITS*UNITS/100) == 0) { printf("%6ld percent completed\n", totalnumber/(UNITS*UNITS*UNITS*UNITS/100)); }
        }
      }

    } }
    fclose(fp);
    fclose(fp1);
    fclose(fp2);
    
    printf(" %10.6f \n", sum);
    
// Free memory
    for (i=0; i<UNITS; i++) { free(pKAAs[i]); }
    free(pKAAs);
    for (i=0; i<UNITS; i++) { free(pLHS[i]); }
    free(pLHS);
    free(pRHS);
    
    free_dvector(pkplus, 1,UNITS);
    free_dvector(pkminus, 1,UNITS);
    free_dmatrix(pMatrIn, 1,n, 1,n);
    free_dmatrix(pLUmatr, 1,n, 1,n);
    free_ivector(pIndx, 1,n);
    free_dmatrix(pMatrRhsRe, 1,UNITS, 1,UNITS);
    free_dmatrix(pMatrRhsIm, 1,UNITS, 1,UNITS);
    free_dvector(pRhs, 1,n);
    free_dvector(pSolutions, 1,n);
    free_dvector(pSolRe, 1, UNITS);
    free_dvector(pSolIm, 1, UNITS);
  }
  else if (fabs(A) > 1) { printf("Strong coupling regime parameter calculation not needed for weak coupling case! \n"); }
  else if (FOURTHORDER == 0) { printf("4-excitation matrix element calculation is switched off!  \n"); }

  printf("higherStrong successfully finished! \n\n");

  return 0;
}

