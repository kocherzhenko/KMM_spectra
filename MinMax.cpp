// Calculate minimum and maximum values of energy, momentum, and intensity in output data file - (c) Aleksey Kocherzhenko, 2014

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Setup.h"

FILE *fp, *fp1;

int getMinMax(FILE *fp, FILE *fp1, int manifold, int number)
{
  int i;
  double E, k, Int;
  double Emax, Emin, kmax, kmin, IntMin, IntMax;

  Emax = -5000.0; Emin = 5000.0; kmax = -5000.0; kmin = 5000.0; IntMin = 5000.0; IntMax = -5000.0; 

  for (i=0; i<number; i++)
  {
    fscanf(fp, " %lf  %lf  %lf ", &E, &k, &Int);
    if (E < Emin) { Emin = E; }
    if (E > Emax) { Emax = E; }
    if (k < kmin) { kmin = k; }
    if (k > kmax) { kmax = k; }
    if (Int < IntMin) { IntMin = Int; }
    if (Int > IntMax) { IntMax = Int; }
  }

  fprintf(fp1, "%d-exc. manifold: \n", manifold);
  fprintf(fp1, " kmin   = %20.15f ; kmax   = %20.15f ; \n", kmin, kmax);
  fprintf(fp1, " Emin   = %20.15f ; Emax   = %20.15f ; \n", Emin, Emax);
  fprintf(fp1, " IntMin = %20.15f ; IntMax = %20.15f ; \n\n", IntMin, IntMax);

  return 0;
}

int main()
{
  int i; // m, n;
  double E, k, Int;
  double Emax, Emin, kmax, kmin, IntMin, IntMax;
  int number;
  int UNITS;
  double COUP, A;
  int FOURTHORDER;

  printf("Starting minMax...\n");
  
  fp = fopen("Input.in", "r");
  fscanf(fp, " %d ", &UNITS);
  fscanf(fp, " %lf ", &COUP);
  fscanf(fp, " %d ", &FOURTHORDER);
  fclose(fp);

  A = EPSILON/2.0/COUP;
  
  if (fabs(A) > 1.0) // Weak coupling, 1- and 3-exc. manifolds
  { 
    printf("Analyzing 1-exc. manifold...\n"); 
    number = UNITS;
    
    fp = fopen("1FinIntens.out", "r");
    fp1 = fopen("EkIrange.out", "w");
    i = getMinMax(fp, fp1, 1, number);
    fclose (fp);
    fclose(fp1);
    
    printf("Analyzing 3-exc. manifold...\n"); 
    number = UNITS*(UNITS-1)*(UNITS-2)/6;

    fp = fopen("3FinIntens.out", "r");
    fp1 = fopen("EkIrange.out", "a");
    i = getMinMax(fp, fp1, 3, number);
    fclose (fp);
    fclose (fp);
  }
  else // Strong coupling, 2- and 4-exc. manifolds
  {    
    printf("Analyzing 2-exc. manifold...\n"); 
    number = UNITS*(UNITS-1)/2;
    
    fp = fopen("2FinIntens.out", "r");
    fp1 = fopen("EkIrange.out", "w");
    i = getMinMax(fp, fp1, 2, number);
    fclose (fp);
    fclose (fp);

    if (FOURTHORDER == 1)
    {    
      printf("Analyzing 4-exc. manifold...\n"); 
      number = UNITS*(UNITS-1)*(UNITS-2)*(UNITS-3)/24;

      fp = fopen("4FinIntens.out", "r");
      fp1 = fopen("EkIrange.out", "a");
      i = getMinMax(fp, fp1, 4, number);
      fclose (fp);
      fclose (fp);
    }
    else { printf ("Calculation for the 4-excitation manifold is turned off! \n"); }
  }
  
  printf("minMax successfully finished!\n\n");
 
  return 0;
}

