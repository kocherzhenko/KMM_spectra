// Program to convert scatter data to grid data - (c) Aleksey Kocherzhenko, 2014

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Setup.h"

#define KMIN -10.0
#define KMAX +10.0
#define EMIN   0.0
#define EMAX   6.0

FILE *fp;

// Get the intensity and number of peaks for every (E, E+dE; k, k+dk) bin
int getGridData(FILE *fp, int **ppCount, double **ppIntens, int number, double stepx, double stepy)
{
  int i, j, m, n, oor = 0;
  double E, k, Int;
  double x, y;

  for (i=0; i<XPOINTS; i++) { for (j=0; j<YPOINTS; j++) { ppCount[i][j] = 0; ppIntens[i][j] = 0.0; } } 

  for (i=0; i<number; i++)
  {
    if ((i % 100000) == 0) printf(" %d of %d lines read from data file... \n", i, number);
    fscanf(fp, " %lf  %lf  %lf ", &E, &k, &Int);
    x = KMIN+stepx; m = 0;
    while (k > x) { x += stepx; m++; }
    y = EMIN+stepy; n = 0;
    while (E > y) { y += stepy; n++; }

    if ((m > XPOINTS) || (n > YPOINTS))
    { 
      oor++;
//      printf("%10.5f %10.5f : %20.10f \n", E, k, Int); 
    } // Point out of range (E or k)
    else
    {
      ppCount[m][n] += 1;
      ppIntens[m][n] += Int; 
    }
  }

  if (oor > 0) { printf(" %d points out of range (E or k)\n", oor); } // Total # of points out of range (E or k)
  
  return 0;
}

int main()
{
  int i, j, m, n;
  double x, y;
  double stepx = (KMAX-KMIN)/XPOINTS;
  double stepy = (EMAX-EMIN)/YPOINTS;
  double E, k, Int;
  int **ppCount;
  double **ppIntens;
  
  double temp, kmin, kmax, Emin, Emax;
  char c;
  int number;
  int UNITS;
  double COUP, A;

  int FOURTHORDER;

  printf("Starting gridConvert...\n"); 

  fp = fopen("Input.in", "r");
  fscanf(fp, " %d ", &UNITS);
  fscanf(fp, " %lf ", &COUP);
  fscanf(fp, " %d ", &FOURTHORDER);
  fclose(fp);

  A = EPSILON/2.0/COUP;

  ppCount = (int **) malloc(XPOINTS*sizeof(int *));
  for (i=0; i<XPOINTS; i++) { ppCount[i] = (int *) malloc(YPOINTS*sizeof(int)); } 
  
  ppIntens = (double **) malloc(XPOINTS*sizeof(double *));
  for (i=0; i<XPOINTS; i++) { ppIntens[i] = (double *) malloc(YPOINTS*sizeof(double)); } 
  
  printf(" %10.5f  %10.5f  \n", stepx, stepy);
 
// Weak coupling regime 
  if (fabs(A) > 1)
  {
  // 1-exc. manifold
    printf("Analyzing the 1-excitation manifold...\n");
    number = UNITS;
    fp = fopen("1FinIntens.out", "r");
    if (fp == NULL) { printf("Cannot open file 1FinIntens.out! Exiting... \n"); return 1; }
    i = getGridData(fp, ppCount, ppIntens, number, stepx, stepy);
    if (i != 0) { printf("Error in getGridData. Exiting...\n"); return 1; }
    fclose (fp);
    
    fp = fopen("ZGrid_1exc.out", "w");
    for (i=0; i<XPOINTS; i++) { for (j=0; j<YPOINTS; j++) {
      if (ppCount[i][j] == 0) { ppIntens[i][j] = 0.0; } 
//      else { ppIntens[i][j] /= ppCount[i][j]; }
      fprintf(fp, " %15.8f  %15.8f  %25.20f \n", KMIN+i*stepx+stepx/2.0, EMIN+j*stepy+stepy/2.0, ppIntens[i][j]/(stepx*stepy));
    } fprintf(fp, "\n"); }
    fclose (fp);

  // 3-exc. manifold
    printf("Analyzing the 3-excitation manifold...\n");
    number = UNITS*(UNITS-1)*(UNITS-2)/6;
    fp = fopen("3FinIntens.out", "r");
    if (fp == NULL) { printf("Cannot open file 3FinIntens.out! Exiting... \n"); return 1; }
    i = getGridData(fp, ppCount, ppIntens, number, stepx, stepy);
    if (i != 0) { printf("Error in getGridData. Exiting...\n"); return 1; }
    fclose (fp);

    fp = fopen("ZGrid_3exc.out", "w");
    for (i=0; i<XPOINTS; i++) { for (j=0; j<YPOINTS; j++) {
      if (ppCount[i][j] == 0) { ppIntens[i][j] = 0.0; } 
//      else { ppIntens[i][j] /= ppCount[i][j]; }
      fprintf(fp, " %15.8f  %15.8f  %25.20f \n", KMIN+i*stepx+stepx/2.0, EMIN+j*stepy+stepy/2.0, ppIntens[i][j]/(stepx*stepy));
    } fprintf(fp, "\n"); }
    fclose (fp);
  }
// Strong coupling regime 
  else
  {
  // 2-exc. manifold
    printf("Analyzing the 2-excitation manifold...\n");
    number = UNITS*(UNITS-1)/2;
    fp = fopen("2FinIntens.out", "r");
    if (fp == NULL) { printf("Cannot open file 2FinIntens.out! Exiting... \n"); return 1; }
    i = getGridData(fp, ppCount, ppIntens, number, stepx, stepy);
    if (i != 0) { printf("Error in getGridData. Exiting...\n"); return 1; }
    fclose (fp);
    
    fp = fopen("ZGrid_2exc.out", "w");
    for (i=0; i<XPOINTS; i++) { for (j=0; j<YPOINTS; j++) {
      if (ppCount[i][j] == 0) { ppIntens[i][j] = 0.0; } 
//      else { ppIntens[i][j] /= ppCount[i][j]; }
      fprintf(fp, " %15.8f  %15.8f  %25.20f \n", KMIN+i*stepx+stepx/2.0, EMIN+j*stepy+stepy/2.0, ppIntens[i][j]/(stepx*stepy));
    } fprintf(fp, "\n"); }
    fclose (fp);

  // 4-exc. manifold
    number = UNITS*(UNITS-1)*(UNITS-2)*(UNITS-3)/24;
    if (FOURTHORDER == 1)
    {
      printf("Analyzing the 4-excitation manifold...\n");
      fp = fopen("4FinIntens.out", "r");
      if (fp == NULL) { printf("Cannot open file 4FinIntens.out! Exiting... \n"); return 1; }
      i = getGridData(fp, ppCount, ppIntens, number, stepx, stepy);
      if (i != 0) { printf("Error in getGridData. Exiting...\n"); return 1; }
      fclose (fp);

      fp = fopen("ZGrid_4exc.out", "w");
      for (i=0; i<XPOINTS; i++) { for (j=0; j<YPOINTS; j++) {
        if (ppCount[i][j] == 0) { ppIntens[i][j] = 0.0; } 
//        else { ppIntens[i][j] /= ppCount[i][j]; }
        fprintf(fp, " %15.8f  %15.8f  %25.20f \n", KMIN+i*stepx+stepx/2.0, EMIN+j*stepy+stepy/2.0, ppIntens[i][j]/(stepx*stepy));
      } fprintf(fp, "\n"); }
      fclose (fp);
    }
    else { printf("Calculation for the 4-excitation manifold is turned off! \n"); }
  }
  
  for (i=0; i<XPOINTS; i++) { free(ppCount[i]); } 
  free(ppCount);
  for (i=0; i<XPOINTS; i++) { free(ppIntens[i]); }
  free(ppIntens); 
  
  printf("gridConvert successfully finished!\n\n");

  return 0;
}

