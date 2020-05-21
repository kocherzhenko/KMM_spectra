// Line Intensity Calculation - (c) Aleksey Kocherzhenko, 2012

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Initialize.h"
#include "Common.h"
#include "Setup.h"

FILE *fp;

// Get vacuum energy
double getVac(double *pk, int UNITS, double COUP)
{ unsigned int i; double Evac = 0.0; double temp1, temp2; double k;
  for (i=1; i<=UNITS; i++) { k = pk[i]; temp1 = EPSILON + 2*COUP*cos(k); temp2 = sqrt( temp1*temp1 + 4*COUP*COUP*sin(k)*sin(k) ); Evac -= (temp2 - temp1)/2.0; }
  return Evac;
}

// Bubble sort function; this is a (SIZE^2) algorithm: can be replaced by a more efficient one!
void bubblesort(double *pe, double *pk, double *pi, unsigned int s)
{
  unsigned int i,j;
  double temp;
    
  for (i=s-1; i>0; i--) { for (j=1; j<=i; j++) {
    if ((*(pe+j-1))>(*(pe+j))) 
    {
      temp=(*(pe+j-1)); *(pe+j-1)=(*(pe+j)); *(pe+j)=temp; 
      temp=(*(pk+j-1)); *(pk+j-1)=(*(pk+j)); *(pk+j)=temp; 
      temp=(*(pi+j-1)); *(pi+j-1)=(*(pi+j)); *(pi+j)=temp; 
    }
  } }
}


// Start of the main program
int main ()
{
  int i, j, k;
  double *pintens;
  double *pk1, *pk2, *pk3, *pk;
  double k1, k2, k3;
  double *penergy;
  double sum;

  double *pkplus, *pkminus;
  double EvacP, EvacM;

  double VacOverlap, InfVacOverlap;
  
  double *pkw, *pew, *piw;
 
  int UNITS; 
  unsigned int dim;
  unsigned int s = 0;
  int SORT;

  complex_t K3V;
  complex_t temp;

  double COUP, A;

  printf("Starting lineInt...\n");

  fp = fopen("Input.in", "r");
  fscanf(fp, " %d ", &UNITS);
  fscanf(fp, " %lf ", &COUP);
  fclose(fp);
  
  A = EPSILON/2.0/COUP;
  dim = UNITS;
  SORT = ((UNITS <= 50) && (UNITS >2)) ? 1 : 0;

//  printf(" %d \n", SORT);

  pintens = (double *) malloc(dim*dim*dim*sizeof(double));
  pk1 = (double *) malloc(dim*dim*dim*sizeof(double));
  pk2 = (double *) malloc(dim*dim*dim*sizeof(double));
  pk3 = (double *) malloc(dim*dim*dim*sizeof(double));
  pk = (double *) malloc(dim*dim*dim*sizeof(double));
  penergy = (double *) malloc(dim*dim*dim*sizeof(double));

  pkplus = dvector(1,UNITS);
  pkminus = dvector(1,UNITS);

// Calculate "+" and "-" vacua
  i = AlphaBetaVals(pkplus, pkminus, UNITS);
//  for (i=1; i<=UNITS; i++) { printf(" %15.8f %15.8f \n", pkplus[i], pkminus[i]); }
  EvacP = getVac(pkminus, UNITS, COUP); EvacM = getVac(pkplus, UNITS, COUP);
  printf("Minus vac. = %25.20f \n", EvacM);
  printf("Plus vac.  = %25.20f \n", EvacP);
 
// Read in |<PHI{-}|PHI{+}>|^2 from "Overlap.out"
  fp = fopen("Overlap.out", "r");
  fscanf(fp, " %lf ", &VacOverlap);
  fclose(fp);
  if (fabs(2*COUP) < EPSILON) { printf("|<PHI{-}|PHI{+}>|^2 = %15.10f \n", VacOverlap); }
  else { printf("|<PHI{-}|sigma(1,x)|PHI{+}>|^2 = %15.10f \n", VacOverlap); }
  InfVacOverlap = VacOverlap;
 
//  printf("Clearing memory...\n");
//  for (i=0; i<dim; i++) { for (j=0; j<dim; j++) { for (k=0; k<dim; k++) { *(pintens+i*dim*dim+j*dim+k) = 0.0; } } }




// Weak coupling regime calculation
  if (fabs(A) > 1)
  {

  // Calculate number of 3-exc. states
    for (i = 0; i<dim; i++) { for (j = i+1; j<dim; j++) { for (k = j+1; k<dim; k++) { s++; } } } //  printf(" %d \n", s);

  // Calculate 0 -> 3 intensities for finite-size systems
    printf("Calculating finite 0 -> 3 intensities...\n");
    fp = fopen("FinK123s.out", "r");
    for (i=0; i<dim; i++) { for (j=0; j<dim; j++) { for (k=0; k<dim; k++) {
          fscanf(fp, "%lf %lf %lf %lf %lf", &k1, &k2, &k3, &(K3V.re), &(K3V.im));
          temp = conj(K3V); temp = prod(temp, K3V); *(pintens+dim*dim*i+dim*j+k) = temp.re/(UNITS*UNITS*UNITS)*VacOverlap;
          *(pk1+i*dim*dim+j*dim+k) = k1; *(pk2+i*dim*dim+j*dim+k) = k2; *(pk3+i*dim*dim+j*dim+k) = k3;
          *(pk+i*dim*dim+j*dim+k) = k1 + k2 + k3;
          *(penergy+i*dim*dim+j*dim+k) = getEn(k1, COUP) + getEn(k2, COUP) + getEn(k3, COUP);
//          printf("%12.8f  %12.8f  %12.8f  :  %12.8f  %12.8f  %12.8f \n", k1, k2, k3, getEn(k1, COUP), getEn(k2, COUP), getEn(k3, COUP));
    } } }
    fclose(fp);
    
    pkw = (double *) malloc(s*sizeof(double));
    pew = (double *) malloc(s*sizeof(double));
    piw = (double *) malloc(s*sizeof(double));

   // Sort 0 -> 3 intensities for finite-size systems by transition energy
   if (SORT != 0) { printf("Sorting finite 0 -> 3 intensities by transition energy...\n"); }
    s = 0;
    for (i=0; i<dim; i++) { for (j=i+1; j<dim; j++) { for (k=j+1; k<dim; k++) {
//      fprintf(fp, " %+15.8lf  %+15.8lf  %+15.8lf   ", *(pk1+i*dim*dim+j*dim+k), *(pk2+i*dim*dim+j*dim+k), *(pk3+i*dim*dim+j*dim+k));
      *(pew+s) = *(penergy+i*dim*dim+j*dim+k)+EvacM-EvacP;
      *(pkw+s) = *(pk1+i*dim*dim+j*dim+k) + *(pk2+i*dim*dim+j*dim+k) + *(pk3+i*dim*dim+j*dim+k); 
      *(piw+s) = *(pintens+i*dim*dim+j*dim+k);
      s++;
    } } }
//    printf(" s = %d \n", s); 
   if (SORT != 0) { bubblesort(pew, pkw, piw, s); }
    
   // Write 0 -> 3 intensities for finite-size systems to file
//    printf("Writing finite 0 -> 3 intensities to file...\n");
    s = 0; sum = 0.0;
    fp = fopen("3FinIntens.out", "w");
    for (i=0; i<dim; i++) { for (j=i+1; j<dim; j++) { for (k=j+1; k<dim; k++) {
      fprintf(fp, " %+25.20lf   ", *(pew+s));
      fprintf(fp, " %+25.20lf   ", *(pkw+s));
      fprintf(fp, " %+25.20lf \n", *(piw+s));
      sum += *(piw+s);
      s++;
    } } } //} fprintf(fp, "\n"); } fprintf(fp, "\n"); }
    fclose(fp);
    printf("3FinIntens.out successfully written...\n");
    printf("Total intensity to 3-excitation manifold (finite systems): %10.6f \n", sum);

    free(pkw); free(pew); free(piw);

   // Calculate 0 -> 1 intensities for finite-size systems
    printf("Calculating finite 0 -> 1 intensities...\n");
    fp = fopen("FinJsVsAlpha.out", "r");
    for (i=0; i<dim; i++)
    {
      fscanf(fp, "%lf %lf %lf", &k1, &(K3V.re), &(K3V.im));
      temp = conj(K3V); temp = prod(temp, K3V); *(pintens+i) = temp.re/UNITS*VacOverlap;
      *(pk1+i) = k1; *(penergy+i) = getEn(k1, COUP); // printf("%10.5f \n", *(penergy+i));
    }
    fclose(fp);

    pkw = (double *) malloc(dim*sizeof(double));
    pew = (double *) malloc(dim*sizeof(double));
    piw = (double *) malloc(dim*sizeof(double));

   // Sort 0 -> 1 intensities for finite-size systems by transition energy
    printf("Sorting finite 0 -> 1 intensities by transition energy...\n");
    for (i=0; i<dim; i++) { *(pew+i) = *(penergy+i)+EvacM-EvacP; *(pkw+i) = *(pk1+i); *(piw+i) = *(pintens+i); }
    bubblesort(pew, pkw, piw, dim);

   // Write 0 -> 1 intensities for finite-size systems to file
//    printf("Writing finite 0 -> 1 intensities to file...\n");
    sum = 0.0;
    fp = fopen("1FinIntens.out", "w");
    for (i=0; i<dim; i++)
    {
      fprintf(fp, "%+25.20lf  %+25.20lf  %+25.20lf \n", *(pew+i), *(pkw+i), *(piw+i)); 
      sum += *(piw+i);
    }
    fclose(fp);
    printf("1FinIntens.out successfully written...\n");
    printf("Total intensity to 1-excitation manifold (finite systems): %10.6f \n", sum);

    free(pkw); free(pew); free(piw);

   // Calculate 0 -> 3 (discrete) intensities for infinite-size systems
    printf("Calculating infinite 0 -> 3 intensities...\n");
    fp = fopen("InfK123s.out", "r");
    for (i=0; i<dim; i++) { for (j=0; j<dim; j++) { for (k=0; k<dim; k++) {
          fscanf(fp, "%lf %lf %lf %lf %lf", &k1, &k2, &k3, &(K3V.re), &(K3V.im));
          temp = conj(K3V); temp = prod(temp, K3V); *(pintens+dim*dim*i+dim*j+k) = temp.re/(UNITS*UNITS*UNITS)*InfVacOverlap;
          *(pk1+i*dim*dim+j*dim+k) = k1; *(pk2+i*dim*dim+j*dim+k) = k2; *(pk3+i*dim*dim+j*dim+k) = k3;
          *(pk+i*dim*dim+j*dim+k) = k1 + k2 + k3;
          *(penergy+i*dim*dim+j*dim+k) = getEn(k1, COUP) + getEn(k2, COUP) + getEn(k3, COUP);
//          printf("%12.8f  %12.8f  %12.8f  :  %12.8f  %12.8f  %12.8f \n", k1, k2, k3, getEn(k1, COUP), getEn(k2, COUP), getEn(k3, COUP));
    } } }
    fclose(fp);

    pkw = (double *) malloc(s*sizeof(double));
    pew = (double *) malloc(s*sizeof(double));
    piw = (double *) malloc(s*sizeof(double));

   // Sort 0 -> 3 (discrete) intensities for infinite-size systems by transition energy
    if (SORT != 0) { printf("Sorting infinite 0 -> 3 intensities by transition energy...\n"); }
    s = 0;
    for (i=0; i<dim; i++) { for (j=i+1; j<dim; j++) { for (k=j+1; k<dim; k++) {
//      fprintf(fp, " %+15.8lf  %+15.8lf  %+15.8lf   ", *(pk1+i*dim*dim+j*dim+k), *(pk2+i*dim*dim+j*dim+k), *(pk3+i*dim*dim+j*dim+k));
      *(pew+s) = *(penergy+i*dim*dim+j*dim+k)+EvacM-EvacP;
      *(pkw+s) = *(pk1+i*dim*dim+j*dim+k) + *(pk2+i*dim*dim+j*dim+k) + *(pk3+i*dim*dim+j*dim+k); 
      *(piw+s) = *(pintens+i*dim*dim+j*dim+k);
      s++;
    } } }
//    printf(" s = %d \n", s); 
    if (SORT != 0) { bubblesort(pew, pkw, piw, s); }
    
   // Write 0 -> 3 (discrete) intensities for infinite-size systems to file
//    printf("Writing infinite 0 -> 3 intensities to file...\n");
    s = 0; sum = 0.0;
    fp = fopen("3InfIntens.out", "w");
    for (i=0; i<dim; i++) { for (j=i+1; j<dim; j++) { for (k=j+1; k<dim; k++) {
      fprintf(fp, " %+25.20lf   ", *(pew+s));
      fprintf(fp, " %+25.20lf   ", *(pkw+s));
      fprintf(fp, " %+25.20lf \n", *(piw+s));
      sum += *(piw+s);
      s++;
    } } } //} fprintf(fp, "\n"); } fprintf(fp, "\n"); }
    fclose(fp);
    printf("3InfIntens.out successfully written...\n");
    printf("Total intensity to 3-excitation manifold (infinite systems): %10.6f \n", sum);

    free(pkw); free(pew); free(piw);

   // Calculate 0 -> 1 (discrete) intensities for infinite-size systems
    printf("Calculating infinite 0 -> 1 intensities...\n");
    fp = fopen("InfJsVsAlpha.out", "r");
    for (i=0; i<dim; i++)
    {
      fscanf(fp, "%lf %lf %lf", &k1, &(K3V.re), &(K3V.im));
      temp = conj(K3V); temp = prod(temp, K3V); *(pintens+i) = temp.re/UNITS*InfVacOverlap;
      *(pk1+i) = k1; *(penergy+i) = getEn(k1, COUP); //printf("%+10.8f \n", *(penergy+i));
    }
    fclose(fp);

    pkw = (double *) malloc(dim*sizeof(double));
    pew = (double *) malloc(dim*sizeof(double));
    piw = (double *) malloc(dim*sizeof(double));

   // Sort 0 -> 1 (discrete) intensities for infinite-size systems by transition energy
    printf("Sorting infinite 0 -> 1 intensities by transition energy...\n");
    for (i=0; i<dim; i++) { *(pew+i) = *(penergy+i)+EvacM-EvacP; *(pkw+i) = *(pk1+i); *(piw+i) = *(pintens+i); }
    bubblesort(pew, pkw, piw, dim);

   // Write 0 -> 1 (discrete) intensities for infinite-size systems to file
//    printf("Writing infinite 0 -> 1 intensities to file...\n")
    sum = 0.0;
    fp = fopen("1InfIntens.out", "w");
    for (i=0; i<dim; i++)
    {
      fprintf(fp, "%+25.20lf  %+25.20lf  %+25.20lf \n", *(pew+i), *(pkw+i), *(piw+i));
      sum += *(piw+i);
    }
    fclose(fp);
    printf("1InfIntens.out successfully written...\n");
    printf("Total intensity to 1-excitation manifold (infinite systems): %10.6f \n", sum);

    free(pkw); free(pew); free(piw);
  }

// Strong coupling regime calculation
  else
  {
  // Calculate number of 2-exc. states
    s = 0;
    for (i = 0; i<dim; i++) { for (j = i+1; j<dim; j++) { s++; } } //  printf(" %d \n", s);

  // Calculate 0 -> 2 intensities for finite-size systems
    printf("Calculating finite 0 -> 2 intensities...\n");
    fp = fopen("FinKsVsAlpha.out", "r");
    for (i=0; i<dim; i++) 
    {
      for (j=0; j<dim; j++) 
      {
        fscanf(fp, "%lf %lf %lf %lf", &k1, &k2, &(K3V.re), &(K3V.im));
        temp = conj(K3V); temp = prod(temp, K3V); *(pintens+dim*i+j) = temp.re/(UNITS*UNITS)*VacOverlap;
        *(pk1+i*dim+j) = k1; *(pk2+i*dim+j) = k2;
        *(pk+i*dim+j) = k1 + k2;
        *(penergy+i*dim+j) = getEn(k1, COUP) + getEn(k2, COUP);
//        printf("%12.8f  %12.8f  :  %12.8f  %12.8f \n", k1, k2, getEn(k1, COUP), getEn(k2, COUP));
      }
    }
    fclose(fp);

    pkw = (double *) malloc(s*sizeof(double));
    pew = (double *) malloc(s*sizeof(double));
    piw = (double *) malloc(s*sizeof(double));

   // Sort 0 -> 2 intensities for finite-size systems by transition energy
    if (SORT != 0) { printf("Sorting finite 0 -> 2 intensities by transition energy...\n"); }
    s = 0;
    for (i=0; i<dim; i++) { for (j=i+1; j<dim; j++) {
//      fprintf(fp, " %+15.8lf  %+15.8lf  %+15.8lf   ", *(pk1+i*dim*dim+j*dim+k), *(pk2+i*dim*dim+j*dim+k), *(pk3+i*dim*dim+j*dim+k));
      *(pew+s) = *(penergy+i*dim+j)+EvacM-EvacP;
      *(pkw+s) = *(pk1+i*dim+j) + *(pk2+i*dim+j); 
      *(piw+s) = *(pintens+i*dim+j);
      s++;
    } } 
//    printf(" s = %d \n", s); 
    if (SORT != 0) { bubblesort(pew, pkw, piw, s); }
    
   // Write 0 -> 2 intensities for finite-size systems to file
//    printf("Writing finite 0 -> 2 intensities to file...\n");
    s = 0; sum = 0.0;
    fp = fopen("2FinIntens.out", "w");
    for (i=0; i<dim; i++) { for (j=i+1; j<dim; j++) {
      fprintf(fp, " %+25.20lf   ", *(pew+s));
      fprintf(fp, " %+25.20lf   ", *(pkw+s));
      fprintf(fp, " %+25.20lf \n", *(piw+s));
      sum += *(piw+s);
      s++;
    } } //} fprintf(fp, "\n"); } fprintf(fp, "\n"); }
    fclose(fp);
    printf("2FinIntens.out successfully written...\n");
    printf("Total intensity to 2-excitation manifold (finite systems): %10.6f \n", sum);

    free(pkw); free(pew); free(piw);

  // Calculate (discrete) 0 -> 2 intensities for infinite-size systems
    printf("Calculating infinite 0 -> 2 intensities...\n");
    fp = fopen("InfKsVsAlpha.out", "r");
    for (i=0; i<dim; i++) 
    {
      for (j=0; j<dim; j++) 
      {
        fscanf(fp, "%lf %lf %lf %lf", &k1, &k2, &(K3V.re), &(K3V.im));
        temp = conj(K3V); temp = prod(temp, K3V); *(pintens+dim*i+j) = temp.re/(UNITS*UNITS)*InfVacOverlap;
        *(pk1+i*dim+j) = k1; *(pk2+i*dim+j) = k2;
        *(pk+i*dim+j) = k1 + k2;
        *(penergy+i*dim+j) = getEn(k1, COUP) + getEn(k2, COUP);
//        printf("%12.8f  %12.8f  :  %12.8f  %12.8f \n", k1, k2, getEn(k1, COUP), getEn(k2, COUP));
      }
    }
    fclose(fp);

    pkw = (double *) malloc(s*sizeof(double));
    pew = (double *) malloc(s*sizeof(double));
    piw = (double *) malloc(s*sizeof(double));

   // Sort 0 -> 2 intensities for finite-size systems by transition energy
    if (SORT != 0) { printf("Sorting infinite 0 -> 2 intensities by transition energy...\n"); }
    s = 0;
    for (i=0; i<dim; i++) { for (j=i+1; j<dim; j++) {
//      fprintf(fp, " %+15.8lf  %+15.8lf  %+15.8lf   ", *(pk1+i*dim*dim+j*dim+k), *(pk2+i*dim*dim+j*dim+k), *(pk3+i*dim*dim+j*dim+k));
      *(pew+s) = *(penergy+i*dim+j)+EvacM-EvacP;
      *(pkw+s) = *(pk1+i*dim+j) + *(pk2+i*dim+j); 
      *(piw+s) = *(pintens+i*dim+j);
      s++;
    } } 
//    printf(" s = %d \n", s); 
    if (SORT != 0) { bubblesort(pew, pkw, piw, s); }

   // Write 0 -> 2 intensities for infinite-size systems to file
//    printf("Writing infinite 0 -> 2 intensities to file...\n");
    s = 0; sum = 0.0;
    fp = fopen("2InfIntens.out", "w");
    for (i=0; i<dim; i++) { for (j=i+1; j<dim; j++) { 
//      fprintf(fp, " %+15.8lf  %+15.8lf  ", *(pk1+i*dim+j), *(pk2+i*dim+j));
      fprintf(fp, " %+25.20lf   ", *(pew+s));
      fprintf(fp, " %+25.20lf   ", *(pkw+s));
      fprintf(fp, " %+25.20lf \n", *(piw+s));
      sum += *(piw+s);
      s++;
    } } //} fprintf(fp, "\n"); } fprintf(fp, "\n"); }
    fclose(fp);
    printf("2InfIntens.out successfully written...\n");
    printf("Total intensity to 2-excitation manifold (infinite systems): %10.6f \n", sum);

    free(pkw); free(pew); free(piw);
  }

/*
  fp = fopen("testeigvals.out", "w");
  for (i=0; i<dim; i++) { for (j=i+1; j<dim; j++) { for (k=j+1; k<dim; k++) {
    fprintf(fp, "%+15.8lf  %+15.8lf  %+15.8lf  ", *(pk1+i*dim*dim+j*dim+k), *(pk2+i*dim*dim+j*dim+k), *(pk3+i*dim*dim+j*dim+k) );
    fprintf(fp, "%+15.8lf  \n", *(penergy+i*dim*dim+j*dim+k)+EvacM);
  }}}//} fprintf(fp, "\n"); } fprintf(fp, "\n"); }
  fclose(fp);
  printf("testeigvals.out successfully written...\n");
*/  

  free(pintens); free(pk1); free(pk2); free(pk3); free(pk); free(penergy);
  free_dvector(pkplus, 1,UNITS); free_dvector(pkminus, 1,UNITS);

  printf("lineInt successfully finished!\n\n");

  return 0;
}

