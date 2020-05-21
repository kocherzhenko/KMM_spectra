// Common Functions for KMM Spectroscopy Programs - (c) Aleksey Kocherzhenko, 2012

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stddef.h>
#include "Setup.h"

typedef struct {double re; double im;} complex_t;

// complex number operations
complex_t prod(complex_t c1, complex_t c2) { complex_t prod; prod.re = c1.re*c2.re-c1.im*c2.im; prod.im = c1.re*c2.im+c1.im*c2.re; return prod; }
complex_t add(complex_t c1, complex_t c2) { complex_t add; add.re = c1.re+c2.re; add.im = c1.im+c2.im; return add; }
complex_t prod(complex_t c, double d) { complex_t prod; prod.re = c.re*d; prod.im = c.im*d; return prod; }
complex_t add(complex_t c, double d) { complex_t add; add.re = c.re+d; add.im = c.im; return add; }
complex_t inv(complex_t c) { complex_t i; i.re = c.re/(c.re*c.re+c.im*c.im); i.im = -c.im/(c.re*c.re+c.im*c.im); return i; }
complex_t neg(complex_t c) { complex_t n; n.re = -c.re; n.im = -c.im; return n; }
complex_t conj(complex_t c) { complex_t c1; c1.re = c.re; c1.im = -c.im; return c1; }
complex_t csqrt(complex_t c) { double r, phi; complex_t cs; 
  r = sqrt(c.re*c.re+c.im*c.im); phi = atan(c.im/c.re); if (c.re < 0.0) { if (c.im >= 0) { phi += PI; }  else { phi -= PI; } } 
  r = sqrt(r); phi = phi/2.0; cs.re = r*cos(phi); cs.im = r*sin(phi); return cs; }
complex_t cexp(complex_t c) { complex_t ce; ce.re = pow(EXP, c.re) * cos(c.im); ce.im = pow(EXP, c.re) * sin(c.im); return ce; }

// parameters in Eq. (14B)

// Calculate THETA(k) = exp[-2*i*theta(k)]
complex_t getTH(double k, double COUP) 
{ 
  complex_t Theta, temp1, temp2; 
  double E0k, Ek;
  double A = EPSILON/2.0/COUP;

  E0k = EPSILON + 2*COUP*cos(k); Ek = sqrt(E0k*E0k+4*COUP*COUP*sin(k)*sin(k));
  Theta.re = E0k/Ek; Theta.im = -2*COUP*sin(k)/Ek;

  if (fabs(A) < 1) { temp1.re = 0.0; temp1.im = k; temp1 = cexp(temp1); Theta = prod(temp1, Theta); }
  
  return Theta; 
}

// Calculate allowed k values for the "+" branch (alpha) and "-" branch (beta)
int AlphaBetaVals(double *pkplus, double *pkminus, int UNITS)
{
  int i, j; int odd; double n = 0.0; double t, *pt = &t;
  
  if ((UNITS % 2) == 0) { odd=0; } else { odd=1; }

  switch (odd)
  {
    case 0 :
    {
      n = -(UNITS-2.0); for (i=1; i<=UNITS; i++) { pkplus[i] = n/UNITS*PI; n = n+2.0; }
      n = -(UNITS-1.0); for (i=1; i<=UNITS; i++) { pkminus[i] = n/UNITS*PI; n = n+2.0; }
      break;
    }
    case 1 :
    {
      n = -(UNITS-1.0); for (i=1; i<=UNITS; i++) { pkplus[i] = n/UNITS*PI; n = n+2.0; }
      n = -(UNITS-2.0); for (i=1; i<=UNITS; i++) { pkminus[i] = n/UNITS*PI; n = n+2.0; }
      break;
    }
  }

  return 0;
}

// Calculate energy of single-exciton state with given k
double getEn(double k, double COUP)
{
   double energy;

   energy = EPSILON + 2*COUP*cos(k);
   energy *= energy;
   energy += 4*COUP*COUP*sin(k)*sin(k);
   energy = sqrt(energy);

   return energy;
}

// Function to print the binary representation of a number
void printbin(unsigned int n)
{
  unsigned int i;
  i=1<<(sizeof(n)*8-1);
  while (i>0)
  {
    if (n & i) {printf("1");} else {printf("0");};
    i>>=1;
  }
  printf("\n");
}

// Function to check whether the number n has a 1 at bit position m
unsigned int checkbin(unsigned int n, unsigned int m)
{
  unsigned int i;
  unsigned int power = 1;

  for (i=1; i<m; i++) { power *= 2; } 
  if (n & power) { i = 1; } else { i = 0; }

  return i;
}

