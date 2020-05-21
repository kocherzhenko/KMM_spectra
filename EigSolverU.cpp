// Eigensolver for X(alpha1, alpha2), a real matrix - (c) Aleksey Kocherzhenko, 2011

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "Setup.h"

// Transformation of Hamiltonian to 3-diagonal form by a sequence of (UNITS2-2) Householder reflections
// Adapted from "Numerical recepies in C" tred2.c
void hhtridiag (double *pmatr, double *pdiag, double *poffdiag, int UNITS)
{
  int l, k, j, i;
  double scale, hh, h, g, f;
  int UNITS2 = 2*UNITS;

  for (i=UNITS2-1; i>=1; i--)
  { 
    l=i-1;
    h=scale=0.0;
    if (l>0)
    {
      for (k=0; k<=l; k++) {scale += fabs(*(pmatr+i*UNITS2+k));}
      if (scale == 0.0) {*(poffdiag+i)=(*(pmatr+i*UNITS2+l));}
      else 
      {
        for (k=0; k<=l; k++) 
        { 
          *(pmatr+i*UNITS2+k) /= scale;
          h += (*(pmatr+i*UNITS2+k))*(*(pmatr+i*UNITS2+k));
        }
        f=*(pmatr+i*UNITS2+l);
        g=(f >= 0.0 ? -sqrt(h) : sqrt(h));
        *(poffdiag+i)=scale*g;
        h -= f*g;
        *(pmatr+i*UNITS2+l) = f-g;
        f = 0.0;
        for (j=0; j<=l; j++)
        {
          *(pmatr+j*UNITS2+i)=(*(pmatr+i*UNITS2+j))/h;
          g = 0.0;
          for (k=0; k<=j; k++) {g += (*(pmatr+j*UNITS2+k))*(*(pmatr+i*UNITS2+k));}
          for (k=j+1; k<=l; k++) {g += (*(pmatr+k*UNITS2+j))*(*(pmatr+i*UNITS2+k));}
          *(poffdiag+j)=g/h;
          f += (*(poffdiag+j))*(*(pmatr+i*UNITS2+j));
        }
        hh=f/(h+h);
        for (j=0; j<=l; j++)
        {
          f=(*(pmatr+i*UNITS2+j));
          *(poffdiag+j)=g=(*(poffdiag+j)-hh*f);
          for (k=0; k<=j; k++) {*(pmatr+j*UNITS2+k) -= (f*(*(poffdiag+k))+g*(*(pmatr+i*UNITS2+k)));}
        }
      }
    }
    else {*(poffdiag+i)=(*(pmatr+i*UNITS2+l));}
    *(pdiag+i)=h;
  }
  *(pdiag)=0.0;
  *(poffdiag)=0.0;
  for (i=0; i<UNITS2; i++)
  {
    l=i-1;
    if (*(pdiag+i)) 
    {
      for (j=0; j<=l; j++)
      {
        g=0.0;
        for (k=0; k<=l; k++) {g += (*(pmatr+i*UNITS2+k))*(*(pmatr+k*UNITS2+j));}
        for (k=0; k<=l; k++) {*(pmatr+k*UNITS2+j) -= g*(*(pmatr+k*UNITS2+i));}
      }
    }
    *(pdiag+i)= (*(pmatr+i*UNITS2+i));
    *(pmatr+i*UNITS2+i)=1.0;
    for (j=0; j<=l; j++) (*(pmatr+j*UNITS2+i))=(*(pmatr+i*UNITS2+j))=0.0;
  }
  
}

// Diagonalization using QL factorization with implicit shifts
// Adapted from "Numerical recepies in C" tqli.c
int qlfact(double* pdiag, double* poffdiag, double *pmatr, int UNITS)
{
  int m,l,iter,i,k = 0;
  double s,r,p,g,f,dd,c,b;
  double temp;
  int UNITS2=2*UNITS;

  for (i=1; i<UNITS2; i++) *(poffdiag+i-1)=*(poffdiag+i);
  *(poffdiag+UNITS2-1)=0.0;

  for (l=0; l<UNITS2; l++)
  {
    iter=0;
    do
    {
      for (m=l; m<UNITS2-1; m++)
      {
        dd=fabs(*(pdiag+m))+fabs(*(pdiag+m+1));
        temp=fabs(*(poffdiag+m))+dd;
        if ((double)(temp == dd)) {break;}
      }
      if (m != l)
      {
        if (iter++ == 30) {printf("Too many iterations in TQLI"); return 1;}
        g=((*(pdiag+l+1))-(*(pdiag+l)))/(2.0*(*(poffdiag+l)));
        r=sqrt((g*g)+1.0);
        g=(*(pdiag+m))-(*(pdiag+l))+(*(poffdiag+l))/(g+SIGN(r,g));
        s=c=1.0;
        p=0.0;

        for (i=m-1;i>=l;i--)
        {
          f=s*(*(poffdiag+i));
          b=c*(*(poffdiag+i));

          *(poffdiag+i+1)=(r=sqrt((f*f)+(g*g)));
          if (r == 0.0)
          {
            *(pdiag+i+1) -= p;
            *(poffdiag+m) = 0.0;
            break;
          }

          s=f/r; c=g/r;
          g=(*(pdiag+i+1))-p;
          r=((*(pdiag+i))-g)*s+2.0*c*b;
          *(pdiag+i+1)=g+(p=s*r);
          g=c*r-b;

          for (k=0; k<UNITS2; k++)
          {
            f=(*(pmatr+k*UNITS2+i+1));
            *(pmatr+k*UNITS2+i+1)=s*(*(pmatr+k*UNITS2+i))+c*f;
            *(pmatr+k*UNITS2+i)=c*(*(pmatr+k*UNITS2+i))-s*f;
          }
        }
 
        if (r == 0 && i>=l) continue;        

        *(pdiag+l)=(*(pdiag+l))-p;
        *(poffdiag+l)=g;
        *(poffdiag+m)=0.0;
      }
    }
    while (m != l); 
  } 
  return 0;
}

