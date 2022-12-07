#include <stdio.h>
#include <stdlib.h>

#include "basic_LinAlg.h"
#include "iter.h"

int main()
{
  long n = 4;
  long m = 4;
  long n1 = 1;
  long m1 = 10;
  double* x = vektor_neu(4);
  double* b = vektor_neu(4);
  double* A = matrix_neu(4,4);
  matrix_laden_ascii("b_ascii.dat", &b,&n,&m);
  int k = matrix_laden_ascii("A_ascii.dat", &A, &n1, &m1);
  for(int i = 0; i<n; i++) {x[i] = 0;}
  for(int i = 1; i<=3;i++) {
    Jacobi_solver(A,x,b,n,20*i);
    printf("\nDas ist x(%i):\n",i*20);
    vektor_ausgeben(A, n, "%10.3e");
    double* b1 = vektor_neu(n);
    vektor_kopieren(b1,b,n);
    for(int i = 0; i<n; i++) {
      b1[i] = -b1[i];
    }
    matrix_mult(1,1,A,x,b1,n,1,n);
    printf("\nDas ist Ax(%i) - b:\n",i*20);
    vektor_ausgeben(b1,n,"%10.3e");
    vektor_freigeben(b1);
  }

  vektor_freigeben(x);
  vektor_freigeben(b);
  matrix_freigeben(A);
  return 0;
}
