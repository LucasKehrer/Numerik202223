#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "basic_LinAlg.h"
#include "iter.h"

/* Optimale Werte für die Parameter. */
#define GAMMA_1 0.01921471959677
#define GAMMA_2 1.98078528040323

int main()
{  
  SparseMatrix_t A;
  long n=225;


  SparseMatrix_laden_ascii("A_sparse_ascii.dat", &A, &n);

  
  //testvektor
  double* tv = vektor_neu(n);
  //setze tv auf 1
  for(int i = 0; i<225; i++) {
    tv[i] = 1;
  }

  double* tv2 = vektor_neu(n);

  Sparse_Matrix_mul(A, tv, tv2, n);
  matrix_ausgeben(tv2, n, 1, " % 10.3e");


  return 0;
}
