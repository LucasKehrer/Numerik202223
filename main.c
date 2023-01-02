#include <stdio.h>
#include <stdlib.h>

#include "basic_LinAlg.h"
#include "cg.h"

#define EPS 1e-12
#define MAX_STEPS 75
#define STEPS 5

int main()
{
  //Die Matrix A(sparse)
  SparseMatrix_t A;
  long n = 225;
  long m = 1;
  SparseMatrix_laden_ascii("A_sparse_ascii.dat", &A, &n);

  //Probevektor x_optimal
  double* x_optimal = vektor_neu(n);
  matrix_laden_ascii("x_ascii.dat", &x_optimal, &n, &m);
  /*Test
  printf("Der optimale Vektor: \n");
  vektor_ausgeben(x_optimal, n, " % 10.3e");
   */

  //Die Lösung des LGS b
  double* b = vektor_neu(n);
  matrix_laden_ascii("b_ascii.dat", &b, &n, &m);
  
  //Nullvektor Test t
  //Nullvektor x_0 
  double* t = vektor_neu(n);
  double* x_0 = vektor_neu(n);
  for (int i = 0; i<n; i++) {
    t[i] = 0;
    x_0[i] = 0;
  }


  //TEST: Sparsematrix Mult
  // SparseMatrix_vektor_mult(1, 0, A, x_optimal, t, n);
  // vektor_ausgeben(t, n, " % 10.3e");
/*
  SparseMatrix_t B = getL(A, n);
  for (int i = 0; i<420; i++) {
    printf("%ld \n", B.columns[i]);
  }
*/

  matrix_ausgeben(SparsetoNormal(A, n),n,n, " % 10.3e");
  return 0;
}

