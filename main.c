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
  double* x = vektor_neu(n);
  double* x_optimal = vektor_neu(n);


  double* t = vektor_neu(n);
  //setze tv auf 1
  for(int i = 0; i<225; i++) {
    x[i] = 0;
    t[i] = 0;
  }

  double* b = vektor_neu(n);
  long m = 1;
  matrix_laden_ascii("b_ascii.dat", &b, &n, &m);
  matrix_laden_ascii("x_exakt_ascii.dat", &b, &n, &m);

  //Sparse_Matrix_mul(A, tv, tv2, n);
  
  

  //TEST Sparse_Matrix_mul(A,x,t,n);

  //matrix_ausgeben(t, n, 1, " % 10.3e");

  
    Sparse_JacobiCheb_solver(A, x, b, n, 1000, GAMMA_1, GAMMA_2);
    for(int j = 0; j<n; j++) {
      t[j] = x_optimal[j] - x[j]; 
    }
    printf("Die norm nach %d Schritten: %f \n", (1)*100, vektor_2Norm(t, n));
  

  Sparse_Matrix_mul(A,x,t,n);
  Sparse_Matrix_mul(A,t,x,n);
  matrix_ausgeben(x, n, 1, " % 10.3e");
  return 0;
}
