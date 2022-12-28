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
  matrix_laden_ascii("x_exakt_ascii.dat", &x_optimal, &n, &m);

  //Sparse_Matrix_mul(A, tv, tv2, n);
  
  

  //TEST Sparse_Matrix_mul(A,x,t,n);

  //matrix_ausgeben(t, n, 1, " % 10.3e");

  //einfach 10 mal laufen lassen, die Funktion signatur passt für das nicht...
  //UND NOCH DEN JAKOBI AN SICH IMPLEMENTIEREN
  
  printf("\n Mit Beschleunigung \n");

  for (int i = 1; i<=10; i++) {
    
    Sparse_JacobiCheb_solver(A, x, b, n, i*100, GAMMA_1, GAMMA_2);
    for(int j = 0; j<n; j++) {
      t[j] = x_optimal[j] - x[j]; 
    }
    printf("Die norm nach %d Schritten: %.10f \n", (i)*100, vektor_2Norm(t, n));

    for(int j = 0; j<225; j++) {
    x[j] = 0;
    t[j] = 0;
  }

  }

  printf("\n Ohne Beschleunigung \n");

  for (int i = 1; i<=10; i++) {
    
    Sparse_Jacobi_solver(A, x, b, n, i*100);
    for(int j = 0; j<n; j++) {
      t[j] = x_optimal[j] - x[j]; 
    }
    printf("Die norm nach %d Schritten: %.10f \n", (i)*100, vektor_2Norm(t, n));

    for(int j = 0; j<225; j++) {
    x[j] = 0;
    t[j] = 0;
  }

  }

  vektor_freigeben(x);
  vektor_freigeben(x_optimal);
  vektor_freigeben(t);
  vektor_freigeben(b);
  SparseMatrix_freigeben(A);

  return 0;
}
