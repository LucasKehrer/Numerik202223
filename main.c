#include <stdio.h>
#include <stdlib.h>

#include "basic_LinAlg.h"
#include "cg.h"

#define EPS 1e-12
#define MAX_STEPS 75
#define STEPS 5

int main()
{
  /*
  * Hinweis: Ich hab helperfunctions in cg.c eingefügt, 
  * deshalb bitte die geänderte Header Datei beachten!
  *
  */
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
    x_0[i] = b[i];
  }

  long r = Sparse_cg(A, x_0, b, n, 1, 75, EPS);
  printf("Abbruch nach %ld iterationen \n", r);
  //vektor_ausgeben(x_0, n, " % 10.10e");


  

  printf("Ohne Vorkoditionierung:--------------------\n");
  for (int i = 1; i<=15; i++) {
    vektor_kopieren(x_0,b,n);
    long r = Sparse_cg(A, x_0, b, n, 0, i*5, EPS);
    if (r+5<=i*5) {
      break;
    }
    printf("Für maxIterationen %d, Abbruch nach %ld iterationen \n", i*5, r);

    for (int i = 0; i<n; i++) {
    t[i] = x_0[i] - x_optimal[i];
    }

    printf("Fehler nach %d Iterationen: %10.10e \n", r, vektor_skalprod(t, t,n) /vektor_skalprod(x_0, x_0,n));
  }


  printf("Mit Vorkoditionierung:--------------------\n");
  for (int i = 1; i<=15; i++) {
    vektor_kopieren(x_0,b,n);
    long r = Sparse_cg(A, x_0, b, n, 1, i*5, EPS);
    if (r+5<=i*5) {
      break;
    }
    printf("Für maxIterationen %d, Abbruch nach %ld iterationen \n", i*5, r);

    for (int i = 0; i<n; i++) {
    t[i] = x_0[i] - x_optimal[i];
    }

    printf("Fehler nach %d Iterationen: %10.10e \n", r, vektor_skalprod(t, t,n) / vektor_skalprod(x_0, x_0,n));
  }

  /*for (int i = 0; i<n; i++) {
    t[i] = 1;
  }
  double *testv= vektor_neu(n);
  SparseMatrix_t Atest = getL(A,n);
  SparseMatrix_vektor_mult(1,0,Atest,t, testv, n);
  vektor_ausgeben(testv, n,  "% 10.3e");
  solve_lowerSparse(Atest, testv, t, n, 1);
  vektor_ausgeben(t, n,  "% 10.3e");
*/


   

  //TEST: Sparsematrix Mult
  // SparseMatrix_vektor_mult(1, 0, A, x_optimal, t, n);
  // vektor_ausgeben(t, n, " % 10.3e");
/*
  SparseMatrix_t B = getL(A, n);
  for (int i = 0; i<420; i++) {
    printf("%ld \n", B.columns[i]);
  }


  SparseMatrix_t testS;
  long nt = 2;
  SparseMatrix_laden_ascii("test_ascii_sparse.dat", &testS, &nt);
  double tx[5];
  double ty[5];
  ty[0] = 3;
  ty[1] = 3;
  ty[2] = 3;
  ty[3] = 3;
  ty[4] = 4;

  printf("%f %f %f \n", testS.values[0], testS.values[1], testS.values[2]);
  printf("%li %li %li \n", testS.rowIndex[0], testS.rowIndex[1], testS.rowIndex[2]);
  printf("%li %li %li \n", testS.columns[0], testS.columns[1], testS.columns[2]);

  printf("Sparse column %li (r[1]-r[0] = %li - %li)\n", Sparse_getRowNumber(testS, 0),A.rowIndex[1],A.rowIndex[0]);
  printf("Sparse column %li\n", Sparse_getRowNumber(testS, 1));
  printf("Colums of i=1 val 0 %li\n", Sparse_getColumn(testS, 1, 0));
  solve_upperSparse(testS, ty, tx, 5, 1);

  printf("%f %f %f %f %f \n", tx[0], tx[1], tx[2], tx[3], tx[4]);
  // matrix_ausgeben(SparsetoNormal(A, n),n,n, " % 10.3e");

*/
  return 0;


}

