#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "basic_LinAlg.h"
#include "QR.h"

int main()
{

  /*
    ...................
    * Die QR zerlegung passt, die hab ich auch im rechner gecheckt
    * Bei der Q^t substituierung kann was schief gehen, idk. obwohl die orthogonale Matrix passt. Habe aber keiner Zeit mehr das genau zu debuggen
    * Müsste mal ein Framework einrichten dafür...
  */
  double* A = matrix_neu(4, 4);
  long n = 4;


  matrix_laden_ascii("A_ascii.dat", &A, &n, &n);

  double* v = vektor_neu(4);
  
  QR_HS_Zerlegung(4,A,v);
  printf("R mit u's auf unterer diagonale und u\n");
  matrix_ausgeben(A, 4 ,4, "% 10.3e");
  printf("\n----\n");
  vektor_ausgeben(v, 4, "% 10.3e");
  
  printf("---------------------\n---------------------\n---------------------\n");


  double * x = vektor_neu(4);
  long m = 1;
  matrix_laden_ascii("b_ascii.dat", &x, &n, &m );

  QR_RueckwSubst(n, A, x);
  printf("nach ruecksubst1: \n");
  vektor_ausgeben(x, 4, " % 10.3e");
  printf("\n----\n");
  QR_HS_QTransX(4, A, v, x);
  printf("Nach Q^t subst: \n");
  vektor_ausgeben(x, 4, " % 10.3e");
  
  matrix_freigeben(A);
  vektor_freigeben(v);
  vektor_freigeben(x);

  printf("---------------------\n---------------------\n---------------------\n");


  return 0;
}
