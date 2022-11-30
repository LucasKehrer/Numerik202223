#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "basic_LinAlg.h"
#include "QR.h"

int main()
{

  /*
    ...................
  */
  double* A = matrix_neu(4, 4);
  long n = 4;


  matrix_laden_ascii("A_ascii.dat", &A, &n, &n);

  double* v = vektor_neu(4);
  
  QR_HS_Zerlegung(4,A,v);

  matrix_ausgeben(A, 4 ,4, "% 10.3e");

  matrix_freigeben(A);
  vektor_freigeben(v);
  return 0;
}
