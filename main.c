#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "basic_LinAlg.h"
#include "CC.h"

int main() {

  /*
    ................................ Tests
  */

  double * SA = symmat_neu(16);
  long n = 10;
  long m = 1;
  matrix_laden_ascii("A_ascii.dat", &SA, &m, &n);
  vektor_ausgeben(SA, 10, "% 10.3e");
  double i = symmat_FrobeniusNorm(SA, 4);
  printf("\n %f \n\n", i);

  symmat_ausgeben(SA, 4,"% 10.3e");
  printf("\n");
  symmat_ausgeben_dreieck(SA, 4, "% 10.3e");

  double * x = vektor_neu(4);
  x[0] = 1;
  x[1] = 2;
  x[2] = 3;
  x[3] = 4;

  double * y = vektor_neu(4);
  y[0] = 1;
  y[1] = 2;
  y[2] = 3;
  y[3] = 4;

  printf("\n");
  symmat_vektor_mult(1,1,SA,x,y,4);
  vektor_ausgeben(y, 4, "% 10.3e");

//DARAN DENKEN
  symmat_freigeben(SA);
  vektor_freigeben(x);
  vektor_freigeben(y);

  return 0;
}
