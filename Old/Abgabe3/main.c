#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "basic_LinAlg.h"
#include "CC.h"

int main() {

  /*
    ................................ Tests
  */

  double * SA = symmat_neu(4);
  double * b = vektor_neu(4);
  long n = 10;
  long m = 1;
  matrix_laden_ascii("A_ascii.dat", &SA, &m, &n);
  //Lade b
  n=4;
  matrix_laden_ascii("b_ascii.dat", &b, &n, &m);
  printf(" C is: \n");
  CC_Zerlegung(4, SA);
  symmat_ausgeben_dreieck(SA, 4, "% 10.3e");

  printf("\n Vorwaertssubst \n");
  
  CC_VorwSubst(4, SA, b);
  vektor_ausgeben(b, 4, "% 10.3e") ;
  printf("\n Rueckwaertssubst \n");

  CC_RueckwSubst(4, SA, b);
  vektor_ausgeben(b, 4, "% 10.3e") ;
printf("\n Probe \n");
  //Probe

  //sry, ich finde die Art wie die deklariert ist seltsam...
  double * x = vektor_neu(4);
  x[0] = 0;
  x[1] = 0;
  x[2] = 0;
  x[3] = 0;

  //Probe
  matrix_laden_ascii("A_ascii.dat", &SA, &m, &n);
  symmat_vektor_mult(1, 0, SA, b, x, 4); //warum nur die cholesky matrix mal b ???
  vektor_ausgeben(x, 4, "% 10.3e");


//DARAN DENKEN :)
  symmat_freigeben(SA);
  vektor_freigeben(b);
  vektor_freigeben(x);

  return 0;
}
