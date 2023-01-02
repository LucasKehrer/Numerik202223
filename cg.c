#include <math.h>
#include <stdio.h>
#include "cg.h"
#include "basic_LinAlg.h"

#define OMEGA 1.581120979521361

/*
 * Berechnet B aus einer Sparse_Matrix und dem Relaxationsparameter omega
*/
double* calculateRelaxationMatrix(SparseMatrix_t A, double omega, long n) {

  //get diagonal of A as vector
  double* dia = vektor_neu(n);
  Spars_getDiag(A, dia, n);

  //get diainv of A as a vector
  double* dia_inv = vektor_neu(n);
  Spars_getDiaginverse(A, dia_inv, n);

  //(1/(2-omega)) = a
  double a = 1/(2-omega);

  //L
  double* L = SparsetoNormal(getL(A,n),n); 

  //wL
  for (int i = 0; i<n; i++) {
    for (int j = 0; j<n; j++) {
      L[i*n+j]=L[i*n+j]*omega;
    }
  }
  
  //L^t
  double* Lt = transposeMatrix(L,n);

  //wL^t
  for (int i = 0; i<n; i++) {
    for (int j = 0; j<n; j++) {
      Lt[i*n+j]=Lt[i*n+j]*omega;
    }
  }

  //(D+wL)
    
  for (int i = 0; i<n; i++) {
    L[i*n+i] = L[i*n+i]+dia[i];
  }

  //a*(D+wL)
  for (int i = 0; i<n; i++) {
    for (int j = 0; j<n; j++) {
      Lt[i*n+j]=Lt[i*n+j]*a;
    }
  }


  // (D+wL^t)
  
  for (int i = 0; i<n; i++) {
    Lt[i*n+i]=Lt[i*n+i]+dia[i];
  }


  //a*(D+WL)*D^-1

  

  //

  return 0;



}

/* Diese Funktion wendet den SSOR-Vorkonditionierer an, w=B^{-1}r.
 *
 * A - Dünnbesetzte Matrix der Dimension nxn (Eingabe)
 * w - Lösung von Bw=r (Ausgabe)
 * r - Rechte Seite (Eingabe)
 * n - Anzahl Zeilen und Spalten von A (Eingabe)
 * omega - Relaxationsparameter aus dem Intervall (0,2) (Eingabe)
 *
 */
static void Sparse_SSOR_prec(SparseMatrix_t A, double *w, double *r, long n, double omega)
{
}

/* Diese Funktion realisiert die Methode der konjugierten Gradienten 
 * (cg-Verfahren) mit Vorkonditionierung zur Approximation der Loesung x eines
 * linearen Gleichungssystems.
 *
 * A - Matrix im sparse Format der Dimension nxn (Eingabe)
 * x - Startvektor der Iteration und anschliessend 
 *     Approximation des Loesungsvektors (Ein- und Ausgabe)
 * b - Rechte Seite des Gleichungssystems (Eingabe)
 * n - Anzahl Zeilen und Spalten von A (Eingabe)
 * prec - Ist prec=0, wird nicht vorkonditioniert. Ist prec!=0, wird
 * der SSOR-Vorkonditionierer verwendet.
 * iter - Anzahl der maximal durchzufuehrenden Iterationen (Eingabe)
 * eps - Abbruchsgenauigkeit (Eingabe)
 *
 * Rueckgabewert: Anzahl der durchgefuehrten Iterationen.
 *
 */
long Sparse_cg(SparseMatrix_t A, double *x, double *b, long n, int prec, long iter, double eps) {
  return 0;
}



