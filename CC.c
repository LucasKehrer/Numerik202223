#include <math.h>
#include <stdio.h>
#include "CC.h"
#include "basic_LinAlg.h"


double min(double a, double b) {
    return a<b ? a : b;
}

double max(double a, double b) {
    return a<b ? b : a;
}

/* Diese Funktion berechnet die Cholesky-Zerlegung einer 
 * symmetrisch positiv definiten Matrix A.
 *
 * n - Anzahl Zeilen und Spalten von A (Eingabe)
 * A - symmetrische Matrix von der lediglich die untere
 *     Dreiecksmatrix abgespeichert ist (Ein- und Ausgabe)
 *
 * Rueckgabewert: Bei Erfolg gibt die Funktion 0 zurueck.
 *                Ist die Matrix A nicht positiv definit,
 *                wird -1 zurueck gegeben.
 *
 * Achtung: Ist das Argument der Wurzelfunktion im 
 *          Cholesky-Algorithmus negativ, so ist die 
 *          Matrix nicht positiv definit und die Funktion
 *          bricht ab.
 *
 * Achtung: Die Matrix A wird mit der Matrix
 *          C ueberschrieben.
 */
long CC_Zerlegung(long n, double *A) {
  double *M = matrix_neu(n,n);
  for(int i= 0; i<n*n; i++) {
    int row = max(i%n,i/n);
    int col = min(i%n,i/n);
    M[i] = A[(row*(row+1)/2) + col];
  }

  for(int i = 0; i<n-1; i++) { //Spalten
    double diag = M[i*n+i];
    if (diag <= 0) {
      return -1;
    }
    A[(i*(i+1)/2)+i] = sqrt(diag);

    for(int j = i+1; j<n; j++) { //Zeilen
      double l = M[j*n+i]/diag;
      A[(j*(j+1)/2)+i] = l*sqrt(diag);

      for (int k = 0; k < n; k++) {
        M[j*n+k] = M[j*n+k] - l * M[i*n+k];
      }
    }
    
    double diago = M[(n*n)-1];
    if (diago <= 0) {
      return -1;
    }
    A[(((n-1)*n)/2)+(n-1)] = sqrt(diago);

  }
  matrix_freigeben(M);
  return 0;
  /*
    ................................
  */
}

/* Diese Funktion realisiert die Vorwaertssubstitution
 * zum Loesen des linearen Gleichungssystems Cy = b,
 * wobei C die untere Dreiecksmatrix der Cholesky-Zerlegung
 * darstellt und in A abgespeichert ist.
 *
 * n - Anzahl Zeilen und Spalten von A (Eingabe)
 * A - Cholesky-Zerlegung der Matrix (Eingabe)
 * b - rechte Seite des Gleichungssystems Cy=b
 *     und Loesungsvektor y nach Beendigung der 
 *     Funktion. (Ein- und Ausgabe)
 *
 * Achtung: Der Vektor b wird mit der Loesung des
 *          linearen Gleichungssystems ueberschrieben.
 */
void CC_VorwSubst(long n, double *A, double *b) {
  double *M = symmat_neu(n);
  symmat_kopieren(M,A,n);
  for (int i = 0; i<n; i++) {
    for (int k = 0; k < i; k++) {
      b[i] = b[i] - M[((i*(i+1))/2)+k];
    }
    b[i] = b[i]/M[((i*(i+1))/2)+i];
    for (int k = i+1;k<n;k++) {
      M[((k*(k+1))/2)+i] = M[((k*(k+1))/2)+i] * b[i];
    }
  }
  symmat_freigeben(M);
  /*
    ................................
  */
}

/* Diese Funktion realisiert die Rueckwaertssubstitution
 * zum Loesen des linearen Gleichungssystems C^T x = b,
 * wobei C die untere Dreiecksmatrix der Cholesky-Zerlegung
 * darstellt in A abgespeiert ist.
 *
 * n - Anzahl Zeilen und Spalten von A (Eingabe)
 * A - Cholesky-Zerlegung der Matrix (Eingabe)
 * b - rechte Seite des Gleichungssystems C^T x=b
 *     und Loesungsvektor x nach Beendigung der 
 *     Funktion. (Ein- und Ausgabe)
 *
 * Achtung: Der Vektor b wird mit der Loesung des
 *          linearen Gleichungssystems ueberschrieben.
 */
void CC_RueckwSubst(long n, double *A, double *b) {
  double *M = symmat_neu(n);
  symmat_kopieren(M,A,n);
  for (int i = n-1; i>=0; i--) {
    for (int k = i+1; k < n; k++) {
      b[i] = b[i] - M[((k*(k+1))/2)+i];
    }
    b[i] = b[i]/M[((i*(i+1))/2)+i];
    for (int k = 0;k<i;k++) {
      M[((i*(i+1))/2)+k] = M[((i*(i+1))/2)+k] * b[i];
    }
  }
  symmat_freigeben(M);
  /*
    ................................
  */
}

