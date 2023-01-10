#include "LR.h"

/* Diese Funktion berechnet die LR-Zerlegung einer 
 * streng regulaeren Matrix A.
 *
 * n - Anzahl Zeilen und Spalten von A (Eingabe)
 * A - Matrix (Ein- und Ausgabe)
 *
 * Achtung: Die Matrix A wird mit den Matrizen
 *          L und R ueberschrieben.
 */
void LR_Zerlegung(long n, double *A) {
  long i,j,k;

  for (i=0; i<n; i++) {
    // Berechnung von R
    for (j=i; j<n; j++) 
      for (k=0; k<=i-1; k++) 
        A[i*n+j] -= A[i*n+k]*A[k*n+j];

    // Berechnung von L
    for (j=i+1; j<n; j++) {
      for (k=0; k<=i-1; k++) 
        A[j*n+i] -= A[j*n+k]*A[k*n+i];
      A[j*n+i] /= A[i*n+i];
    }
  }
}

/* Diese Funktion realisiert die Vorwaertssubstitution
 * zum Loesen des linearen Gleichungssystems Ly = b,
 * wobei L die untere Dreiecksmatrix der LR-Zerlegung
 * ist und zusammen mit R kombiniert in A abgespeiert
 * ist.
 *
 * n - Anzahl Zeilen und Spalten von A (Eingabe)
 * A - LR-Zerlegung einer Matrix (Eingabe)
 * b - rechte Seite des Gleichungssystems Ly=b
 *     und Loesungsvektor y nach Beendigung der 
 *     Funktion. (Ein- und Ausgabe)
 *
 * Achtung: Der Vektor b wird mit der Loesung des
 *          linearen Gleichungssystems ueberschrieben.
 */
void LR_VorwSubst(long n, double *A, double *b) {
  long i,j;

  for (i=1; i<n; i++) 
    for (j=0; j<i; j++)
      b[i] -= A[i*n+j]*b[j];
}

/* Diese Funktion realisiert die Rueckwaertssubstitution
 * zum Loesen des linearen Gleichungssystems Rx = b,
 * wobei R die obere Dreiecksmatrix der LR-Zerlegung
 * ist und zusammen mit L kombiniert in A abgespeiert
 * ist.
 *
 * n - Anzahl Zeilen und Spalten von A (Eingabe)
 * A - LR-Zerlegung einer Matrix (Eingabe)
 * b - rechte Seite des Gleichungssystems Rx=b
 *     und Loesungsvektor x nach Beendigung der 
 *     Funktion. (Ein- und Ausgabe)
 *
 * Achtung: Der Vektor b wird mit der Loesung des
 *          linearen Gleichungssystems ueberschrieben.
 */
void LR_RueckwSubst(long n, double *A, double *b) {
  long i,j;

  for (i=n-1; i>=0; i--) {
    for (j=i+1; j<n; j++)
      b[i] -= A[i*n+j]*b[j];

    b[i] /= A[i*n+i];
  }
}

