#include "iter.h"
#include "basic_LinAlg.h"

/* Diese Funktion realisiert das Jacobi-Verfahren zur
 * Approximation der Loesung x eines linearen Gleichungssystems.
 *
 * A - Matrix der Dimension nxn (Eingabe)
 * x - Startvektor der Itteration und anschliessend
 *     Approximation des Loesungsvektors (Ein- und Ausgabe)
 * b - Rechte Seite des Gleichungssystems (Eingabe)
 * n - Anzahl Zeilen und Spalten von A (Eingabe)
 * iter - Anzahl der durchzufuehrenden Iterationen (Eingabe)
 *
 */
void Jacobi_solver(double *A, double *x, double *b, long n, long iter)
{
  // Vektor Ausgabe
  // TODO muss in main!
  /*
  if(iter == 60 || iter == 40 || iter == 20) {
    vektor_ausgeben(x, n, "%10.3e");
    double *C = vektor_neu(n);
    vektor_kopieren(C,b,n);
    for(int i = 0; i<n; i++) {
      C[i] = -C[i];
    }
    matrix_mult(1,1,A,x,C,n,1,n);
    vektor_ausgeben(C,n,"%10.3e");
  }
  */
  
  if(iter == 0) {
    return;
  } else {
    double *xkop = vektor_neu(n);
    vektor_kopieren(xkop, x, n);
    for(int i = 0; i<n; i++) {
      double aii = A[i*n+i];
      double sum = 0;
      for (int j = 0; j<n; j++) {
        if (j==i) {
          continue;
        } else {
          sum += A[i*n+j] * xkop[j] - b[i];
        }
      }
      x[i] = -sum/aii;
    }

    vektor_freigeben(xkop);
    Jacobi_solver(A,x,b,n,iter - 1);
  }
}

/* Diese Funktion realisiert das Gauss-Seidel-Verfahren zur
 * Approximation der Loesung x eines linearen Gleichungssystems.
 *
 * A - Matrix der Dimension nxn (Eingabe)
 * x - Startvektor der Itteration und anschliessend
 *     Approximation des Loesungsvektors (Ein- und Ausgabe)
 * b - Rechte Seite des Gleichungssystems (Eingabe)
 * n - Anzahl Zeilen und Spalten von A (Eingabe)
 * iter - Anzahl der durchzufuehrenden Iterationen (Eingabe)
 *
 * Hinweis: Spezialfall des Relaxations-Verfahren (SOR)
 */
void GaussSeidel_solver(double *A, double *x, double *b, long n, long iter)
{
  SQR_solver(A,x,b,n,iter,1);
}

/* Diese Funktion realisiert das Relaxations-Verfahren zur
 * Approximation der Loesung x eines linearen Gleichungssystems.
 *
 * A - Matrix der Dimension nxn (Eingabe)
 * x - Startvektor der Itteration und anschliessend
 *     Approximation des Loesungsvektors (Ein- und Ausgabe)
 * b - Rechte Seite des Gleichungssystems (Eingabe)
 * n - Anzahl Zeilen und Spalten von A (Eingabe)
 * iter - Anzahl der durchzufuehrenden Iterationen (Eingabe)
 * omega - Relaxationsparameter, i.A. aus dem Intervall (0,1]
 *         fuer symmetrisch positiv definite Matrizen kann man
 *         omega aus dem Intervall (0,2) waehlen (Eingabe)
 *
 * Hinweis: Spezialfall des Relaxations-Verfahren (SOR)
 */
void SOR_solver(double *A, double *x, double *b, long n, long iter,
                double omega)
{
  if(iter == 0) {
    return;
  } else {
    double *xkop = vektor_neu(n);
    vektor_kopieren(xkop, x, n);
    for(int i = 0; i<n; i++) {
      double aii = A[i*n+i];
      double sum1 = 0;
      for (int j = 0; j<i; j++) {
        sum1 += A[i*n+j] * x[j];
      }
      double sum2 = 0;
      for (int j = i+1; j<n; j++) {
        sum2 += A[i*n+j] * xkop[j];
      }
      x[i] = (1-omega)*xkop[i] - (omega/aii)*sum1+sum2-b[i];
    }

    vektor_freigeben(xkop);
    SQR_solver(A,x,b,n,iter - 1, omega);
}

