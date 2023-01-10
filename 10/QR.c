#include <math.h>

#include "QR.h"
#include "basic_LinAlg.h"

/* Diese Funktion berechnet die QR-Zerlegung einer regulaeren 
 * Matrix A mit Hilfe von Householder-Spiegelungen.
 *
 * n - Anzahl Zeilen und Spalten von A (Eingabe)
 * A - Matrix (Ein- und Ausgabe)
 * u - Vektor (Ausgabe)
 *
 * Achtung: Die Matrix A wird mit den Vektoren u_i aus der
 *          Householder-Spiegelung und der Matrix R
 *          ueberschrieben. Die Elemente u_i^(i) werden
 *          seperat im Vektor u abgespeichert.
 *
 * Achtung: Im Vektor u muss genuegend Platz vorhanden sein.
 *
 * Achtung: Da A regulaer ist, gilt im Algorithmus stets mu_i > 0.
 */
void QR_HS_Zerlegung(long n, double *A, double *u) {
  long i,j,k;
  double beta, sigma, lambda, mu;

  for (i=0; i<n-1; i++) {
    lambda = (A[i*n+i]<0) ? 1 : -1;

    mu = 0;
    for (k=i; k<n; k++) 
      mu += A[k*n+i]*A[k*n+i];
    mu = sqrt(mu);

    sigma = sqrt( 2*mu*(mu+fabs(A[i*n+i])) );
    u[i] = ( A[i*n+i]-lambda*mu ) / sigma;
    for (k=i+1; k<n; k++) // speichere u_i in A ab
      A[k*n+i] /= sigma;

    A[i*n+i] = lambda*mu;
    for (j=i+1; j<n; j++) {
      // berechne beta = (A_j,u_i)
      beta = A[i*n+j]*u[i];
      for (k=i+1; k<n; k++)
        beta += A[k*n+j]*A[k*n+i];

      // berechne neue Matrixeintraege a_kj
      A[i*n+j] -= 2*beta*u[i];
      for (k=i+1; k<n; k++)
        A[k*n+j] -= 2*beta*A[k*n+i];
    }
  }
}

/* Diese Funktion berechnet Q^T x, wobei Q durch die sukzessive
 * Multiplikation mit den Householder-Matrizen realisiert wird.
 * Diese Matrizen sind gegeben durch die Vektoren u_i, die in
 * der unteren Dreiecksmatrix von A abgespeichert sind, und den
 * Elementen u_i^(i), die im Vektor u gegeben sind.
 *
 * n - Anzahl Zeilen und Spalten von A (Eingabe)
 * A - QR-Zerlegung erzeugt von QR_HS_Zerlegung (Eingabe)
 * u - Vektor mit den Elementen u_i^(i) (Eingabe)
 * x - Vektor, an den Q^T multipliziert wird
 *     und Ergebnisvektor nach Beendigung der 
 *     Funktion. (Ein- und Ausgabe)
 *
 * Achtung: Der Vektor x wird mit dem Ergebnis ueberschrieben.
 */
void QR_HS_QTransX(long n, double *A, double *u, double *x) {
  long i, k;
  double beta;

  for (i=0; i<n-1; i++) {
    // berechne beta = (u_i,x)
    beta = u[i]*x[i];
    for (k=i+1; k<n; k++)
      beta += A[k*n+i]*x[k];

    // realisiere x = Q_i x
    x[i] -= 2*beta*u[i];
    for (k=i+1; k<n; k++)
      x[k] -= 2*beta*A[k*n+i];
  }
}

/* Diese Funktion realisiert die Rueckwaertssubstitution
 * zum Loesen des linearen Gleichungssystems Rx = b,
 * wobei R die obere Dreiecksmatrix der QR-Zerlegung
 * ist und zusammen mit Q kombiniert in A abgespeiert
 * ist.
 *
 * Diese Funktion ist identisch mit der Rueckwaertssubstitution
 * einer LR-Zerlegung!
 *
 * n - Anzahl Zeilen und Spalten von A (Eingabe)
 * A - QR-Zerlegung einer Matrix (Eingabe)
 * b - rechte Seite des Gleichungssystems Rx=b
 *     und Loesungsvektor x nach Beendigung der 
 *     Funktion. (Ein- und Ausgabe)
 *
 * Achtung: Der Vektor b wird mit der Loesung des
 *          linearen Gleichungssystems ueberschrieben.
 */
void QR_RueckwSubst(long n, double *A, double *b) {
  long i,j;

  for (i=n-1; i>=0; i--) {
    for (j=i+1; j<n; j++)
      b[i] -= A[i*n+j]*b[j];

    b[i] /= A[i*n+i];
  }
}

/* Diese Funktion berechnet die QR-Transformation RQ, wobei Q durch 
 * die sukzessive Multiplikation mit den Householder-Matrizen realisiert 
 * wird. Diese Matrizen sind gegeben durch die Vektoren u_i, die in
 * der unteren Dreiecksmatrix von A abgespeichert sind, und den
 * Elementen u_i^(i), die im Vektor u gegeben sind.
 *
 * n  - Anzahl Zeilen und Spalten von A (Eingabe)
 * A  - QR-Zerlegung erzeugt von QR_HS_Zerlegung (Eingabe)
 * u  - Vektor mit den Elementen u_i^(i) (Eingabe)
 * RQ - Vektor, an den Q^T multipliziert wird
 *      und Ergebnisvektor nach Beendigung der 
 *      Funktion. (Ausgabe)
 *
 * Achtung: Die Matrix RQ wird ueberschrieben, in ihr muss genuegend
 *          Speicherplatz vorhanden sein.
 */
void QR_HS_Transformation(long n, double *A, double *u, double *RQ) {
  long i,j,k;
  double beta;

  // kopiere R (aus A) nach RQ
  for (j=0; j<n; j++)
    for (k=0; k<n; k++)
      if (j<=k) 
        RQ[j*n+k] = A[j*n+k];
      else
        RQ[j*n+k] = 0;

  // Multipliziere sukzessive mit Q_i von rechts
  // hierbei bezeichne RQ_i^j die j-te Zeile aus RQ_i,
  // wobei RQ_i = R Q_0 Q_1 ... Q_{i-1}, i>1 und RQ_0 = R.
  // Somit ist RQ_n die fertige QR-Transformation
  for (i=0; i<n-1; i++) {
    for (j=0; j<n; j++) {
      // Berechne ( RQ_i^j , u_i )
      beta = RQ[j*n+i]*u[i];
      for (k=i+1; k<n; k++)
        beta += RQ[j*n+k]*A[k*n+i];

      // Berechne j-te Zeile von RQ_i Q_i
      RQ[j*n+i] -= 2*beta*u[i];
      for (k=i+1; k<n; k++) 
        RQ[j*n+k] -= 2*beta*A[k*n+i];
    }
  }
}
