#include <math.h>
#include <stdio.h>

#include "QR.h"


/*
* Helper function for sig(x)
*/
int signum(double x) {
  if (x >= 0) {
    return 1;
  } else {
    return -1;
  }
}



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
void QR_HS_Zerlegung(long n, double *A, double *u)
{
/*
  

    // Für i = 1 bis i = n-1
 for (int i = 0; i < n-1; i++) {
    int myi_sum = 0;
    for (int k = i; k<n ; k++) {
        myi_sum += A[k*n+i]*A[k*n+i];
    }
    double myi = sqrt(myi_sum);

    double lami = -signum(A[i*n+i]);

  if (myi == 0) {
    u[i] = 1;
  } else {

    double sigi = sqrt(2*myi*(myi + abs(A[i*n+i])));
    //süß, ne
    u[i] = (A[i*n+i] - lami*myi)/sigi;
    
    // U_k^(i)
    for (int k = i+1; k<n; k++) {
      u[k] = A[k*n+1]/sigi;
    }

    //a_ii
    A[i*n+i] = lami*myi;

    

    for (int j = i; j<n; j++) {
      double betai = 0;
      for (int k = i; k<n; k++ ) {
        betai += A[k*n+j]*u[k];
      }

      for (int k = i; k<n; k++) {
        A[k*n+j] = A[k*n+j]- 2*betai*u[k];
      }

    }


    //Speiche u_k^(i)
    for (int k= i+1; k<n; k++) {
      printf("Hallo %d %f \n", i, u[k]);
      A[k*n+i] = u[k];
    }
    



  }   
 } 
*/

 for (int i = 0; i<n-1; i++) {
    double betai = 0;
    for (int k = i; k<n; k++) {
      //myi
      betai += A[k*n+i]*A[k*n+i];
    }
      double alpha = sqrt(betai);
      if (alpha == 0) {
        u[i] = 0; //nicht 1?
      } else {
        double c_i = 1/(betai + alpha * abs(A[i*n+i]));
        if (A[i*n+i] < 0) {
          alpha = alpha * (-1);
        }
        u[i] = A[i*n+i]+alpha;
        printf("Test 2 = %f \n", (-1)*alpha);
        A[i*n+i] = (-1) * alpha;
        for (int k = i+1; k<n; k++) {
          double sigma = u[i]*A[i*n+k];
          for (int j = i+1; j<n; j++) {
            sigma+= A[j*n+i]*A[j*n+k];
          } 
          sigma = sigma * c_i;
          A[i*n+k] -= u[i] * sigma;

          for (int j = i+1; j<n; j++) {
            A[j*n+k] -= A[j*n+i]*sigma;
          } 
        }


    
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
void QR_HS_QTransX(long n, double *A, double *u, double *x)
{

  /*
    ...................
  */


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
void QR_RueckwSubst(long n, double *A, double *b)
{

  /*
    ...................
  */

}
