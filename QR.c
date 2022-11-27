#include <math.h>

#include "QR.h"

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
    ...................
  */

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
