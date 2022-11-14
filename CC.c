#include <math.h>
#include "CC.h"

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
  /*
    ................................
  */
}

/* Diese Funktion realisiert die Vorwaertssubstitution
 * zum Loesen des linearen Gleichungssystems Cy = b,
 * wobei C die untere Dreiecksmatrix der Cholesky-Zerlegung
 * darstellt und in A abgespeiert ist.
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
  /*
    ................................
  */
}

