#include <math.h>
#include <stdio.h>
#include "cg.h"
#include "basic_LinAlg.h"

#define OMEGA 1.581120979521361

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

