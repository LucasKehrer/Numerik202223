#include <math.h>
#include <stdio.h>
#include "cg.h"
#include "basic_LinAlg.h"

#define OMEGA 1.581120979521361

void solve_diagonal(SparseMatrix_t A, double *x, long n){
    for (int i = 0; i < n; i++) {
    for (int j = 0; j < Sparse_Zeilenindex(A, i); j++) {
            if(Sparse_Spaltenindex(A, i, j)==i){
                        x[i]=x[i]*A.values[A.rowIndex[i]+j];
            }
    }
  }
}

void solve_lowerSparse(SparseMatrix_t A, double* x, double* y, long n, double omega){
}
void solve_upperSparse(SparseMatrix_t A, double* x, double* y, long n, double omega){
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
    double* temp=vektor_neu(n);
    solve_lowerSparse(A, r, temp, n, omega);
    solve_diagonal(A, temp, n);
    solve_upperSparse(A, temp, w, n, omega);
    vektor_skalieren(2-omega, w, n);

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
    double* r = vektor_neu(n);
    double* s = vektor_neu(n);
    vektor_kopieren(r, x, n);
    SparseMatrix_vektor_mult(1, -1, A, r, b, n);
    int cur_iter=0;
    if(prec!=0){
        double* w = vektor_neu(n);
        Sparse_SSOR_prec(A, w, r, n, OMEGA);
        vektor_kopieren(s, w, n);
        double alpha;
        double beta;
        double* w_temp=vektor_neu(n);
        double* r_temp=vektor_neu(n);
        double* temp=vektor_neu(n);
        vektor_kopieren(w_temp, w, n);
        vektor_kopieren(r_temp, r, n);
        while(cur_iter<iter && vektor_skalprod(r, r, n) < eps*vektor_skalprod(b, b, n)){
            vektor_kopieren(temp, s, n);
            SparseMatrix_vektor_mult(1, 0, A, temp, temp, n);
            alpha=vektor_skalprod(w, r, n)/vektor_skalprod(temp, s, n);
            vektor_addieren(-alpha, s, x, n);
            vektor_addieren(-alpha, temp, r_temp, n);
            Sparse_SSOR_prec(A, w, r, n, OMEGA);
            beta=vektor_skalprod(w_temp, r_temp, n)/vektor_skalprod(w, r, n);
            vektor_skalieren(beta, s, n);
            vektor_addieren(1, w_temp, s, n);
            vektor_kopieren(w, w_temp, n);
            vektor_kopieren(r, r_temp, n);
            cur_iter++;
        }
        vektor_freigeben(w_temp);
        vektor_freigeben(r_temp);
        vektor_freigeben(w);
        vektor_freigeben(temp);
        return cur_iter;
    }
    else{
        vektor_kopieren(s, r, n);
        double alpha;
        double beta;
        double* r_temp=vektor_neu(n);
        double* temp=vektor_neu(n);
        vektor_kopieren(r_temp, r, n);
        vektor_kopieren(temp, s, n);
        while(cur_iter<iter && vektor_skalprod(r, r, n) < eps*vektor_skalprod(b, b, n)){
            SparseMatrix_vektor_mult(1, 0, A, temp, temp, n);
            alpha=vektor_skalprod(r, s, n)/vektor_skalprod(temp, s, n);
            vektor_addieren(-alpha, s, x, n);
            vektor_addieren(-alpha, temp, r_temp, n);
            beta=vektor_skalprod(r_temp, r_temp, n)/vektor_skalprod(r, r, n);
            vektor_skalieren(beta, s, n);
            vektor_addieren(1, r_temp, s, n);
            vektor_kopieren(r, r_temp, n);
            cur_iter++;
        }
        vektor_freigeben(r_temp);
        vektor_freigeben(temp);
        return cur_iter;
    }
    vektor_freigeben(r);
    vektor_freigeben(s);
}

