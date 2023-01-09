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
  double* T = vektor_neu(n);
  solve_lowerSparse(A,r,T,n,omega);
  //D*(D+wL)^-1
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < Sparse_getRowNumber(A, i); j++) {
      if(Sparse_getColumn(A, i, j)==i){
        T[i]=T[i]*A.values[A.rowIndex[i]+j];
      }
    }
  }
  solve_upperSparse(A,T,w,n,omega);
  for (int i = 0; i < n; i++) {
    w[i] = (2-omega) * w[i];
  }
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
  if(prec != 0) {
    return Sparse_cg_prec(A,x,b,n,prec,iter,eps);
  }
  long iterN = 0;
  double* r = vektor_neu(n);
  vektor_kopieren(r,b,n);
  SparseMatrix_vektor_mult(1,-1, A, x, r, n);
  double* s = vektor_neu(n);
  vektor_kopieren(s, r, n);
  
  double alpha = 0;
  double beta = 0;
  double* z = vektor_neu(n);
  double* r_temp = vektor_neu(n);
  while ((iterN < iter) && (vektor_skalprod(r, r, n) >= eps*vektor_skalprod(b, b, n))) {
    SparseMatrix_vektor_mult(1,0,A,s,z,n);
    alpha = vektor_skalprod(r,s,n)/vektor_skalprod(z, s, n);
    vektor_kopieren(r_temp, r, n);
    for (int i = 0; i < n; i++) {
      x[i] = x[i] - alpha*s[i];
      r[i] = r[i] - alpha*z[i];
    }
    beta = vektor_skalprod(r,r,n) / vektor_skalprod(r_temp,r_temp,n);
    for (int i = 0; i<n; i++) {
      s[i] = r[i] + beta * s[i]; 
    }
    iterN+=1;
  }

  vektor_freigeben(r);
  vektor_freigeben(s);
  vektor_freigeben(z);
  vektor_freigeben(r_temp);
  return iterN;
}

//mit preconditioner
long Sparse_cg_prec(SparseMatrix_t A, double *x, double *b, long n, int prec, long iter, double eps) {
  long iterN = 0;
  double* r = vektor_neu(n);
  vektor_kopieren(r,b,n);
  SparseMatrix_vektor_mult(1,-1, A, x, r, n);
  double* s = vektor_neu(n);
  
  
  double alpha = 0;
  double beta = 0;
  double* z = vektor_neu(n);
  double* r_temp = vektor_neu(n);
  double* w_temp = vektor_neu(n);

  double* w = vektor_neu(n);
  Sparse_SSOR_prec(A, w, r, n, OMEGA);
  vektor_kopieren(s, w, n);

  while ((iterN < iter) && (vektor_skalprod(r, r, n) >= eps*vektor_skalprod(b, b, n))) {
    SparseMatrix_vektor_mult(1,0,A,s,z,n);
    alpha = vektor_skalprod(w,r,n)/vektor_skalprod(z, s, n);
    vektor_kopieren(r_temp, r, n);
    for (int i = 0; i < n; i++) {
      x[i] = x[i] - alpha*s[i];
      r[i] = r[i] - alpha*z[i];
    }
    vektor_kopieren(w_temp, w, n);
    Sparse_SSOR_prec(A, w, r, n, OMEGA);
    beta = vektor_skalprod(w,r,n) / vektor_skalprod(w_temp,r_temp,n);
    for (int i = 0; i<n; i++) {
      s[i] = w[i] + beta * s[i]; 
    }
    iterN+=1;
  }

  vektor_freigeben(r);
  vektor_freigeben(s);
  vektor_freigeben(z);
  vektor_freigeben(r_temp);
  return iterN;
}



/*
 * Inverse untere Dreiecksmatrix x = L^-1 y 
 * y = Lösung Lx = y
 * L = untere Dreiecksmatrix
 * x = Da wird das Ergebnis gespeichert
*/
void solve_lowerSparse(SparseMatrix_t A, double *y, double *x, long n, double omega) {

    //Lösche x
    for (int i = 0; i < n; i++) {
      x[i] = 0;
    }

    for (int i = 0; i < n ; i++) {
      x[i] = y[i];
       

      for (int j = 0; j < Sparse_getRowNumber(A, i); j++) {
        if (Sparse_getColumn(A, i, j) < i) {
          x[i]-= omega * A.values[A.rowIndex[i]+j] * x[Sparse_getColumn(A, i, j)];
        }
      }
    //Hole diaeintrag
    for (int j2 = 0; j2 < Sparse_getRowNumber(A, i); j2++) {
        if (Sparse_getColumn(A, i, j2) == i) {
          x[i] *= 1/(A.values[A.rowIndex[i]+j2]);
        }
    }
    
  }
}


void solve_upperSparse(SparseMatrix_t A, double *y, double *x, long n, double omega) {

    //Lösche x
    for (int i = n-1; i >= 0; i--) {
      x[i] = 0;
    }
    
    for (int i = n-1; i >= 0; i--) {
      x[i] += y[i];
      for (int j = 0; j < Sparse_getRowNumber(A, i); j++) {
        if (Sparse_getColumn(A,i,j)>i) {
          x[i] -= (omega * A.values[A.rowIndex[i]+j] * x[Sparse_getColumn(A,i,j)]);
        }
      }
      

      for (int j = 0; j < Sparse_getRowNumber(A, i); j++) {
        if (Sparse_getColumn(A,i,j) == i) {
          x[i] *= 1/(A.values[A.rowIndex[i]+j]);
        }
      }
    }
    
  }

