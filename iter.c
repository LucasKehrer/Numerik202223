#include "iter.h"
#include "basic_LinAlg.h"
//testing 
#include "stdio.h"

/* Diese Funktion wertet Tschebyscheff-Polynome T_n. 
 *
 * n - Ordnung des Tschebyscheff-Polynoms T_n.
 * x - Auswertungspunkt.
 * Die Funktion gibt T_n(x) zurück.
 */
static double ChebPoly (long n, double x)
{
  double t_k1 = x, t_k = 1.;
  double t_k2 = 2. * x * t_k1 - t_k;

  if (n == 0)
    {
      return t_k;
    }
  else if (n == 1)
    {
      return t_k1;
    }
  else if (n == 2)
    {
      return t_k2;
    }
  
  for (long k = 2; k < n; k++)
    {
      t_k = t_k1;
      t_k1 = t_k2;
      t_k2 = 2. * x * t_k1 - t_k;
    }

  return t_k2;
}

/*
* A bit functionaliity for CSR-Matrix
*/

/*
 * Returns die Anzahl an Elementen in einer Zeile
 * ACHTUNG: Angepasst, das Row 0 = die erste ist
 * ACHTUNG: Schmeißt keinen Fehler, wenn i = n ist (bzw gibt nen outofbound halt..)
 */

long Sparse_getRowNumber(SparseMatrix_t A, long i) {
  if (i < 0) {
    return 0;
  }
  return A.rowIndex[i+1]-A.rowIndex[i];
}

/*
 * Returned die aktuelle Spalte für einen Eintrag
 */
long Sparse_getColumn(SparseMatrix_t A, long i, long j) {
  return A.columns[A.rowIndex[i]+j];
}  


/*
 * Die Nachfolgende Funktion berechnet Ax, wobei A eine CSR-Matrix ist
 *
 * A - CSR-Matrix
 * x - Vektor zum dranmultiplizieren
 * y - Vektor, in den das Ergebnis gespeichert werden soll
 */
void Sparse_Matrix_mul(SparseMatrix_t A, double *x, double *y, long n) {
  for (int i = 0; i < n; i++) {
    y[i] = 0;
    for (int j = 0; j < Sparse_getRowNumber(A, i); j++) {
      
      y[i] += A.values[A.rowIndex[i]+j] * x[Sparse_getColumn(A, i, j)];

      /* Test
       *if (i<100) {
       *printf("y %d = %f mit A.values = %d \n", i, y[i], i+Sparse_getColumn(A, i, j));
       *}
       */
    }
  }
}
/*Speichert die Inverse der Diagonale in einen Vektor*/
void Spars_getDiaginverse(SparseMatrix_t A, double *x, long n) {
  for (int i = 0; i < n; i++) {
    x[i] = 0;
    for (int j = 0; j < Sparse_getRowNumber(A, i); j++) {
      if ((i == Sparse_getColumn(A, i, j)) && (A.values[A.rowIndex[i]+j] != 0)) {
        x[i] = 1/(A.values[A.rowIndex[i]+j]);
      }
    }
  }
}

void multiply_Diavektor(double *B_inv, double *x, long n) {
  for (int i = 0; i < n; i++) {
    x[i] = x[i]*B_inv[i];  
    }
}


void Sparse_JacobiChebHelper_solver(SparseMatrix_t A, double gamma, double twodivGamma1Gamma2, double *diaInv_vektor, double *x_k, double *xk_1, double *b, long k, long n) {
//HIER
  double omega = 1 ;
  if (k != 0) {
  omega = 2 * gamma * ((ChebPoly(k,gamma)) / (ChebPoly(k+1, gamma)));
  }

  //printf("omega %d: %f \n",k+1, omega);
  double *temp = vektor_neu(n);

  //Ax_k
  Sparse_Matrix_mul(A, x_k, temp, n);
  //printf("erste Zeile Ax_k: %f \n", temp[0]);

  //b-Ax_k = a1
  for (int i = 0; i<n; i++) {
    temp[i] = b[i] - temp[i];
  }
  //printf("erste Zeile b-Ax_k: %f \n", temp[0]);

  //B^-1 * a1
  for (int i = 0; i<n; i++) {
    //passt das ? 
    temp[i] = diaInv_vektor[i] * temp[i];
  }
  //printf("erste Zeile B^-1*(b-Ax_k) mit diaelim %f: %f \n",diaInv_vektor[0], temp[0]);

  // 2/g2+g1
  for (int i = 0; i<n; i++) {
    temp[i] = twodivGamma1Gamma2 * temp[i];
  }
  //printf("erste Zeile 2/g2+g1: %f \n", temp[0]);
  // +xk -xk-1
  for (int i = 0; i<n; i++) {
    temp[i] = temp[i]+x_k[i]-xk_1[i];
  }
  //printf("erste Zeile x_k - x_k-1 (mit %f %f): %f \n",x_k[0], xk_1[0], temp[0]);

  //*omega
  for (int i = 0; i<n; i++) {
    temp[i] = omega*temp[i];
  }
  //printf("erste Zeile omega * ..: %f \n", temp[0]);
  //+xk-1
  for (int i = 0; i<n; i++) {
    temp[i] = temp[i] + xk_1[i];
  }
  //printf("+ x_k-1: %f \n", temp[0]);

  vektor_kopieren(xk_1, x_k, n);
  vektor_kopieren(x_k, temp, n);

  vektor_freigeben(temp);
  //printf("passt das mit dem Kopieren?: %f \n", xk_1[0]);
}


/* Diese Funktion realisiert das Jacobi-Verfahren mit
 * Tschebyscheff-Beschleunigung zur Approximation der Loesung x eines linearen
 * Gleichungssystems.
 *
 * A - Matrix der Dimension nxn (Eingabe)
 * x - Startvektor der Itteration und anschliessend
 *     Approximation des Loesungsvektors (Ein- und Ausgabe)
 * b - Rechte Seite des Gleichungssystems (Eingabe)
 * n - Anzahl Zeilen und Spalten von A (Eingabe)
 * iter - Anzahl der durchzufuehrenden Iterationen (Eingabe)
 * gamma_1 - Parameter aus Aufgabe 6.3 
 * gamma_2 - Parameter aus Aufgabe 6.3
 */
void Sparse_JacobiCheb_solver(SparseMatrix_t A, double *x, double *b, long n,
			      long iter, double gamma_1, double gamma_2)
{
  double gamma = (-1)*((gamma_2+gamma_1)/(gamma_2-gamma_1));
  //printf("gamma: %f \n", gamma);
  double twodivGamm1Gamma2 = 2/(gamma_2+gamma_1);
  //printf("1/g2+g1: %f \n", twodivGamm1Gamma2);
  double *B_inv = vektor_neu(n);
  Spars_getDiaginverse(A, B_inv, n);

  double *x_km1 = vektor_neu(n);

  for (int i = 0; i < n; i++) {
    x[i] = 0;
    x_km1[i] = 0; 
    }

//minimum eine Iteration
  Sparse_JacobiChebHelper_solver(A, gamma, twodivGamm1Gamma2, B_inv, x, x_km1, b, 0, n); //letztes passt? k ist aktueller step
  long k = iter -1;
  while (k >= 1) {
    Sparse_JacobiChebHelper_solver(A, gamma, twodivGamm1Gamma2, B_inv, x, x_km1, b, iter - k, n);
    k -= 1;
  }
  //printf("passt kopie?: %f \n", x[0]);
  vektor_freigeben(B_inv);
  vektor_freigeben(x_km1);

}

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
void Sparse_Jacobi_solver(SparseMatrix_t A, double *x, double *b, long n, long iter)
{

  double *B_inv = vektor_neu(n);
  Spars_getDiaginverse(A, B_inv, n);

  double *x_km1 = vektor_neu(n);
  for (int i = 0; i < n; i++) {
    x[i] = 0;
    x_km1[i] = 0; 
    }

  Sparse_JacobiChebHelper_solver(A, 1, 1, B_inv, x, x_km1, b, 0, n); //letztes passt? k ist aktueller step
  long counter = iter -1;
  while (counter >= 1) {
    Sparse_JacobiChebHelper_solver(A, 1, 1, B_inv, x, x_km1, b, 0, n);
    counter -= 1;
  }

  vektor_freigeben(B_inv);
  vektor_freigeben(x_km1);
}


