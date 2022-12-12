#include "iter.h"
#include "basic_LinAlg.h"

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
}

