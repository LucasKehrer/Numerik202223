#include <math.h>
#include <stdio.h>

#include "QR.h"
#include "basic_LinAlg.h"


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
for (int j = 0; j < n-1; j++){

  double beta = 0;
  for (int k = j; k<n; k++) {
    beta+=A[k*n+j]*A[k*n+j];
  }

  double alpha = sqrt(beta);

  if (alpha == 0) {
    u[j] = 0;
  } else {
    double c_j = 1/(beta+alpha*fabs(A[j*n+j]));
    if (A[j*n+j] < 0) {
      alpha = -alpha;
    } 
    u[j] = A[j*n+j];
    A[j*n+j]= (-alpha);

    for (int k = j+1; k<n; k++) {
      double sum = u[j] * A[j*n+k];
      for (int i = j+1; i<n; i++) {
        sum+=sum*c_j;
      }
      sum = sum * c_j;
      A[j*n+k] = A[j*n+k]-u[j]*sum;
      for (int i = j+1; i<n; i++) {
        A[i*n+k] = A[i*n+k] - A[i*n+j] * sum;
      }
    }
  }

printf("DEBUG: \n");
matrix_ausgeben(A, 4 ,4, "% 10.3e");
printf("\n");
vektor_ausgeben(u, 4, "% 10.3e");

}

*/

  //für alle HAuptminore außer der letzte
  for (int i = 0; i < (n-1); i++) {
    
    //berechne Norm erster Spalte
    double alpha = 0;
    for (int k = i; k < n; k++ ) {
      alpha+=A[k*n+i]*A[k*n+i];
    }
    alpha = sqrt(alpha);

    //printf("Alpha %d=%f\n", i, alpha);

    
    //hole u
    //u_1^(1) (zwischeschritt) = aii - alpha
    u[i] = A[i*n+i] - alpha;

    //printf("A %f\n", A[i*n+i] );
    //printf("u %f\n", u[i] );
    //hole norm splate 1 mit verändertem a
    double phi = u[i]*u[i];
    //printf("phiTemp %d=%f\n", i, phi);
    for (int k = i+1; k < n; k++ ) {
      phi += A[k*n+i]*A[k*n+i];
      //printf("phiTemp %d=%f\n", i, phi);
    } 

    phi = sqrt(phi);
    //printf("phi %d=%f\n", i, phi);

    //speichere u's sonstige
    u[i] = u[i]/phi;
    //printf("in %d u_%d = %f\n",i,i, u[i] );

    for (int k = i+1; k < n; k++ ) {
      // printf("A %f phi %f div %f\n", A[k*n+i], phi, A[k*n+i]/ phi );
      double t = A[k*n+i]/ phi;
      u[k] = t;
     
     // printf("in %d u_%d = %f\n",i,k, u[k] );
    }

    //berechne a's
    double * A2 = matrix_neu(4,4); 
    matrix_kopieren(A2, A, 4, 4);
     for (int k = i; k<n; k++) {
        for (int j = i; j<n; j++) { 
        
        double beta = 0;
        for (int g = i; g< n; g++) {
          beta += u[k]*u[g]*A2[g*n+j];
       //   printf("vv^t_i = k=%d, g=%d; %f %f = %f, A2[g*n+j] = %f produkt = %f\n", k, g, u[k],u[g],u[k]*u[g] ,A2[g*n+j], u[k]*u[g]*A2[g*n+j]);
        }
       // printf("beta= %f \n", beta);
        double t = A2[k*n+j] - 2*beta;
       // printf("old A (%d,%d)= %f newA= %f \n",k,j,A2[k*n+j], t);
        A[k*n+j] = t; 
       // printf("-----------------\n");
        }
        }

      for (int k = i+1; k < n; k++ ) {
        A[k*n+i] = u[k];
      }
  	  matrix_freigeben(A2);
      //matrix_ausgeben(A, 4, 4, " % 10.3e");

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
  double * x2 = vektor_neu(n); 
  double * A2 = matrix_neu(n,n);
  for (int j = 0; j<n; j++){
      A2[j*n+j] = 1;      
  }

  for (int i = n-2; i>=0; i--) {

    double * A3 = matrix_neu(4,4);
   for (int j = 0; j<n; j++){
      A3[j*n+j] = 1;      
    }
  
  for (int k = i; k<n; k++) {
      if (k==i) {
        x2[k]=u[k];
      } else {
        x2[k] = A[k*n+i];
      }  
    }
  

    for (int k = i; k<n; k++) {
        for (int j = i; j<n; j++) { 
        
        double beta = 0;
        for (int g = i; g< n; g++) {
          beta += x2[k]*x2[g]*A2[g*n+j];
      // printf("vv^t_i = k=%d, g=%d; %f %f = %f, A2[g*n+j] = %f produkt = %f\n", k, g, u[k],u[g],u[k]*u[g] ,A2[g*n+j], u[k]*u[g]*A2[g*n+j]);
        }
      //printf("beta= %f \n", beta);
        double t = A2[k*n+j] - 2*beta;
       //printf("old A (%d,%d)= %f newA= %f \n",k,j,A2[k*n+j], t);
        A3[k*n+j] = t; 
       //printf("-----------------\n");
        }
        }
    matrix_kopieren(A2, A3, n, n);
    matrix_freigeben(A3);

  }  
  printf("\n PROBE: Das ist Q:\n");
  matrix_ausgeben(A2, n, n, " % 10.3e");

  matrix_vektor_mult(1, 0, A2, x, x2, 4, 4);
  vektor_kopieren(x,x2, n);
  vektor_freigeben(x2);
  matrix_freigeben(A2);
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

 for (int i = n-1; i>=0; i--) {
  double t = 0;
  for (int j = i+1; j<n; j++) {
    t -= b[j]*A[i*n+j];
    //printf ("t = %f\n", t);
  }
  b[i] = (b[i] + t) /A[i*n+i];
  //printf("Le b = %f\n", b[i]);
 }

}
