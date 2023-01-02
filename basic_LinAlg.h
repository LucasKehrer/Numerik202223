/*
 *  Praktische Mathematik
 *  im Sommersemester 2015
 *  von Steffen Weisser
 *
 *  29. Mai 2015
 *
 *  Version 3.0
 *
 *  Dieses Paket zum Rechnen mit Matrizen und Vektoren
 *  basiert auf der Datei matrix2.c aus der Vorlesung
 *  Modellieren und Programmieren vom WS 2014/15.
 *
 */

#ifndef BASIC_LINALG_H
#define BASIC_LINALG_H

/********************************************************************************/
/*********** Definitionen von Datenstrukturen ***********************************/
/********************************************************************************/

typedef struct SparseMatrix_s {
  double *values;  // Eintr�ge
  long *columns;   // Spaltenindices
  long *rowIndex;  // Index des ersten Eintrags pro Zeile
} SparseMatrix_t;


/********************************************************************************/
/*********** Arbeiten mit Vektoren **********************************************/
/********************************************************************************/

/* Diese Funktion legt einen Vektor der Laenge n dynamisch an. 
 * Sie gibt einen Zeiger auf den Vektor zurueck. 
 */
double *vektor_neu(long n);   // Laenge des Vektors

/* Es wird der Vektor x, der zuvor mit vektor_neu() angelegt
 * wurde, wieder frei gegeben.
 */
void vektor_freigeben(double *x);  // Vektor zum freigeben

/* Es wird der Vektor y nach x kopiert. Hierbei muss in x 
 * genuegend Speicher vorhanden sein.
 */
void vektor_kopieren(double *x,         // Ziel (Rueckgabe)
                     const double *y,   // Quelle
                     long n);           // Laenge

/* Der Vektor x wird skalliert:  x = alpha * x
 */
void vektor_skalieren(double alpha,  // Skallierungsparameter
                      double *x,     // Vektor x (Ein- und Ausgabe)
                      long n);       // Laenge

/* Die Funktion berechnet die 2-Norm des Vektors x. */
double vektor_2Norm(const double *x,   // Matrix
                    long n);           // Laenge

/* Es wird das Skalarprodukt von den Vektoren x und y 
 * berechnet und als Rueckgabewert zurueck gegeben.
 */
double vektor_skalprod(const double *x,  // Vektor x
                       const double *y,  // Vektor y
                       long n);          // Laenge

/* Es werden zwei Vektoren addiert und das Ergebnis wird in y 
 * zurueck gegeben. Genauer gilt:   y = alpha * x + y
 */
void vektor_addieren(double alpha,      // reelle Zahl mit der x skalliert wird
                     const double *x,   // Vektor x
                     double *y,         // Vektor y (Eingabe und Rueckgabe)
                     long n);           // Laenge

/********************************************************************************/
/*********** Arbeiten mit Matrizen **********************************************/
/********************************************************************************/

/* Diese Funktion legt eine Matrix mit n Zeilen und m Spalten 
 * dynamisch an. Sie gibt einen Zeiger auf die Matrix zurueck. 
 */
double *matrix_neu(long m,    // Anzahl der Zeilen
                   long n);   // Anzahl der Spalten

/* Es wird die Matrix A, die zuvor mit matrix_neu() angelegt
 * wurde, wieder frei gegeben.
 */
void matrix_freigeben(double *A);  // Matrix zum freigeben

/* Es wird die Matrix B nach A kopiert. Hierbei muss in A 
 * genuegend Speicher vorhanden sein.
 */
void matrix_kopieren(double *A,         // Ziel (Rueckgabe)
                     const double *B,   // Quelle
                     long m,            // Anzahl Zeilen
                     long n);           // Anzahl Spalten

/* Die Funktion berechnet die Frobeniusnorm von A. */
double matrix_FrobeniusNorm(const double *A,   // Matrix
                            long m,            // Anzahl Zeilen
                            long n);           // Anzahl Spalten

/* Es werden zwei Matrizen addiert und das Ergebnis wird in B 
 * zurueck gegeben. Genauer gilt:   B = alpha * A + B
 */
void matrix_addieren(double alpha,      // reelle Zahl mit der A skalliert wird
                     const double *A,   // Matrix A
                     double *B,         // Matrix B (Eingabe und Rueckgabe)
                     long m,            // Anzahl Zeilen
                     long n);           // Anzahl Spalten

/* Multiplikation und Addition von Matrizen. Es gilt

   C = alpha * A*B + beta * C

   Hierbei haben die Matrizen die folgenden Anzahlen von Zeilen und Spalten:

     A  -  m x l
     B  -  l x n
     C  -  m x n

   ACHTUNG: Die Matrix C muss bereits angelegt sein!
*/
void matrix_mult(double alpha,     // reelle Zahl mit der A*B skalliert wird
                 double beta,      // reelle Zahl mit der C skalliert wird
                 const double *A,  // Matrix A
                 const double *B,  // Matrix B 
                 double *C,        // Matrix C (Eingabe und Rueckgabe)
                 long m,           // Anzahl Zeilen
                 long n,           // Anzahl Zeilen
                 long l);          // Anzahl Spalten

/* Matrix-Vektor-Multiplikation. Es gilt

   y = alpha * A*x + beta * y

   Hierbei haben die Matrizen die folgenden Anzahlen von Zeilen und Spalten:

     A  -  m x n
     x  -  n
     y  -  m
*/
void matrix_vektor_mult(double alpha,     // reelle Zahl mit der A*x skalliert wird
                        double beta,      // reelle Zahl mit der y skalliert wird
                        const double *A,  // Matrix A
                        const double *x,  // Vektor x 
                        double *y,        // Vektor y (Eingabe und Rueckgabe)
                        long m,           // Anzahl Zeilen von A
                        long n);          // Anzahl Spalten von A

/********************************************************************************/
/*********** Lesen und Schreiben von Matrizen und Vektoren **********************/
/********************************************************************************/

/* Diese Funktion gibt den Vektor a auf dem Bildschirm aus.
 * Das Ausgabeformat kann z.B. mit format=" % 10.3e" veraendert
 * werden. */
void vektor_ausgeben(const double *a,      // Matrix
                     long n,               // Laenge
                     const char *format);  // Formatierung fuer printf

/* Diese Funktion gibt die Matrix A auf dem Bildschirm aus.
 * Das Ausgabeformat kann z.B. mit format=" % 10.3e" veraendert
 * werden. */
void matrix_ausgeben(const double *A,      // Matrix
                     long m,               // Anzahl Zeilen
                     long n,               // Anzahl Spalten
                     const char *format);  // Formatierung fuer printf

/* Diese Funktion gibt die symmetrische Matrix A auf dem Bildschirm aus. 
 * Hierbei wird beruecksichtigt, dass lediglich die untere Dreiecksmatrix
 * von A abgespeichert ist. Das Ausgabeformat kann z.B. mit format=" % 10.3e" 
 * veraendert werden. */
void symmat_ausgeben(const double *A,      // Matrix
                     long n,               // Anzahl Zeilen und Spalten
                     const char *format);  // Formatierung fuer printf

/* Diese Funktion gibt die symmetrische Matrix A (oder untere Dreiecksmatrix)
 * auf dem Bildschirm aus, wobei lediglich die Eintraege auf und unterhald 
 * der Diagonale ber"ucksichtigt werden. Das Ausgabeformat kann z.B. mit 
 * format=" % 10.3e" veraendert werden. */
void symmat_ausgeben_dreieck(const double *A,      // Matrix
                             long n,               // Anzahl Zeilen und Spalten
                             const char *format);  // Formatierung fuer printf

/* Die Funktion liest eine Matrix von der Standardeingabe. */
void matrix_einlesen(double **A,   // Rueckgabe: Matrix (Call by Reference)
                     long *m,      // Rueckgabe: Anzahl Zeilen (Call by Reference)
                     long *n);     // Rueckgabe: Anzahl Spalten (Call by Reference)

/* Die Funktion laed eine Matrix aus einer Datei im ascii-Format. */
int matrix_laden_ascii(char* dateiname, // Datei aus der die Matrix geladen wird
                       double **A,  // Rueckgabe: Matrix (Call by Reference)
                       long *m,     // Rueckgabe: Anzahl Zeilen (Call by Reference)
                       long *n);    // Rueckgabe: Anzahl Spalten (Call by Reference)

/* Die Funktion laed eine Matrix aus einer Datei im bin-Format. */
int matrix_laden_bin(char* dateiname, // Datei aus der die Matrix geladen wird
                     double **A,   // Rueckgabe: Matrix (Call by Reference)
                     long *m,      // Rueckgabe: Anzahl Zeilen (Call by Reference)
                     long *n);     // Rueckgabe: Anzahl Spalten (Call by Reference)

/* Die Funktion speichert eine Matrix in einer Datei im ascii-Format. */
int matrix_speichern_ascii(char* dateiname,      // Datei zum speichern der Matrix
                           double *A,            // Matrix
                           long m,               // Anzahl Zeilen
                           long n,               // Anzahl Spalten
                           const char *format);  // Formatierung fuer fprintf

/* Die Funktion speichert eine Matrix in einer Datei im bin-Format. */
int matrix_speichern_bin(char* dateiname, // Datei aus der die Matrix geladen wird
                         double *A,       // Matrix
                         long m,          // Anzahl Zeilen
                         long n);         // Anzahl Spalten

/* Die Funktion laed eine sparse Matrix aus einer Datei im ascii-Format. */
int SparseMatrix_laden_ascii(char* dateiname,   // Datei aus der die Matrix geladen wird
                             SparseMatrix_t *A, // Rueckgabe: Matrix (Call by Reference)
                             long *n);          // Rueckgabe: Anzahl der Zeilen (Call by Reference)

/* Die Funktion laed eine sparse Matrix aus einer Datei im bin-Format. */
int SparseMatrix_laden_bin(char* dateiname,   // Datei aus der die Matrix geladen wird
                           SparseMatrix_t *A, // Rueckgabe: Matrix (Call by Reference)
                           long *n);          // Rueckgabe: Anzahl Zeilen (Call by Reference)

/* Die Funktion speichert eine sparse Matrix in einer Datei im ascii-Format. */
int SparseMatrix_speichern_ascii(char* dateiname,      // Datei zum speichern der Matrix
                                 SparseMatrix_t A,     // Matrix
                                 long n,               // Anzahl Zeilen
                                 const char *format);  // Formatierung fuer fprintf

/* Die Funktion speichert eine sparse Matrix in einer Datei im bin-Format. */
int SparseMatrix_speichern_bin(char* dateiname,  // Datei aus der die Matrix geladen wird
                               SparseMatrix_t A, // Matrix
                               long n);          // Anzahl der Zeilen

/********************************************************************************/
/*********** Arbeiten mit symmetrischen Matrizen ********************************/
/********************************************************************************/

/* Diese Funktion legt eine symmetrische Matrix mit n Zeilen und Spalten 
 * dynamisch an. Hierbei wird nur Speicherplatz fuer die untere Dreiecksmatrix
 * angelegt. Die Funktion gibt einen Zeiger auf die Matrix zurueck. 
 */
double *symmat_neu(long n);

/* Es wird die Matrix A, die zuvor mit symmat_neu() angelegt
 * wurde, wieder frei gegeben.
 */
void symmat_freigeben(double *A);

/* Es wird die symmetrische Matrix B nach A kopiert. Hierbei muss in A 
 * genuegend Speicher vorhanden sein.
 */
void symmat_kopieren(double *A,         // Ziel (Rueckgabe)
                     const double *B,   // Quelle
                     long n);           // Anzahl Zeilen und Spalten

/* Die Funktion berechnet die Frobeniusnorm von der symmetrischen Matrix A. */
double symmat_FrobeniusNorm(const double *A,   // Matrix
                            long n);           // Anzahl Zeilen und Spalten

/* Es werden zwei symmetrische Matrizen addiert und das Ergebnis wird in B 
 * zurueck gegeben. Genauer gilt:   B = alpha * A + B
 */
void symmat_addieren(double alpha,      // reelle Zahl mit der A skalliert wird
                     const double *A,   // Matrix A
                     double *B,         // Matrix B (Eingabe und Rueckgabe)
                     long n);           // Anzahl Zeilen und Spalten

/* Matrix-Vektor-Multiplikation fuer symmetrische Matrizen. Es gilt

   y = alpha * A*x + beta * y

   Hierbei ist A eine n x n symmetrisch Matrix, von der nur die untere 
   Dreiecksmatrix abgespeichert ist.
*/
void symmat_vektor_mult(double alpha,     // reelle Zahl mit der A*x skalliert wird
                        double beta,      // reelle Zahl mit der y skalliert wird
                        const double *A,  // Matrix A
                        const double *x,  // Vektor x 
                        double *y,        // Vektor y (Eingabe und Rueckgabe)
                        long n);          // Anzahl Zeilen und Spalten von A

/********************************************************************************/
/*********** Funktionen fuer duennbesetzte Matrizen (sparse) ********************/
/********************************************************************************/

//Meine Funktionen von Blatt 7 (Die Vektormult ist leicht verändert für das skalar an y)

long Sparse_getRowNumber(SparseMatrix_t A, long i);
long Sparse_getColumn(SparseMatrix_t A, long i, long j);

//Transformiert matrix in sparsematrix
SparseMatrix_t transformSparse (double* A, long n);

void Spars_getDiag(SparseMatrix_t A, double *x, long n);
double* SparsetoNormal(SparseMatrix_t sA, long n);

//Helper, den ich später vllt brauche
void Spars_getDiaginverse(SparseMatrix_t A, double *x, long n);

double* transposeMatrix(double* A, long n);

SparseMatrix_t getL (SparseMatrix_t A, long n);
/* Matrix-Vektor-Multiplikation fuer eine quadratische sparse Matrizen. Es gilt

   y = alpha * A*x + beta * y

   Hierbei ist A eine duennbesetzte n x n Matrix im sparse Format.
*/
void SparseMatrix_vektor_mult(double alpha,      // reelle Zahl mit der A*x skalliert wird
                              double beta,       // reelle Zahl mit der y skalliert wird
                              SparseMatrix_t A,  // Matrix A
                              const double *x,   // Vektor x 
                              double *y,         // Vektor y (Eingabe und Rueckgabe)
                              long n);           // Anzahl Zeilen von A

/* Diese Funktion gibt den Speicher innerhalb einer sparse Matrix Struktur frei.
*/
void SparseMatrix_freigeben(SparseMatrix_t A);

#endif
