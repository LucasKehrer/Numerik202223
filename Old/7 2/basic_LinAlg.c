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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

// eigene Header-Datei einbinden!!
#include "basic_LinAlg.h"

/********************************************************************************/
/*********** Funktionen fuer Matrizen *******************************************/
/********************************************************************************/

double *matrix_neu(long m,    // Anzahl der Zeilen
                   long n) {  // Anzahl der Spalten

  // WICHTIG: sizeof(double)
  double *A = (double *) calloc( m*n, sizeof(double) );

  if (A == NULL) {  // Fehlerbehandlung
    fprintf(stderr, "Fehler beim anlegen einer Matrix.\n");
    exit(1); // Programm mit Rueckgabewert 1 beenden!!
  }

  /* Variable A zurueckgeben (Kopie des lokalen Zeigers)
     Dies ist unproblematisch, da der Speicher, auf den A zeigt, 
     erst dann freigegeben wird, wenn wir dies explizit veranlassen. */
  return A;
}

void matrix_freigeben(double *A) { // Matrix zum freigeben
  free(A);
  // A = NULL; // haette keinen Effekt, da Call by Value
}

void matrix_ausgeben(const double *A,      // Matrix
                     long m,               // Anzahl Zeilen
                     long n,               // Anzahl Spalten
                     const char *format) { // Formatierung fuer printf
  long i;
  static char form[32] = " % 7.3f";

  if (format != NULL)
    strcpy(form, format);

  for (i=0; i<m*n; i++)  
    if ( (i+1)%n )  // Matrixeintraege ausgeben
      printf(form, A[i]);
    else {          // letzten Eintrag in der Zeile ausgeben
      printf(form, A[i]);
      printf("\n");
    }
}

void matrix_kopieren(double *A,         // Ziel (Rueckgabe)
                     const double *B,   // Quelle
                     long m,            // Anzahl Zeilen
                     long n) {          // Anzahl Spalten
  long i;

  for (i=0; i<m*n; i++)    // Schleife ueber alle Eintraege
    A[i] = B[i];           // Matrixeintraege kopieren
}

double matrix_FrobeniusNorm(const double *A,   // Matrix
                            long m,            // Anzahl Zeilen
                            long n) {          // Anzahl Spalten
  long i;
  double norm = 0;  // WICHTIG: Initialisieren wegen Summation

  for (i=0; i<m*n; i++)       // Schleife ueber alle Elemente
    // Matrixeintraege quadrieren und aufsummieren
    norm = norm + A[i] * A[i];

  return sqrt(norm);
}

// B = alpha * A + B
void matrix_addieren(double alpha,      // reelle Zahl mit der A skalliert wird
                     const double *A,   // Matrix A
                     double *B,         // Matrix B (Eingabe und Rueckgabe)
                     long m,            // Anzahl Zeilen
                     long n) {          // Anzahl Spalten
  long i;

  for (i=0; i<m*n; i++)   // Schleife ueber alle Matrixeintraege
      // einzelne Matrixeintraege skallieren und addieren
      B[i] = B[i] + alpha * A[i]; 
}

/*
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
                 long l) {         // Anzahl Spalten
  long i, j, k;

  /*
    Hierin wird das Skalarprodukt von einer Zeile aus A und 
    einer Spalte aus B zwischengespeichert.
  */
  double summe; 

  for (i=0; i<m; i++) {   // Schleife ueber alle Zeilen
    for (j=0; j<n; j++) { // Schleife ueber alle Spalten
      /*
        Skalarprodukt von der i-ten Zeile von A und
        der j-ten Spalte von B berechnen.
      */
      summe = 0;  // WICHTIG: Initialisieren wegen Summation
      for (k=0; k<l; k++)
        summe = summe + A[i*l+k] * B[k*n+j];

      /*
        einzelne Matrixeintraege von C berechnen
      */
      C[i*n+j] = alpha * summe + beta * C[i*n+j]; 
    }
  }
}

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
                        long n) {         // Anzahl Spalten von A
  matrix_mult(alpha,beta,A,x,y,m,1,n);
}

void matrix_einlesen(double **A,   // Rueckgabe: Matrix (Call by Reference)
                     long *m,      // Rueckgabe: Anzahl Zeilen (Call by Reference)
                     long *n) {    // Rueckgabe: Anzahl Spalten (Call by Reference)
  long i, j;

  /* 
    Bei scanf wird der Adressoperator & fuer die Argumente m und n
    nicht benoetigt. Es handelt sich hierbei bereits um Zeiger auf
    die Variablen wegen Call by Reference.
  */
  printf("Anzahl Zeilen und Spalten (\'m n\'): ");
  if ( scanf("%ld %ld", m, n) != 2 || *m<1 || *n<1 ) {
    fprintf(stderr, "Fehler beim Einlesen der Anzahl der Zeilen und Spalten"
                    "fur die anzulegende Matrix. Es muss gelten: m>0 und n>0\n");
    exit(1); 
  }

  /* 
    Soll die Matrix A bzw die Laengen m, n verwendet werden, muessen sie
    dereferenziert werden. Wegen Call by Reference handelt es sich hierbei
    jeweils um einen Zeiger auf die entsprechende Variable.
  */
  *A = matrix_neu(*m, *n);

  for (i=0; i<*m; i++) {     // Schleife ueber alle Zeilen
    for (j=0; j<*n; j++) {   // Schleife ueber alle Spalten
      // Matrixeintraege einlesen
      printf("%3ld. Zeile, %3ld. Spalte: ", i, j);
      scanf("%lf", &(*A)[i*(*n)+j]);
    }
    printf("\n");
  }
}

int matrix_laden_ascii(char* dateiname, // Datei aus der die Matrix geladen wird
                       double **A,  // Rueckgabe: Matrix (Call by Reference)
                       long *m,     // Rueckgabe: Anzahl Zeilen (Call by Reference)
                       long *n) {   // Rueckgabe: Anzahl Spalten (Call by Reference)
  long i, j;
  char c;
  FILE *datei;  // Dateizeiger deklarieren

  datei = fopen(dateiname, "r");  // Datei zum lesen oeffnen
  
  // Fehlerbehandlung, falls Oeffnen fehlgeschlagen
  if (datei == NULL) 
    return -1; 

  while ( (c=fgetc(datei)) == '#' ) // Kommentare ueberspringen
    while ( fgetc(datei) != '\n' );

  ungetc(c, datei); // ein Zeichen wurde zu viel gelesen. Zurueck damit!

  /* 
    Anzahl Zeilen und Spalten lesen

    Bei fscanf wird der Adressoperator & fuer die Argumente m und n
    nicht benoetigt. Es handelt sich hierbei bereits um Zeiger auf
    die Variablen wegen Call by Reference.
  */
  if (fscanf(datei, "%ld %ld\n", m, n) != 2) {
    fclose(datei);  // Datei schliessen
    return -2; 
  }

  /* 
    Soll die Matrix A bzw die Laengen m, n verwendet werden, muessen sie
    dereferenziert werden. Wegen Call by Reference handelt es sich hierbei
    jeweils um einen Zeiger auf die entsprechende Variable.
  */
  *A = matrix_neu(*m, *n);

  // Matrixeintraege einlesen
  for (i=0; i<*m; i++) {     // Schleife ueber alle Zeilen
    for (j=0; j<*n; j++) {   // Schleife ueber alle Spalten
      if (fscanf(datei, "%lf", &(*A)[i*(*n)+j]) != 1) {
        fclose(datei);              // Datei schliessen
        matrix_freigeben(*A);      // angelegte Matrix freigeben
        *A = NULL;                  // zur Sicherheit, da A nicht gelesen
        return -2;
      }
    }
  }

  fclose(datei);  // Datei schliessen

  return 0;
}

int matrix_laden_bin(char* dateiname, // Datei aus der die Matrix geladen wird
                     double **A,   // Rueckgabe: Matrix (Call by Reference)
                     long *m,      // Rueckgabe: Anzahl Zeilen (Call by Reference)
                     long *n) {    // Rueckgabe: Anzahl Spalten (Call by Reference)
  char c;
  FILE *datei;  // Dateizeiger deklarieren

  datei = fopen(dateiname, "r");  // Datei zum lesen oeffnen
  
  // Fehlerbehandlung, falls Oeffnen fehlgeschlagen
  if (datei == NULL) 
    return -1; 

  while ( (c=fgetc(datei)) == '#' ) // Kommentare ueberspringen
    while ( fgetc(datei) != '\n' );

  ungetc(c, datei); // ein Zeichen wurde zu viel gelesen. Zurueck damit!

  /* 
    Anzahl Zeilen und Spalten lesen

    Bei fscanf wird der Adressoperator & fuer die Argumente m und n
    nicht benoetigt. Es handelt sich hierbei bereits um Zeiger auf
    die Variablen wegen Call by Reference.
  */
  if (fscanf(datei, "%ld %ld\n", m, n) != 2) {
    fclose(datei);  // Datei schliessen
    return -2; 
  }

  /* 
    Soll die Matrix A bzw die Laengen m, n verwendet werden, muessen sie
    dereferenziert werden. Wegen Call by Reference handelt es sich hierbei
    jeweils um einen Zeiger auf die entsprechende Variable.
  */
  *A = matrix_neu(*m, *n);

  // Matrixeintraege einlesen
  if (fread(*A, sizeof(double), (*m)*(*n), datei) != (*m)*(*n)) {
    fclose(datei);              // Datei schliessen
    matrix_freigeben(*A);      // angelegte Matrix freigeben
    *A = NULL;                  // zur Sicherheit, da A nicht gelesen
    return -2;
  }

  fclose(datei);  // Datei schliessen

  return 0;
}

int matrix_speichern_ascii(char* dateiname,      // Datei zum speichern der Matrix
                           double *A,            // Matrix
                           long m,               // Anzahl Zeilen
                           long n,               // Anzahl Spalten
                           const char *format) { // Formatierung fuer fprintf
  long i, j;
  static char form[32] = " % .15f";
  FILE *datei;  // Dateizeiger deklarieren

  if (format != NULL)
    strcpy(form, format);

  datei = fopen(dateiname, "w");  // Datei zum schreiben oeffnen
  
  // Fehlerbehandlung, falls Oeffnen fehlgeschlagen
  if (datei == NULL) 
    return -1; 

  fprintf(datei, "# Erstellt durch matrix_speichern_ascii()\n");
  fprintf(datei, "# Dateiname: %s\n", dateiname);

  // Anzahl Zeilen und Spalten schreiben
  fprintf(datei, "%ld %ld\n", m, n);

  // Matrixeintraege schreiben
  for (i=0; i<m; i++) {     // Schleife ueber alle Zeilen
    for (j=0; j<n; j++) {   // Schleife ueber alle Spalten
      fprintf(datei, form, A[i*n+j]);
    }
    fprintf(datei, "\n");
  }

  fclose(datei);  // Datei schliessen

  return 0;
}

int matrix_speichern_bin(char* dateiname, // Datei aus der die Matrix geladen wird
                         double *A,       // Matrix
                         long m,          // Anzahl Zeilen
                         long n) {        // Anzahl Spalten
  FILE *datei;  // Dateizeiger deklarieren

  datei = fopen(dateiname, "w");  // Datei zum lesen oeffnen
  
  // Fehlerbehandlung, falls Oeffnen fehlgeschlagen
  if (datei == NULL) 
    return -1; 

  fprintf(datei, "# Erstellt durch matrix_speichern_bin()\n");
  fprintf(datei, "# Dateiname: %s\n", dateiname);

  // Anzahl Zeilen und Spalten schreiben
  fprintf(datei, "%ld %ld\n", m, n);

  // Matrix schreiben
  if (fwrite(A, sizeof(double), m*n, datei) != m*n) {
    fclose(datei);  // Datei schliessen
    return -2;
  }

  fclose(datei);  // Datei schliessen

  return 0;
}

int SparseMatrix_laden_ascii(char* dateiname,   // Datei aus der die Matrix geladen wird
                             SparseMatrix_t *A, // Rueckgabe: Matrix (Call by Reference)
                             long *n) {         // Rueckgabe: Anzahl der Zeilen (Call by Reference)
  long i;
  long anzahlEintraege;
  char c;
  FILE *datei;  // Dateizeiger deklarieren

  datei = fopen(dateiname, "r");  // Datei zum lesen oeffnen
  
  // Fehlerbehandlung, falls Oeffnen fehlgeschlagen
  if (datei == NULL) 
    return -1; 

  while ( (c=fgetc(datei)) == '#' ) // Kommentare ueberspringen
    while ( fgetc(datei) != '\n' );

  ungetc(c, datei); // ein Zeichen wurde zu viel gelesen. Zurueck damit!

  /* 
    Anzahl Zeilen und Spalten lesen sowie die Anzahl der nicht null Eintraege

    Bei fscanf wird der Adressoperator & fuer das Argument n
    nicht benoetigt. Es handelt sich hierbei bereits um Zeiger auf
    die Variablen wegen Call by Reference.
  */
  if (fscanf(datei, "%ld %ld\n", n, &anzahlEintraege) != 2) {
    fclose(datei);  // Datei schliessen
    return -2; 
  }

  /* 
    Soll die Matrix A bzw die Laengen n verwendet werden, muessen sie
    dereferenziert werden. Wegen Call by Reference handelt es sich hierbei
    jeweils um einen Zeiger auf die entsprechende Variable.
  */
  A->values = matrix_neu(anzahlEintraege, 1);

  // Matrixeintraege einlesen
  for (i=0; i<anzahlEintraege; i++) {     // Schleife ueber alle Eintraege
    if (fscanf(datei, "%lf", &(A->values[i])) != 1) {
      fclose(datei);               // Datei schliessen
      matrix_freigeben(A->values); // angelegte Matrix freigeben
      A->values = NULL;            // zur Sicherheit, da A nicht gelesen
      return -2;
    }
  }

  A->columns = (long *) calloc( anzahlEintraege, sizeof(long) );

  if (A->columns == NULL) {  // Fehlerbehandlung
    fprintf(stderr, "Fehler beim anlegen einer Matrix.\n");
    exit(1); // Programm mit Rueckgabewert 1 beenden!!
  }

  // Spaltenindizes einlesen
  for (i=0; i<anzahlEintraege; i++) {     // Schleife ueber alle Eintraege
    if (fscanf(datei, "%ld", &(A->columns[i])) != 1) {
      fclose(datei);                // Datei schliessen
      matrix_freigeben(A->values);  // angelegte Matrix freigeben
      free(A->columns);             // angelegte Matrix freigeben
      A->values = NULL;             // zur Sicherheit, da A nicht gelesen
      A->columns = NULL;            // zur Sicherheit, da A nicht gelesen
      return -2;
    }
  }

  A->rowIndex = (long *) calloc( *n+1, sizeof(long) );

  if (A->rowIndex == NULL) {  // Fehlerbehandlung
    fprintf(stderr, "Fehler beim anlegen einer Matrix.\n");
    exit(1); // Programm mit Rueckgabewert 1 beenden!!
  }

  // rowIndex einlesen
  for (i=0; i<*n+1; i++) {     // Schleife ueber alle Zeilen+1
    if (fscanf(datei, "%ld", &(A->rowIndex[i])) != 1) {
      fclose(datei);                 // Datei schliessen
      matrix_freigeben(A->values);   // angelegte Matrix freigeben
      free(A->columns);              // angelegte Matrix freigeben
      free(A->rowIndex);             // angelegte Matrix freigeben
      A->values = NULL;              // zur Sicherheit, da A nicht gelesen
      A->columns = NULL;             // zur Sicherheit, da A nicht gelesen
      A->rowIndex = NULL;            // zur Sicherheit, da A nicht gelesen
      return -2;
    }
  }

  fclose(datei);  // Datei schliessen

  return 0;
}

int SparseMatrix_laden_bin(char* dateiname,   // Datei aus der die Matrix geladen wird
                           SparseMatrix_t *A, // Rueckgabe: Matrix (Call by Reference)
                           long *n) {         // Rueckgabe: Anzahl Zeilen (Call by Reference)
  char c;
  long anzahlEintraege;
  FILE *datei;  // Dateizeiger deklarieren

  datei = fopen(dateiname, "r");  // Datei zum lesen oeffnen
  
  // Fehlerbehandlung, falls Oeffnen fehlgeschlagen
  if (datei == NULL) 
    return -1; 

  while ( (c=fgetc(datei)) == '#' ) // Kommentare ueberspringen
    while ( fgetc(datei) != '\n' );

  ungetc(c, datei); // ein Zeichen wurde zu viel gelesen. Zurueck damit!

  /* 
    Anzahl Zeilen und nicht null Eintraege lesen

    Bei fscanf wird der Adressoperator & fuer die Argumente m und n
    nicht benoetigt. Es handelt sich hierbei bereits um Zeiger auf
    die Variablen wegen Call by Reference.
  */
  if (fscanf(datei, "%ld %ld\n", n, &anzahlEintraege) != 2) {
    fclose(datei);  // Datei schliessen
    return -2; 
  }

  /* 
    Soll die Matrix A bzw die Laengen m, n verwendet werden, muessen sie
    dereferenziert werden. Wegen Call by Reference handelt es sich hierbei
    jeweils um einen Zeiger auf die entsprechende Variable.
  */
  A->values = matrix_neu(anzahlEintraege, 1);

  // Matrixeintraege einlesen
  if (fread(A->values, sizeof(double), anzahlEintraege, datei) != anzahlEintraege) {
    fclose(datei);                 // Datei schliessen
    matrix_freigeben(A->values);   // angelegte Matrix freigeben
    A->values = NULL;              // zur Sicherheit, da A nicht gelesen
    return -2;
  }

  A->columns = (long *) calloc( anzahlEintraege, sizeof(long) );

  if (A->columns == NULL) {  // Fehlerbehandlung
    fprintf(stderr, "Fehler beim anlegen einer Matrix.\n");
    exit(1); // Programm mit Rueckgabewert 1 beenden!!
  }

  // Spaltenindizes einlesen
  if (fread(A->columns, sizeof(long), anzahlEintraege, datei) != anzahlEintraege) {
    fclose(datei);                // Datei schliessen
    matrix_freigeben(A->values);  // angelegte Matrix freigeben
    free(A->columns);             // angelegte Matrix freigeben
    A->values = NULL;             // zur Sicherheit, da A nicht gelesen
    A->columns = NULL;            // zur Sicherheit, da A nicht gelesen
    return -2;
  }

  A->rowIndex = (long *) calloc( *n+1, sizeof(long) );

  if (A->rowIndex == NULL) {  // Fehlerbehandlung
    fprintf(stderr, "Fehler beim anlegen einer Matrix.\n");
    exit(1); // Programm mit Rueckgabewert 1 beenden!!
  }

  // rowIndex einlesen
  if (fread(A->rowIndex, sizeof(long), *n+1, datei) != *n+1) {
    fclose(datei);                 // Datei schliessen
    matrix_freigeben(A->values);   // angelegte Matrix freigeben
    free(A->columns);              // angelegte Matrix freigeben
    free(A->rowIndex);             // angelegte Matrix freigeben
    A->values = NULL;              // zur Sicherheit, da A nicht gelesen
    A->columns = NULL;             // zur Sicherheit, da A nicht gelesen
    A->rowIndex = NULL;            // zur Sicherheit, da A nicht gelesen
    return -2;
  }

  fclose(datei);  // Datei schliessen

  return 0;
}

int SparseMatrix_speichern_ascii(char* dateiname,      // Datei zum speichern der Matrix
                                 SparseMatrix_t A,     // Matrix
                                 long n,               // Anzahl Zeilen
                                 const char *format) { // Formatierung fuer fprintf
  long i;
  static char form[32] = " % 7.3f";
  FILE *datei;  // Dateizeiger deklarieren

  if (format != NULL)
    strcpy(form, format);

  datei = fopen(dateiname, "w");  // Datei zum schreiben oeffnen
  
  // Fehlerbehandlung, falls Oeffnen fehlgeschlagen
  if (datei == NULL) 
    return -1; 

  fprintf(datei, "# Erstellt durch SparseMatrix_speichern_ascii()\n");
  fprintf(datei, "# Dateiname: %s\n", dateiname);

  // Anzahl Zeilen und Spalten schreiben
  fprintf(datei, "%ld %ld\n", n, A.rowIndex[n]);

  // Matrixeintraege schreiben
  for (i=0; i<A.rowIndex[n]; i++) {
    fprintf(datei, form, A.values[i]);
  }
  fprintf(datei, "\n");

  // Spaltenindizes schreiben
  for (i=0; i<A.rowIndex[n]; i++) {
    fprintf(datei, " %6ld", A.columns[i]);
  }
  fprintf(datei, "\n");

  // rowIndex schreiben
  for (i=0; i<n+1; i++) {
    fprintf(datei, " %6ld", A.rowIndex[i]);
  }
  fprintf(datei, "\n");

  fclose(datei);  // Datei schliessen

  return 0;
}

int SparseMatrix_speichern_bin(char* dateiname,  // Datei aus der die Matrix geladen wird
                               SparseMatrix_t A, // Matrix
                               long n) {         // Anzahl der Zeilen
  FILE *datei;  // Dateizeiger deklarieren

  datei = fopen(dateiname, "w");  // Datei zum lesen oeffnen
  
  // Fehlerbehandlung, falls Oeffnen fehlgeschlagen
  if (datei == NULL) 
    return -1; 

  fprintf(datei, "# Erstellt durch SparseMatrix_speichern_bin()\n");
  fprintf(datei, "# Dateiname: %s\n", dateiname);

  // Anzahl Zeilen und nicht null Eintraege schreiben
  fprintf(datei, "%ld %ld\n", n, A.rowIndex[n]);

  // Matrixeintraege schreiben
  if (fwrite(A.values, sizeof(double), A.rowIndex[n], datei) != A.rowIndex[n]) {
    fclose(datei);  // Datei schliessen
    return -2;
  }

  // Spaltenindizes schreiben
  if (fwrite(A.columns, sizeof(long), A.rowIndex[n], datei) != A.rowIndex[n]) {
    fclose(datei);  // Datei schliessen
    return -2;
  }

  // rowIndex schreiben
  if (fwrite(A.rowIndex, sizeof(long), n+1, datei) != n+1) {
    fclose(datei);  // Datei schliessen
    return -2;
  }

  fclose(datei);  // Datei schliessen

  return 0;
}
/********************************************************************************/
/*********** Funktionen fuer Vektoren *******************************************/
/********************************************************************************/

void vektor_ausgeben(const double *x,      // Matrix
                     long n,               // Laenge
                     const char *format) { // Formatierung fuer printf
  matrix_ausgeben(x, n, 1, format);
}

/* Diese Funktion legt einen Vektor der Laenge n dynamisch an. 
 * Sie gibt einen Zeiger auf den Vektor zurueck. 
 */
double *vektor_neu(long n) {   // Laenge des Vektors
  return matrix_neu(n,1);
}

/* Es wird der Vektor x, der zuvor mit vektor_neu() angelegt
 * wurde, wieder frei gegeben.
 */
void vektor_freigeben(double *x) {  // Vektor zum freigeben
  matrix_freigeben(x);
}

/* Es wird der Vektor y nach x kopiert. Hierbei muss in x 
 * genuegend Speicher vorhanden sein.
 */
void vektor_kopieren(double *x,         // Ziel (Rueckgabe)
                     const double *y,   // Quelle
                     long n) {          // Laenge
  matrix_kopieren(x,y,n,1);
}

/* Der Vektor x wird skalliert:  x = alpha * x
 */
void vektor_skalieren(double alpha,  // Skallierungsparameter
                      double *x,     // Vektor x (Ein- und Ausgabe)
                      long n) {      // Laenge
  long i;
  for (i=0; i<n; i++)
    x[i] *= alpha;
}

/* Die Funktion berechnet die 2-Norm des Vektors x. */
double vektor_2Norm(const double *x,   // Matrix
                    long n) {          // Laenge
  return matrix_FrobeniusNorm(x,n,1);
}

/* Es wird das Skalarprodukt von den Vektoren x und y 
 * berechnet und als Rueckgabewert zurueck gegeben.
 */
double vektor_skalprod(const double *x,  // Vektor x
                       const double *y,  // Vektor y
                       long n) {         // Laenge
  long i;
  double erg = 0;

  for (i=0; i<n; i++)
    erg += x[i] * y[i];

  return erg;
}

/* Es werden zwei Vektoren addiert und das Ergebnis wird in y 
 * zurueck gegeben. Genauer gilt:   y = alpha * x + y
 */
void vektor_addieren(double alpha,      // reelle Zahl mit der x skalliert wird
                     const double *x,   // Vektor x
                     double *y,         // Vektor y (Eingabe und Rueckgabe)
                     long n) {          // Laenge
  matrix_addieren(alpha,x,y,n,1);
}

/********************************************************************************/
/*********** Funktionen fuer symmetrische Matrizen ******************************/
/********************************************************************************/

double *symmat_neu(long n) {  // Anzahl der Zeilen und Spalten
  return vektor_neu(n*(n+1)/2);
}

void symmat_freigeben(double *A) { // Matrix zum freigeben
  free(A);
  // A = NULL; // haette keinen Effekt, da Call by Value
}

void symmat_ausgeben(const double *A,      // Matrix
                     long n,               // Anzahl Zeilen und Spalten
                     const char *format) { // Formatierung fuer printf
  long i,j;
  static char form[32] = " % 7.3f";

  if (format != NULL)
    strcpy(form, format);

  for (i=0; i<n; i++) {
    for (j=0; j<n; j++)  
      if ( i>=j )     // unterhalb oder auf der Diagonalen
        printf(form, A[i*(i+1)/2+j]);
      else            // ueberhalb der Diagonalen
        printf(form, A[j*(j+1)/2+i]);

    printf("\n");
  }
}

void symmat_ausgeben_dreieck(const double *A,      // Matrix
                             long n,               // Anzahl Zeilen und Spalten
                             const char *format) { // Formatierung fuer printf
  long i,j;
  static char form[32] = " % 7.3f";

  if (format != NULL)
    strcpy(form, format);

  for (i=0; i<n; i++) {
    for (j=0; j<=i; j++)  
      printf(form, A[i*(i+1)/2+j]); // unterhalb oder auf der Diagonalen

    printf("\n");
  }
}

void symmat_kopieren(double *A,         // Ziel (Rueckgabe)
                     const double *B,   // Quelle
                     long n) {          // Anzahl Zeilen und Spalten
  vektor_kopieren(A,B,n*(n+1)/2);
}

double symmat_FrobeniusNorm(const double *A,   // Matrix
                            long n) {          // Anzahl Zeilen und Spalten
  long i,j;
  double norm = 0;  // WICHTIG: Initialisieren wegen Summation

  for (i=0; i<n; i++) {     // Schleife ueber alle Elemente
    // Matrixeintraege quadrieren und aufsummieren
    norm = norm + A[i*(i+1)/2+i] * A[i*(i+1)/2+i];   // Diagonalelemente
    for (j=0; j<i; j++)
      norm = norm + 2 * A[i*(i+1)/2+j] * A[i*(i+1)/2+j]; // nicht Diagonalelemente
  }

  return sqrt(norm);
}

// B = alpha * A + B
void symmat_addieren(double alpha,      // reelle Zahl mit der A skalliert wird
                     const double *A,   // Matrix A
                     double *B,         // Matrix B (Eingabe und Rueckgabe)
                     long n) {          // Anzahl Zeilen und Spalten
  vektor_addieren(alpha,A,B,n*(n+1)/2);
}

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
                        long n) {         // Anzahl Zeilen und Spalten von A
  long i,j;

  for (i=0; i<n; i++) {
    y[i] *= beta;

    for (j=0; j<n; j++)
      if (i>=j)
        y[i] += alpha * A[i*(i+1)/2+j]*x[j];
      else
        y[i] += alpha * A[j*(j+1)/2+i]*x[j];
  }
}

/********************************************************************************/
/*********** Funktionen fuer duennbesetzte Matrizen (sparse) ********************/
/********************************************************************************/

/* Diese Funktion gibt den Speicher innerhalb einer sparse Matrix Struktur frei.
*/
void SparseMatrix_freigeben(SparseMatrix_t A) {
  free(A.values);
  free(A.columns);
  free(A.rowIndex);
  //A->values = NULL;   // haette keinen Effekt, da Call by Value und
  //A->columns = NULL;  // somit eine Kopie der Struktur uebergeben wird
  //A->rowIndex = NULL; 
}
