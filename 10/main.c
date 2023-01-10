#include <math.h>
#include <stdio.h>
#include "LR.h"
#include "QR.h"
#include "basic_LinAlg.h"



/*
 * The following function implements the direct Iteration
*/
void directIter(double* A, double* y, long n, long iter, double ewopt) {
    
    //filepointer für output
    FILE *fp = NULL;
    //öffne oder erstelle datei in fb
    fp = fopen("Teil1A.txt", "w");
    double* x = vektor_neu(n);
    double* xtemp = vektor_neu(n);

    for (int i = 0; i < iter; i++) {
        
        vektor_kopieren(x,y,n);
        vektor_skalieren(1/vektor_2Norm(x, n), x, n);
        matrix_vektor_mult(1, 0, A, x, y, n, n);

        vektor_kopieren(x,y,n);

        vektor_kopieren(xtemp,x,n);
        matrix_vektor_mult(1,0,A,x, xtemp, n,n);
        //matrix_ausgeben(A,n,n," % 10.3e");
        //vektor_ausgeben(xtemp, n, " % 10.3e");
        double pred = vektor_skalprod(y,xtemp,n)/vektor_skalprod(y,y,n);
        fprintf(fp, "%d\t%f\n", i+1, fabs(ewopt-pred)/ (fabs(ewopt)));

       // printf("%d\t%f\n", i+1, fabs(ewopt-pred)/(fabs(ewopt)));
    }

    vektor_freigeben(xtemp);
    fclose(fp);
    vektor_freigeben(x);
}

/*
 * Implemts the inverse method 
*/
void inverseIter(double* A, double* y, double alpha, long n, long iter, double ewopt) {

    //filepointer für output
    FILE *fp = NULL;
    //öffne oder erstelle datei in fb
    fp = fopen("Teil1B.txt", "w");


    double* At = matrix_neu(n,n);
    matrix_kopieren(At, A, n, n);
    for (int i = 0; i<n; i++) {
        At[i*n+i]-= alpha;
    }
    //matrix_ausgeben(At, n,n," % 10.3e");
    LR_Zerlegung(n, At);

    //matrix_ausgeben(At, n,n," % 10.3e");

    double* x = vektor_neu(n);
    double* xtemp = vektor_neu(n);
    
    for (int i = 0; i < iter; i++) {
        
        vektor_kopieren(x,y,n);
        vektor_skalieren(1/vektor_2Norm(x, n), x, n);
        //vektor_ausgeben(x,n," % 10.3e");
        LR_VorwSubst(n, At, x);
        LR_RueckwSubst(n, At, x);
        vektor_kopieren(y,x,n);

        
        
        vektor_kopieren(xtemp,x,n);
        matrix_vektor_mult(1,0,A,x, xtemp, n,n);
        //matrix_ausgeben(A,n,n," % 10.3e");
        //vektor_ausgeben(xtemp, n, " % 10.3e");
        double pred = vektor_skalprod(y,xtemp,n)/vektor_skalprod(y,y,n);
        fprintf(fp, "%d\t%f\n", i+1, fabs(ewopt-pred)/ (fabs(ewopt)));

       // printf("%d\t%f\n", i+1, fabs(ewopt-pred)/(fabs(ewopt)));
        
    }
    vektor_freigeben(xtemp);
    
    vektor_freigeben(x);
    matrix_freigeben(At);
    fclose(fp);
}


void QRAlgorithmus(double* A, long n, long iter, double* EWoptimal) {
    //filepointer für output
    FILE *fp1 = NULL;
    FILE *fp2 = NULL;
    FILE *fp3 = NULL;
    FILE *fp4 = NULL;
    //öffne oder erstelle datei in fb
    fp1 = fopen("Teil2EW1.txt", "w");
    fp2 = fopen("Teil2EW2.txt", "w");
    fp3 = fopen("Teil2EW3.txt", "w");
    fp4 = fopen("Teil2EW4.txt", "w");

    double* u = vektor_neu(n);
    double* At = matrix_neu(n,n);
    for(int i = 0; i<iter; i++) {

    QR_HS_Zerlegung(n,A,u);
    QR_HS_Transformation(n,A,u,At);

    matrix_kopieren(A,At,n,n);
    fprintf(fp1, "%d\t%f\n", i+1, fabs(EWoptimal[0]-A[0])/ (EWoptimal[0]));
    fprintf(fp2, "%d\t%f\n", i+1, fabs(EWoptimal[1]-A[n+1])/ (EWoptimal[1]));
    fprintf(fp3, "%d\t%f\n", i+1, fabs(EWoptimal[2]-A[2*n+2])/ (EWoptimal[2]));
    fprintf(fp4, "%d\t%f\n", i+1, fabs(EWoptimal[3]-A[3*n+3])/ (EWoptimal[3]));
    }
    
    matrix_freigeben(At);

    fclose(fp1);
    fclose(fp2);
    fclose(fp3);
    fclose(fp4);
}



int main() {
    long n = 4; 
    double* x_start = vektor_neu(n);
    x_start[0]=1;
    for (int i = 1; i<n; i++) {
        x_start[i] = 0;
    }

    double* A = matrix_neu(n,n);
    matrix_laden_ascii("A_ascii.dat", &A, &n, &n);

    directIter(A,x_start,n,20,6);

    //Ich will ein schlechts beispiel, aber alles konvergiert nach 2 iterationschritten
    x_start[0]=1;
    x_start[0]=0;
    x_start[0]=0;
    x_start[0]=0;
    

    inverseIter(A,x_start, 1.90, n, 20, 2);

    x_start[0]=6;
    x_start[1]=4;
    x_start[2]=4;
    x_start[3]=2;

    QRAlgorithmus(A,n,20,x_start);

    vektor_freigeben(x_start);
    matrix_freigeben(A);
}



