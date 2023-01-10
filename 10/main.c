#include <math.h>
#include <stdio.h>
#include "LR.h"
#include "QR.h"
#include "basic_LinAlg.h"



/*
 * The following function implements the direct Iteration
*/
void directIter(double* A, double* y, long n, long iter, double* x_optimal) {
    
    //filepointer für output
    FILE *fp = NULL;
    //öffne oder erstelle datei in fb
    fp = fopen("Teil1A.txt", "w");
    double* x = vektor_neu(n);

    for (int i = 0; i < iter; i++) {
        
        vektor_kopieren(x,y,n);
        vektor_skalieren(1/vektor_2Norm(x, n), x, n);
        matrix_vektor_mult(1, 0, A, x, y, n, n);

        vektor_kopieren(x,y,n);
        vektor_addieren(-1, x_optimal, x, n);
        fprintf(fp, "%d\t%f\n", i+1, vektor_2Norm(x, n)/vektor_2Norm(x_optimal, n));

        printf("%d\t%f\n", i+1, vektor_2Norm(x, n)/vektor_2Norm(x_optimal, n));
    }

    fclose(fp);
    vektor_freigeben(x);
}

/*
 * Implemts the inverse method 
*/
void inverseIter(double* A, double* y, double alpha, long n, long iter, double* x_optimal) {

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
    
    for (int i = 0; i < iter; i++) {
        
        vektor_kopieren(x,y,n);
        vektor_skalieren(1/vektor_2Norm(x, n), x, n);
        //vektor_ausgeben(x,n," % 10.3e");
        LR_VorwSubst(n, At, x);
        LR_RueckwSubst(n, At, x);
        vektor_ausgeben(x,n," % 10.3e");
        vektor_kopieren(y,x,n);


        vektor_addieren(-1, x_optimal, x, n);
        fprintf(fp, "%d\t%f\n", i+1, vektor_2Norm(x, n)/vektor_2Norm(x_optimal, n));

        printf("%d\t%f\n", i+1, vektor_2Norm(x, n)/vektor_2Norm(x_optimal, n));
    }
    vektor_ausgeben(x,n," % 10.3e");
    vektor_freigeben(x);
    matrix_freigeben(At);
    fclose(fp);
}




int main() {
    long n = 4; 
    double* x_start = vektor_neu(n);
    double* x_EV6 = vektor_neu(n);
    x_start[0]=1;
    for (int i = 1; i<n; i++) {
        x_start[i] = 0;
    }
    x_EV6[0]=3;
    x_EV6[1]=-3;
    x_EV6[2]=-3;
    x_EV6[3]=3;


    double* A = matrix_neu(n,n);
    matrix_laden_ascii("A_ascii.dat", &A, &n, &n);

    //directIter(A,x_start,n,20,x_EV6);

    double* x_EV2 = vektor_neu(n);
    x_start[0]=1;
    for (int i = 1; i<n; i++) {
        x_start[i] = 0;
    }
    x_EV2[0]=1;
    x_EV2[1]=1;
    x_EV2[2]=1;
    x_EV2[3]=1;

    inverseIter(A,x_start, 1.90, n, 20, x_EV2);

    for(int i = 0; i<n; i++) {
        A[i*n+i] -= 1.9;
    }
    LR_Zerlegung(n,A);
    LR_VorwSubst(n, A, x_start);
    LR_RueckwSubst(n,A,x_start);
    vektor_ausgeben(x_start, n, " % 10.3e");

    vektor_freigeben(x_start);
    vektor_freigeben(x_EV6);
    matrix_freigeben(A);
}



