#ifndef QR_H
#define QR_H

void QR_HS_Zerlegung(long n, double *A, double *u);
void QR_HS_QTransX(long n, double *A, double *u, double *x);
void QR_HS_Transformation(long n, double *A, double *u, double *RQ);

void QR_RueckwSubst(long n, double *A, double *b);

#endif
