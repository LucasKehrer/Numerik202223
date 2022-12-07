#ifndef ITER_H
#define ITER_H

void Jacobi_solver(double *A, double *x, double *b, long n, long iter);
void SOR_solver(double *A, double *x, double *b, long n, long iter,
                double omega);
void GaussSeidel_solver(double *A, double *x, double *b, long n, long iter);

#endif
