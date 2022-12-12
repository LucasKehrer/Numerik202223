#ifndef ITER_H
#define ITER_H

#include "basic_LinAlg.h"

void Sparse_JacobiCheb_solver(SparseMatrix_t A, double *x, double *b, long n, long iter, double gamma_1, double gamma_2);
void Sparse_Jacobi_solver(SparseMatrix_t A, double *x, double *b, long n, long iter);
void Sparse_Matrix_mul(SparseMatrix_t A, double *x, double *y, long n);

#endif
