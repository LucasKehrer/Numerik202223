#ifndef CG_H
#define CG_H

#include "basic_LinAlg.h"

long Sparse_cg(SparseMatrix_t A, double *x, double *b, long n, int prec, long iter,
               double eps);
long Sparse_cg_prec(SparseMatrix_t A, double *x, double *b, long n, int prec, long iter, double eps);
//tests
void solve_lowerSparse(SparseMatrix_t A, double *y, double *x, long n, double omega);
void solve_upperSparse(SparseMatrix_t A, double *y, double *x, long n, double omega);
#endif
