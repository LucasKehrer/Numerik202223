#ifndef CG_H
#define CG_H

#include "basic_LinAlg.h"

long Sparse_cg(SparseMatrix_t A, double *x, double *b, long n, int prec, long iter,
               double eps);

#endif
