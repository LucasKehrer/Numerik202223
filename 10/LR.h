#ifndef LR_H
#define LR_H

void LR_Zerlegung(long n, double *A);
void LR_VorwSubst(long n, double *A, double *b);
void LR_RueckwSubst(long n, double *A, double *b);

#endif
