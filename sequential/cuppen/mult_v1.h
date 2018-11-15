#ifndef MULT_H
#define MULT_H

#ifdef __cplusplus
extern "C" {
#endif

void mult_blas( double *c , double *a , double *b, int m, int k, int n );

#ifdef __cplusplus
}
#endif

#endif /* MULT_H */

