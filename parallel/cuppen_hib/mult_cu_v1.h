#ifndef MULT_CUDA_H
#define MULT_CUDA_H

#ifdef __cplusplus
extern "C" {
#endif

void mult_cuda( double *c , double *a , double *b, int m, int k, int n );

#ifdef __cplusplus
}
#endif

#endif /* MULT_CUDA_H */

