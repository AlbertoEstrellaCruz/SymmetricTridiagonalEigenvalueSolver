#ifndef RANK_ONE_CU_H
#define RANK_ONE_CU_H

#ifdef __cplusplus
extern "C" {
#endif

void rank_one_cuda( double *q, double *lambda, double *d, double *v, int n, double eps );

#ifdef __cplusplus
}
#endif

#endif /* RANK_ONE_CU_H */

