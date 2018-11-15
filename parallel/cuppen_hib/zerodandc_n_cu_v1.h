#ifndef ZERODANC_N_CU_H
#define ZERODANC_N_CU_H

#ifdef __cplusplus
extern "C" {
#endif

__global__ void solve_secular_cu( double *delta, double *lambda, double *dmlambda, double *d, double *v2, int n, double eps );

#ifdef __cplusplus
}
#endif

#endif /* ZERODANC_N_CU_H */

