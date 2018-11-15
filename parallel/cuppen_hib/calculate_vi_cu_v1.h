#ifndef CALCULATE_VI_CU_H
#define CALCULATE_VI_CU_H

#ifdef __cplusplus
extern "C" {
#endif

__global__ void correct_v_cu( double *v_corr, double *d, double *dfmlc, double *v, int n);

#ifdef __cplusplus
}
#endif

#endif /* CALCULATE_VI_CU_H */

