#ifndef CALCULATE_NORMQ_CU_H
#define CALCULATE_NORMQ_CU_H

#ifdef __cplusplus
extern "C" {
#endif

__global__ void calculate_normq_cu( double *normq, double *dfmlc, double *v, int n );

#ifdef __cplusplus
}
#endif

#endif /* CALCULATE_NORMQ_CU_H */

