#ifndef CALCULATE_Q_CU_H
#define CALCULATE_Q_CU_H

#ifdef __cplusplus
extern "C" {
#endif

__global__ void calculate_q_cu( double *q, double *v, double *dfmlc, double *normq, int n );

#ifdef __cplusplus
}
#endif

#endif /* CALCULATE_Q_CU_H */

