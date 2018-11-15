#ifndef HIB_TIME_H
#define HIB_TIME_H

#ifdef __cplusplus
extern "C" {
#endif

#include <cuda_runtime.h>

typedef struct _hib_cuda_timer
{
	cudaEvent_t start;
	cudaEvent_t stop;
} hib_cuda_timer_t;

double get_posix_time();

double get_cycles_time();

double get_openmp_time();

double get_mpi_time();

void start_cuda_time();

double finish_cuda_time();

#ifdef __cplusplus
}
#endif

#endif /* HIB_TIME_H */

