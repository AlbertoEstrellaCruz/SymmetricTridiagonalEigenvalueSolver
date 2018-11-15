#include <stdlib.h>
#include <stdio.h>

#include "cublas_v2.h"

#include "mult_cu_v1.h"

#ifdef __cplusplus
extern "C" 
#endif
void mult_cuda( double *c , double *a , double *b, int m, int k, int n )
{
	cublasStatus_t stat;
	cublasHandle_t handle;

	double alpha;
	double beta;

	double *d_c;
	double *d_a;
	double *d_b;

	alpha = 1.0;
	beta = 0.0;

	cudaMalloc( &d_a, m * k * sizeof( double ) );
	cudaMalloc( &d_b, k * n * sizeof( double ) );
	cudaMalloc( &d_c, m * n * sizeof( double ) );

	cudaMemcpy( d_a, a, m * k * sizeof( double ), cudaMemcpyHostToDevice );
	cudaMemcpy( d_b, b, k * n * sizeof( double ), cudaMemcpyHostToDevice );

	stat = cublasCreate( &handle );

	if ( stat != CUBLAS_STATUS_SUCCESS ) {
        fprintf ( stderr, "CUBLAS initialization failed\n" );
		return;
    }

	stat = cublasDgemm( handle, CUBLAS_OP_N, CUBLAS_OP_N, m, n, k, &alpha, d_a, m, d_b, k, &beta, d_c, m );

	if ( stat != CUBLAS_STATUS_SUCCESS ) {
        fprintf ( stderr, "CUBLAS cublasDgemm failed\n" );
		return;
    }

	cudaMemcpy( c, d_c, m * n * sizeof( double ), cudaMemcpyDeviceToHost );

	cudaFree( d_c );
	cudaFree( d_b );
	cudaFree( d_a );

	cublasDestroy(handle);
}

