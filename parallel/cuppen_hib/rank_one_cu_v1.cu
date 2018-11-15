#include <stdlib.h>
#include <stdio.h>

#include "zerodandc_n_cu_v1.h"
#include "calculate_vi_cu_v1.h"
#include "calculate_normq_cu_v1.h"
#include "calculate_q_cu_v1.h"

#include "rank_one_cu_v1.h"

__device__ void calculate_v2_cu_r( double *v2, double *v, int r )
{
	v2[r] = v[r] * v[r];
}

__global__ void calculate_v2_cu( double *v2, double *v, int n )
{
	int r;

	r = threadIdx.x + blockDim.x * blockIdx.x;

	while ( r < n ) //for ( r = 0; r < n; ++r )
	{
		calculate_v2_cu_r( v2, v, r );

		r += blockDim.x * gridDim.x;
	}
}

#ifdef __cplusplus
extern "C" 
#endif
void rank_one_cuda( double *q, double *lambda, double *d, double *v, int n, double eps )
{
	int n_size;
	int nxn_size;

	double *q_d;
	double *lambda_d;
	double *d_d;
	double *v_d;

	double *delta_d;
	double *dfmlc_d;
	double *v2_d;
	double *v_corr_d;
	double *normq_d;

	dim3 n_grid_dim;
	dim3 n_block_dim;
	dim3 nxn_grid_dim;
	dim3 nxn_block_dim;

	n_size = n * sizeof( double );
	nxn_size = n * n_size;

	n_block_dim = 32;
	n_grid_dim = n / 32 + 1;

	nxn_block_dim.x = n_block_dim.x;
	nxn_block_dim.y = n_block_dim.x;
	nxn_grid_dim.x = n_grid_dim.x;
	nxn_grid_dim.y = n_grid_dim.x;

	cudaMalloc( &q_d, nxn_size );

	cudaMalloc( &delta_d, nxn_size );
	cudaMalloc( &dfmlc_d, nxn_size );

	cudaMalloc( &lambda_d, n_size );
	cudaMalloc( &d_d, n_size );
	cudaMalloc( &v_d, n_size );

	cudaMalloc( &v2_d, n_size );
	cudaMalloc( &v_corr_d, n_size );
	cudaMalloc( &normq_d, n_size );

	cudaMemcpy( v_d, v, n_size, cudaMemcpyHostToDevice );
	cudaMemcpy( d_d, d, n_size, cudaMemcpyHostToDevice );

  // Calculando v al cuadrado

	calculate_v2_cu<<< n_grid_dim, n_block_dim >>>( v2_d, v_d, n );

  // Calculando los valores propios de la ecuacion secular

	solve_secular_cu<<< n_grid_dim, n_block_dim >>>( delta_d, lambda_d, dfmlc_d, d_d, v2_d, n, eps );

  // Calculando la correccion de v 

	correct_v_cu<<< n_grid_dim, n_block_dim >>>( v_corr_d, d_d, dfmlc_d, v_d, n );

  // Calculando los vectores propios en Q

	calculate_normq_cu<<< n_grid_dim, n_block_dim >>>( normq_d, dfmlc_d, v_corr_d, n);

	calculate_q_cu<<< nxn_grid_dim, nxn_block_dim >>>( q_d, v_corr_d, dfmlc_d, normq_d, n );

	cudaMemcpy( lambda, lambda_d, n_size, cudaMemcpyDeviceToHost );
	cudaMemcpy( q, q_d, nxn_size, cudaMemcpyDeviceToHost );

	cudaFree( normq_d );
	cudaFree( v_corr_d );
	cudaFree( v2_d );
	cudaFree( v_d );
	cudaFree( d_d );
	cudaFree( lambda_d );
	cudaFree( dfmlc_d );
	cudaFree( delta_d );
	cudaFree( q_d );
}

