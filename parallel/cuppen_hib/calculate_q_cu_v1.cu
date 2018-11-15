#include "calculate_q_cu_v1.h"

__device__ void calculate_q_cu_fc( double *q, double *v, double *dfmlc,  double *normq, int f, int c, int n )
{
	q[ f + n * c ] = ( v[f] / dfmlc[ f + n * c ] ) / normq[c];
}

__global__ void calculate_q_cu( double *q, double *v, double *dfmlc, double *normq, int n )
{
	int c;
	int f;

	c = threadIdx.y + blockDim.y * blockIdx.y;

	while ( c < n ) //for ( c = 0; c < n; ++c )
    {
		f = threadIdx.x + blockDim.x * blockIdx.x;

    	while ( f < n ) //for ( f = 0; f < n; ++f)
     	{
    	    calculate_q_cu_fc( q, v, dfmlc,  normq, f, c, n );

			f += blockDim.x * gridDim.x;
    	}

		c += blockDim.y * gridDim.y;
	}
}

