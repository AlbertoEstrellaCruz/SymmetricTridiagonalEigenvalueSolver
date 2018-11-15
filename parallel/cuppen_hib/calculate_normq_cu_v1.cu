#include <math.h>

#include "calculate_normq_cu_v1.h"

__device__ void calculate_normq_cu_c( double *normq, double *dfmlc, double *v, int c, int n )
{
	int f;

	normq[c] = 0.0;

    for ( f = 0; f < n; ++f )
	{
   	    normq[c] += ( ( v[f] * v[f] ) / ( dfmlc[ f + n * c ] * dfmlc[ f + n * c ] ) );
    }

    normq[c]  = sqrt( normq[c] );
}

__global__ void calculate_normq_cu( double *normq, double *dfmlc, double *v, int n )
{
	int c;

	c = threadIdx.x + blockDim.x * blockIdx.x;

	while ( c < n ) //for ( c = 0; c < n; ++c )
    {
		calculate_normq_cu_c( normq, dfmlc, v, c, n );

		c += blockDim.x * gridDim.x;
	}
}

