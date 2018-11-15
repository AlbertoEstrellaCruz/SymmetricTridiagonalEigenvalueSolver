#include <math.h>

#include "calculate_vi_cu_v1.h"

__device__ void calculate_vi_cu( double *pvi, double *d, double *dfmlc, double *v, int i, int n )
{
	int j;

	double vi;

	vi = dfmlc[ i + n * i ];

	for ( j = 0; j < i ; ++j )
	{
    	vi *= ( dfmlc[ i + n * j ] / ( d[j] - d[i] ) );
	}

	for ( j = i + 1; j < n; ++j )
	{
    	vi *= ( dfmlc[ i + n * j ] / ( d[j] - d[i] ) );
	}

	if ( n % 2 == 1 )
	{
    	vi = -vi;
	}
/*
	if ( vi < 0.0 )
	{
    	printf( "i = %d\n", i );
    	printf( "vi = %.15e\n", vi );
	}
*/
	vi = sqrt(vi);

	if ( v[i] < 0.0 )
	{
    	vi = -vi;
	}

	*pvi = vi;
}

#ifdef __cplusplus
extern "C" 
#endif
__global__ void correct_v_cu( double *v_corr, double *d, double *dfmlc, double *v, int n)
{
	int i;

	i = threadIdx.x + blockDim.x * blockIdx.x;

	while ( i < n ) //for ( i = 0; i < n; ++i )
	{
		calculate_vi_cu( v_corr + i, d, dfmlc, v, i, n );

		i += blockDim.x * gridDim.x;
	}
}

