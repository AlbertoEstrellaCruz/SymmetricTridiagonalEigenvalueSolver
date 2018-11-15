#include <stdio.h>
#include <math.h>

#include "calculate_vi_v1.h"

void calculate_vi( double *pvi, double *d, double *dfmlc, double *v, int i, int n )
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

	if ( vi < 0.0 )
	{
    	printf( "i = %d\n", i );
    	printf( "vi = %.15e\n", vi );
	}

	vi = sqrt(vi);

	if ( v[i] < 0.0 )
	{
    	vi = -vi;
	}

	*pvi = vi;

}

