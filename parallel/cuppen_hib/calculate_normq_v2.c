#include <math.h>

#include "hib_defs_v1.h"
#include "calculate_normq_v1.h"

void calculate_normq_c( double *normq, double *dfmlc, double *v, int c, int n )
{
	int f;

	normq[c] = 0.0;

    for ( f = 0; f < n; ++f )
	{
   	    normq[c] += ( ( v[f] * v[f] ) / ( dfmlc[ f + n * c ] * dfmlc[ f + n * c ] ) );
    }

    normq[c]  = sqrt( normq[c] );
}

void calculate_normq( double *normq, double *dfmlc, double *v, int n )
{
	int c;

	#pragma omp parallel for num_threads( HIB_OMP_NT ) default( none ) \
	 private( c ) shared( normq, dfmlc, v, n )
	for ( c = 0; c < n; ++c )
    {
		calculate_normq_c( normq, dfmlc, v, c, n );
	}
}

