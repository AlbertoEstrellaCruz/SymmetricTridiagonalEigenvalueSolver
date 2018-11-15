#include <stdio.h>

#include "hib_defs_v1.h"
#include "calculate_q_v1.h"

void calculate_q_fc( double *q, double *v, double *dfmlc,  double *normq, int f, int c, int n )
{
	q[ f + n * c ] = ( v[f] / dfmlc[ f + n * c ] ) / normq[c];
}

void calculate_q( double *q, double *v, double *dfmlc, double *normq, int n )
{
	int c;
	int f;

	#pragma omp parallel for num_threads( HIB_OMP_NT ) default( none ) \
	 private( c,f ) shared( q, v, dfmlc, normq, n )
	for ( c = 0; c < n; ++c )
    {
    	for ( f = 0; f < n; ++f)
     	{
    	    calculate_q_fc( q, v, dfmlc,  normq, f, c, n );
    	}
	}
}

