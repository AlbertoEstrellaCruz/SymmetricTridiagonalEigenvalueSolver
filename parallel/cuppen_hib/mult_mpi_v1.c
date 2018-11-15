#include <stdlib.h>
#include <stdio.h>

#include "hib_defs_v1.h"
#include "commfunc_mpi_v1.h"
#include "mult_v1.h"
#include "mult_cu_v1.h"

#include "mult_mpi_v1.h"

void mult_hib_0( double *c, double *a, double *b, int m, int k, int n )
{
	int n0;
	int n1;

	n0 = n / 2;
	n1 = n - n0;

	send_int( k, 1 );

	send_vector( a, m * k, 1 );	

	send_vector( b + k * n0, k * n1, 1 );

	if ( n0 < MIN_N_MULT_CU)
	{
		mult_blas( c, a, b, m, k, n0 );
	}
	else
	{
		mult_cuda( c, a, b, m, k, n0 );
	}

	receive_vector( c + m * n0, m * n1, 1 );
}

void mult_hib_1()
{
	int m;
	int k;
	int n1;
	int mxk;
	int kxn1;

	double *a;
	double *b;
	double *c;

	receive_int( &k, 0 );

	receive_new_vector( &a , &mxk, 0 );

	m = mxk / k;

	receive_new_vector( &b , &kxn1, 0 );	

	n1 = kxn1 / k;

	c = ( double* ) malloc( m * n1 * sizeof( double ) );

	if ( n1 < MIN_N_MULT_CU)
	{
		mult_blas( c, a, b, m, k, n1 );
	}
	else
	{
		mult_cuda( c, a, b, m, k, n1 );
	}

	send_vector( c, m * n1, 0 );

	free( c );
	free( b );
	free( a );
}

