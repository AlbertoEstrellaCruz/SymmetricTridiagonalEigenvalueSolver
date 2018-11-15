#include <stdlib.h>
#include <stdio.h>

#include "hib_defs_v1.h"
#include "zerodandc_n_v1.h"
#include "calculate_vi_v1.h"
#include "calculate_normq_v1.h"
#include "calculate_q_v1.h"

#include "rank_one_v1.h"

void calculate_v2_r( double *v2, double *v, int r )
{
	v2[r] = v[r] * v[r];
}

void calculate_v2( double *v2, double *v, int n )
{
	int r;

	#pragma omp parallel for num_threads( HIB_OMP_NT ) default( none ) \
	 private( r ) shared( v2, v, n )
	for ( r = 0; r < n; ++r )
	{
		calculate_v2_r( v2, v, r );
	}
}

void solve_secular( double *delta, double *lambda, double *dmlambda, double *d, double *v2, int n, double eps )
{
	int i;

	#pragma omp parallel for num_threads( HIB_OMP_NT ) default( none ) \
	 private( i ) shared( delta, lambda, dmlambda, d, v2, n, eps )
	for ( i = 0; i < n; ++i )
	{
		zerodandc_n( delta + i * n, lambda + i, dmlambda + i * n, d, v2, i, n, eps );
	}
}

void correct_v( double *v_corr, double *d, double *dfmlc, double *v, int n)
{
	int i;

	#pragma omp parallel for num_threads( HIB_OMP_NT ) default( none ) \
	 private( i ) shared( v_corr, d, dfmlc, v, n )
	for ( i = 0; i < n; ++i )
	{
		calculate_vi( v_corr + i, d, dfmlc, v, i, n );
	}
}

void rank_one( double *q, double *lambda, double *d, double *v, int n, double eps )
{
	double *delta;
	double *dfmlc;
	double *v2;
	double *v_corr;	
	double *normq;

	delta = ( double* ) malloc( n * n * sizeof( double ) );
	dfmlc = ( double* ) malloc( n * n * sizeof( double ) );
	v2 = ( double* ) malloc( n * sizeof( double ) );

	v_corr = ( double* ) malloc( n * sizeof( double ) );
	
	normq = ( double* ) malloc( n * sizeof( double ) );
	
	calculate_v2( v2, v, n );

  // Calculando los valores propios de la ecuacion secular

	solve_secular( delta, lambda, dfmlc, d, v2, n, eps );

	free(delta);
	free(v2);

  // Calculando la correccion de v 

	correct_v( v_corr, d, dfmlc, v, n );

  // Calculando los vectores propios en Q

	calculate_normq( normq, dfmlc, v_corr, n);

	calculate_q( q, v_corr, dfmlc, normq, n );

	free(normq);
	free(v_corr);
	free(dfmlc);
}

