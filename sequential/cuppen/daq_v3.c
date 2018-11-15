#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "io_v2.h"
#include "mag_v2.h"
#include "merge_v1.h"
#include "deflate_v1.h"
#include "rank_one_v1.h"
#include "mult_v1.h"
 
#include "daq_v1.h"

void apply_givens_rot( double *x, double *y, double c, double s, int n )
{
	int i;

	double xi;
	double yi;

	for ( i = 0; i < n; ++i )
	{
    	xi = x[i];
    	yi = y[i];
    
    	x[i] = c * xi - s * yi;
    	y[i] = s * xi + c * yi;

	}
}

void unrotate( double *q, double *g_s, double *g_c, int *g_i,  int *g_j , int n_giv, int n)
{
	int k;
	int i;
	int j;
	double c;
	double s;

	for ( k = 0; k < n_giv; ++k )
	{    
    	i = g_i[k];
    	j = g_j[k];
    	c = g_c[k];
    	s = g_s[k];

		apply_givens_rot( q + i * n, q + j * n, c, s, n );

	}
}

void daq( double *q_t, double *lambda_r1, double *a, double *b, int n, double eps )
{
	int m;
	int mm1;
	int k;
	int c;
	int r;
	int n_not_def;
	int n_giv;

	double rho_ori;
	double rho_scale;

	int *p_sort;
	int *p_def;
	int *g_i;
	int *g_j;

	double *a_1;
	double *b_1;
	double *a_2;
	double *b_2;
	double *q_1;
	double *q_2;
	double *d_1;
	double *d_2;
	double *v;
	double *q_12;
	double *v_sort;
	double *d_sort;
	double *v_def_ori;
	double *d_def_ori;
	double *g_s;
	double *g_c;
	double *v_def;
	double *d_def;
	double *lambda_r1_ori;
	double *q_r1;

	m = n / 2;
	mm1 = m - 1;

	rho_ori = b[ mm1 ];

	a_1 = ( double* ) malloc( m * sizeof( double ) );
	b_1 = ( double* ) malloc( mm1 * sizeof( double ) );

  // Scale a and b to get rho = 1/2 and split problem

	for ( k = 0; k < m; ++k )
	{
    	a_1[k] = ( 1.0 / (2.0 * rho_ori) ) * a[k];
	}

	for ( k = 0; k < mm1; ++k )
	{
    	b_1[k] = ( 1.0 / (2.0 * rho_ori) ) * b[k];
	}

	a_2 = ( double* ) malloc( m * sizeof( double ) );
	b_2 = ( double* ) malloc( mm1 * sizeof( double ) );

	for ( k = 0; k < m; ++k )
	{
    	a_2[k] =  ( 1.0 / (2.0 * rho_ori) ) * a[ m+k ];
	}

	for ( k = 0; k < mm1; ++k )
	{
    	b_2[k] =  ( 1.0 / (2.0 * rho_ori) ) * b[ m+k ];
	}

	rho_scale = 0.5;

	a_1[mm1] = a_1[mm1] - rho_scale;
	a_2[0] = a_2[0] - rho_scale;

  // Solve subproblems

	q_1 = ( double* ) malloc( m * m * sizeof( double ) );
	q_2 = ( double* ) malloc( m * m * sizeof( double ) );

	d_1 = ( double* ) malloc( m * sizeof( double ) );
	d_2 = ( double* ) malloc( m * sizeof( double ) );

	call_dstedc( d_1, q_1, a_1, b_1, m );
	call_dstedc( d_2, q_2, a_2, b_2, m );

	free( b_2 );
	free( a_2 );
	free( b_1 );
	free( a_1 );

  // Scale v to get norm one

	v = ( double* ) malloc( n * sizeof( double ) );

	for ( k = 0; k < m; ++k )
	{
    	v[k] = ( 1.0 / sqrt(2.0) ) * q_1[ mm1 + m * k ];
	}

	for ( k = 0; k < m; ++k )
	{
    	v[ m+k ] = ( 1.0 / sqrt(2.0) ) * q_2[ m * k ];
	}

	d_sort = ( double* ) malloc( n * sizeof ( double ) );
	p_sort = ( int* ) malloc( n * sizeof( int ) );

  // Sort d_1 and d_2 equivalen to PDD'

	merge_d( d_sort, p_sort, d_1, d_2, m, m, n );

	free( d_2 );
	free( d_1 );

  // Rearrenge v with Pv where (Pv)' = v'P' so v_sort*v_sort' is eq (Pv)(Pv)'= P(v*v')P'

	v_sort = ( double* ) malloc( n * sizeof( double ) );

	for ( k = 0; k < n; ++k )
	{
    	v_sort[k] = v[ p_sort[k] ];
	}

	free( v );

	v_def_ori = ( double* ) malloc( n * sizeof( double ) );
	d_def_ori = ( double* ) malloc( n * sizeof( double ) );

	p_def = ( int* ) malloc( n * sizeof( int ) ) ;

	g_s = ( double* ) malloc( n * sizeof( double ) );
	g_c = ( double* ) malloc( n * sizeof( double ) );
	g_i = ( int* ) malloc( n * sizeof( int ) );
	g_j = ( int* ) malloc( n * sizeof( int ) );

	//write_vd( "/home/alberto/Dropbox/Tesis/matrices/mtvd2", v_sort, d_sort, n );
	//read_vd( &v_sort, &d_sort, &n, "/home/alberto/Dropbox/Tesis/matrices/mtvd1" );

  // Deflate problem to reduce size of problem and reach condition to solve secular equation

	deflate( g_s, g_c, g_i, g_j, p_def, d_def_ori, v_def_ori, &n_not_def, &n_giv, d_sort, v_sort, n, eps );

	free( v_sort );
	free( d_sort );

	d_def = ( double* ) malloc( n * sizeof( double ) );

	for (  k = 0; k < n; ++k )
	{
    	d_def[k] = d_def_ori[ p_def[k] ];
	}

	free( d_def_ori );

	v_def = ( double* ) malloc( n * sizeof( double ) );

	for ( k = 0; k < n; ++k )
	{
    	v_def[k] = v_def_ori[ p_def[k] ];
	}

	free( v_def_ori );
	
	lambda_r1_ori = ( double* ) malloc( n_not_def * sizeof( double ) );
	q_r1 = ( double* ) malloc( n_not_def * n_not_def * sizeof( double ) );

  // Solve rank-one update

	rank_one( q_r1, lambda_r1_ori, d_def, v_def, n_not_def, eps );

	free( v_def );

  // Fill complete array of eigenvalues

	for ( k = 0; k < n_not_def; ++k )
	{
    	lambda_r1[k] = lambda_r1_ori[k];
	}

	free( lambda_r1_ori );

	for ( k = n_not_def; k < n; ++k )
	{
    	lambda_r1[k] = d_def[k];
	}

	free( d_def );

  // Reset i j givens coordinates

	for ( k = 0; k < n_giv; ++k )
	{
    	g_i[k] = p_sort[ g_i[k] ];
    	g_j[k] = p_sort[ g_j[k] ];
	}

  // Build q_12 matrix wit q_1 and q_2

	q_12 = ( double* ) calloc( n * n, sizeof( double ) );

	for ( c = 0; c < m; ++c )
	{
    	for ( r  = 0; r < m; ++r )
		{
    	    q_12[ r + n * c ] = q_1[ r + m * c ];
	    }
	}

	free( q_1 );

	for ( c = 0; c < m; ++c )
	{
    	for ( r = 0; r < m; ++r  )
		{
    	    q_12[ ( r + m ) + n * ( c + m ) ] = q_2[ r + m * c ];
    	}
	}

	free( q_2 );

  // Unrotate matrix Q_12 = Q_12 * G_def'

	unrotate( q_12, g_s, g_c, g_i,  g_j , n_giv, n);

	free( g_j );
	free( g_i );
	free( g_c );
	free( g_s );

  // Permute columns AP'

	for ( c  = 0; c < n; ++c )
	{
    	k = p_sort[ p_def[c] ];
	    for ( r = 0; r < n; ++r )
		{
    	    q_t[ r + n * c ] = q_12[ r + n * k ];
	    }
	}

	free( p_def );
	free( p_sort );

  // Copy unchange multiplication result Q_12(:, 0:n_not_def-1)

	for ( c = 0; c < n_not_def; ++c )
	{
    	for ( r = 0; r < n; ++r )
		{
    	    q_12[ r + n * c ] = q_t[ r + n * c ];
		}
	}

  // Compute reduced multiplication Q_12(:, 0:n_not_def-1) * Q_r1

	mult_blas( q_t , q_12 , q_r1, n, n_not_def, n_not_def );

	free( q_12 );
	free( q_r1 );

  // Unscale lambda to get eigenvalues of the original matrix

	for ( r  = 0; r < n; ++r )
	{
    	lambda_r1[ r ] = (2.0 * rho_ori ) * lambda_r1[ r ];
	}
}

