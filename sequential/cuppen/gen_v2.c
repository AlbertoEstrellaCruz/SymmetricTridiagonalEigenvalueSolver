#include <stdlib.h>
#include <time.h>
#include <math.h>

#include "magma.h"
#include "magma_lapack.h"

#include "gen_v2.h"


int compare_doubles (const void *a, const void *b)
{
	double ad;
	double bd;

	ad = *( double* ) a;
	bd = *( double* ) b;

	if ( ad <  bd ) return -1;
	else if ( ad >  bd ) return 1;
	else if ( a < b ) return -1;
	else if ( a > b) return 1;
	else return 0;
}

void norm2(double *norm, double *v, double n )
{
	int r;
	double norm2;
	double abs_vr;

	norm2 = 0.0;

	for ( r = 0; r < n; ++r )
	{
		abs_vr = fabs( v[r] );
		norm2 += abs_vr * abs_vr;
	}
	
	*norm = sqrt( norm2 );
	
}


void gen_ab_rand( double *a, double *b, int n ) 
{
	int nm1;
	int r;
	double rn;
	double lim_inf;
	double lim_sup;	

	nm1 = n - 1;

	srand( time( NULL ) );

	lim_inf = 0.0; //-3.0 * n;
	lim_sup = 1.0; //3.0 * n;

	for (r = 0; r < n; ++r)
	{ 
		rn =  ( double ) rand()  / RAND_MAX ;
		a[r] = lim_inf + ( lim_sup - lim_inf ) * rn;
	}

	lim_inf = 0.0; //-30.0 * n;
	lim_sup = 1.0; //30.0 * n;

	for (r = 0; r < nm1; ++r)
	{
		rn = ( double ) rand() / RAND_MAX;
		b[r] = lim_inf + ( lim_sup - lim_inf ) * rn;
	}
}

void gen_ab_mag( double *a, double *b, int n )
{
	int r;
	double lim_inf;
	double lim_sup;

	magma_int_t nm1;
	magma_int_t idist;
	magma_int_t iseed[4]; // [0, 4095] = %4096,  4th must be odd

	nm1 = n - 1;
	idist = 1;

	srand( time( NULL ) );

	iseed[0] = rand() % 4096;
	iseed[1] = rand() % 4096;
	iseed[2] = rand() % 4096;
	iseed[3] = rand() % 4096 / 2 + 1;

	lapackf77_dlarnv( &idist , iseed, &n, a );
	lapackf77_dlarnv( &idist , iseed, &nm1, b );

	lim_inf = 0.0; //-3.0 * n;
	lim_sup = 1.0; //3.0 * n;

	for (r = 0; r < n; ++r)
	{ 
		a[r] = lim_inf + ( lim_sup - lim_inf ) * a[r];
	}

	lim_inf = 0.0; //-30.0 * n;
	lim_sup = 1.0; //30.0 * n;

	for (r = 0; r < nm1; ++r)
	{
		b[r] = lim_inf + ( lim_sup - lim_inf ) * b[r];
	}
}

void gen_vd_rand( double *v, double *d, int n ) 
{
	int r;
	double rn;
	double lim_inf;
	double lim_sup;
	double norm;

	srand( time( NULL ) );

	lim_inf = 0.0; //-10.0 * n;
	lim_sup = 1.0; //10.0 * n;

	for (r = 0; r < n; ++r)
	{ 
		rn =  ( double ) rand()  / RAND_MAX ;
		v[r] = lim_inf + ( lim_sup - lim_inf ) * rn;
	}

	norm2( &norm, v, n );

	for (r = 0; r < n; ++r)
	{
		v[r] = v[r] / norm;
	}	

	lim_inf = 0.0; //-100.0 * n;
	lim_sup = 1.0; //100.0 * n;

	for (r = 0; r < n; ++r)
	{
		rn = ( double ) rand() / RAND_MAX;
		d[r] = lim_inf + ( lim_sup - lim_inf ) * rn;
	}

	qsort( d, n, sizeof( double ), compare_doubles );
}

void gen_vd_mag( double *v, double *d, int n ) 
{
	int r;
	double lim_inf;
	double lim_sup;
	double norm;

	magma_int_t idist;
	magma_int_t iseed[4]; // [0, 4095] = %4096,  4th must be odd

	idist = 1;

	srand( time( NULL ) );

	iseed[0] = rand() % 4096;
	iseed[1] = rand() % 4096;
	iseed[2] = rand() % 4096;
	iseed[3] = rand() % 4096 / 2 + 1;

	lapackf77_dlarnv( &idist , iseed, &n, v );	
	lapackf77_dlarnv( &idist , iseed, &n, d );

	lim_inf = 0.0; //-10.0 * n;
	lim_sup = 1.0; //10.0 * n;

	for (r = 0; r < n; ++r)
	{ 
		v[r] = lim_inf + ( lim_sup - lim_inf ) * v[r];
	}

	norm2( &norm, v, n );

	for (r = 0; r < n; ++r)
	{
		v[r] = v[r] / norm;
	}	

	lim_inf = 0.0; //-100.0 * n;
	lim_sup = 1.0; //100.0 * n;

	for (r = 0; r < n; ++r)
	{
		d[r] = lim_inf + ( lim_sup - lim_inf ) * d[r];
	}

	qsort( d, n, sizeof( double ), compare_doubles );
}

