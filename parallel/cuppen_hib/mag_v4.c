#include <stdlib.h>
#include <stdio.h>

#include "hib_time_v2.h"

#include "magma.h"
#include "magma_lapack.h"

#include "mag_v2.h"

void call_dstedc( double *lambda, double *q, double *a, double *b, int n )
{
  // Declaration of variables

	char compz;
	double *e;
	double *work;
	magma_int_t lwork;
	magma_int_t *iwork;
	magma_int_t liwork;
	magma_int_t info;

	double opti_lwork;
	magma_int_t opti_liwork;

	magma_int_t ione;
	magma_int_t nm1;
/*
	real_Double_t start;
	real_Double_t end;
	double cpu_time;
*/
  // Initialization of variables

	compz = 'I';

  // Query for workspace sizes	

	lwork = -1;
	liwork = -1;

	lapackf77_dstedc( &compz, &n, NULL, NULL, NULL, &n, &opti_lwork, &lwork, &opti_liwork, &liwork, &info );

  // Alloc worspace

	lwork = ( magma_int_t ) opti_lwork;
	liwork = opti_liwork;

	//printf("lwork=%d, liwork=%d\n", lwork, liwork);

	work = ( double* ) malloc( lwork * sizeof( double ) );
	iwork = ( magma_int_t* ) malloc( liwork * sizeof( magma_int_t ) );

  // Copy values of a and b to avoid overwrite or destruction

	e = ( double* ) malloc( (n-1) * sizeof( double ) );

	ione = 1;
	nm1 = n - 1;

	blasf77_dcopy( &n, a, &ione, lambda, &ione );
	blasf77_dcopy( &nm1, b, &ione, e, &ione );

	//start =  get_posix_time();

  // Compute eigensystem using lapack

	lapackf77_dstedc( &compz, &n, lambda, e, q, &n, work, &lwork, iwork, &liwork, &info );
/*
	end = get_posix_time() - start;
	cpu_time = end;

	printf( " call_dstedc cpu time : %7.5f sec.\n" , cpu_time ); // Lapack time
*/
  // Free memory

	free( e );
	free( iwork );
	free( work );
}

void call_dstedx( double *lambda, double *q, double *a, double *b, int n ) 
{
  // Declaration of variables

	magma_range_t range;
	double *e;
	double *work;
	magma_int_t lwork;
	magma_int_t *iwork;
	magma_int_t liwork;
	double *dwork;
	magma_int_t info;

	double opti_lwork;
	magma_int_t opti_liwork;

	magma_int_t ione;
	magma_int_t nm1;
/*
	real_Double_t start;
	real_Double_t end;
	double gpu_time;
*/
  // Initialization of variables

	range = MagmaRangeAll;
	//printf("range = %d \n", range);

  // Query for workspace sizes

	magma_dstedx( range, n, 0.0, 0.0, 0, 0, NULL, NULL, NULL, n, &opti_lwork, -1, &opti_liwork, -1, NULL, &info );

  // Alloc worspace

	lwork = ( magma_int_t ) opti_lwork;
	liwork = opti_liwork;

	work = ( double* ) malloc( lwork * sizeof( double ) );
	iwork = ( magma_int_t* ) malloc( liwork * sizeof( magma_int_t ) );

	magma_dmalloc( &dwork, 3*n*(n/2 + 1) );

  // Copy values of a and b to avoid overwrite or destruction

	e = ( double* ) malloc( (n-1) * sizeof( double ) );

	ione = 1;
	nm1 = n - 1;

	blasf77_dcopy( &n, a, &ione, lambda, &ione );
	blasf77_dcopy( &nm1, b, &ione, e, &ione );

	//start =  get_posix_time();

  // Compute eigensystem using magma

	magma_dstedx( range, n, 0.0, 0.0, 0, 0, lambda, e, q, n, work, lwork, iwork, liwork, dwork, &info );
/*
	end = get_posix_time() - start;
	gpu_time = end;

	printf( " call_dstedx gpu time : %7.5f sec.\n" , gpu_time );	// Magma version time
*/
  // Free memory

	magma_free( dwork );

	free( e );
	free( iwork );
	free( work );
}

