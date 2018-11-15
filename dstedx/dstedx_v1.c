#include <stdlib.h>
#include <stdio.h>

#include "cublas_v2.h"
#include "magma.h"

#include "gen_v3.h"
#include "io_v2.h"
#include "err_v2.h"
#include "hib_time_v2.h"
#include "mag_v2.h"

int main( int argc, char **argv )
{
	int n;
	//int nm1;
/*
	double err_l;
	double err_q;
*/
	double start_time;
	double final_time;

	double *a;
	double *b;
	double *q;
	double *l;

    if (argc != 4) 
	{
        fprintf(stderr, "Usage : %s <matrix_size> <path_file_matrix> <path_file_solution> \n", argv[0]);
        exit(EXIT_FAILURE);
    }

	n = atoi( argv[1] );
	//nm1 = n - 1;

	//a = ( double* ) malloc( n * sizeof( double ) );
	//b = ( double* ) malloc( nm1 * sizeof( double ) );
	l = ( double* ) malloc( n * sizeof( double ) );
	q = ( double* ) malloc( n * n * sizeof( double ) );

	magma_init();

  // Generate random problem

	//gen_ab_mag( a, b, n );

	read_ab( &a, &b, &n, argv[2] );
	//write_ab( "/home/alberto/Dropbox/Tesis/matrices/mtab2", a, b, n );

	start_time = get_posix_time();

	call_dstedx( l, q, a, b, n );

	final_time = get_posix_time();

	printf( "%d\t%.7f\n", n, final_time - start_time );

	write_lq( argv[3], l, q, n );
	//read_lq( &l, &q, &n, "/home/alberto/Dropbox/Tesis/matrices/slq2" );
/*
	err_t_l( &err_l, a, b, q, l, n );
	err_t_q( &err_q, q, n );
	printf("  eigenpairs error: %.16e \n", err_l);
	printf("  othogonality error : %.16e \n", err_q);
*/
	magma_finalize();

	free( b );
	free( a );
	free( q );
	free( l );
	//free( b );
	//free( a );

	return EXIT_SUCCESS;
}

