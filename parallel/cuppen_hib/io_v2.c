#include <stdlib.h>
#include <stdio.h>

#include "io_v2.h"

void write_vd( char *file_path, double *v, double *d, int n )
{
	int64_t ni;
	FILE *stream;

	ni = ( int64_t ) n;

	stream = fopen( file_path, "wb" );

	fwrite( &ni, sizeof( int64_t ), 1, stream );

	fwrite( v, sizeof( double ), n, stream );

	fwrite( d, sizeof( double ), n, stream );

	fclose(stream);
}

void read_vd( double **v, double **d, int *n, char *file_path )
{
	int64_t ni;
	FILE *stream;

	stream = fopen( file_path,  "rb" );	

	fread( &ni, sizeof(int64_t), 1, stream );
	
	*n = ( int ) ni;

	*v = ( double* ) malloc( ni * sizeof( double ) );
	*d = ( double* ) malloc( ni * sizeof( double ) ); 

	fread( *v, sizeof( double ), ni,  stream );
	fread( *d, sizeof( double ), ni, stream );

	fclose(stream);
}

void write_lq( char *file_path, double *l , double *q, int n )
{
	int64_t ni;
	FILE *stream;

	ni = ( int64_t ) n;

	stream = fopen( file_path, "wb" );

	fwrite( &ni, sizeof( int64_t ), 1, stream );

	fwrite( l, sizeof( double ), n, stream );

	fwrite( q, sizeof( double ), n*n, stream );

	fclose(stream);
}

void read_lq( double **l, double **q, int *n, char *file_path)
{
	int64_t ni;
	int nxn;
	FILE *stream;

	stream = fopen( file_path,  "rb" );	

	fread( &ni, sizeof(int64_t), 1, stream );

	*n = ( int ) ni;
	nxn = ( int ) ( ni * ni );

	*l = ( double* ) malloc( ni * sizeof( double ) );
	*q = ( double* ) malloc( nxn * sizeof( double ) ); 

	fread( *l, sizeof( double ), ni,  stream );
	fread( *q, sizeof( double ), nxn, stream );

	fclose(stream);	
}

void write_ab( char *file_path, double *a, double *b, int n ) 
{
	int64_t ni;
	FILE *stream;

	ni = ( int64_t ) n;

	stream = fopen( file_path, "wb" );

	fwrite( &ni, sizeof( int64_t ), 1, stream );

	fwrite( a, sizeof( double ), n, stream );

	fwrite( b, sizeof( double ), n-1, stream );

	fclose(stream);
}

void read_ab( double **a, double **b, int *n, char *file_path ) 
{	
	int64_t ni;
	int nm1;
	FILE *stream;

	stream = fopen( file_path,  "rb" );	

	fread( &ni, sizeof(int64_t), 1, stream );
	
	*n = ( int ) ni;
	nm1 = ( int ) ni - 1;

	*a = ( double* ) malloc( ni * sizeof( double ) );
	*b = ( double* ) malloc( nm1 * sizeof( double ) ); 

	fread( *a, sizeof( double ), ni,  stream );
	fread( *b, sizeof( double ), nm1, stream );

	fclose(stream);
}

