#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>

#include "commfunc_mpi_v1.h"

void send_vector( double *v , int n, int rank_to )
{	
	MPI_Send( v, n, MPI_DOUBLE, rank_to, 0, MPI_COMM_WORLD );
}

void receive_vector( double *v , int n, int rank_from )
{
	MPI_Status status;
	
	MPI_Recv( v, n, MPI_DOUBLE, rank_from, MPI_ANY_TAG, MPI_COMM_WORLD, &status );
}

void receive_new_vector( double **v , int *n, int rank_from )
{
	MPI_Status status;
	
	MPI_Probe(rank_from, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

	MPI_Get_count(&status, MPI_DOUBLE, n);

	*v = ( double* ) malloc( (*n) * sizeof(double) );
	
	MPI_Recv( *v, *n, MPI_DOUBLE, rank_from, MPI_ANY_TAG, MPI_COMM_WORLD, &status );

}

void send_int( int i, int rank_to )
{
	MPI_Send( &i, 1, MPI_INT, rank_to, 0, MPI_COMM_WORLD );
}

void receive_int( int *i, int rank_from )
{
	MPI_Status status;
	
	MPI_Recv( i, 1, MPI_INT, rank_from, MPI_ANY_TAG, MPI_COMM_WORLD, &status );
}

