#include <stdlib.h>
#include <stdio.h>

#include <mpi.h>

#include "cublas_v2.h"

#include "cuppen_hib_v1.h"

int main( int argc, char** argv )
{
	int rank;
	int size;

	MPI_Init( &argc, &argv );

	MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    MPI_Comm_size( MPI_COMM_WORLD, &size );

	if ( 0 == rank )
	{
		node_0( rank, size, argc, argv );
	}
	else if ( 1 == rank )
	{
		node_1( rank, size );
	}

    MPI_Finalize();

    return EXIT_SUCCESS;
}

