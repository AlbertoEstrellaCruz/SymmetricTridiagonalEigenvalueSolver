#include <time.h>
#include <mpi.h>
#include <omp.h>

#include "hib_time_v2.h"

double get_posix_time()
{
	struct timespec now;

	clock_gettime( CLOCK_REALTIME, &now );

	return now.tv_sec + now.tv_nsec * 1e-9;
}

double get_cycles_time()
{
	clock_t ticks;

	ticks = clock();

	return ( ( double ) ticks ) / CLOCKS_PER_SEC;
}

double get_openmp_time()
{
	return omp_get_wtime();
}

double get_mpi_time()
{
	return MPI_Wtime();
}

void start_cuda_time( hib_cuda_timer_t *t )
{
	cudaEventCreate( &t->start );
	cudaEventCreate( &t->stop );

	cudaEventRecord( t->start, 0 );
}

double finish_cuda_time( hib_cuda_timer_t *t )
{
	float elapsed_time;

	cudaEventRecord( t->stop, 0 );
	cudaEventSynchronize( t->stop );

	cudaEventElapsedTime( &elapsed_time, t->start, t->stop );

	cudaEventDestroy( t->start );
	cudaEventDestroy( t->stop );

	return ( double ) elapsed_time * 1e-3;
}

