#include "magma.h"
#include "magma_lapack.h"

#include "mult_v1.h"

void mult_blas( double *c , double *a , double *b, int m, int k, int n )
{
	double alpha;
	double beta;

	alpha = 1.0;
	beta = 0.0;

	blasf77_dgemm(  "N", "N", &m, &n, &k, &alpha, a, &m, b, &k, &beta, c, &m );
}

