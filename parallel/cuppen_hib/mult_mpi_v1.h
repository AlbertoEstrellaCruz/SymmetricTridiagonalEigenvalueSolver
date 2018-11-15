#ifndef MULT_MPI_H
#define MULT_MPI_H

#ifdef __cplusplus
extern "C" {
#endif

void mult_hib_0( double *c , double *a , double *b, int m, int k, int n );
void mult_hib_1();

#ifdef __cplusplus
}
#endif

#endif /* MULT_MPI_H */

