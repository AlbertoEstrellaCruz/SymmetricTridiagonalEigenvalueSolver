#ifndef COMMFUNC_MPI_H
#define COMMFUNC_MPI_H

#ifdef __cplusplus
extern "C" {
#endif

void send_vector( double *v, int n, int rank_to );
void receive_vector( double *v, int n, int rank_from );
void receive_new_vector( double **v, int *n, int rank_from );

void send_int( int i, int rank_to );
void receive_int( int *i, int rank_from );

#ifdef __cplusplus
}
#endif

#endif /* COMMFUNC_MPI_H */

