#ifndef IO_H
#define IO_H

#ifdef __cplusplus
extern "C" {
#endif

//All reading/writing values must have 8 bytes, doubles and ints
void write_vd( char *file_path, double *v, double *d, int n );
void read_vd( double **v, double **d, int *n, char *file_path );
void write_lq( char *file_path, double *l , double *q, int n );
void read_lq( double **l, double **q, int *n, char *file_path );
void write_ab( char *file_path, double *a, double *b, int n );
void read_ab( double **a, double **b, int *n, char *file_path );

#ifdef __cplusplus
}
#endif

#endif /* IO_H */

