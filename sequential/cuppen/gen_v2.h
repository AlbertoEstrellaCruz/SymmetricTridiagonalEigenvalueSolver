#ifndef GEN_H
#define GEN_H

#ifdef __cplusplus
extern "C" {
#endif

void gen_ab_rand( double *a, double *b, int n );
void gen_ab_mag( double *a, double *b, int n );
void gen_vd_rand( double *v, double *d, int n );
void gen_vd_mag( double *v, double *d, int n );

#ifdef __cplusplus
}
#endif

#endif /* GEN_H */

