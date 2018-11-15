#ifndef MAG_H
#define MAG_H

#ifdef __cplusplus
extern "C" {
#endif

void call_dstedc( double *lambda, double *q, double *a, double *b, int n );
void call_dstedx( double *lambda, double *q, double *a, double *b, int n );

#ifdef __cplusplus
}
#endif

#endif /* MAG_H */

