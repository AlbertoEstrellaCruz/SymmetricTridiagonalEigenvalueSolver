#ifndef ERR_H
#define ERR_H

#ifdef __cplusplus
extern "C" {
#endif

void err_mr1_l( double *err_l, double *d, double *v, double *q, double *lambda, int n );
void err_t_l( double *err_l, double *a, double *b, double *q, double *lambda, int n );
void err_t_q( double *err_q, double *q, int n );

#ifdef __cplusplus
}
#endif

#endif /* ERR_H */

