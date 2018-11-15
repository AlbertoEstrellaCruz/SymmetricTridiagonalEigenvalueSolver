#ifndef DEFLATE_H
#define DEFLATE_H

#ifdef __cplusplus
extern "C" {
#endif

void deflate( double *g_s, double *g_c, int *g_i, int *g_j, int *p_def, double *d_def, double *v_def, int *pn_not_def, int *pn_giv, double *d_sort, double *v_sort, int n, double eps );

#ifdef __cplusplus
}
#endif

#endif /* DEFLATE_H */

