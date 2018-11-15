#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "deflate_v1.h"

double get_max_d( double *d, int n )
{
	double abs_d_0;
	double abs_d_nm1;
	double max_d;

	abs_d_0 = fabs( d[0] );
	abs_d_nm1 = fabs( d[ n-1 ] );

	if ( abs_d_0  > abs_d_nm1 )
	{
    	max_d = abs_d_0;
	}
	else
	{
    	max_d = abs_d_nm1;
	}

	return max_d;
}

double get_max_v( double *v, int n )
{
	int r;

	double max_v;
	double abs_v_r;

	max_v = fabs( v[0] );
	abs_v_r = 0.0;

	for ( r = 1; r < n; ++r )
	{
	    abs_v_r = fabs( v[r] );

	    if ( abs_v_r > max_v )
		{
	        max_v = abs_v_r;
	    }
	}

	return max_v;
}

double get_dlapy2( double a, double b)
{
	double a_abs;
	double b_abs;
	double max_ab;
	double min_ab;
	double min_over_max;
	double r;

	a_abs = fabs( a );
	b_abs = fabs( b );

	if ( a > b )
	{
    	max_ab = a_abs;
    	min_ab = b_abs;
	}
	else
	{
    	max_ab = b_abs;
    	min_ab = a_abs;
	}

	if ( 0.0 == min_ab )
	{
    	r = max_ab;
	}
	else
	{
    	min_over_max = min_ab / max_ab;
    	r = max_ab * sqrt( 1.0 + min_over_max * min_over_max );
	}

	return r;
}

void get_givens_csr( double *pc, double *ps, double *proot, double a, double b )
{
	double c;
	double s;
	double root;

	if ( 0.0 == b )
    {
    	c = 1.0;
    	s = 0.0;
    	root = a;
    }
	else
    {
    	root = get_dlapy2( a, b );
    
    	s = -b / root;
    
    	c = a / root;
	}

	*pc = c;
	*ps = s;
	*proot = root;
}

void deflate( double *g_s, double *g_c, int *g_i, int *g_j, int *p_def, double *d_def, double *v_def, int *pn_not_def, int *pn_giv, double *d_sort, double *v_sort, int n, double eps )
{
	int row;
	int n_not_def;
	int n_giv;
	int indx_def;
	int prev_row;
	int new_row;

	double max_d;
	double max_v;
	double max_dv;
	double tol_v;
	double tol_d;
	double a;
	double b;
	double dif;
	double c;
	double s;
	double root;

	int *deflated;

  // Use two loops for cache in c code

	for ( row = 0; row < n; ++row )
	{
	    d_def[row] = d_sort[row];
	}

	for ( row = 0; row < n; ++row )
	{
	    v_def[row] = v_sort[row];
	}

	max_d = get_max_d( d_sort, n );
	max_v = get_max_v( v_sort, n );

	if ( max_d > max_v )
	{
    	max_dv = max_d;
	}
	else
	{
    	max_dv = max_v;
	}

	tol_v = 8.0 * eps * max_dv;
	tol_d = 8.0 * eps * max_dv;

	deflated = ( int* ) calloc( n, sizeof( int ) );

	n_giv = -1;

	indx_def = n;

	for ( row = 0; row < n; ++row )
    {
    	if ( fabs( v_def[row] ) <= tol_v )
        {
			//v_def[row] = 0.0; // not necessary

    	    deflated[row] = 1;

    	    --indx_def;

    	    p_def[indx_def] = row;
    	}
	}

	// get prev and new  bottom top

	prev_row = n - 1;

	while ( prev_row >= 0 && deflated[prev_row] != 0 )
    {
	    --prev_row;
	}

	new_row = prev_row - 1;

	while ( new_row >= 0 )
    {
	    while ( new_row >= 0 && deflated[new_row] != 0 )
        {
	        --new_row;    
	    }

	    if ( new_row >= 0 )
        {
	      // detect deflation
        
	        a = v_def[new_row];
        
	        b = v_def[prev_row];
        
	        dif = d_def[new_row] - d_def[prev_row];
        
	        get_givens_csr( &c, &s, &root, a, b );
        
	        if ( fabs( dif *  c * s ) <= tol_d )
            {
	            //printf( "Givens Rotation\n" );

	          // save rotation

	            ++n_giv;
            
	            g_c[n_giv] = c;
            
	            g_s[n_giv] = s;
            
	            g_i[n_giv] = new_row;
            
	            g_j[n_giv] = prev_row;
            
	          // do deflation
            
	            v_def[new_row] = root; // c*a - s*b;
            
	            //v_def[prev_row] = 0.0; // c*a - s*b; // not necessary
            
	          // with d very necessary

				a = d_def[new_row];

				b = d_def[prev_row];
            
	            d_def[new_row] = a * ( c*c ) + b * ( s*s );
            
	            d_def[prev_row] = a * ( s*s ) + b * ( c*c );
            
	          // permutation

	            deflated[prev_row] = 2;

	            --indx_def;

	            p_def[indx_def] = prev_row;
	        }
        
	        prev_row = new_row;
        
	        new_row = prev_row - 1;
	    }
	}

  // complete permutation

	n_not_def = -1;

	for ( row = 0; row < n; ++row )
    {
	    if ( 0 == deflated[row] )
        {
	        ++n_not_def;

	        p_def[n_not_def] = row;
	    }
	}

	free( deflated );

	++n_not_def;
	++n_giv;
/*
	if (n_not_def == indx_def)
	{
	    printf( "n_not_def = %d\n", n_not_def );
	    printf( "n_giv = %d\n", n_giv );
	}
*/
	*pn_not_def = n_not_def;
	*pn_giv = n_giv;
}

