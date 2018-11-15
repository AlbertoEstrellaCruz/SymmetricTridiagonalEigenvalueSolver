#include <stdio.h>
#include <math.h>

#include "err_v2.h"

void err_mr1_l( double *err_l, double *d, double *v, double *q, double *lambda, int n )
{
	int c;
	int r;
	int k;
	int nxc;

	double e;
	double e_rc;

	e = 0.0;

	for ( c = 0;  c < n; ++c )
	{
		nxc = n * c;

    	for ( r = 0; r < n; ++r )
        {
    	    e_rc = 0.0;

    	    for ( k = 0; k < r; ++k )
			{
    	        e_rc += v[r] * v[k] * q[ k + nxc ];
    	    }

    	    e_rc += ( d[r] + v[r] * v[r] ) * q[ r + nxc ];

    	    for ( k = r + 1; k < n; ++k )
			{
    	        e_rc += v[r] * v[k] * q[ k+ nxc ];
    	    }

    	    e_rc = fabs( e_rc - q[ r + nxc ] * lambda[c] );
        
    	    e += e_rc *e_rc;
    	}
	}

	*err_l = sqrt(e);
}

void err_t_l( double *err_l, double *a, double *b, double *q, double *lambda, int n )
{
	int c;
	int r;
	int nxc;

	double e;
	double e_rc;

	e = 0.0;

	for ( c = 0; c < n; ++c )
    {
		nxc = n * c;

    	e_rc = a[0] * q[ nxc ] + b[0] * q[ 1 + nxc ];

    	e_rc = fabs( e_rc - q[ nxc ] * lambda[c] );
    	
    	e = e + e_rc * e_rc;    
	}

	for ( c = 0; c < n; ++c )
    {
		nxc = n * c;

    	e_rc = b[ n-2 ] * q[ n-2 + nxc ] + a[ n-1 ] * q[ n-1 + nxc ];
    
    	e_rc = fabs( e_rc - q[ n-1 + nxc ] * lambda[c] );

    	e = e + e_rc * e_rc;
	}

	for ( c = 0; c < n; ++c )
	{
		nxc = n * c;

    	for ( r = 1; r < (n-1); ++r )
        {
        	e_rc = b[ r-1 ] * q[ r-1 + nxc ] + a[r] * q[ r + nxc ] + b[r] * q[ r+1  + nxc ];

        	e_rc = fabs( e_rc - q[ r + nxc ] * lambda[c] );
        
        	e = e + e_rc * e_rc;

    	}
	}

	*err_l = sqrt(e);
}

void err_t_q(double *err_q, double *q, int n )
{
	int c;
	int r;
	int k;

	double e;
	double e_rc;

	e = 0.0;

	for ( c = 0; c < n; ++c )
	{

		for ( r = 0; r < n; ++r )
		{
        
        	e_rc = 0.0;
        
        	for (k = 0; k < n; ++k)
			{
            	e_rc += q[ r + n * k ] * q[ c + n * k ];
        	}
        
        	if ( r == c )
			{
            	e_rc -= 1.0;
        	}
        
        	e_rc = fabs( e_rc );
        
        	e += e_rc * e_rc;
        
    	}
	}

	*err_q = sqrt(e);
}

