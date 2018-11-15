#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "zerodandc_n_v1.h"

void zerodandc_n( double *delta, double *plambda, double *dmlambda, double *d, double *v2, int i, int n, double eps )
{
	int max_iter;
	int n_iter;
	int j;

	double di;
	double eta;
	double g;
	double dip1;
	double lambda;
	double psi;
	double phi;
	double f;
	double mu;
	double dpsi;
	double dphi;
	double Deltai;
	double Deltaip1;
	double a;
	double b;
	double c;

	max_iter = 60;

	di = d[i];

	eta = 1.0;

	n_iter = 0;

	g = 1.0;

	if ( i < ( n - 1 ) )
    {
    	dip1 = d[ i+1 ];
    
	  // Determining on which half is the root
    
    	lambda = ( di + dip1 ) / 2;
        
    	psi = 0.0;
    	for ( j = 0; j <= i; ++j )
		{
        	psi += v2[j] / ( d[j] - lambda );
    	}
    
    	phi = 0.0;
    	for ( j = n-1; j > i; --j ) // zero is on the left half of the interval
		{
        	phi += v2[j] / ( d[j] - lambda );
    	}
    
    	f = 1.0 + phi + psi;
    
    	if ( f > 0.0 )
        {

	  // Translating the origin to di
        
        	for ( j = 0; j < n; ++j )
			{
            	delta[j] = d[j] - d[i];
        	}
        
        	mu = lambda - di;
        	dip1 -= di; //delta[i+1]
       
          // Finding root
        
	        while ( fabs(g) > eps && fabs( eta / mu ) > 8.0 * eps && max_iter > n_iter )
            {
    	        psi = 0.0;
    	        dpsi = 0.0;
    	        for ( j = 0; j <= i; ++j )
				{
    	            psi += v2[j] / ( delta[j] - mu );
    	            dpsi += v2[j] / ( ( delta[j] - mu ) * ( delta[j] - mu ) );
    	        }
            
    	        phi = 0.0;
    	        dphi = 0.0;
    	        for ( j = n-1; j > i; --j )
				{
    	            phi += v2[j] / ( delta[j] - mu );
    	            dphi += v2[j] / ( ( delta[j] - mu ) * ( delta[j] - mu ) );
    	        }
            
    	        g = 1.0 + phi + psi;

    	        // Solve for zero
            
    	        Deltai = -mu; // delta[i]-mu because delta[i] = d[i] -d[i] = 0
    	        Deltaip1 = dip1 - mu; // delta[i+1] - mu
    	        
    	        a = ( Deltai + Deltaip1 ) * g - Deltai * Deltaip1 * ( dpsi + dphi );
	            b = Deltai * Deltaip1 * g;
    	        c = g - Deltai * dpsi - Deltaip1 * dphi;
            
    	        if ( a > 0.0 )
				{
    	            eta = ( 2.0 * b ) / ( a + sqrt( a * a - 4.0 * b * c ) );
				}
    	        else
				{
    	            eta = ( a - sqrt( a * a - 4.0 * b * c ) ) / ( 2.0 * c );
    	        }
            
    	        mu += eta;
            
    	        if ( 0.0 == mu )
				{
    	            //printf( "mu = %.15lf\n", mu );
    	            mu = eps;
    	        }
            
    	        ++n_iter;
            
    	    }
        
    	    lambda = mu + di; 
        }
	    else // zero is on the right half of the interval
        {

		  // translating the origin to dip1

        	for ( j = 0; j < n; ++j )
			{
    	        delta[j] = d[j] - dip1;
    	    }

    	    mu = lambda - dip1;
    	    di -= dip1; // delta[i]
        
    	  // finding root
        
    	    while ( fabs(g) > eps && fabs( eta / mu ) > 8.0 * eps && max_iter > n_iter ) 
            {
    	        psi = 0.0;
    	        dpsi = 0.0;
    	        for ( j = 0 ; j <= i; ++j )
				{
    	            psi += v2[j] / ( delta[j] - mu );
    	            dpsi += v2[j] / ( ( delta[j] - mu ) * ( delta[j] - mu ) );
    	        }
            
    	        phi = 0.0;
	            dphi = 0.0;
    	        for ( j = n-1; j > i; --j )
				{
    	            phi += v2[j] / ( delta[j] - mu );
    	            dphi += v2[j] / ( ( delta[j] - mu ) * (delta[j] - mu) );
    	        }
            
    	        g = 1.0 + phi + psi;

    	      // Solve for zero
            
    	        Deltai = di - mu; // delta[i] - mu
    	        Deltaip1 = -mu; // delta[i+1] - mu because delta[i+1] = d[i+1] - d[i+1] = 0
            
    	        a = ( Deltai + Deltaip1 ) * g - Deltai * Deltaip1 * ( dpsi + dphi );
    	        b = Deltai * Deltaip1 * g;
    	        c = g - Deltai * dpsi - Deltaip1 * dphi;
            
    	        if ( a > 0.0 )
				{
    	            eta = ( 2.0 * b ) / ( a + sqrt( a * a - 4 * b * c ) );
				}
    	        else
				{
    	            eta = ( a - sqrt( a * a - 4 * b * c ) ) / ( 2 * c );
    	        }
            
    	        mu += eta;
            
    	        if ( 0.0 == mu )
				{
    	            //printf( "mu = %.15lf\n", mu );
    	            mu = -eps;
    	        }
            
    	        ++n_iter;
            
    	    }
        
    	    lambda = mu + dip1;
        
    	}
    
    }
	else // i == n-1
    {
	    dip1 = d[i] + 1.0;
    
	  // Determining on which half is the root

	    lambda = d[i] + 0.5;
    
	    psi = 0.0;
	    for ( j = 0; j <= i-1; ++j )
		{
    	    psi += v2[j] / ( d[j] - lambda );
	    }

	    phi = v2[i] / ( d[i] - lambda );
    
	    f = 1.0 + phi + psi;
    
	    if  ( f > 0.0 ) // zero is on the left half of the interval
        {
	      // Translating the origin to di
        
    	    for ( j = 0; j < n; ++j )
			{
    	        delta[j] = d[j] - di;
			}
        
    	    mu = lambda - di;
    	    dip1  = delta[i-1]; // delta[n-1]
        
    	  // fiding root
        
    	    while ( fabs(g) > eps && fabs( eta / mu ) > 8.0 * eps && max_iter > n_iter )
            {
	            psi = 0.0;
	            dpsi = 0.0;
	            for ( j = 0; j < i; ++j )
				{
	                psi += v2[j] / ( delta[j] - mu );
	                dpsi += v2[j] / ( ( delta[j] - mu ) * ( delta[j] - mu ) );
	            }
            
	            phi = v2[i] / ( delta[i] - mu );
	            dphi = v2[i] / ( ( delta[i] - mu ) * ( delta[i] - mu ) );
	            
	            g = 1.0 + phi + psi;
            
	          // Solve for zero
            
	            Deltai = -mu; // delta[ n-1 ] - mu because delta[ n-1 ] = d[n] - d[n] = 0
	            Deltaip1 = dip1 - mu; // delta[ n-2 ] - mu
            
	            a = ( Deltaip1 + Deltai ) * g - Deltaip1 * Deltai * ( dpsi + dphi );
	            b = Deltaip1 * Deltai * g;
	            c = g - Deltaip1 * dpsi - Deltai * dphi;
            
	            if ( a < 0.0 )
				{
	                eta = ( 2.0 * b ) / ( a - sqrt( a * a - 4.0 * b * c ) );
				}
	            else
				{
	                eta = ( a + sqrt( a * a - 4.0 * b * c ) ) / ( 2.0 * c );
	            }
            
	            mu += eta;
            
	            if ( 0.0 == mu )
				{
	                //printf( "mu = %.15lf\n", mu );
	                mu = eps;
	            }
            
	            ++n_iter;            
	        }
        
	        lambda = mu + di;
        }
	    else // zero is on the right half of the interval
        {
	      // translating the origin to dip1
        
	        for ( j = 0; j < n; ++j )
			{
	            delta[j] = d[j] - dip1;
	        }
        
	        mu = lambda - dip1;
        
	      // finding root
        
	        while ( fabs(g) > eps && fabs( eta / mu ) > 8.0 * eps && max_iter > n_iter )
			{            
	            psi = 0.0;
	            dpsi = 0.0;
	            for ( j = 0; j < i; ++j )
				{
	                psi += v2[j] / ( delta[j] - mu );
	                dpsi += v2[j] / ( ( delta[j] - mu ) * ( delta[j] - mu ) );
	            }
            
	            phi = v2[i] / ( delta[i] - mu );
	            dphi = v2[i] / ( ( delta[i] - mu ) * ( delta[i] - mu ) );
            
	            g = 1.0 + phi + psi;
            
	          // Solve for zero
            
	            Deltai = delta[i] - mu; // delta[ n-1 ] - mu
	            Deltaip1 = delta[ i-1 ] - mu; // delta[ n-2 ] -mu
            
	            a = ( Deltaip1 + Deltai ) * g - Deltaip1 * Deltai * ( dpsi + dphi );
	            b = Deltaip1 * Deltai * g;
	            c = g - Deltaip1 * dpsi - Deltai * dphi;
            
	            if ( a < 0.0 )
				{
	                eta = ( 2.0 * b ) / ( a - sqrt(a * a - 4.0 * b * c ) );
				}
	            else
				{
	                eta = ( a + sqrt(a * a - 4.0 * b * c ) ) / ( 2.0 * c );
	            }
            
	            mu += eta;
            
	            if ( 0.0 == mu )
				{
	                //printf( "mu = %.15lf\n", mu );
	                mu = -eps;
				}
            
	            ++n_iter;
            
	        }
        
	        lambda = mu + dip1;
        
	    }
	}

	for ( j = 0; j < n; ++j )
	{
	    dmlambda[j]  = delta[j] - mu;
	}

	*plambda = lambda;

	//printf( "i = %d\n", i );
	//printf( "n_iter = %d\n", n_iter );
	//printf( "g = %.15e\n", g );
}

