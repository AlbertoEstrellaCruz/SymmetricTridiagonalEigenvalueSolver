#include "merge_v1.h"

void merge_d( double *s, int *p, double *l, double *r,  int nl, int nr, int n )
{
	int i;
	int j;
	int k;

	i = 0;
	j = 0;
	k = 0;

	for ( k = 0; k < n; ++k )
    {
    	if ( i == nl )
        {
    	    s[k] = r[j];
    	    p[k] = nl + j;
    	    ++j;
        }
    	else if ( j == nr )
        {
    	    s[k] = l[i];
    	    p[k] = i;
    	    ++i;
        }
    	else if ( r[j] < l[i] )
        {
    	    s[k] = r[j];
    	    p[k] = nl + j;
    	    ++j;
        }
    	else
        {
    	    s[k] = l[i];
    	    p[k] = i;
    	    ++i;        
    	}
	}
}

