#ifndef __GLOBAL_H__
#define __GLOBAL_H__


typedef long long INT64;
#define TRACE(a)	{fprintf(stderr,"TRACE(" __FUNCTION__ "):" a); \
			fprintf(stderr,"\n");}

void *do_alloc(long a, long b);

#define TESTSTATUS( status ) \
  { if ( (status)->statusCode ) { \
  	fprintf(stderr,"** LAL status encountered in file \"%s\" function \"%s\" line %d\n", __FILE__, __FUNCTION__, __LINE__);\
  	REPORTSTATUS((status)); exit(-1); \
	}  }

#include "grid.h"

static float inline AM_response(int segment, SKY_GRID *grid, int point, float *coeff)
{
float a;
int ii;
a=0.0;
for(ii=0;ii<GRID_FIT_COUNT;ii++)
	a+=coeff[segment*GRID_FIT_COUNT+ii]*grid->e[ii+GRID_FIT_START][point];
#if 0 /* just for testing */
return a;
#else
return (a*a);
#endif
}

static float inline sqr_f(float a)
{
return (a*a);
}

#endif
