#ifndef __GLOBAL_H__
#define __GLOBAL_H__

/* These affect compilation of entire program. 
   Their are #define's not variables for efficiency */

#define WEIGHTED_SUM



typedef long long INT64;
typedef float SUM_TYPE;
typedef short COUNT_TYPE;



#define TRACE(a)	{fprintf(stderr,"TRACE(" __FUNCTION__ "):" a); \
			fprintf(stderr,"\n");}

void *do_alloc(long a, long b);

#define TESTSTATUS( status ) \
  { if ( (status)->statusCode ) { \
  	fprintf(stderr,"** LAL status encountered in file \"%s\" function \"%s\" line %d\n", __FILE__, __FUNCTION__, __LINE__);\
  	REPORTSTATUS((status)); exit(-1); \
	}  }

#include "grid.h"

static float inline sqr_f(float a)
{
return (a*a);
}

#endif
