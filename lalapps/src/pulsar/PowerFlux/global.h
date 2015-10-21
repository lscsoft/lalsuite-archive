/*
*  Copyright (C) 2007 Vladimir Dergachev
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with with program; see the file COPYING. If not, write to the
*  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
*  MA  02111-1307  USA
*/

#ifndef __GLOBAL_H__
#define __GLOBAL_H__

/* These affect compilation of entire program. 
   Their are #define's not variables for efficiency */

#define WEIGHTED_SUM



typedef long long INT64;
typedef float SUM_TYPE;
typedef short COUNT_TYPE;

#define M_1_3 0.333333333333333333333333333333
#define M_1_6 0.166666666666666666666666666666666

#define MEMUSAGE	((long)sbrk(0))

#define TRACE(a)	{fprintf(stderr,"TRACE(__FUNCTION__):" a); \
			fprintf(stderr,"\n");}

void *do_alloc(long a, long b);

#define aligned_alloca(a) ((void *)(((unsigned long)(alloca(a+63))+63) & ~63))

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

#define CHECKPOINT   {\
	fprintf(stderr, "CHECKPOINT function %s line %d file %s\n", __FUNCTION__, __LINE__, __FILE__); \
	}

#endif
