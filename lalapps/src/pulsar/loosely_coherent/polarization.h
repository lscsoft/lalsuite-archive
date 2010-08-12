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

#ifndef __POLARIZATION_H__
#define __POLARIZATION_H__

#include <math.h>
#include "grid.h"
#include "global.h"

typedef struct S_POLARIZATION {
	float orientation;
	float  plus_proj;
	float  cross_proj;
	struct S_POLARIZATION *conjugate; /* plus for cross and cross for plus */
	/* For linear polarizations plus_factor is always 1.0 and
	   cross factor is always 0.0 */
	float plus_factor;
	float cross_factor;
	
	SKY_GRID_TYPE *AM_coeffs;
	
	char *name;

	} POLARIZATION;

extern int no_am_response;

inline static void precompute_am_constants(float *e, float ra, float dec)
{
float e2, e3, e4, e5;

e2=cos(M_PI_2-dec);
e3=sin(M_PI_2-dec);
e4=cos(ra);
e5=sin(ra);
/* unit vector */
e[0]=e3*e4;
e[1]=e3*e5;
e[2]=e2;
/* other useful values */
e[3]=e3;
e[4]=e4;
e[5]=e5;
/* these values are needed for regression of plus and cross */
e[6]=e4*e5;
e[7]=e3*e4;
e[8]=e3*e5;

e[9]=e2*e2*e4;
e[10]=e2*e2*e5;
e[11]=e3*e3*e4;
e[12]=e3*e3*e5;
e[13]=e2*e3*e4;
e[14]=e2*e3*e5;

e[15]=e2*e4*e4;
e[16]=e2*e5*e5;
e[17]=e3*e4*e4;
e[18]=e3*e5*e5;
e[19]=e2*e4*e5;

e[20]=e2*e2*e4*e4;
e[21]=e2*e2*e5*e5;
e[22]=e3*e3*e4*e4;
e[23]=e3*e3*e5*e5;

e[24]=e2*e2*e4*e5;
e[25]=e3*e3*e4*e5;
}

static float inline F_plus_coeff(int segment, float *e, float *coeff)
{
float a;
int ii;
if(no_am_response)return 1.0;
a=0.0;
for(ii=0;ii<GRID_FIT_COUNT;ii++)
	a+=coeff[segment*GRID_FIT_COUNT+ii]*e[ii+GRID_FIT_START];
return a;
}

static float inline F_plus(int segment, SKY_GRID *grid, int point, float *coeff)
{
float a;
int ii;
if(no_am_response)return 1.0;
a=0.0;
for(ii=0;ii<GRID_FIT_COUNT;ii++)
	a+=coeff[segment*GRID_FIT_COUNT+ii]*grid->e[ii+GRID_FIT_START][point];
return a;
}

static float inline AM_response(int segment, SKY_GRID *grid, int point, float *coeff)
{
float a;
a=F_plus(segment, grid, point, coeff);
return (a*a);
}

#endif
