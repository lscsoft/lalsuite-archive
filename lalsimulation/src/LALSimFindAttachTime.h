/*
 * Copyright (C) 2015 S. Babak
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

#ifndef _LALSIMFINDATTACHTIME_H
#define _LALSIMFINDATTACHTIME_H

#include <math.h>

#include <lal/LALGSL.h>

#include <lal/LALDatatypes.h>
//#include <lal/LALSimInspiral.h>
#include <lal/TimeSeries.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/TimeSeries.h>
#include <lal/Units.h>


#include <lal/LALSimIMR.h>
#include <lal/Date.h>
#include <lal/Units.h>
#include <lal/SeqFactories.h>

//#include "LALSimIMREOBNRv2.h"
#include "LALSimIMRSpinEOB.h"


#include <gsl/gsl_linalg.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>




#if defined(__cplusplus)
extern "C" {
#elif 0
} /* so that editors will match preceding brace */
#endif 
 
double  XLALSimLocateOmegaTime(
    REAL8Array *dynamicsHi,
    unsigned int numdynvars,
    unsigned int retLenHi,
    SpinEOBParams   seobParams,
    SpinEOBHCoeffs  seobCoeffs,
    REAL8 m1,
    REAL8 m2,
    REAL8 *radiusVec,
    int *found
);


double XLALSimLocateAmplTime(
    REAL8Vector *timeHi, 
    COMPLEX16Vector *hP22,
    REAL8 *radiusVec,
    int *found
);
    
    
#endif    