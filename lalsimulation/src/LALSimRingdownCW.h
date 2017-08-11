#ifndef _LALSIM_RINGDOWN_CW_W
#define _LALSIM_RINGDOWN_CW_H
/* ************************************************************  */
/*
 * Copyright (C) 2016 Lionel London
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

 #if defined(__cplusplus)
 extern "C" {
 #elif 0
 } /* so that editors will match preceding brace */
 #endif

/* Include the desired Libs */
#include <stdbool.h>
#include <math.h>
#include <complex.h>
/* LAL specific libs  */

#include <lal/LALDatatypes.h>
#include <lal/LALStdlib.h>
#include <lal/LALSimInspiral.h>
#include <lal/LALSimIMR.h>
#include <lal/LALConfig.h>
#include <lal/LALConstants.h>
#include <lal/Date.h>
#include <lal/FrequencySeries.h>
#include <lal/Units.h>

/*
* Basic NOTE(s): Prototypes for LALSimRingdownCW.c
*/

/*
* Domain mapping for dimnesionless BH spin
*/
static double KAPPA( double jf, int l, int m );

/*
* Dimensionless QNM Frequencies: Note that name encodes date of writing
*/
static complex double CW07102016( double kappa,  /* Domain mapping for  remnant BH's spin (Dimensionless) */
                          int l,        /* Polar eigenvalue */
                          int input_m,  /* Azimuthal eigenvalue*/
                          int n );      /* Overtone Number*/

/* ************************************************************  */
#endif	/* of #ifndef _LALSIM_RINGDOWN_MMRDNS_H */
