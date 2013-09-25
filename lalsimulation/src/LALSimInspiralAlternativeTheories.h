/* Copyright (C) 2013 Walter Del Pozzo
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
 
#ifndef _LALSIMINSPIRALALTERNATIVETHEORIES_H
#define _LALSIMINSPIRALALTERNATIVETHEORIES_H

#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h> 
#include <string.h>
#include <math.h>
#include <lal/LALConstants.h>
#include <lal/LALMalloc.h>
#include <lal/LALError.h>
#include <lal/LALSimInspiralTestGRParams.h>

/************************************************************* 
 * The following set of functions compute the PPE parameters given 
 * a certain class of alternative gravity. The list is given in 
 * http://arxiv.org/pdf/1105.2088v2.pdf
 *************************************************************/
/*************************************************************
 * Function to compute the value of betaPPE from distance
 * Compton wavelength of the graviton, redshift (if available, 
 * default=0) and chirp mass. All quantities should be in meters
 * (http://arxiv.org/pdf/gr-qc/0504017v2.pdf )
 *************************************************************/
 
void XLALSimInspiralComputePPEparametersForMassiveGravity(LALSimInspiralTestGRParam **parameter, /* the test parameter structure output */ 
														REAL8 wavelength,					  /* Compton wavelength of the graviton in meters */ 
														REAL8 distance,						  /* distance to the source in meters, use the luminosity distance for now, but there should be an additional 1+z */
														REAL8 chirpmass,					  /* chirp mass of the system in meters */
														REAL8 z);						  /* redshift to the source */

/*************************************************************
 * Function to compute the value of the PPE parameters 
 * for the Brans-Dicke theory (http://arxiv.org/pdf/gr-qc/0504017v2.pdf )
 *************************************************************/
 
void XLALSimInspiralComputePPEparametersForBransDicke(LALSimInspiralTestGRParam **parameter, /* the test parameter structure output */ 
														REAL8 alpha1,					    /* scalar charge for body 1 */ 
														REAL8 alpha2,						/* scalar charge for body 2 */
														REAL8 omegaBD,					    /* BD parameter */
														REAL8 eta);						/* symmetric mass ratio */
														
/* http://arxiv.org/abs/0912.2724*/
void XLALSimInspiralComputePPEparametersForVariableG(LALSimInspiralTestGRParam **parameter, /* the test parameter structure output */ 
														REAL8 Gdot,					    /* time variation of G */
														REAL8 chirpmass,				/* chirp mass */
														REAL8 z);						/* redshift */

/* http://arxiv.org/pdf/1101.2921v2.pdf */

void XLALSimInspiralComputePPEparametersForQuadraticGravity(LALSimInspiralTestGRParam **parameter, /* the test parameter structure output */ 
														REAL8 eta,					/* chirp mass */
														REAL8 zeta);					/* coupling constant */
														
#endif /* _LALSIMINSPIRALALTERNATIVETHEORIES_H */