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
 
#include  <lal/LALSimInspiralAlternativeTheories.h>
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
														REAL8 z)							  /* redshift to the source */
{
    REAL8 lambdaG=pow(10.0,wavelength);
	REAL8 coefficient=-LAL_PI*LAL_PI*distance*chirpmass/(lambdaG*lambdaG*(1.0+z));
	if (!XLALSimInspiralTestGRParamExists(*parameter,"betaPPE")) XLALSimInspiralAddTestGRParam(parameter,"betaPPE",coefficient);
	else XLALSimInspiralSetTestGRParam(*parameter,"betaPPE",coefficient);
	if (!XLALSimInspiralTestGRParamExists(*parameter,"bPPE")) XLALSimInspiralAddTestGRParam(parameter,"bPPE",-1.0);
	else XLALSimInspiralSetTestGRParam(*parameter,"bPPE",-1.0);
}

/*************************************************************
 * Function to compute the value of the PPE parameters 
 * for the Brans-Dicke theory (http://arxiv.org/pdf/gr-qc/0504017v2.pdf )
 *************************************************************/
 
void XLALSimInspiralComputePPEparametersForBransDicke(LALSimInspiralTestGRParam **parameter, /* the test parameter structure output */ 
														REAL8 alpha1,					    /* scalar charge for body 1 */ 
														REAL8 alpha2,						/* scalar charge for body 2 */
														REAL8 omegaBD,					    /* BD parameter */
														REAL8 eta)							/* symmetric mass ratio */
														
{
    REAL8 coefficient=5.0*pow(eta,2.0/5.0)*(alpha1-alpha2)*(alpha1-alpha2)/(43008.0*omegaBD);
	if (!XLALSimInspiralTestGRParamExists(*parameter,"betaPPE")) XLALSimInspiralAddTestGRParam(parameter,"betaPPE",coefficient);
	else XLALSimInspiralSetTestGRParam(*parameter,"betaPPE",coefficient);
	if (!XLALSimInspiralTestGRParamExists(*parameter,"bPPE")) XLALSimInspiralAddTestGRParam(parameter,"bPPE",-7.0/3.0);
	else XLALSimInspiralSetTestGRParam(*parameter,"bPPE",-7.0/3.0);
}

/* http://arxiv.org/abs/0912.2724*/
void XLALSimInspiralComputePPEparametersForVariableG(LALSimInspiralTestGRParam **parameter, /* the test parameter structure output */ 
														REAL8 Gdot,					    /* time variation of G */
														REAL8 chirpmass,				/* chirp mass */
														REAL8 z)						/* redshift */
														
{
	/* order of magnitude for Gdot~8×10^−7 * G *yr^-1 */
	REAL8 Mz = chirpmass*(1.0+z);
	if (!XLALSimInspiralTestGRParamExists(*parameter,"betaPPE")) XLALSimInspiralAddTestGRParam(parameter,"betaPPE",(-25.0/65536.0)*Gdot*Mz);
	else XLALSimInspiralSetTestGRParam(*parameter,"betaPPE",(-25.0/65536.0)*Gdot*Mz);
	if (!XLALSimInspiralTestGRParamExists(*parameter,"bPPE")) XLALSimInspiralAddTestGRParam(parameter,"bPPE",-13.0/3.0);
	else XLALSimInspiralSetTestGRParam(*parameter,"bPPE",-13.0/3.0);
	if (!XLALSimInspiralTestGRParamExists(*parameter,"alphaPPE")) XLALSimInspiralAddTestGRParam(parameter,"alphaPPE",-(5.0/512.0)*Gdot*Mz);
	else XLALSimInspiralSetTestGRParam(*parameter,"alphaPPE",-(5.0/512.0)*Gdot*Mz);
	if (!XLALSimInspiralTestGRParamExists(*parameter,"aPPE")) XLALSimInspiralAddTestGRParam(parameter,"aPPE",-8.0/3.0);
	else XLALSimInspiralSetTestGRParam(*parameter,"aPPE",-8.0/3.0);
}

/* http://arxiv.org/pdf/1101.2921v2.pdf */

void XLALSimInspiralComputePPEparametersForQuadraticGravity(LALSimInspiralTestGRParam **parameter, /* the test parameter structure output */ 
														REAL8 eta,					/* chirp mass */
														REAL8 zeta)					/* coupling constant */
														
{	
	/* zeta constrained to be < 10^7 */
	if (!XLALSimInspiralTestGRParamExists(*parameter,"betaPPE")) XLALSimInspiralAddTestGRParam(parameter,"betaPPE",(50.0/3.0)*zeta);
	else XLALSimInspiralSetTestGRParam(*parameter,"betaPPE",(50.0/3.0)*zeta);
	if (!XLALSimInspiralTestGRParamExists(*parameter,"bPPE")) XLALSimInspiralAddTestGRParam(parameter,"bPPE",4.0/3.0);
	else XLALSimInspiralSetTestGRParam(*parameter,"bPPE",4.0/3.0);
	if (!XLALSimInspiralTestGRParamExists(*parameter,"alphaPPE")) XLALSimInspiralAddTestGRParam(parameter,"alphaPPE",(5.0/6.0)*zeta*pow(eta,-4.0/5.0));
	else XLALSimInspiralSetTestGRParam(*parameter,"alphaPPE",(5.0/6.0)*zeta*pow(eta,-4.0/5.0));
	if (!XLALSimInspiralTestGRParamExists(*parameter,"aPPE")) XLALSimInspiralAddTestGRParam(parameter,"aPPE",4.0/3.0);
	else XLALSimInspiralSetTestGRParam(*parameter,"aPPE",4.0/3.0);
}