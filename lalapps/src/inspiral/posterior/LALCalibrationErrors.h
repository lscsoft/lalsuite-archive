/*
*  Copyright (C) 2011 Walter Del Pozzo, Salvatore Vitale, Tjonnie Li
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

/* <lalVerbatim> */
#ifndef _LALCALIBRATIONERRORS_H  /* Double-include protection. */
#define _LALCALIBRATIONERRORS_H
#include <lal/LALInspiral.h>
#include "LALInspiralMCMC.h"
#include "LALInspiralMCMCUser.h"
#include <lal/LALDatatypes.h>
#include <math.h>
#include <lal/FrequencySeries.h>

// GSL PACKAGES
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#ifdef  __cplusplus   /* C++ protection. */
extern "C" {
#endif


NRCSID( LALCALIBRATIONERRORSH, "$Id$" );

#define LALCALIBRATIONERRORSH_ENULL 1
#define LALCALIBRATIONERRORSH_EDIV0 2
#define LALCALIBRATIONERRORSH_MSGENULL "Null pointer"
#define LALCALIBRATIONERRORSH_MSGEDIV0 "Division by zero"



void SampleCalibrationErrorsAmplitude(REAL8 *logF, INT4 length, INT4 IFO, INT4 seed, REAL8 *errors);
void SampleCalibrationErrorsPhase(REAL8 *logF, INT4 length, INT4 IFO, INT4 seed, REAL8 *errors);

void FitNoiseRealisation(
														LALStatus						*status,
														INT4								R,		/* Polynomial fitter order */
														INT4								N,		/* total amount of samples (NB: needs to be odd) */
														REAL8								*y,		/* Noise sample to be fitted */
														REAL8								dlogf, /* stepsize in logf */
														REAL8								*D			/* Fitting coefficients */
														);
														
void InvertMatrixSVD (	gsl_matrix *A,gsl_matrix *InvA,	int	N);		
											
REAL8 ConvertCoefficientsToFunction(REAL8 *coeff, REAL8 f, REAL8 cen);																				
/*
REAL8 Amp_H1(REAL8 f);
REAL8 Amp_L1(REAL8 f);
REAL8 Amp_V1(REAL8 f);
REAL8 Ph_H1(REAL8 f);
REAL8 Ph_L1(REAL8 f);
REAL8 Ph_V1(REAL8 f);
*/
void CreateErrorStreams(LALMCMCInput *inputMCMC,CHAR *IFOname, int i,int seed);
void ApplyCalibrationErrorsToData(LALMCMCInput inputMCMC,COMPLEX16FrequencySeries *CalibNoise,CHAR *IFOname,int i, int seed);
void ApplyCalibrationErrorsToWaveform(COMPLEX16FrequencySeries *injF,COMPLEX16FrequencySeries *CalibInj,CHAR *IFOname, int i, int seed);
void CalibPolar(COMPLEX16FrequencySeries *injF, COMPLEX16FrequencySeries *calibInjF, CHAR *IFOname,int seed);
REAL8 GenerateFrequencySamples(REAL8 f_min, REAL8 f_max, UINT4 Length);
#ifdef  __cplusplus
}
#endif  /* C++ protection. */
#endif  /* Double-include protection. */
/* </lalVerbatim> */
