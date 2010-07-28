/*
*  Copyright (C) 2009 Tjonnie Li, Chris Van Den Broeck
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

/*------------------------ Included Libraries ------------------------*/
#ifndef _LALINSPIRALCOMPUTEFISHER_H /* Double-include protection*/
#define _LALINSPIRALCOMPUTEFISHER_H

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <lal/LALStdio.h>
#include <lal/LALDatatypes.h>
#include <lal/LALConstants.h>
#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>
#include <lal/Sequence.h>
#include <lal/LALInspiralBank.h>
#include <lal/LALStatusMacros.h>
#include <lal/LALMalloc.h>
#include <lal/LALNoiseModels.h>
#include <lal/MatrixUtils.h>
#include <lal/FindChirpPTF.h>
#include <lal/SimulateCoherentGW.h>
#include <lal/GeneratePPNInspiral.h>
#include <lal/LALInspiral.h>
#include <lal/PrintFTSeries.h>
#include <lal/TimeFreqFFT.h>
#include <lal/TimeSeries.h>
#include <lal/FrequencySeries.h>
#include <lal/Units.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#ifdef __cplusplus
extern "C" {
#endif

NRCSID( LALINSPIRALCOMPUTEFISHERH, "$Id: LALInspiralComputeFisher.h, v1.0 2 February 2010 $");

/*----------------------- Define error messages ----------------------*/

#define LALINSPIRALCOMPUTEFISHERC_EINPUT 1
#define LALINSPIRALCOMPUTEFISHERC_EINPUT_MSG "Invalid Input"

#define LALINSPIRALCOMPUTEFISHERC_EMBAD 2
#define LALINSPIRALCOMPUTEFISHERC_EMBAD_MSG "Invalid Mass Parameter"

#define LALINSPIRALCOMPUTEFISHERC_EMEM 3
#define LALINSPIRALCOMPUTEFISHERC_EMEM_MSG "Out of Memory"

#define LALINSPIRALCOMPUTEFISHERC_ETSCREATE 4
#define LALINSPIRALCOMPUTEFISHERC_ETSCREATE_MSG "Error Creating TimeSeries"

#define LALINSPIRALCOMPUTEFISHERC_EPAR 5
#define LALINSPIRALCOMPUTEFISHERC_EPAR_MSG "Error Creating PPNparams"

#define LALINSPIRALCOMPUTEFISHERC_EGW 6
#define LALINSPIRALCOMPUTEFISHERC_EGW_MSG "Error Creating CoherentGW"

#define LALINSPIRALCOMPUTEFISHERC_EINPUTPARAMS 7
#define LALINSPIRALCOMPUTEFISHERC_EINPUTPARAMS_MSG "Invalid input PPNConsistencyParamStruc"

#define LALINSPIRALCOMPUTEFISHERC_EDUPLPARAMS 8
#define LALINSPIRALCOMPUTEFISHERC_EDUPLPARAMS_MSG "Invalid Duplicate PPNConsistencyParamStruc"

/***********************************************************************
 *
 *  Declaring structures
 *
 **********************************************************************/

/* ALGORITHM CONTROL STRUCTURE */



typedef struct tagPPNConsistencyParamStruc{
  /* Passed parameters. */
  SkyPosition position; /* location of source on sky */
  REAL4 psi;            /* polarization angle (radians) */
  LIGOTimeGPS epoch;    /* start time of output time series */

  /* Input parameters. */
  REAL8 mTot_real8; /* total system mass (Msun) */
  REAL8 eta_real8;  /* mass ratio */
  REAL8 delta;      /* sqrt(1-4eta) */
  REAL4 mTot;       /* total system mass (Msun) */
  REAL4 eta;        /* mass ratio */
  REAL4 d;          /* distance (metres) */
  REAL4 inc;        /* inclination angle (radians) */
  REAL4 cosI;				/* cosine of inclination angle */
  REAL4 sinI;				/* sine of inclination angle */
  REAL4 phi;        /* coalescence phase (radians) */
  REAL8 deltaT;     /* requested sampling interval (s) */
  REAL4 fStartIn;   /* requested start frequency (Hz) */
  REAL4 fStopIn;    /* requested stop frequency (Hz) */
  UINT4 lengthIn;   /* maximum length of waveform */
  REAL4Vector *ppn; /* post-Newtonian selection parameters */
  INT4 ampOrder;    /* PN amplitude selection 0-5 */
  /* PN phasing coefficients for use in AmpCorConsistency */
  REAL4 phi0, phi2, phi3, phi4, phi5, phi5l, phi6, phi6l, phi7;
  

  /* Output parameters. */
  REAL8 tc;         /* time to coalescence from start of waveform */
  REAL4 dfdt;       /* maximum value of df*dt over any timestep */
  REAL4 fStart;     /* actual start frequency (Hz) */
  REAL4 fStop;      /* actual stop frequency (Hz) */
  UINT4 length;     /* length of signal generated */
  INT4 termCode;    /* termination code */
  const CHAR *termDescription; /* description of termination code */
}PPNConsistencyParamStruc;

typedef struct tagFisherACS
{
	/* GENERAL CONTROL PARAMETERS */
	INT4 verbose_switch;
	INT4 printall;
	char folder[128];
	LIGOTimeGPS epoch;
	INT4 N_inputparams;
	REAL8 dynRange;
	
	/* TIME/FREQUENCY SERIES */
	INT4 N_ts;
	INT4 N_fs;
	REAL8 deltaT_ts;
	REAL8 deltaF_fs;
	
	/* DIFFERENTIATION PARAMETERS */
	INT4 xsamples_tot;
	INT4 linReg_countMax;
	INT4 linReg_countMin;
	INT4 linReg_loopMax;
	REAL8 deltaPhi_max;
	REAL8 reldelta;
	REAL8 absdelta_adjustUp;
	REAL8 absdelta_adjustDown;
	INT4 xsamples_adjustUp;
	INT4 xsamples_adjustDown;
	
	/* DETECTOR FREQUENCY WINDOW */
	REAL8 fstart;
	REAL8 fstop;
	
	/* OUTPUT FILE */
  FILE *psd_out;
  FILE *derivatives_out; 
  FILE *fourierderivs_out;  
  FILE *fisher_out;
  FILE *cov_out; 
  FILE *hx_out;
  
}FisherACS;

/***********************************************************************
 *
 *  Declaring function
 *
 **********************************************************************/

void
LALGeneratePPNAmpCorConsistency( LALStatus     *,
			CoherentGW    *output,
			PPNConsistencyParamStruc *params );

void LALInspiralComputeFisherMatrix (
    LALStatus                              *status,  
    REAL8FrequencySeries                   *psd,
    PPNConsistencyParamStruc               *params,
    FisherACS															 ACS	 );
    
void LALInspiralComputeFisherComponents (
    LALStatus															 *status,
    InspiralMetric                         *Comp,
    REAL4TimeSeries                        *hderiv1,
    REAL4TimeSeries                        *hderiv2,
    REAL8FrequencySeries                   *psd,
    UINT4                                  compid,
    FisherACS															 ACS);
    
void LALInspiralComputeDerivatives (
    LALStatus                          *status,
    REAL4TimeSeries                    *hderiv,
    PPNConsistencyParamStruc                      *PPNparams,
    INT4                              paramid,
    FisherACS														ACS);
    
void LALInspiralComputeSNR (
		LALStatus															 *status,
		REAL8FrequencySeries                   *psd,
    PPNConsistencyParamStruc                          *PPNparams,
    FisherACS															 ACS); 

void LALInspiralComputeDerivatives_linReg(
		LALStatus																*status,			// Status pointer
		PPNConsistencyParamStruc														*InputParams,	// Input parameters
		REAL8																		deltax,				// total width of differentiation
		INT4																		paramid,			// parameter choice
		INT4																		icur,				// start waveform element
		INT4																		iend,					// end waveform element (length)
		REAL4TimeSeries													*hderiv,			// derivative calculated with linReg algorithm
    INT4                                    xsamples,
		FisherACS																ACS);    
    
void LALInspiralComputeDerivatives_5point(
		LALStatus																*status,			// Status pointer
		PPNConsistencyParamStruc														*InputParams,	// Input parameters
		INT4																		paramid,			// parameter choice
		REAL8																		deltax,				// differentiation width
		INT4																		istart,				// starting time
		INT4																		iend,					// end time
		REAL4TimeSeries													*hderiv);			// TimeSerie holding derivative
    
void XLALDestroyCoherentGW (
    CoherentGW  *GravWav );
    
void XLALCopyPPNConsistencyParamStruc (
      PPNConsistencyParamStruc             *inputParams,
      PPNConsistencyParamStruc             *duplicateParams );
    
void LALInspiralInvertMatrix (
    LALStatus															*status,
    REAL4                                 *Input,
    REAL4                                 *Inverse,
    FisherACS															ACS);
    
void LALInspiralCombinePlusCross(
		LALStatus															*status,			// Status pointer
		CoherentGW														GravWav,			// Input waveform, data[2i] = hp, data[2i+1] = hc
		REAL4TimeSeries												*ht);					// combined h(t) from hp and hc    

#ifdef __cplusplus
}
#endif /* C++ protection*/
#endif /* Double-include protection */

