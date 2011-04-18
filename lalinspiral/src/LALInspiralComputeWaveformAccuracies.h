/*
*  Copyright (C) 2011 Tjonnie Li
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
#ifndef _LALINSPIRALCOMPUTEWAVEFORMACCURACIES_H  /* Double-include protection. */
#define _LALINSPIRALCOMPUTEWAVEFORMACCURACIES_H 

#include <lal/LALDatatypes.h>
#include <lal/LALStdio.h>
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
#include <lal/LALDetectors.h>
#include <lal/Date.h>
#include <lal/DetResponse.h>
#include <lal/LALError.h>
#include <lal/GenerateInspiral.h>
#include <lal/TimeDelay.h>
#include <lal/GeneratePPNAmpCorConsistency.h>

#ifdef  __cplusplus   /* C++ protection. */
extern "C" {
#endif


NRCSID( LALINSPIRALCOMPUTEWAVEFORMACCURACIESH, "$Id$" );

#define LALINSPIRALCOMPUTEWAVEFORMACCURACIESH_ENULL 1
#define LALINSPIRALCOMPUTEWAVEFORMACCURACIESH_ENULL_MSG "Null pointer"

void LALInspiralComputeDhInDh (
														LALStatus															 *status,
														PPNConsistencyParamStruc               *exactParams,
														PPNConsistencyParamStruc               *modelParams,
														REAL8																	 fStart,
														REAL8																	 fStop,
														REAL8FrequencySeries									 *psd);

void LALInspiralComputeSNR (
														LALStatus															 *status,
														PPNConsistencyParamStruc               *PPNparams,
														REAL8																	 fStart,
														REAL8																	 fStop,
														REAL8FrequencySeries									 *psd,
														REAL8																	 *SNR);
														
void LALInspiralCombinePlusCross(
																 LALStatus									*status,			// Status pointer
																 CoherentGW									*GravWav,			// Input waveform, data[2i] = hp, data[2i+1] = hc
																 PPNConsistencyParamStruc		*params,				// input parameters
																 REAL4TimeSeries						*ht);					// combined h(t) from hp and hc		
																 
void XLALDestroyCoherentGW (CoherentGW  *GravWav   )	;															 											

#ifdef  __cplusplus
}
#endif  /* C++ protection. */
#endif  /* Double-include protection. */
/* </lalVerbatim> */
