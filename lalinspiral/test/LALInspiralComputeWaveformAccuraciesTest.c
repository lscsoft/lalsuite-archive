/*
*  Copyright (C) 2007 Tjonnie Li
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
#include <stdlib.h>
#include <lal/LALStdlib.h>
#include <lal/LALInspiralComputeWaveformAccuracies.h>

NRCSID( LALINSPIRALCOMPUTEWAVEFORMACCURACIESTESTC, "$Id$" );

extern int lalDebugLevel;

int main( int argc, char **argv )
{
  // LAL ERROR REPORTING SYSTEM
  static LALStatus stat;
  lalDebugLevel = 0;

  /*********************************************************************
   *
   *  CREATE PARAMETER STRUCTURE
   * 
   ********************************************************************/
   
  PPNConsistencyParamStruc params;
  params.position.longitude = longitude; /*acos(gsl_ran_flat (p, -1.0, 1.0));*/
  params.position.latitude 	= latitude;/*gsl_ran_flat (p, 0, LAL_TWOPI);*/
  params.position.system 		= COORDINATESYSTEM_EQUATORIAL;
  params.epoch.gpsSeconds		= gpsSeconds;
  params.epoch.gpsNanoSeconds = gpsNanoSeconds;
  params.psi								= psi; /*acos(gsl_ran_flat(p,-1.0,1.0));*/
  params.deltaT 						= testACS.deltaT_ts;
  params.mTot_real8 				= m1 + m2;
  params.eta_real8  				= m1 * m2 / pow(m1 + m2, 2.);
  params.inc 								= inc;
  params.tc  								= tc;
  params.phi 								= phic;
  params.d   								= dist;
  params.fStartIn 					= fStartIn;
  params.fStopIn  					= fStopIn;
  params.ampOrder 					= ampOrder;
	
  /*********************************************************************
   *
   *  Create PSD
   * 
   ********************************************************************/    
  
  /* PSD */
  REAL8FrequencySeries *psd;
  
  // CREATE FREQUENCY SERIES FOR PSD COMPUTING
  psd = XLALCreateREAL8FrequencySeries("PSD_LIGO",  &(params.epoch), 0.0, testACS.deltaF_fs, &lalStrainUnit, testACS.N_fs);
  
  // Computing LIGO PSD
  LALNoiseSpectralDensity( &status, psd->data, noisemodel, psd->deltaF );
  
  /*********************************************************************
   *
   *  
   * 
   ********************************************************************/
	
	// CLOSING PROGRAM
  REPORTSTATUS( &stat );
  return stat.statusCode;
}
/* </lalVerbatim> */
