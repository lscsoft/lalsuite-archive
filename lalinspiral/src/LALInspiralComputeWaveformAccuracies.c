/*
*  Copyright (C) 2007 Jolien Creighton
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
#include <lal/LALStdlib.h>
#include <lal/LALInspiralComputeWaveformAccuracies.h>

NRCSID( LALINSPIRALCOMPUTEWAVEFORMACCURACIESC, "$Id$" );


/* <lalVerbatim file="LALInspiralComputeSNRCP">  */
void LALInspiralComputeDhInDh (
														LALStatus															 *status,
														PPNConsistencyParamStruc               *exactParams,
														PPNConsistencyParamStruc               *modelParams,
														REAL8																	 fStart,
														REAL8																	 fStop,
														REAL8FrequencySeries									 *psd)
/* </lalVerbatim> */
{
	
  /*********************************************************************
   *
   *  Error handling
   * 
   ********************************************************************/ 
  INITSTATUS( status, "LALInspiralComputeSNR", LALINSPIRALCOMPUTEWAVEFORMACCURACIESC);
  ATTATCHSTATUSPTR(status);
	
  /*********************************************************************
   *
   *  Temporary Variables
   * 
   ********************************************************************/ 
  
	//Counters
	INT4 i,j,k;
	
	//TIME/FREQUENCY SERIES SIZE
	INT4 N_ts;
	INT4 N_fs;
	REAL8 deltaT = 1.0/16384.0;
	REAL8 deltaF = 0.0;
	
	// DYNAMICAL RANGE
	REAL8 dynRange = 1E+27;
	REAL8 invpsd = 0.0;
	REAL8 dhindh = 0.0;
	
  // PARAMETERS AND WAVEFORMS
  CoherentGW    							exactWave;		// Holding waveform, loaded onto ht
  CoherentGW    							modelWave;		// Holding waveform, loaded onto ht
	
	// TIME/FREQUENCY SERIES
	REAL4TimeSeries 				*exact_ht;				// Initial Waveform
	REAL4TimeSeries 				*model_ht;				// Initial Waveform
	REAL4TimeSeries 				*dht;				// Initial Waveform
	REAL4TimeSeries 				*invFT;			// Store inverse Fourier Transform
  COMPLEX8FrequencySeries	*dHf;				// FT of ht
  COMPLEX8FrequencySeries	*inprod_integrand;				// FT of ht
	
	// FOURIER TRANSFORMS
	RealFFTPlan *fwdRealPlan = NULL;	// FT plan
  RealFFTPlan *revRealPlan = NULL;	// RFT plan
  
  /*********************************************************************
   *
   *  COMPUTE COHERENTGW FROM INPUT PARAMETERS
   * 
   ********************************************************************/  
  
  /* CLEARING COHERENTGW */
	memset(&exactWave, 0, sizeof(CoherentGW) );  
	memset(&modelWave, 0, sizeof(CoherentGW) );  
	
  // GENERATE WAVEFORM AND COMBINE HP AND HC
  LALGeneratePPNAmpCorConsistency(status->statusPtr, &exactWave, exactParams);
  LALGeneratePPNAmpCorConsistency(status->statusPtr, &modelWave, modelParams);
  
  // SET SIZE TIME/FREQUENCY SERIES
  N_ts = modelWave.f->data->length;
  N_fs = N_ts/2+1;
  
  // SET STEPSIZE TIME/FREQUENCY SERIES
  deltaF = 1.0/((REAL8) N_ts*deltaT);
  
  /*********************************************************************
   *
   *  CREATING TIME/FREQUENCY SERIES
   * 
   ********************************************************************/
  
	exact_ht 				= XLALCreateREAL4TimeSeries("exact_ht", &(exactParams->epoch), 0.0, deltaT, &lalStrainUnit, N_ts);
	model_ht 				= XLALCreateREAL4TimeSeries("exact_ht", &(modelParams->epoch), 0.0, deltaT, &lalStrainUnit, N_ts);
	dht 				= XLALCreateREAL4TimeSeries("dht", &(modelParams->epoch), 0.0, deltaT, &lalStrainUnit, N_ts);
	invFT			=	XLALCreateREAL4TimeSeries("invFT", &(exactParams->epoch), 0.0, deltaT,&lalStrainUnit, N_ts);
	dHf				= XLALCreateCOMPLEX8FrequencySeries("dHf", &(exactParams->epoch), 0.0, deltaF, &lalStrainUnit, N_fs);
	inprod_integrand	= XLALCreateCOMPLEX8FrequencySeries("inprod_integrand", &(exactParams->epoch), 0.0, deltaF, &lalStrainUnit, N_fs);

	// CLEAR TIME SERIES	
	for(j=0;j<N_ts;j++){
		exact_ht->data->data[j] 		= 0.0; 
		model_ht->data->data[j] 		= 0.0; 
		dht->data->data[j] 		= 0.0; 
		invFT->data->data[j]	= 0.0;
	}  
	
	for(j=0; j<N_fs; j++){
		dHf->data->data[j].re				= 0.0; 
		dHf->data->data[j].im				= 0.0; 
		inprod_integrand->data->data[j].re	= 0.0;
		inprod_integrand->data->data[j].im	= 0.0;
	}
  
  // FFT
  fwdRealPlan = XLALCreateForwardREAL4FFTPlan(N_ts,0);
  revRealPlan = XLALCreateReverseREAL4FFTPlan(N_ts,0);

	/*********************************************************************
   *
   *  Compute dh
   * 
   ********************************************************************/
	
	LALInspiralCombinePlusCross(status->statusPtr, &exactWave, exactParams, exact_ht);  
	LALInspiralCombinePlusCross(status->statusPtr, &modelWave, modelParams, model_ht);
	
	for(i=0;i<N_ts;i++){
		dht->data->data[i] = model_ht->data->data[i] - exact_ht->data->data[i];
	}  
	
	/*********************************************************************
   *
   *  Compute INNER PRODUCT
   * 
   ********************************************************************/
  
  // PERFORM FFT ON HT
  LALTimeFreqRealFFT( status->statusPtr, dHf, dht, fwdRealPlan );
  
  // Scale using dynamical range
  for (i = 0; i < N_fs; ++i)
  {
    dHf->data->data[i].re *= dynRange; 
    dHf->data->data[i].im *= dynRange; 
  } 
  
  // Calculate Optimal SNR
	for (k = 0; k < N_fs; ++k)
  {
    if (psd->data->data[k] == 0.0 || ((REAL8) k)*deltaF < fStart || ((REAL8) k)*deltaF > fStop )
      invpsd = 0.0;
    else
      invpsd = 1.0 / ( dynRange * dynRange * psd->data->data[k]*2.0);
		
    inprod_integrand->data->data[k].re = (dHf->data->data[k].re * dHf->data->data[k].re + dHf->data->data[k].im * dHf->data->data[k].im) * invpsd;
    inprod_integrand->data->data[k].im = 0.0;
  }
  
  // Reserve FT, FIRST ELEMENT IS THE INTEGRAL!!
	LALFreqTimeRealFFT( status->statusPtr, invFT, inprod_integrand, revRealPlan );
	
	// Compute Total SNR
	dhindh = sqrt(4.0*invFT->data->data[0]); 
  
  /*********************************************************************
   *
   *  FREE MEMORY
   * 
   ********************************************************************/
   
  XLALDestroyREAL4TimeSeries(exact_ht);					exact_ht = NULL;
  XLALDestroyREAL4TimeSeries(model_ht);					model_ht = NULL;
  XLALDestroyREAL4TimeSeries(dht);					dht = NULL;
  XLALDestroyREAL4TimeSeries(invFT);					invFT = NULL;
	XLALDestroyCOMPLEX8FrequencySeries(dHf);			dHf = NULL;
	XLALDestroyCOMPLEX8FrequencySeries(inprod_integrand);			inprod_integrand = NULL;
	
	XLALDestroyREAL4FFTPlan(fwdRealPlan);
  XLALDestroyREAL4FFTPlan(revRealPlan);
   
  /*********************************************************************
   *
   *  Detach Error handling
   * 
   ********************************************************************/
  
  
  DETATCHSTATUSPTR( status );
  RETURN(status);
}


/* <lalVerbatim file="LALInspiralComputeSNRCP">  */
void LALInspiralComputeSNR (
														LALStatus															 *status,
														PPNConsistencyParamStruc               *PPNparams,
														REAL8																	 fStart,
														REAL8																	 fStop,
														REAL8FrequencySeries									 *psd,
														REAL8																	 *SNR)
/* </lalVerbatim> */
{
	
  /*********************************************************************
   *
   *  Error handling
   * 
   ********************************************************************/ 
  INITSTATUS( status, "LALInspiralComputeSNR", LALINSPIRALCOMPUTEWAVEFORMACCURACIESC);
  ATTATCHSTATUSPTR(status);
	
  /*********************************************************************
   *
   *  Temporary Variables
   * 
   ********************************************************************/ 
  
	//Counters
	INT4 i,j,k;
	
	//TIME/FREQUENCY SERIES SIZE
	INT4 N_ts;
	INT4 N_fs;
	REAL8 deltaT = 1.0/16384.0;
	REAL8 deltaF = 0.0;
	
	// DYNAMICAL RANGE
	REAL8 dynRange = 1E+27;
	
	// REAL NUMBERS
  REAL8 invpsd;										// Inverse PSD		
  REAL8 SNRtot;										// Total SNR	
  
  // PARAMETERS AND WAVEFORMS
  CoherentGW    	SNRwaveform;		// Holding waveform, loaded onto ht
	
	// TIME/FREQUENCY SERIES
	REAL4TimeSeries 				*ht;				// Initial Waveform
	REAL4TimeSeries 				*invFT;			// Store inverse Fourier Transform
  COMPLEX8FrequencySeries	*Hf;				// FT of ht
  COMPLEX8FrequencySeries	*snrpower;	// PSD of SNR	
	
	// FOURIER TRANSFORMS
	RealFFTPlan *fwdRealPlan = NULL;	// FT plan
  RealFFTPlan *revRealPlan = NULL;	// RFT plan
  
  
  /*********************************************************************
   *
   *  COMPUTE COHERENTGW FROM INPUT PARAMETERS
   * 
   ********************************************************************/  
  
  /* CLEARING COHERENTGW */
	memset(&SNRwaveform, 0, sizeof(CoherentGW) );  
	
  // GENERATE WAVEFORM AND COMBINE HP AND HC
  LALGeneratePPNAmpCorConsistency(status->statusPtr, &SNRwaveform, PPNparams);
  
  // SET SIZE TIME/FREQUENCY SERIES
  N_ts = SNRwaveform.f->data->length;
  N_fs = N_ts/2+1;
  
  // SET STEPSIZE TIME/FREQUENCY SERIES
  deltaF = 1.0/((REAL8) N_ts*deltaT);
  
  //LALStrainInGeocentricTime(status->statusPtr,ACS->detector, &SNRparams, SNRwaveform, ht);	   
	
  /*********************************************************************
   *
   *  CREATING TIME/FREQUENCY SERIES
   * 
   ********************************************************************/
  
	ht 				= XLALCreateREAL4TimeSeries("ht", &(PPNparams->epoch), 0.0, deltaT, &lalStrainUnit, N_ts);
	invFT			=	XLALCreateREAL4TimeSeries("invFT", &(PPNparams->epoch), 0.0, deltaT,&lalStrainUnit, N_ts);
	Hf				= XLALCreateCOMPLEX8FrequencySeries("Hf", &(PPNparams->epoch), 0.0, deltaF, &lalStrainUnit, N_fs);
	snrpower	= XLALCreateCOMPLEX8FrequencySeries("snrpower", &(PPNparams->epoch), 0.0, deltaF, &lalStrainUnit, N_fs);

	// CLEAR TIME SERIES	
	for(j=0;j<N_ts;j++){
		ht->data->data[j] 		= 0.0; 
		invFT->data->data[j]	= 0.0;
	}  
	
	for(j=0; j<N_fs; j++){
		Hf->data->data[j].re				= 0.0; 
		Hf->data->data[j].im				= 0.0; 
		snrpower->data->data[j].re	= 0.0;
		snrpower->data->data[j].im	= 0.0;
	}
  
  // FFT
  fwdRealPlan = XLALCreateForwardREAL4FFTPlan(N_ts,0);
  revRealPlan = XLALCreateReverseREAL4FFTPlan(N_ts,0);
	
	LALInspiralCombinePlusCross(status->statusPtr, &SNRwaveform, PPNparams, ht);
  
  /*********************************************************************
   *
   *  Compute SNR
   * 
   ********************************************************************/
  
  // PERFORM FFT ON HT
  LALTimeFreqRealFFT( status->statusPtr, Hf, ht, fwdRealPlan );
  
  // Scale using dynamical range
  for (i = 0; i < N_fs; ++i)
  {
    Hf->data->data[i].re *= dynRange; 
    Hf->data->data[i].im *= dynRange; 
  }  
	
	// Calculate Optimal SNR
	for (k = 0; k < N_fs; ++k)
  {
    if (psd->data->data[k] == 0.0 || ((REAL8) k)*deltaF < fStart || ((REAL8) k)*deltaF > fStop )
      invpsd = 0.0;
    else
      invpsd = 1.0 / ( dynRange * dynRange * psd->data->data[k]*2.0);
		
    snrpower->data->data[k].re = (Hf->data->data[k].re * Hf->data->data[k].re + Hf->data->data[k].im * Hf->data->data[k].im) * invpsd;
    snrpower->data->data[k].im = 0.0;
  }
  
  // Reserve FT, FIRST ELEMENT IS THE INTEGRAL!!
	LALFreqTimeRealFFT( status->statusPtr, invFT, snrpower, revRealPlan );
	
	// Compute Total SNR
	SNRtot = sqrt(4.0*invFT->data->data[0]); 
  *SNR = SNRtot;
  
  /*********************************************************************
   *
   *  Compute SNR
   * 
   ********************************************************************/
  
	XLALDestroyREAL4TimeSeries(ht);							ht = NULL;
	XLALDestroyREAL4TimeSeries(invFT);					invFT = NULL;
	XLALDestroyCOMPLEX8FrequencySeries(Hf);			Hf = NULL;
	XLALDestroyCOMPLEX8FrequencySeries(snrpower); snrpower =NULL;
  XLALDestroyCoherentGW(&SNRwaveform);
  XLALDestroyREAL4FFTPlan(fwdRealPlan);
  XLALDestroyREAL4FFTPlan(revRealPlan);
	
  /*********************************************************************
   *
   *  Detach Error handling
   * 
   ********************************************************************/
  
  DETATCHSTATUSPTR( status );
  RETURN(status);
}


void LALInspiralCombinePlusCross(
																 LALStatus									*status,			// Status pointer
																 CoherentGW									*GravWav,			// Input waveform, data[2i] = hp, data[2i+1] = hc
																 PPNConsistencyParamStruc		*params,				// input parameters
																 REAL4TimeSeries						*ht)					// combined h(t) from hp and hc
{
	/*********************************************************************
	 * 
	 * This function combines the output of a CoherentGW structure that has
	 * the data in the form of hp = data[2i], hc = data[2i+1], where i
	 * is the ith element of the array.
	 * 
	 * It uses a beam pattern function appropriate for an L-shaped detector
	 * such as LIGO or VIRGO. The use different setups (such as ET) 
	 * has not been implemented yet.
	 * 
	 *********************************************************************/	
	
	
	/*********************************************************************
	 * 
	 * LAL error handling
	 * 
	 *********************************************************************/	
  
  INITSTATUS( status, "LALInpiralCombinePlusCross", LALINSPIRALCOMPUTEWAVEFORMACCURACIESC);
  ATTATCHSTATUSPTR(status);	
  
	/*********************************************************************
	 * 
	 * Check input parameters
	 * 
	 *********************************************************************/	  
  

	
	/*********************************************************************
	 * 
	 * Create Temporary variables
	 * 
	 *********************************************************************/		
	
	/* COUNTERS */
	UINT4 i;
	
	/* BINARY SYSTEM PARAMETERS */
	REAL8		theta;	// Zenith Angle
	REAL8		phi;		// Azimuthal Angle
	REAL4		psi;		// Polarisation Angle
	
	/* BEAM PATTERN FUNCTIONS */
	REAL4		Fp = 0.0;
	REAL4		Fc = 0.0;
	
	/* TimeSeries for storing hplus and hcross from Coherent GW */
	REAL4TimeSeries			*hp;
	REAL4TimeSeries			*hc;

	hp = XLALCreateREAL4TimeSeries("hp", &(ht->epoch), 0.0, ht->deltaT, &lalStrainUnit, ht->data->length);
	hc = XLALCreateREAL4TimeSeries("hc", &(ht->epoch), 0.0, ht->deltaT, &lalStrainUnit, ht->data->length);
	
	/* CLEARING OUT TIMESERIES */
	for(i=0;i<ht->data->length;i++)
	{
		hp->data->data[i] = 0.0;
		hc->data->data[i] = 0.0;
	}
	
	/*********************************************************************
	 * 
	 * Read the two polariations in CoherentGW into hp and hc
	 * 
	 *********************************************************************/	
  
  if(ht->data->length > GravWav->f->data->length)
  {
		for(i=0;i<GravWav->f->data->length;i++)
		{
			hp->data->data[ht->data->length-GravWav->f->data->length+i] = GravWav->h->data->data[2*i];
			hc->data->data[ht->data->length-GravWav->f->data->length+i] = GravWav->h->data->data[2*i+1];
		}
	}
	else
	{
		for(i=0;i<ht->data->length;i++)
		{
			hp->data->data[i] = GravWav->h->data->data[2*i];
			hc->data->data[i] = GravWav->h->data->data[2*i+1];
		}
	}
	
	/*********************************************************************
	 * 
	 * Calculate Fp and Fc from the input parameters 
	 * - Polarisation angle psi
	 * - Azimuthal angle phi
	 * - Zenith angle, theta
	 * 
	 *********************************************************************/								
  
	theta = params->position.latitude;
	phi		= params->position.longitude;
	psi		= params->psi;
	
	Fp = 1.0/2.0 * ( 1 + pow(cos(theta),2.0) ) * cos(2.0*phi) * cos(2.0*psi)
	- cos(theta) * sin(2.0*phi) * sin(2.0*psi);
	
	
	//fprintf(stdout, "Fp = %e \n", Fp);
	
	Fc =  1.0/2.0 * (1 + pow(cos(theta),2.0) ) * cos(2.0*phi) * sin(2.0*psi)
	+ cos(theta) * sin(2.0*phi) * cos(2.0*psi);
	
	//fprintf(stdout, "Fc = %e \n", Fc);		
  
	/*********************************************************************
	 * 
	 * Calculate Detector Response for equitorial coordinates
	 * 
	 *********************************************************************/		
  
	
	/*********************************************************************
	 * 
	 * Combining hp and hc together with Fp and Fc to get h(t)
	 * 
	 *********************************************************************/
  
  for(i=0;i<ht->data->length;i++)
  {
		ht->data->data[i] = Fp * hp->data->data[i] + Fc * hc->data->data[i];
		//fprintf(stdout, "%e | %e | %e | %e | %e \n",Fp, Fc, hp->data->data[i],hc->data->data[i],ht->data->data[i]);
	}
	
	/*********************************************************************
	 * 
	 * Clean up
	 * 
	 *********************************************************************/	
	
	// TIME SERIES
	TRY( XLALDestroyREAL4TimeSeries(hp),	status);	
	TRY( XLALDestroyREAL4TimeSeries(hc),	status);
	
	/*********************************************************************
	 * 
	 * Detach LAL error handling
	 * 
	 *********************************************************************/	  
  
  DETATCHSTATUSPTR(status);
  RETURN(status);	
}

void XLALDestroyCoherentGW (
														CoherentGW  *GravWav   )
{
	
  if( GravWav )
  {
    if(GravWav->h     && GravWav->h->data)      XLALDestroyVectorSequence(GravWav->h->data);
    if(GravWav->a     && GravWav->a->data)      XLALDestroyVectorSequence(GravWav->a->data);
    if(GravWav->f     && GravWav->f->data)      XLALDestroyREAL4TimeSeries(GravWav->f); 
    if(GravWav->phi   && GravWav->phi->data)    XLALDestroyREAL8TimeSeries(GravWav->phi); 
    if(GravWav->shift && GravWav->shift->data)  XLALDestroyREAL4TimeSeries(GravWav->shift);
  }
  return;
}

/* </lalVerbatim> */
