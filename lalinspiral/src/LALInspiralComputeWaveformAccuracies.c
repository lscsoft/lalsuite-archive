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
														REAL8																	 *dhindh)
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
  //printf("exact_fstartin = %e | model_fstartin = %e \n", exactParams->fStopIn, modelParams->fStopIn);
  LALGeneratePPNAmpCorConsistency(status->statusPtr, &exactWave, exactParams);
  //REPORTSTATUS(status->statusPtr);
  LALGeneratePPNAmpCorConsistency(status->statusPtr, &modelWave, modelParams);
  //REPORTSTATUS(status->statusPtr);
  
  // SET SIZE TIME/FREQUENCY SERIES
  N_ts = modelWave.f->data->length;
  N_fs = N_ts/2+1;
  
  // SET STEPSIZE TIME/FREQUENCY SERIES
  deltaF = 1.0/((REAL8) N_ts*deltaT);
  
  /*********************************************************************
   *
   *  Create PSD
   * 
   ********************************************************************/    
  
  /* PSD */
  REAL8FrequencySeries *psd;
  
  // PSD MODEL
  void (*noisemodel)(LALStatus*,REAL8*,REAL8) = LALAdvLIGOPsd;
  
  // CREATE FREQUENCY SERIES FOR PSD COMPUTING
  psd = XLALCreateREAL8FrequencySeries("PSD_LIGO",  &(exactParams->epoch), 0.0, deltaF, &lalStrainUnit, N_fs);
  
  // Computing LIGO PSD
  LALNoiseSpectralDensity( status->statusPtr, psd->data, noisemodel, psd->deltaF );
  /*
  FILE *psd_out;
  psd_out =fopen("psd.dat","w");
  for(i=0; i<N_fs; i++){
  	fprintf(psd_out, "%e\t%e\n", (REAL8) i, psd->data->data[i]);
  }
  fclose(psd_out);
  */
  
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
	/*
	FILE *exact_out;
	exact_out = fopen("exact.dat", "w");
	FILE *model_out;
	model_out = fopen("model.dat", "w");
	
	for(i=0; i<N_ts; i++){
		fprintf(exact_out, "%e\t%e\t%e\t%e\t%e\n", i*deltaT, exactWave.h->data->data[2*i], exactWave.h->data->data[2*i+1], exactWave.phi->data->data[i], cos(exactWave.phi->data->data[i]));
		fprintf(model_out, "%e\t%e\t%e\t%e\t%e\n", i*deltaT, modelWave.h->data->data[2*i], modelWave.h->data->data[2*i+1], modelWave.phi->data->data[i], cos(modelWave.phi->data->data[i]));
	}
	fclose(exact_out);
	fclose(model_out);
		
	FILE *dht_out;
	dht_out=fopen("dht.dat","w");
	*/
	for(i=0;i<N_ts;i++){
		dht->data->data[i] = model_ht->data->data[i] - exact_ht->data->data[i];
		//fprintf(dht_out, "%e\t%e\t%e\t%e\n", (REAL8) i, model_ht->data->data[i], exact_ht->data->data[i], dht->data->data[i]);
	}  
	//fclose(dht_out);
	
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
	*dhindh = sqrt(4.0*invFT->data->data[0]); 
  
  /*********************************************************************
   *
   *  FREE MEMORY
   * 
   ********************************************************************/
   
  XLALDestroyREAL4TimeSeries(exact_ht);					exact_ht = NULL;
  XLALDestroyREAL4TimeSeries(model_ht);					model_ht = NULL;
  XLALDestroyREAL4TimeSeries(dht);					dht = NULL;
  XLALDestroyREAL4TimeSeries(invFT);					invFT = NULL;
  XLALDestroyREAL8FrequencySeries(psd);					psd = NULL;
	XLALDestroyCOMPLEX8FrequencySeries(dHf);			dHf = NULL;
	XLALDestroyCOMPLEX8FrequencySeries(inprod_integrand);			inprod_integrand = NULL;
	
	XLALDestroyCoherentGW(&exactWave);
	XLALDestroyCoherentGW(&modelWave);
	
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
														LALStatus												*status,
														InspiralTemplate               	*inputParams,
														REAL8														*dphis,
														REAL8														*SNR)
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
  //printf("CREATING VARIABLES \n");
	//Counters
	INT4 i,j,k;
	
	//TIME/FREQUENCY SERIES SIZE
	INT4 N_ts = 0;
	INT4 N_fs = 0;
	REAL8 deltaT = 1.0/inputParams->tSampling;
	REAL8 deltaF = 0.0;
	
  // SET SIZE TIME/FREQUENCY SERIES - TO BE REVISED
  //N_ts = 1000;
  //N_fs = N_ts/2+1;
  
  // SET STEPSIZE TIME/FREQUENCY SERIES
  //deltaF = 1.0/((REAL8) N_ts*deltaT);  	
	
	// DYNAMICAL RANGE
	REAL8 dynRange = 1E+27;
	
	// REAL NUMBERS
  REAL8 invpsd;										// Inverse PSD		
  REAL8 SNRtot;										// Total SNR	
  
  // PARAMETERS AND WAVEFORMS
  CoherentGW    	SNRwaveform;		// Holding waveform, loaded onto ht
	
	// TIME/FREQUENCY SERIES
	REAL4TimeSeries 				*ht = NULL;				// Initial Waveform
	REAL4TimeSeries 				*invFT = NULL;			// Store inverse Fourier Transform
  COMPLEX8FrequencySeries	*Hf = NULL;				// FT of ht
  REAL4FrequencySeries	*Hf_F2 = NULL;				// FT of ht
  COMPLEX8FrequencySeries	*snrpower = NULL;	// PSD of SNR	
	
	// FOURIER TRANSFORMS
	RealFFTPlan *fwdRealPlan = NULL;	// FT plan
  RealFFTPlan *revRealPlan = NULL;	// RFT plan
  
  /*********************************************************************
   *
   *  CREATING TIME/FREQUENCY SERIES
   * 
   ********************************************************************/
  //printf("MAKING TIME/FREQ SERIES \n");
  // END TIME USED AS START TIME - NEEDS TO BE FIXED
	/*
	invFT			=	XLALCreateREAL4TimeSeries("invFT", &(inputParams->end_time), 0.0, deltaT,&lalStrainUnit, N_ts);
	Hf				= XLALCreateCOMPLEX8FrequencySeries("Hf", &(inputParams->end_time), 0.0, deltaF, &lalStrainUnit, N_fs);
	Hf_F2			= XLALCreateREAL4FrequencySeries("Hf_F2", &(inputParams->end_time), 0.0, deltaF/2.0, &lalStrainUnit, 2*N_fs);
	snrpower	= XLALCreateCOMPLEX8FrequencySeries("snrpower", &(inputParams->end_time), 0.0, deltaF, &lalStrainUnit, N_fs);

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
  */
  // FFT
  
  
  /*********************************************************************
   *
   *  COMPUTE WAVEFORM, FFT IF WAVEFORM IS IN TIMEDOMAIN
   * 
   ********************************************************************/  
	
	//printf("COMPUTING WAVEFORMS \n");
	
	switch( inputParams->approximant ){
		case AmpCorPPNTest:
			/* CLEARING COHERENTGW */
			memset(&SNRwaveform, 0, sizeof(CoherentGW) ); 
			PPNConsistencyParamStruc *inputParamsConsistency; 
			inputParamsConsistency = (PPNConsistencyParamStruc *) LALMalloc( sizeof(PPNConsistencyParamStruc) );
			memset(inputParamsConsistency, 0, sizeof(PPNConsistencyParamStruc));
			XLALCopyPPNConsistencyFromInspiralTemplate(inputParams, inputParamsConsistency, dphis);
			LALGeneratePPNAmpCorConsistency(status->statusPtr, &SNRwaveform, inputParamsConsistency);
			
			// SET NTS/NFS DT/DF
			N_ts = SNRwaveform.f->data->length;
			N_fs = N_ts/2+1;
			deltaT = 1.0/inputParams->tSampling;
			deltaF = 1.0/(N_ts*deltaT);
			
			// CREATE NECESSARY SERIES
			ht 				= XLALCreateREAL4TimeSeries("ht", &(inputParams->end_time), 0.0, deltaT, &lalStrainUnit, N_ts);
			for(j=0;j<N_ts;j++){ ht->data->data[j] = 0.0;}  
			
			LALInspiralCombinePlusCross(status->statusPtr, &SNRwaveform, inputParamsConsistency, ht);
			
			FILE *ht_out;
			ht_out = fopen("ht.dat", "w");
			for(j=0;j<(INT4)SNRwaveform.h->data->length; j++){
				fprintf(ht_out, "%e\t%e\t%e\t%e\n", j*ht->deltaT, SNRwaveform.h->data->data[2*j], SNRwaveform.h->data->data[2*j+1], ht->data->data[j]);
			}
			fclose(ht_out);
			
			Hf = XLALCreateCOMPLEX8FrequencySeries("Hf", &(inputParams->end_time), 0.0, deltaF, &lalStrainUnit, N_fs);
			for(j=0; j<N_fs; j++){Hf->data->data[j].re = 0.0; Hf->data->data[j].im = 0.0;}
			
			fwdRealPlan = XLALCreateForwardREAL4FFTPlan(N_ts,0);
			LALTimeFreqRealFFT( status->statusPtr, Hf, ht, fwdRealPlan );
			
			XLALDestroyREAL4TimeSeries(ht);	ht = NULL;
			XLALDestroyREAL4FFTPlan(fwdRealPlan);			fwdRealPlan=NULL;
			
			REPORTSTATUS( status->statusPtr );
			break;
			
		case TaylorF2Test:
			deltaT = 1.0/inputParams->tSampling;
			deltaF = 1.0/inputParams->tC;
			N_ts = (INT4) floor(inputParams->tC/deltaT)+1;
			N_fs = N_ts/2+1;
			
			Hf = XLALCreateCOMPLEX8FrequencySeries("Hf", &(inputParams->end_time), 0.0, deltaF, &lalStrainUnit, N_fs);
			Hf_F2 = XLALCreateREAL4FrequencySeries("Hf_F2", &(inputParams->end_time), 0.0, deltaF/2.0, &lalStrainUnit, 2*N_fs);
			for(j=0; j<2*N_fs; j++){Hf_F2->data->data[j] = 0.0;}
			for(j=0; j<N_fs; j++){Hf->data->data[j].re = 0.0; Hf->data->data[j].im = 0.0;}
		
			LALInspiralParameterCalc(status->statusPtr,inputParams);
			LALInspiralRestrictedAmplitude(status->statusPtr,inputParams);
			
			expnCoeffs ak;
			expnFunc expnFunction;
			
			memset(&ak,0,sizeof(expnCoeffs));
			
			LALInspiralSetup(status->statusPtr, &ak, inputParams);
			LALInspiralChooseModel(status->statusPtr, &expnFunction, &ak, inputParams);
			
			LALInspiralStationaryPhaseApprox2(status->statusPtr,Hf_F2->data,inputParams);
			
			for(j=0; j<N_fs; j++){
				Hf->data->data[j].re = Hf_F2->data->data[j];
				Hf->data->data[j].im = Hf_F2->data->data[2*N_fs-j];
			}
			
			XLALDestroyREAL4FrequencySeries(Hf_F2);					Hf_F2 = NULL;
			
			break;		
		
		default:
			printf("ERROR: UNKNOWN APPROXIMANT\n");
			//exit(-1);
  }
  
  /*********************************************************************
   *
   *  Create PSD
   * 
   ********************************************************************/    
  //printf("COMPUTING PSD \n");
  /* PSD */
  REAL8FrequencySeries *psd;
  
  // PSD MODEL
  void (*noisemodel)(LALStatus*,REAL8*,REAL8) = LALAdvLIGOPsd;
  
  // CREATE FREQUENCY SERIES FOR PSD COMPUTING
  psd = XLALCreateREAL8FrequencySeries("PSD_LIGO",  &(inputParams->end_time), 0.0, deltaF, &lalStrainUnit, N_fs);
  
  // Computing LIGO PSD
  LALNoiseSpectralDensity( status->statusPtr, psd->data, noisemodel, psd->deltaF );   	
  
  /*********************************************************************
   *
   *  Compute SNR
   * 
   ********************************************************************/
	
  // Scale using dynamical range
  
  FILE *Hf_file;
	Hf_file = fopen("Hf.dat", "w");
			
  for (i = 0; i < N_fs; ++i)
  {
    Hf->data->data[i].re *= dynRange; 
    Hf->data->data[i].im *= dynRange; 
    
		fprintf(Hf_file, "%e\t%e\n", i*Hf->deltaF, sqrt(Hf->data->data[i].re*Hf->data->data[i].re+Hf->data->data[i].im*Hf->data->data[i].im)/dynRange);
  }  
	
	fclose(Hf_file);
	
	// Calculate Optimal SNR
	snrpower = XLALCreateCOMPLEX8FrequencySeries("snrpower", &(inputParams->end_time), 0.0, deltaF, &lalStrainUnit, N_fs);
	for(j=0; j<N_fs; j++){snrpower->data->data[j].re = 0.0; snrpower->data->data[j].im = 0.0;}	
	
	for (k = 0; k < N_fs; ++k)
  {
    if (psd->data->data[k] == 0.0 || ((REAL8) k)*deltaF < inputParams->fLower || ((REAL8) k)*deltaF > inputParams->fCutoff )
      invpsd = 0.0;
    else
      invpsd = 1.0 / ( dynRange * dynRange * psd->data->data[k]*2.0);
		
    snrpower->data->data[k].re = (Hf->data->data[k].re * Hf->data->data[k].re + Hf->data->data[k].im * Hf->data->data[k].im) * invpsd;
    snrpower->data->data[k].im = 0.0;
  }
  
  // Reserve FT, FIRST ELEMENT IS THE INTEGRAL!!
  invFT	=	XLALCreateREAL4TimeSeries("invFT", &(inputParams->end_time), 0.0, deltaT,&lalStrainUnit, N_ts);
	for(j=0;j<N_ts;j++){invFT->data->data[j] = 0.0;}
	
	revRealPlan = XLALCreateReverseREAL4FFTPlan(N_ts,0);
	LALFreqTimeRealFFT( status->statusPtr, invFT, snrpower, revRealPlan );
	
	// Compute Total SNR
	SNRtot = sqrt(4.0*invFT->data->data[0]); 
  *SNR = SNRtot;
  
  /*********************************************************************
   *
   *  Compute SNR
   * 
   ********************************************************************/
  //printf("CLEANING \n");
	XLALDestroyREAL4TimeSeries(invFT);					invFT = NULL;
	XLALDestroyCOMPLEX8FrequencySeries(Hf);			Hf = NULL;
	XLALDestroyCOMPLEX8FrequencySeries(snrpower); snrpower =NULL;
	XLALDestroyREAL8FrequencySeries(psd);					psd = NULL;
  if(inputParams->approximant==AmpCorPPNTest) XLALDestroyCoherentGW(&SNRwaveform);
  XLALDestroyREAL4FFTPlan(revRealPlan);			revRealPlan=NULL;
	
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
