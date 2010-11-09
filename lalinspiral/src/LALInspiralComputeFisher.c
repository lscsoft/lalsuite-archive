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

/*
 
 TO DO LIST:
 
 1. include derivatives for phi_T
 2. take derivative in time by using SG filter
 3. Make testACS a pointer that is passed to save space
 
 */

#if 0
<lalVerbatim file="LALInspiralComputeFisherCV">
Author: Tjonnie Li, Chris Van Den Broeck
$Id$
</lalVerbatim>

<lalLaTeX>
\subsection{Module \texttt{LALInspiralComputeFisher.c}}

Module to compute the Fisher matrix for a variety of time domain waveforms.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{XLALInspiralComputeFisherMatrixCP}
\idx{XLALInspiralComputeFisherMatrix()}
\begin{itemize}
\item \texttt{fisher,} Output, the Fisher matrix at the point defined by \texttt{params}
\item \texttt{psd,} Input, the power spectral density of the data
\item \texttt{params,} Input, the parameters where Fisher must be computed.
\end{itemize}

\input{XLALInspiralComputeTDWDerivCP}
\idx{XLALInspiralComputeTDWDeriv()}
\begin{itemize}
\item \texttt{Wderiv,} Output, the time derivative of waveform at the point
point defined by \texttt{params}
\item \texttt{psd,}  Input, the power spectral density of the data
\item \texttt{params,} Input, the parameters where Fisher must be computed
\item \texttt{paramid,} Input, id of the parameter to take derivative on
\item \texttt{initdelta,} Input, initial difference in parameters
\item \texttt{tolerance,} Input, stop iteration when difference between two
bisections is smaller than tolerance.
\end{itemize}


\subsubsection*{Description}
We calculate the usual Fisher matrix.

\subsubsection*{Algorithm}

\subsubsection*{Uses}

\begin{verbatim}
LALMalloc
LALFree
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{LALInspiralComputeFisherCV}}

</lalLaTeX>
#endif

#include <lal/LALStdlib.h>
#include <lal/LALInspiralComputeFisher.h>

NRCSID(LALINSPIRALCOMPUTEFISHERC, "$Id: LALInspiralComputeFisher.c, v1.0 2 February 2010 $" );

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

void XLALCopyPPNConsistencyParamStruc (
																			 PPNConsistencyParamStruc             *inputParams,
																			 PPNConsistencyParamStruc             *duplicateParams )
{
  // Counter
  UINT4 i = 0;
  
  if ( inputParams && inputParams->ppn && duplicateParams)
  {
		/* Passed parameters. */
		duplicateParams->position = inputParams->position; /* location of source on sky */
		duplicateParams->psi = inputParams->psi;            /* polarization angle (radians) */
		duplicateParams->epoch = inputParams->epoch;    /* start time of output time series */
		
		/* Input parameters. */
		duplicateParams->mTot_real8 = inputParams->mTot_real8; /* total system mass (Msun) */
		duplicateParams->eta_real8 = inputParams->eta_real8;  /* mass ratio */
		duplicateParams->delta = inputParams->delta;      /* sqrt(1-4eta) */
		duplicateParams->mTot = inputParams->mTot;       /* total system mass (Msun) */
		duplicateParams->eta = inputParams->eta;        /* mass ratio */
		duplicateParams->d = inputParams->d;          /* distance (metres) */
		duplicateParams->inc = inputParams->inc;        /* inclination angle (radians) */
		duplicateParams->cosI = inputParams->cosI;				/* cosine of inclination angle */
		duplicateParams->sinI = inputParams->sinI;				/* sine of inclination angle */
		duplicateParams->phi = inputParams->phi;        /* coalescence phase (radians) */
		duplicateParams->deltaT = inputParams->deltaT;     /* requested sampling interval (s) */
		duplicateParams->fStartIn = inputParams->fStartIn;   /* requested start frequency (Hz) */
		duplicateParams->fStopIn = inputParams->fStopIn;    /* requested stop frequency (Hz) */
		duplicateParams->lengthIn = inputParams->lengthIn;   /* maximum length of waveform */
		duplicateParams->ampOrder = inputParams->ampOrder;    /* PN amplitude selection 0-5 */
		
		/* PN phasing coefficients for use in AmpCorConsistency */
		for (i=0; i<9; i++) {
			duplicateParams->phasePNparams[i] = inputParams->phasePNparams[i];
		}
		
		/*
		 duplicateParams->phi0 = inputParams->phi0;
		 duplicateParams->phi2 = inputParams->phi2;
		 duplicateParams->phi3 = inputParams->phi3;
		 duplicateParams->phi4 = inputParams->phi4;
		 duplicateParams->phi5 = inputParams->phi5;
		 duplicateParams->phi5l = inputParams->phi5l;
		 duplicateParams->phi6 = inputParams->phi6;
		 duplicateParams->phi6l = inputParams->phi6l;
		 duplicateParams->phi7 = inputParams->phi7;
		 */
    
    // Check if *.ppn vector is created
    if ( !(duplicateParams->ppn) )
    { // if not, create vector with the same length as input.ppn
      duplicateParams->ppn = XLALCreateREAL4Vector(inputParams->ppn->length);
      for (i=0; i<inputParams->ppn->length; i++)
      {
        //fprintf(stdout, "Copy PPNConsistencyParamStruc -> created new vector \n");
        duplicateParams->ppn->data[i] = inputParams->ppn->data[i];
      }      
    }
    else if( duplicateParams->ppn->length == inputParams->ppn->length)
    { // check if the lengths are the same before filling
      for (i=0; i<inputParams->ppn->length; i++)
      {
        //fprintf(stdout, "Copy PPNConsistencyParamStruc -> vector exists \n");
        duplicateParams->ppn->data[i] = inputParams->ppn->data[i];
      }
    }
    else
    {
			fprintf(stdout, "Could not copy PPNConsistencyParamStruc \n");
    }
  }
  else
  {
    fprintf(stdout, "Empty input parameter \n");
  }
	
  return;
}

/* <lalVerbatim file="LALInspiralComputeFisherMatrixCP">  */
void LALInspiralComputeFisherMatrix (
																		 LALStatus                              *status,  
																		 REAL8FrequencySeries                   *psd,
																		 PPNConsistencyParamStruc                          *params,
																		 REAL4																	 *Fisher,
																		 FisherACS															 *ACS 
																		 )
/* </lalVerbatim> */
{ 
  if(ACS->verbose_switch == 1) fprintf(stdout, "Running LALInspiralComputeFisherMatrix \n");
  
  /*********************************************************************
   *
   *  Error handling
   * 
   ********************************************************************/ 
  
  INITSTATUS( status, "LALInspiralComputeDerivatives", LALINSPIRALCOMPUTEFISHERC);
  ATTATCHSTATUSPTR(status);
	
  /*********************************************************************
   *
   *  CREATING TEMPORARY VARIABLES
   * 
   ********************************************************************/ 
  
  if(ACS->verbose_switch == 1) fprintf(stdout, "... Initialising Variables \n");
  
  /* COUNTERS */
  INT4 i, j, k;         
  
  /* STRUCTURES */
  REAL4TimeSeries *hderiv[ACS->N_inputparams]; // Waveform derivatives
  
  /* FISHER MATRIX */
  for(i=0;i<(ACS->N_inputparams*ACS->N_inputparams);i++)
  {
		Fisher[i] = 0.0;
	}
  
  /* LOCAL COPY OF INPUT PARAMETERS */
  PPNConsistencyParamStruc PPNparams;
	
  // Fourier Transform of derivatives
  COMPLEX8FrequencySeries *Hfderiv[ACS->N_inputparams];  // FT of derivatives
  
  /*********************************************************************
   *
   *  OUTPUT FILES
   * 
   ********************************************************************/   
  
  //FILE *fisher_out;
  //FILE *cov_out;
  FILE *derivatives_out;
  FILE *fourierderivs_out;
  
  //char fisher_name[128];
  //char cov_name[128];
  char derivatives_name[128];
  char fourierderivatives_name[128];
	
	//char folder[] = ACS->folder;
  
  //sprintf(fisher_name, "%s%s%d%s", ACS->folder, "fisher", ACS->seed, ".dat");
  //sprintf(cov_name, "%s%s%d%s", ACS->folder, "covariance", ACS->seed, ".dat");
  sprintf(derivatives_name, "%s%s%d%s", ACS->folder, "derivatives", ACS->seed, ".dat");
  sprintf(fourierderivatives_name, "%s%s%d%s", ACS->folder, "fourierderivatives", ACS->seed, ".dat");
  
  //sprintf(traject_name,"%s%s%d%s", ACS->folder, traject_base, paramid, ".dat");  
  //sprintf(hx_name, "%s%s%d%s", ACS->folder, hx_base, paramid, ".dat");
  //sprintf(linReg_name, "%s%s%d%s", ACS->folder, linReg_base, paramid, ".dat");
  
  //fisher_out = fopen(fisher_name, "w");
  //cov_out		= fopen(cov_name, "w");
  derivatives_out = fopen(derivatives_name ,"w");
  fourierderivs_out = fopen(fourierderivatives_name ,"w");
  
  /*********************************************************************
   *
   *  SETTING VARIABLES
   * 
   ********************************************************************/
  
  if(ACS->verbose_switch == 1) fprintf(stdout, "... Setting Variables \n");
  
  /* COPY INPUT PARAMETERS */
  PPNparams.ppn = XLALCreateREAL4Vector( params->ppn->length );
  TRY( XLALCopyPPNConsistencyParamStruc(params, &PPNparams), status);
  
  // Creating Time/Frequency Series for derivatives
  for (i = 0; i < ACS->N_inputparams; ++i) 
  {
    if( (hderiv[i] = (REAL4TimeSeries *) LALMalloc( sizeof(REAL4TimeSeries) ) ) == NULL  &&
			 (Hfderiv[i] = (COMPLEX8FrequencySeries *) LALMalloc( sizeof(COMPLEX8FrequencySeries) ) ) == NULL )
    {
      ABORT(status, LALINSPIRALCOMPUTEFISHERC_EMEM, LALINSPIRALCOMPUTEFISHERC_EMEM_MSG);
    }
    else
    { 
      /* CREATE TIME/FREQUENCY SERIES */
      hderiv[i] = XLALCreateREAL4TimeSeries("hderiv", &(ACS->epoch), 0.0, ACS->deltaT_ts, &lalStrainUnit, ACS->N_ts);
      Hfderiv[i] = XLALCreateCOMPLEX8FrequencySeries("Hfderiv", &(ACS->epoch), 0.0, ACS->deltaF_fs, &lalStrainUnit, ACS->N_ts/2+1);
      
      /* CLEARING TIME/FREQUENCY SERIES */
      for (j=0; j<ACS->N_ts; j++)
      {
        if (j<ACS->N_ts)       {hderiv[i]->data->data[j] = 0.0;}
        if (j<(ACS->N_fs)) {Hfderiv[i]->data->data[j].re = 0.0; Hfderiv[i]->data->data[j].im = 0.0;}
      }
      
    }
  }  
  
  /*********************************************************************
   *
   *  COMPUTING SNR
   * 
   ********************************************************************/
  
  if(ACS->verbose_switch == 1) fprintf(stdout, "... Computing Optimal SNR \n");
  
  /* COMPUTING SNR */
  REAL8 SNR;
  TRY( LALInspiralComputeSNR(status->statusPtr, psd, &PPNparams, &SNR, ACS), status);
  ACS->SNR = SNR;
  
	
  /*********************************************************************
   *
   *  COMPUTING DERIVATIVES
   * 
   ********************************************************************/
  
  if(ACS->verbose_switch == 1) fprintf(stdout, "... Computing derivatives \n");
  
  /* VARIALBES FOR SWITCHING BETWEEN MASSES AND PHI_0/PHI_2 DERIVATIVES */
	
	/* TEMPORARY TIME SERIES */
	REAL4TimeSeries *tmphderiv[2];
	
	tmphderiv[0] = (REAL4TimeSeries *) LALMalloc( sizeof( REAL4TimeSeries ) );
	tmphderiv[1] = (REAL4TimeSeries *) LALMalloc( sizeof( REAL4TimeSeries ) );
	
	tmphderiv[0] = XLALCreateREAL4TimeSeries("tmphderiv0", &(ACS->epoch), 0.0, ACS->deltaT_ts, &lalStrainUnit, ACS->N_ts);
	tmphderiv[1] = XLALCreateREAL4TimeSeries("tmphderiv1", &(ACS->epoch), 0.0, ACS->deltaT_ts, &lalStrainUnit, ACS->N_ts);
	
	for (j=0; j<((INT4)tmphderiv[0]->data->length); j++)
	{
		tmphderiv[0]->data->data[j] = 0.0;
		tmphderiv[1]->data->data[j] = 0.0;
	}
	
	//derivatives of mtot and eta for calculating new hderiv 0 and 1
	REAL4 dmtotdphi0 = pow(params->eta_real8,3.0/8.0)	*	(-3715.0+2772.0*params->eta_real8)	*	(params->mTot_real8)	*	pow(5.0*LAL_MTSUN_SI*params->mTot_real8,5.0/8.0)
	/(-1486.0+2772.0*params->eta_real8);
	REAL4 dmtotdphi2 = 12096.0	*	pow(params->eta_real8,5.0/8.0)	*	(params->mTot_real8)	*	pow(LAL_MTSUN_SI*params->mTot_real8,3.0/8.0)
	/(pow(5.0,5.0/8.0)	*	(-743.0+1386.0*params->eta_real8));
	REAL4 detadphi0 = 3.0	*	pow(5.0,5.0/8.0)	*	pow(params->eta_real8,11.0/8.0)	*	(743.0+924.0*params->eta_real8)	*	pow(LAL_MTSUN_SI*params->mTot_real8,5.0/8.0)
	/(-1486.0+2772.0*params->eta_real8);
	REAL4 detadphi2 = 4032.0	*	pow(5.0,3.0/8.0)	*	pow(params->eta_real8,13.0/8.0)	*	pow(LAL_MTSUN_SI*params->mTot_real8,3.0/8.0)
	/(743.0-1386.0*params->eta_real8);
	
  for (i = 0; i<ACS->N_inputparams; i++)
  {
    LALInspiralComputeDerivatives(status->statusPtr, hderiv[i], params, i, ACS);
  }
	
	/* CHANGE TO PHI_0 AND PHI_2 DERIVATIVES */
	
	if(ACS->mass2phi_switch == 1)
	{
		for(j=0; j<((INT4)hderiv[0]->data->length);j++)
		{
			tmphderiv[0]->data->data[j] = hderiv[0]->data->data[j];
			tmphderiv[1]->data->data[j] = hderiv[1]->data->data[j];
		}
		
		for(j=0; j<((INT4)hderiv[0]->data->length);j++)
		{
			hderiv[0]->data->data[j] = tmphderiv[0]->data->data[j]*dmtotdphi0+tmphderiv[1]->data->data[j]*detadphi0;
			hderiv[1]->data->data[j] = tmphderiv[0]->data->data[j]*dmtotdphi2+tmphderiv[1]->data->data[j]*detadphi2;
		}
	}
  
  for (i=0; i<ACS->N_inputparams; i++) 
	{
		LALInspiralWaveTaper(status->statusPtr, hderiv[i]->data, 3);
		if(ACS->SavitskyGolay_switch==1){for(k=0;k<1E3;k++)
		{
			if(ACS->verbose_switch == 1)printf("\r... Smoothing - %d percent done", k/10);
			LALSavitskyGolayFilter(status->statusPtr, ACS->SG_R, ACS->SG_M, hderiv[i], hderiv[i], 0);
		}}
	}
  
  if(ACS->verbose_switch == 1) fprintf(stdout, "\n");
	
	if(ACS->verbose_switch == 1)printf("... Writing to file \n");
	if(ACS->printall == 1)
	{
		for (i = 0; i < ACS->N_ts; ++i)
		{
			fprintf(derivatives_out, "%e\t", i*ACS->deltaT_ts);
			for (j = 0; j < ACS->N_inputparams; j++)
			{
				fprintf(derivatives_out, "%e", hderiv[j]->data->data[i]);
				if(j==(ACS->N_inputparams-1)) fprintf(derivatives_out, "\n");
				else fprintf(derivatives_out, "\t");
			}
		}
	}
	
	XLALDestroyREAL4TimeSeries(tmphderiv[0]);
	XLALDestroyREAL4TimeSeries(tmphderiv[1]);
	
	
  /*********************************************************************
   *
   *  COMPUTE FOURIER TRANSFORMS OF DERIVATIVES
   * 
   ********************************************************************/
	if(ACS->verbose_switch == 1)printf("... Computing Fourier Transforms \n");
	
	/* PERFORM FFT PLAN */
	
  for (i = 0; i < ACS->N_inputparams; ++i)
  {
    RealFFTPlan *fwdRealPlanDeriv;         // FFT plan for derivs
    fwdRealPlanDeriv = NULL;
    fwdRealPlanDeriv  = XLALCreateForwardREAL4FFTPlan(ACS->N_ts,0); 
    LALTimeFreqRealFFT( status->statusPtr, Hfderiv[i], hderiv[i], fwdRealPlanDeriv );
    TRY( LALDestroyRealFFTPlan( status->statusPtr, &fwdRealPlanDeriv), status); fwdRealPlanDeriv = NULL;
  } 
  
  /* RESCALE WITH DYNAMICAL RANGE */
  REAL8 tempRe[ACS->N_inputparams];
  REAL8 tempIm[ACS->N_inputparams];
  for(i=0; i<ACS->N_inputparams; i++)
  {
		tempRe[i] = 0.0;
		tempIm[i] = 0.0;
	}
	
  
  for (i = 0; i < ACS->N_fs; ++i)
  {
    //temp0Re = Hfderiv[0]->data->data[j].re*ACS->dynRange;
    //temp0Im = Hfderiv[0]->data->data[j].im*ACS->dynRange;
    //temp1Re = Hfderiv[1]->data->data[j].re*ACS->dynRange;
    //temp1Im = Hfderiv[1]->data->data[j].im*ACS->dynRange;
    //temp2Re = Hfderiv[2]->data->data[j].re*ACS->dynRange;
    //temp2Im = Hfderiv[2]->data->data[j].im*ACS->dynRange;
    //temp3Re = Hfderiv[3]->data->data[j].re*ACS->dynRange;
    //temp3Im = Hfderiv[3]->data->data[j].im*ACS->dynRange;
    //temp4Re = Hfderiv[4]->data->data[j].re*ACS->dynRange;
    //temp4Im = Hfderiv[4]->data->data[j].im*ACS->dynRange;
    
    if(ACS->printall == 1)fprintf(fourierderivs_out, "%e\t", i*ACS->deltaF_fs); 
    
    for(j=0; j<ACS->N_inputparams; j++)
    {
			
			tempRe[j] = 0.0;
			tempIm[j] = 0.0;
			tempRe[j] = Hfderiv[j]->data->data[i].re*ACS->dynRange;
			tempIm[j] = Hfderiv[j]->data->data[i].im*ACS->dynRange;
			
			if(ACS->printall == 1)
			{
				fprintf(fourierderivs_out, "%e", sqrt(tempRe[j]*tempRe[j] + tempIm[j]*tempIm[j])/ACS->dynRange); 
				if(j == (ACS->N_inputparams -1)){ fprintf(fourierderivs_out, "\n"); }
				else{ fprintf(fourierderivs_out, "\t"); }
			}
		}
  }
  
	
  /*********************************************************************
   *
   *  COMPUTE FISHER MATRIX
   * 
   ********************************************************************/
	
  if(ACS->verbose_switch == 1) fprintf(stdout, "... Computing Fisher matrix components \n");
  
  /* COMPUTE FISHER COMPONENTS */
  for (i = 0; i < ACS->N_inputparams; ++i)
  {
    for (j = 0; j < i + 1; ++j) 
    {
      LALInspiralComputeFisherComponents( status->statusPtr, Fisher, hderiv[i], hderiv[j], psd, i*ACS->N_inputparams + j, ACS);
    }
  }
	
  /* FILLING FISHER MATRIX */
  for (i = 0; i < ACS->N_inputparams; ++i)
  {
    for (j = i + 1; j < ACS->N_inputparams; ++j)
      Fisher[i*ACS->N_inputparams + j] = Fisher[j*ACS->N_inputparams + i];
  }
	
	//if( (i%ACS->N_inputparams) == (ACS->N_inputparams-1) )
	//{
	//fprintf(invert_out, "%e\n", Input[i]); 
	//fprintf(stdout, "%e\n", Input[i]); 
	//}
	//else
	//{
	//fprintf(invert_out, "%e\t", Input[i]);
	//fprintf(stdout, "%e\t", Input[i]);
	//}
  
  /*********************************************************************
   *
   *  CLEANING UP
   * 
   ********************************************************************/
  if(ACS->verbose_switch == 1) fprintf(stdout, "... Cleaning Up \n");
  
  /* DESTROY LOCAL COPY INPUT PARAMETERS */
  TRY( XLALDestroyREAL4Vector( PPNparams.ppn ), status); 
  
  /* DESTROY TIME/FREQUENCY SERIES */
  for(i=0; i<ACS->N_inputparams; i++)
  {
    if(ACS->verbose_switch == 1)fprintf(stdout, "paramid = %d \n", i);
    XLALDestroyREAL4TimeSeries( hderiv[i] ); hderiv[i] = NULL;
    XLALDestroyCOMPLEX8FrequencySeries( Hfderiv[i] ); Hfderiv[i] = NULL;
  }
  
  /* DESTROY FFT PLAN */
	
	/* CLOSE OUTPUT FILES */
	//fclose(fisher_out);
	//fclose(cov_out);
	fclose(derivatives_out);
	fclose(fourierderivs_out);
	
  /*********************************************************************
   *
   *  DETACH LAL ERROR HANDLING
   * 
   ********************************************************************/
  DETATCHSTATUSPTR( status );
  if(ACS->verbose_switch == 1)fprintf(stdout, "Running XLALInspiralComputeFisherMatrixCP - done \n");
  RETURN(status);
	
}

/* <lalVerbatim file="XLALInspiralComputeDerivativesCP">  */
void LALInspiralComputeDerivatives (
																		LALStatus                          *status,
																		REAL4TimeSeries                    *hderiv,
																		PPNConsistencyParamStruc           *PPNparams,
																		INT4                               paramid,
																		FisherACS														*ACS)
/* </lalVerbatim> */
{
  if(ACS->verbose_switch == 1)fprintf(stdout, "Running LALInspiralComputeDerivatives: paramid = %d \n", paramid);
  
  /*********************************************************************
   *
   *  Error handling
   * 
   ********************************************************************/ 
  
  INITSTATUS( status, "LALInspiralComputeDerivatives", LALINSPIRALCOMPUTEFISHERC);
  ATTATCHSTATUSPTR(status);
  
  /*********************************************************************
   *
   *  Check input parameters
   * 
   ********************************************************************/ 
  
  ASSERT( hderiv, status, LALINSPIRALCOMPUTEFISHERC_EINPUT, LALINSPIRALCOMPUTEFISHERC_EINPUT_MSG );
  ASSERT( PPNparams, status, LALINSPIRALCOMPUTEFISHERC_EINPUT, LALINSPIRALCOMPUTEFISHERC_EINPUT_MSG );
  
  /*********************************************************************
   *
   *  Output Files
   * 
   ********************************************************************/ 
  
  if(ACS->verbose_switch==1){fprintf(stdout, "... Creating Output Files \n");}
	FILE *phi_out;
	phi_out = fopen("/home/tjonnie/Programs/lscsoft/ampCorPPN_filter/src/lalsuite/lalinspiral/test/phi.dat","w");
  
  
  /*********************************************************************
   *
   *  Temporary Variables
   * 
   ********************************************************************/ 
  
  if(ACS->verbose_switch==1){fprintf(stdout, "... Initialising Variables \n");}
  
  // Counters
  INT4 i,j;
  
  // Initialise variables for general derivative analysis
  REAL8 absdelta          = 0.0;    // deltaX
  //INT4 ht_len							= 0;
  
  // test
  INT4 icur = 0;
  INT4 iend = 0;
  INT4 deltax_dir = 0;
	REAL8 deltax_max = 0.0;
	REAL8 deltax_min = 0.0;
	INT4 loopcount = 0;
	INT4 deltaL = 0;				// Length difference in differentiation region
	INT4 waveformlength_min = ACS->N_ts;
	INT4 xsamples_tot = ACS->xsamples_tot;
	REAL8 deltaPhi_max = ACS->deltaPhi_max;
	INT4 linReg_countMax = ACS->linReg_countMax;
	
  // Initialise TimeSeries
  REAL4TimeSeries *hderiv_tmp;         // Storing derivative
  REAL4TimeSeries *phiStart;      // start phi
  REAL4TimeSeries *phiEnd;        // end phi
  REAL4TimeSeries *AStart;        // start amplitude
  REAL4TimeSeries *AEnd;        // end amplitude
  
  
	// Initialise parameter structure
	PPNConsistencyParamStruc htPPNparams;
	htPPNparams.ppn = XLALCreateREAL4Vector(PPNparams->ppn->length);
	
	// Initialise waveforms (hxwaveform declared/destroyed in analysis loop)
	CoherentGW htwaveform;
	memset(&htwaveform, 0, sizeof(CoherentGW));
  
  /*********************************************************************
   *
   *  Allocating Memory
   * 
   ********************************************************************/
  
  if(ACS->verbose_switch==1){fprintf(stdout, "... Allocating memory \n");}
  
  // Allocating Memory for TimeSeries
  if(     (     hderiv_tmp           = (REAL4TimeSeries *) LALMalloc( sizeof( REAL4TimeSeries) ) ) == NULL 
		 &&  (phiStart         = (REAL4TimeSeries *) LALMalloc( sizeof( REAL4TimeSeries) ) ) == NULL 
		 &&  (phiEnd         = (REAL4TimeSeries *) LALMalloc( sizeof( REAL4TimeSeries) ) ) == NULL 
		 &&  (AStart         = (REAL4TimeSeries *) LALMalloc( sizeof( REAL4TimeSeries) ) ) == NULL 
		 &&  (AEnd           = (REAL4TimeSeries *) LALMalloc( sizeof( REAL4TimeSeries) ) ) == NULL)
  {
    TRY( LALFree(hderiv_tmp), status);             hderiv_tmp              = NULL;
    TRY( LALFree(phiStart), status);               phiStart              = NULL;
    TRY( LALFree(phiEnd), status);                 phiEnd              = NULL;
    TRY( LALFree(AStart), status);                 AStart              = NULL;
    TRY( LALFree(AEnd), status);                 	 AEnd              = NULL;
    ABORT( status, LALINSPIRALCOMPUTEFISHERC_EMEM, LALINSPIRALCOMPUTEFISHERC_EMEM_MSG);
  }
  else
  {
    // Creating TimeSeries
    if ( XLALCreateREAL4TimeSeries("test", &(ACS->epoch), 0.0, ACS->deltaT_ts, &lalStrainUnit, ACS->N_ts) == NULL)
    {
      if(ACS->verbose_switch == 1)fprintf(stdout, "... Setting variables and allocating memory - Creating TimeSeries -> Error, destroy TimeSeries and Abort \n");
      ABORT( status, LALINSPIRALCOMPUTEFISHERC_ETSCREATE, LALINSPIRALCOMPUTEFISHERC_ETSCREATE_MSG );
    }
    else
    {
      hderiv_tmp              = XLALCreateREAL4TimeSeries("hderiv_tmp", &(ACS->epoch), 0.0, ACS->deltaT_ts, &lalStrainUnit, ACS->N_ts);
      phiStart            = XLALCreateREAL4TimeSeries("phiStart", &(ACS->epoch), 0.0, ACS->deltaT_ts, &lalStrainUnit, ACS->N_ts);
      phiEnd              = XLALCreateREAL4TimeSeries("phiEnd", &(ACS->epoch), 0.0, ACS->deltaT_ts, &lalStrainUnit, ACS->N_ts);
      AStart              = XLALCreateREAL4TimeSeries("AStart", &(ACS->epoch), 0.0, ACS->deltaT_ts, &lalStrainUnit, ACS->N_ts);
      AEnd              = XLALCreateREAL4TimeSeries("AEnd", &(ACS->epoch), 0.0, ACS->deltaT_ts, &lalStrainUnit, ACS->N_ts);
      
			for (j = 0; j < ((INT4) ACS->N_ts); j++)
      {
        hderiv_tmp->data->data[j]              = 0.0;
      }  
    }
    
    // Clearing waveform holding time series
    //if (    hderiv_tmp->data == NULL      &&
		//phiStart->data == NULL    &&
		//phiEnd->data == NULL    &&
		//AStart->data == NULL    &&
		//AEnd->data == NULL          )
    //{
		//fprintf(stdout, "... Setting variables and allocating memory - Creating TimeSeries -> Error, destroy TimeSeries and Abort \n");
		//TRY( XLALDestroyREAL4TimeSeries(hderiv_tmp), status);
		//TRY( XLALDestroyREAL4TimeSeries(phiStart), status);
		//TRY( XLALDestroyREAL4TimeSeries(phiEnd), status);
		//TRY( XLALDestroyREAL4TimeSeries(AStart), status);
		//TRY( XLALDestroyREAL4TimeSeries(AEnd), status);
		//ABORT( status, LALINSPIRALCOMPUTEFISHERC_ETSCREATE, LALINSPIRALCOMPUTEFISHERC_ETSCREATE_MSG);
    //}
    //else
    //{
		
    //} 
  }
	
  if(ACS->verbose_switch==1){fprintf(stdout, "... Allocating memory - done \n");}
	
  /*********************************************************************
   *
   *  Traverse in time to calculate derivative
   * 
   ********************************************************************/
	
	if(ACS->verbose_switch==1){fprintf(stdout, "... Calculating Derivative \n");   }
	
  // Load relevant parameter depending on paramid
  switch(paramid)
  {
    case 0:
      absdelta = PPNparams->mTot_real8*ACS->reldelta;
      break;
    case 1:
      absdelta = PPNparams->eta_real8*ACS->reldelta;    
      break;
    case 2:
      if(PPNparams->phi==0.0){absdelta = 1.0*ACS->reldelta*500.0;}
      else{absdelta = PPNparams->phi*ACS->reldelta*500.0;}
      break;
    case 3:
      absdelta = ACS->deltaT_ts;
      break;
		case 4:
			if(PPNparams->cosI == 0.0) absdelta = 1.0*ACS->reldelta*10.0;
			else {absdelta = PPNparams->cosI*ACS->reldelta*10.0;}
			break;
		case 5:
			absdelta = PPNparams->phasePNparams[ACS->testPhaseParam]*ACS->reldelta;
			break;
    default:
      ABORT( status, LALINSPIRALCOMPUTEFISHERC_EINPUT, LALINSPIRALCOMPUTEFISHERC_EINPUT_MSG );
      break;
  }	
	
  /* COMPUTE DERIVATIVES */
  if(paramid==3)
  {
		printf("--> Computing time derivative\n");
		/* COPY INPUT PARAMETERS */
		PPNConsistencyParamStruc hxCentralPPNparams;
		hxCentralPPNparams.ppn = XLALCreateREAL4Vector(PPNparams->ppn->length);
		TRY( XLALCopyPPNConsistencyParamStruc(PPNparams, &hxCentralPPNparams), status);
		
		/* CREATE WAVEFORM */
		CoherentGW hxCentralWaveform;
		memset(&hxCentralWaveform,0,sizeof(CoherentGW));
		TRY( LALGeneratePPNAmpCorConsistency( status->statusPtr, &hxCentralWaveform, &hxCentralPPNparams), status);
		waveformlength_min = ((INT4) hxCentralWaveform.f->data->length);
		
		
		REAL4TimeSeries *ht_tmp;
		ht_tmp = XLALCreateREAL4TimeSeries("ht_tmp", &(ACS->epoch), 0.0, ACS->deltaT_ts, &lalStrainUnit, ACS->N_ts);
		
		for(j=0;j<ACS->N_ts;j++)
		{
			ht_tmp->data->data[j] = 0.0;
			hderiv_tmp->data->data[j] = 0.0;
		}
		
		LALInspiralCombinePlusCross(status->statusPtr, hxCentralWaveform, hxCentralPPNparams, ht_tmp, ACS);
		
		LALSavitskyGolayFilter(status->statusPtr, 8, 4, ht_tmp, hderiv_tmp, 1);
	}
	else
	{
		while(icur < (ACS->N_ts-1))
		{
		loopcount++;
		if(loopcount == 1){deltax_min = absdelta; deltax_max=absdelta;}
		
		/* RESTART BOUNDS AFTER 25 ITERATIONS */
		if((loopcount%25) == 0)
		{
			
			linReg_countMax += linReg_countMax/10;
			deltaPhi_max *= 1.1;
			
			/* RESET DELTAX BOUNDS */
			deltax_min = absdelta; 
			deltax_max = absdelta;				
		}
		
		// Initialise parameter structure
		PPNConsistencyParamStruc hxStartPPNparams;               // params start differentiation region
		PPNConsistencyParamStruc hxEndPPNparams;                 // params end differentiation region
		hxStartPPNparams.ppn = XLALCreateREAL4Vector(PPNparams->ppn->length);
		hxEndPPNparams.ppn = XLALCreateREAL4Vector(PPNparams->ppn->length);
		TRY( XLALCopyPPNConsistencyParamStruc(PPNparams, &hxStartPPNparams), status);
		TRY( XLALCopyPPNConsistencyParamStruc(PPNparams, &hxEndPPNparams), status);
		
		// Initialise waveforms (hxwaveform declared/destroyed in analysis loop)
		CoherentGW hxStartWaveform;                			// waveform start of differentiation region
		CoherentGW hxEndWaveform;                				// waveform end of differentiation region
		memset(&hxStartWaveform,0,sizeof(CoherentGW));
		memset(&hxEndWaveform, 0, sizeof(CoherentGW));
		
		/* CLEAR TIME SERIES */
		for (j = 0; j < ((INT4) ACS->N_ts); j++)
		{
			//hderiv_tmp->data->data[j]              = 0.0;
			phiStart->data->data[j]            = 0.0;
			phiEnd->data->data[j]              = 0.0;
			AStart->data->data[j]              = 0.0;
			AEnd->data->data[j]              = 0.0;
		}      
		
		/* USE EDGES OF DIFFERENTIATION REGION */
		switch(paramid)
		{
			case 0:
				hxStartPPNparams.mTot_real8 	-= xsamples_tot/2 * absdelta;
				hxEndPPNparams.mTot_real8 		+= xsamples_tot/2 * absdelta;
				
				/* REFILL PHI_I FOR Start */
				if(ACS->testPhaseParam!=0) hxStartPPNparams.phasePNparams[0] = -pow(hxStartPPNparams.eta_real8,-3.0/8.0)*pow(5.0*LAL_MTSUN_SI*hxStartPPNparams.mTot_real8,-5.0/8.0); 
				if(ACS->testPhaseParam!=1) hxStartPPNparams.phasePNparams[1] = -(3715.0/8064.0 + 55.0/96.0*hxStartPPNparams.eta_real8)*pow(hxStartPPNparams.eta_real8,-5.0/8.0)*pow(5.0*LAL_MTSUN_SI*hxStartPPNparams.mTot_real8,-3.0/8.0);
				if(ACS->testPhaseParam!=2) hxStartPPNparams.phasePNparams[2] = 3.0/4.0*LAL_PI*pow(hxStartPPNparams.eta_real8,-0.75)*pow(5.0*LAL_MTSUN_SI*hxStartPPNparams.mTot_real8,-0.25); 
				if(ACS->testPhaseParam!=3) hxStartPPNparams.phasePNparams[3] = -(9275495.0/14450688.0 + 284875.0/258048.0*hxStartPPNparams.eta_real8 + 1855.0/2048.0*pow(hxStartPPNparams.eta_real8,2.0))*pow(hxStartPPNparams.eta_real8,-7.0/8.0)*pow(5.0*LAL_MTSUN_SI*hxStartPPNparams.mTot_real8,-1.0/8.0);
				if(ACS->testPhaseParam!=4) hxStartPPNparams.phasePNparams[4] = -1.0/hxStartPPNparams.eta_real8*(-38645.0/172032.0 + 65.0/2048.0*hxStartPPNparams.eta_real8)*LAL_PI*log(hxStartPPNparams.eta_real8/(5.0*LAL_MTSUN_SI*hxStartPPNparams.mTot_real8));
				if(ACS->testPhaseParam!=5) hxStartPPNparams.phasePNparams[5] = -1.0/hxStartPPNparams.eta_real8*(-38645.0/172032.0 + 65.0/2048.0*hxStartPPNparams.eta_real8)*LAL_PI; 
				if(ACS->testPhaseParam!=6) hxStartPPNparams.phasePNparams[6] = -(831032450749357.0/57682522275840.0 - 53.0/40.0*LAL_PI*LAL_PI - 107.0/56.0*LAL_GAMMA + 107.0/448.0*log(hxStartPPNparams.eta_real8/(256*5.0*LAL_MTSUN_SI*hxStartPPNparams.mTot_real8)) + (-123292747421.0/4161798144.0 + 2255.0/2048.0*LAL_PI*LAL_PI + 385.0/48.0*(-1987.0/3080.0) - 55.0/16.0*(-11831.0/9240.0))*hxStartPPNparams.eta_real8 + 154565.0/1835008.0*pow(hxStartPPNparams.eta_real8,2.0) - 1179625.0/1769472.0*pow(hxStartPPNparams.eta_real8,3.0))*pow(hxStartPPNparams.eta_real8,-9.0/8.0)*pow(5.0*LAL_MTSUN_SI*hxStartPPNparams.mTot_real8,1.0/8.0);
				if(ACS->testPhaseParam!=7) hxStartPPNparams.phasePNparams[7] = -107.0/448.0*pow(hxStartPPNparams.eta_real8,-9.0/8.0)*pow(5.0*LAL_MTSUN_SI*hxStartPPNparams.mTot_real8,1.0/8.0);
				if(ACS->testPhaseParam!=8) hxStartPPNparams.phasePNparams[8] = -(188516689.0/173408256.0 + 488825.0/516096.0*hxStartPPNparams.eta_real8 - 141769.0/516096.0*pow(hxStartPPNparams.eta_real8,2.0))*LAL_PI*pow(hxStartPPNparams.eta_real8,-5.0/4.0)*pow(5.0*LAL_MTSUN_SI*hxStartPPNparams.mTot_real8,1.0/4.0); 		
				
				/* REFILL PHI_I FOR End */
				if(ACS->testPhaseParam!=0) hxEndPPNparams.phasePNparams[0] = -pow(hxEndPPNparams.eta_real8,-3.0/8.0)*pow(5.0*LAL_MTSUN_SI*hxEndPPNparams.mTot_real8,-5.0/8.0); 
				if(ACS->testPhaseParam!=1) hxEndPPNparams.phasePNparams[1] = -(3715.0/8064.0 + 55.0/96.0*hxEndPPNparams.eta_real8)*pow(hxEndPPNparams.eta_real8,-5.0/8.0)*pow(5.0*LAL_MTSUN_SI*hxEndPPNparams.mTot_real8,-3.0/8.0);
				if(ACS->testPhaseParam!=2) hxEndPPNparams.phasePNparams[2] = 3.0/4.0*LAL_PI*pow(hxEndPPNparams.eta_real8,-0.75)*pow(5.0*LAL_MTSUN_SI*hxEndPPNparams.mTot_real8,-0.25); 
				if(ACS->testPhaseParam!=3) hxEndPPNparams.phasePNparams[3] = -(9275495.0/14450688.0 + 284875.0/258048.0*hxEndPPNparams.eta_real8 + 1855.0/2048.0*pow(hxEndPPNparams.eta_real8,2.0))*pow(hxEndPPNparams.eta_real8,-7.0/8.0)*pow(5.0*LAL_MTSUN_SI*hxEndPPNparams.mTot_real8,-1.0/8.0);
				if(ACS->testPhaseParam!=4) hxEndPPNparams.phasePNparams[4] = -1.0/hxEndPPNparams.eta_real8*(-38645.0/172032.0 + 65.0/2048.0*hxEndPPNparams.eta_real8)*LAL_PI*log(hxEndPPNparams.eta_real8/(5.0*LAL_MTSUN_SI*hxEndPPNparams.mTot_real8));
				if(ACS->testPhaseParam!=5) hxEndPPNparams.phasePNparams[5] = -1.0/hxEndPPNparams.eta_real8*(-38645.0/172032.0 + 65.0/2048.0*hxEndPPNparams.eta_real8)*LAL_PI; 
				if(ACS->testPhaseParam!=6) hxEndPPNparams.phasePNparams[6] = -(831032450749357.0/57682522275840.0 - 53.0/40.0*LAL_PI*LAL_PI - 107.0/56.0*LAL_GAMMA + 107.0/448.0*log(hxEndPPNparams.eta_real8/(256*5.0*LAL_MTSUN_SI*hxEndPPNparams.mTot_real8)) + (-123292747421.0/4161798144.0 + 2255.0/2048.0*LAL_PI*LAL_PI + 385.0/48.0*(-1987.0/3080.0) - 55.0/16.0*(-11831.0/9240.0))*hxEndPPNparams.eta_real8 + 154565.0/1835008.0*pow(hxEndPPNparams.eta_real8,2.0) - 1179625.0/1769472.0*pow(hxEndPPNparams.eta_real8,3.0))*pow(hxEndPPNparams.eta_real8,-9.0/8.0)*pow(5.0*LAL_MTSUN_SI*hxEndPPNparams.mTot_real8,1.0/8.0);
				if(ACS->testPhaseParam!=7) hxEndPPNparams.phasePNparams[7] = -107.0/448.0*pow(hxEndPPNparams.eta_real8,-9.0/8.0)*pow(5.0*LAL_MTSUN_SI*hxEndPPNparams.mTot_real8,1.0/8.0);
				if(ACS->testPhaseParam!=8) hxEndPPNparams.phasePNparams[8] = -(188516689.0/173408256.0 + 488825.0/516096.0*hxEndPPNparams.eta_real8 - 141769.0/516096.0*pow(hxEndPPNparams.eta_real8,2.0))*LAL_PI*pow(hxEndPPNparams.eta_real8,-5.0/4.0)*pow(5.0*LAL_MTSUN_SI*hxEndPPNparams.mTot_real8,1.0/4.0); 					
				break;
				
			case 1:
				hxStartPPNparams.eta_real8 	-= xsamples_tot/2 * absdelta;
				hxEndPPNparams.eta_real8 		+= xsamples_tot/2 * absdelta;
				
				/* REFILL PHI_I FOR Start */
				if(ACS->testPhaseParam!=0) hxStartPPNparams.phasePNparams[0] = -pow(hxStartPPNparams.eta_real8,-3.0/8.0)*pow(5.0*LAL_MTSUN_SI*hxStartPPNparams.mTot_real8,-5.0/8.0); 
				if(ACS->testPhaseParam!=1) hxStartPPNparams.phasePNparams[1] = -(3715.0/8064.0 + 55.0/96.0*hxStartPPNparams.eta_real8)*pow(hxStartPPNparams.eta_real8,-5.0/8.0)*pow(5.0*LAL_MTSUN_SI*hxStartPPNparams.mTot_real8,-3.0/8.0);
				if(ACS->testPhaseParam!=2) hxStartPPNparams.phasePNparams[2] = 3.0/4.0*LAL_PI*pow(hxStartPPNparams.eta_real8,-0.75)*pow(5.0*LAL_MTSUN_SI*hxStartPPNparams.mTot_real8,-0.25); 
				if(ACS->testPhaseParam!=3) hxStartPPNparams.phasePNparams[3] = -(9275495.0/14450688.0 + 284875.0/258048.0*hxStartPPNparams.eta_real8 + 1855.0/2048.0*pow(hxStartPPNparams.eta_real8,2.0))*pow(hxStartPPNparams.eta_real8,-7.0/8.0)*pow(5.0*LAL_MTSUN_SI*hxStartPPNparams.mTot_real8,-1.0/8.0);
				if(ACS->testPhaseParam!=4) hxStartPPNparams.phasePNparams[4] = -1.0/hxStartPPNparams.eta_real8*(-38645.0/172032.0 + 65.0/2048.0*hxStartPPNparams.eta_real8)*LAL_PI*log(hxStartPPNparams.eta_real8/(5.0*LAL_MTSUN_SI*hxStartPPNparams.mTot_real8));
				if(ACS->testPhaseParam!=5) hxStartPPNparams.phasePNparams[5] = -1.0/hxStartPPNparams.eta_real8*(-38645.0/172032.0 + 65.0/2048.0*hxStartPPNparams.eta_real8)*LAL_PI; 
				if(ACS->testPhaseParam!=6) hxStartPPNparams.phasePNparams[6] = -(831032450749357.0/57682522275840.0 - 53.0/40.0*LAL_PI*LAL_PI - 107.0/56.0*LAL_GAMMA + 107.0/448.0*log(hxStartPPNparams.eta_real8/(256*5.0*LAL_MTSUN_SI*hxStartPPNparams.mTot_real8)) + (-123292747421.0/4161798144.0 + 2255.0/2048.0*LAL_PI*LAL_PI + 385.0/48.0*(-1987.0/3080.0) - 55.0/16.0*(-11831.0/9240.0))*hxStartPPNparams.eta_real8 + 154565.0/1835008.0*pow(hxStartPPNparams.eta_real8,2.0) - 1179625.0/1769472.0*pow(hxStartPPNparams.eta_real8,3.0))*pow(hxStartPPNparams.eta_real8,-9.0/8.0)*pow(5.0*LAL_MTSUN_SI*hxStartPPNparams.mTot_real8,1.0/8.0);
				if(ACS->testPhaseParam!=7) hxStartPPNparams.phasePNparams[7] = -107.0/448.0*pow(hxStartPPNparams.eta_real8,-9.0/8.0)*pow(5.0*LAL_MTSUN_SI*hxStartPPNparams.mTot_real8,1.0/8.0);
				if(ACS->testPhaseParam!=8) hxStartPPNparams.phasePNparams[8] = -(188516689.0/173408256.0 + 488825.0/516096.0*hxStartPPNparams.eta_real8 - 141769.0/516096.0*pow(hxStartPPNparams.eta_real8,2.0))*LAL_PI*pow(hxStartPPNparams.eta_real8,-5.0/4.0)*pow(5.0*LAL_MTSUN_SI*hxStartPPNparams.mTot_real8,1.0/4.0); 		
				
				/* REFILL PHI_I FOR End */
				if(ACS->testPhaseParam!=0) hxEndPPNparams.phasePNparams[0] = -pow(hxEndPPNparams.eta_real8,-3.0/8.0)*pow(5.0*LAL_MTSUN_SI*hxEndPPNparams.mTot_real8,-5.0/8.0); 
				if(ACS->testPhaseParam!=1) hxEndPPNparams.phasePNparams[1] = -(3715.0/8064.0 + 55.0/96.0*hxEndPPNparams.eta_real8)*pow(hxEndPPNparams.eta_real8,-5.0/8.0)*pow(5.0*LAL_MTSUN_SI*hxEndPPNparams.mTot_real8,-3.0/8.0);
				if(ACS->testPhaseParam!=2) hxEndPPNparams.phasePNparams[2] = 3.0/4.0*LAL_PI*pow(hxEndPPNparams.eta_real8,-0.75)*pow(5.0*LAL_MTSUN_SI*hxEndPPNparams.mTot_real8,-0.25); 
				if(ACS->testPhaseParam!=3) hxEndPPNparams.phasePNparams[3] = -(9275495.0/14450688.0 + 284875.0/258048.0*hxEndPPNparams.eta_real8 + 1855.0/2048.0*pow(hxEndPPNparams.eta_real8,2.0))*pow(hxEndPPNparams.eta_real8,-7.0/8.0)*pow(5.0*LAL_MTSUN_SI*hxEndPPNparams.mTot_real8,-1.0/8.0);
				if(ACS->testPhaseParam!=4) hxEndPPNparams.phasePNparams[4] = -1.0/hxEndPPNparams.eta_real8*(-38645.0/172032.0 + 65.0/2048.0*hxEndPPNparams.eta_real8)*LAL_PI*log(hxEndPPNparams.eta_real8/(5.0*LAL_MTSUN_SI*hxEndPPNparams.mTot_real8));
				if(ACS->testPhaseParam!=5) hxEndPPNparams.phasePNparams[5] = -1.0/hxEndPPNparams.eta_real8*(-38645.0/172032.0 + 65.0/2048.0*hxEndPPNparams.eta_real8)*LAL_PI; 
				if(ACS->testPhaseParam!=6) hxEndPPNparams.phasePNparams[6] = -(831032450749357.0/57682522275840.0 - 53.0/40.0*LAL_PI*LAL_PI - 107.0/56.0*LAL_GAMMA + 107.0/448.0*log(hxEndPPNparams.eta_real8/(256*5.0*LAL_MTSUN_SI*hxEndPPNparams.mTot_real8)) + (-123292747421.0/4161798144.0 + 2255.0/2048.0*LAL_PI*LAL_PI + 385.0/48.0*(-1987.0/3080.0) - 55.0/16.0*(-11831.0/9240.0))*hxEndPPNparams.eta_real8 + 154565.0/1835008.0*pow(hxEndPPNparams.eta_real8,2.0) - 1179625.0/1769472.0*pow(hxEndPPNparams.eta_real8,3.0))*pow(hxEndPPNparams.eta_real8,-9.0/8.0)*pow(5.0*LAL_MTSUN_SI*hxEndPPNparams.mTot_real8,1.0/8.0);
				if(ACS->testPhaseParam!=7) hxEndPPNparams.phasePNparams[7] = -107.0/448.0*pow(hxEndPPNparams.eta_real8,-9.0/8.0)*pow(5.0*LAL_MTSUN_SI*hxEndPPNparams.mTot_real8,1.0/8.0);
				if(ACS->testPhaseParam!=8) hxEndPPNparams.phasePNparams[8] = -(188516689.0/173408256.0 + 488825.0/516096.0*hxEndPPNparams.eta_real8 - 141769.0/516096.0*pow(hxEndPPNparams.eta_real8,2.0))*LAL_PI*pow(hxEndPPNparams.eta_real8,-5.0/4.0)*pow(5.0*LAL_MTSUN_SI*hxEndPPNparams.mTot_real8,1.0/4.0); 					
				break;
				
				//case 2:
				//break;
			case 2:
				hxStartPPNparams.phi 	-= xsamples_tot/2 * absdelta;
				hxEndPPNparams.phi 		+= xsamples_tot/2 * absdelta;				
				break;
			case 3:
				break;
			case 4:
				hxStartPPNparams.cosI 	-= xsamples_tot/2 * absdelta;
				hxEndPPNparams.cosI 		+= xsamples_tot/2 * absdelta;
				
				/* MIRROR AROUND -1 IF CROSSING -1 */
				if(hxStartPPNparams.cosI < -1.0 ){ hxStartPPNparams.cosI = -2.0 - hxStartPPNparams.cosI; hxStartPPNparams.sinI = -hxStartPPNparams.sinI;}
				
				/* MIRROR AROUND +1 IF CROSSING +1 */
				if(hxEndPPNparams.cosI > 1.0 ){ hxEndPPNparams.cosI = 2.0 - hxEndPPNparams.cosI; hxEndPPNparams.sinI = -hxEndPPNparams.sinI;}
				
				break;	
				
			case 5:
				hxStartPPNparams.phasePNparams[ACS->testPhaseParam] 	-= xsamples_tot/2 * absdelta;
				hxEndPPNparams.phasePNparams[ACS->testPhaseParam] 		+= xsamples_tot/2 * absdelta;
				break;
				
			default:
				break;
		}
		
		/* COMPUTE WAVEFORMS, NO NEED TO COMBINE HP AND HC AS ONLY PHASE IS REQUIRED */
		TRY( LALGeneratePPNAmpCorConsistency( status->statusPtr, &hxStartWaveform, &hxStartPPNparams), status);
		TRY( LALGeneratePPNAmpCorConsistency( status->statusPtr, &hxEndWaveform, &hxEndPPNparams), status);
		
		/* LOAD PHASE INTO TIMESERIES */
		ASSERT( (((INT4)hxStartWaveform.phi->data->length) < ACS->N_ts), status, LALINSPIRALCOMPUTEFISHERC_ETSSIZE, LALINSPIRALCOMPUTEFISHERC_ETSSIZE_MSG );
		//printf("%d | %d \n", ACS->N_ts, hxStartWaveform.phi->data->length);
		for(j=0;j<((INT4)hxStartWaveform.phi->data->length);j++)
		{
			phiStart->data->data[ACS->N_ts-1-j] = hxStartWaveform.phi->data->data[hxStartWaveform.phi->data->length-1-j];
			AStart->data->data[ACS->N_ts-1-j] = hxStartWaveform.h->data->data[2*(hxStartWaveform.phi->data->length-1-j)];
		}
		if(paramid==3){for(j=0;j<((INT4)hxEndWaveform.phi->data->length)-xsamples_tot;j++)
		{
			phiEnd->data->data[(ACS->N_ts-hxEndWaveform.phi->data->length)+j+xsamples_tot] = hxEndWaveform.phi->data->data[j];
			AEnd->data->data[(ACS->N_ts-hxEndWaveform.phi->data->length)+j+xsamples_tot] = hxEndWaveform.h->data->data[2*j];
		}}
		else{for(j=0;j<((INT4)hxEndWaveform.phi->data->length);j++)
		{
			phiEnd->data->data[ACS->N_ts-1-j] = hxEndWaveform.phi->data->data[hxEndWaveform.phi->data->length-1-j];
			AEnd->data->data[ACS->N_ts-1-j] = hxEndWaveform.h->data->data[2*(hxEndWaveform.phi->data->length-1-j)];
		}}
		
		/* COMPUTE LENGTH DIFFERENCE */
		deltaL = abs(hxStartWaveform.f->data->length-hxEndWaveform.f->data->length);
		
		
		/*if(nprints<3){
		 for(j=0;j<ACS->N_ts;j++) {
		 fprintf(phi_out, "%d\t%e\n", j, phiEnd->data->data[j]-phiStart->data->data[j] ); nprints++;
		 }}*/
		
		/* Test region in time for which current deltax is suitable */
		for(i=icur;i<ACS->N_ts;i++)
		{
			/* PRINT CURRENT STATUS */
			if(ACS->verbose_switch==1)fprintf(stdout, "\r i = %d | icur = %d | deltax = %e |  deltaA = %e | deltaL = %d | dPhi = %e | dPhi_max = %e ", i, icur,  absdelta, fabs(1.0-AStart->data->data[i]/AEnd->data->data[i]),deltaL, fabs(phiEnd->data->data[i]-phiStart->data->data[i]), deltaPhi_max);
			
			//fprintf(phi_out, "%d\t%e\n", i, phiEnd->data->data[i]-phiStart->data->data[i] );
			
			/* CHECK IF AMOUNT OF EDGES IS HIGHER THAN MINIMUM) */
			/* TOTAL LENGTH OF THE WAVEFORM SHOULD BE 2*ACS->N_TS INSTEAD OF ACS->N_TS -> CHANGE INTO ASTART/AEND LATER */
			if( 	(( (fabs(1.0-AStart->data->data[i]/AEnd->data->data[i]) < 1.0E-7) && (AStart->data->data[i] != 0.0 && AEnd->data->data[i] != 0.0))
						 && ( (fabs(phiStart->data->data[i]-phiEnd->data->data[i]) == 0.0) && (phiStart->data->data[i]!=0.0) && (phiEnd->data->data[i]!=0.0) ) )
				 || ((deltaL < ACS->linReg_countMin) && (paramid != 2 && paramid != 3 && paramid != 4 && paramid != 5) ) 
				 
				 )
			{
				if (ACS->verbose_switch==1)
				{
					if ( (fabs(1.0-AStart->data->data[i]/AEnd->data->data[i]) < 1.0E-7) && (AStart->data->data[i] != 0.0 && AEnd->data->data[i] != 0.0)) if(ACS->verbose_switch == 1)fprintf(stdout, " delta_Amp = 0");
					if ( (fabs(phiStart->data->data[i]-phiEnd->data->data[i]) == 0.0) && (phiStart->data->data[i]!=0.0) && (phiEnd->data->data[i]!=0.0) ) if(ACS->verbose_switch == 1)fprintf(stdout, "| delta_Phi = 0 \n");
					if ( (deltaL < ACS->linReg_countMin) && (paramid != 2 && paramid != 3 && paramid != 4 && paramid != 5) ) if(ACS->verbose_switch == 1)fprintf(stdout, "| deltaL < deltaL_min \n");
				}
				deltax_dir = 1;
				if( ((i-icur)>(ACS->N_ts/10)) && (icur<ACS->N_ts-ACS->N_ts/10)) iend = i;
				break;
			}				
			
			/* CHECK IF DELTA_PHI IS SMALL ENOUGH AND IF AMOUNT OF EDGES IS SMALLER THAN MAXIMUM */
			if( ((deltaL > linReg_countMax) && (paramid != 2 && paramid != 3 && paramid != 4 && paramid != 5)) || ((fabs(phiStart->data->data[i]-phiEnd->data->data[i]) > deltaPhi_max) && (phiStart->data->data[i]!=0.0) && (phiEnd->data->data[i]!=0.0) ) )
			{
				if((fabs(phiStart->data->data[i]-phiEnd->data->data[i]) > deltaPhi_max) && (phiStart->data->data[i]!=0.0) && (phiEnd->data->data[i]!=0.0)) if(ACS->verbose_switch == 1)fprintf(stdout, "| deltaPhi > deltaPhi_Max \n");
				if( (deltaL > linReg_countMax) && (paramid != 2 && paramid != 3 && paramid != 4 && paramid != 5) )if(ACS->verbose_switch==1)if(ACS->verbose_switch == 1)fprintf(stdout, "| deltaL > deltaL_Max \n");
				deltax_dir = 2;
				if( ((i-icur)>(ACS->N_ts/10)) && (icur<ACS->N_ts-ACS->N_ts/10)) iend = i;
				break;
			}
			if(i==ACS->N_ts-1) if(ACS->verbose_switch==1)fprintf(stdout, "| end reached \n");
		}
		
		/* Calculate derivative for suitable region */
		//if(icur != i && (((i-icur)>(ACS->N_ts/10)) || (icur>(ACS->N_ts-ACS->N_ts/10))) )
		if(i>=ACS->N_ts-1 )
		{
			/* KEEP TRACK OF MINIMAL WAVELENGTH */
			if(waveformlength_min > ((INT4)hxEndWaveform.f->data->length)) waveformlength_min = hxEndWaveform.f->data->length;
			if(waveformlength_min > ((INT4)hxStartWaveform.f->data->length)) waveformlength_min = hxStartWaveform.f->data->length;
			
			/* CALCULATE DERIVATIVE */
			LALInspiralComputeDerivatives_3PolyChi2(status->statusPtr, PPNparams, absdelta, paramid, icur, i, hderiv_tmp, xsamples_tot, ACS);
			icur = i;
			loopcount = 0;
			
			/* RESET DELTAX BOUNDS */
			deltax_min = absdelta; 
			deltax_max = absdelta;	
		}
		else if(icur != i) {
			if(ACS->verbose_switch == 1) fprintf(stdout, "... Differentiable region too small, retry \n");
			//if(i > 600)for(j=0;j<ACS->N_ts;j++) fprintf(phi_out,"%d\t%e\t%e\t%e \n", j, phiStart->data->data[j], phiEnd->data->data[j], fabs(phiEnd->data->data[j]-phiStart->data->data[j]));
			//if(i > 600)break;
		}
		//fprintf(traject_out, "%d\t%e\t%e\t%e \n",loopcount, absdelta, deltax_min, deltax_max);
		
		/* Update deltax */
		switch(deltax_dir)
		{
			case 1:
				if(paramid == 3){xsamples_tot = xsamples_tot*ACS->xsamples_adjustUp;}
				else
				{
					if(deltax_max == absdelta) {deltax_min = absdelta; deltax_max = absdelta * ACS->absdelta_adjustUp; absdelta *= ACS->absdelta_adjustUp;}
					else{deltax_min = absdelta; absdelta = (deltax_max + deltax_min)/2.0;}
				}
				break;
				
			case 2:
				if(paramid == 3)
				{
					xsamples_tot -= xsamples_tot/ACS->xsamples_adjustDown;
					if(xsamples_tot<10){ LALInspiralComputeDerivatives_5point(status->statusPtr, PPNparams, paramid, absdelta, icur, ((INT4)hxEndWaveform.phi->data->length)-xsamples_tot, hderiv_tmp, ACS);}
				}
				else 
				{
					if(deltax_min == absdelta){deltax_max = absdelta; deltax_min = absdelta*ACS->absdelta_adjustDown; absdelta *= ACS->absdelta_adjustDown;}
					else{deltax_max = absdelta; absdelta = (deltax_max+deltax_min)/2.0;}
				}
				break;
				
			default:
				if(ACS->verbose_switch == 1)fprintf(stdout, "... no directional change \n");
				break;
		}
		
		if(xsamples_tot<10) break;
		
		if(loopcount == 5E3)
		{
			LALInspiralComputeDerivatives_5point(status->statusPtr, PPNparams, paramid, absdelta, icur, (INT4)hxEndWaveform.phi->data->length, hderiv_tmp, ACS); 
			break;
		}
		
		/* CLEANING UP */
		TRY( XLALDestroyREAL4Vector( hxStartPPNparams.ppn ), status);
		TRY( XLALDestroyREAL4Vector( hxEndPPNparams.ppn ), status);
		TRY( XLALDestroyCoherentGW(&hxStartWaveform), status); 
		TRY( XLALDestroyCoherentGW(&hxEndWaveform), status); 
		
	}
	}
	/* FOR LOG DL, DERIVATIVE IS H(T) AGAIN */
  //else if(paramid == 2)
  //{
	//TRY( XLALCopyPPNConsistencyParamStruc(PPNparams, &htPPNparams), status);
	//TRY( LALGeneratePPNAmpCorConsistency( status->statusPtr, &htwaveform, &htPPNparams), status);
	//LALInspiralCombinePlusCross(status->statusPtr,htwaveform, htPPNparams, hderiv_tmp, ACS);
	////ht_len = htwaveform.f->data->length;
	////for(j=0;j<((INT4)htwaveform.f->data->length);j++){hderiv_tmp->data->data[(ACS->N_ts-htwaveform.f->data->length)+j] = htwaveform.h->data->data[2*j];}		
	//}
  
  for(i=1;i<(waveformlength_min+1);i++)
  {
		switch(paramid)
		{
			case 0:
				hderiv->data->data[ACS->N_ts-i] = hderiv_tmp->data->data[ACS->N_ts-i];
				break;
			case 1:
				hderiv->data->data[ACS->N_ts-i] = hderiv_tmp->data->data[ACS->N_ts-i];
				break;
				//case 2:
				//hderiv->data->data[ACS->N_ts-i] = -1.0*hderiv_tmp->data->data[ACS->N_ts-i];
				//break;
			case 2:
				hderiv->data->data[ACS->N_ts-i] = hderiv_tmp->data->data[ACS->N_ts-i];
				break;
			case 3:
				hderiv->data->data[ACS->N_ts-i] = hderiv_tmp->data->data[ACS->N_ts-i];
        break;
			case 4:
				hderiv->data->data[ACS->N_ts-i] = hderiv_tmp->data->data[ACS->N_ts-i];
				break;
				
			case 5:
				hderiv->data->data[ACS->N_ts-i] = hderiv_tmp->data->data[ACS->N_ts-i];
				break;
				
			default:
				if(ACS->verbose_switch == 1)fprintf(stdout, "XLALInspiralComputeDerivatives failed: invalid derivative variable \n");
				break;
		}    
	}
  
  if(ACS->verbose_switch==1){fprintf(stdout, "... Calculating Derivative - done \n");   }
  
  /*********************************************************************
   *
   *  Clean up
   * 
   ********************************************************************/
  if(ACS->verbose_switch==1){fprintf(stdout, "... Cleaning Up \n");}
  
  // OUTPUT FILES
  //fclose(hx_out);
  //fclose(traject_out);
  //fclose(linReg_out);
	fclose(phi_out);
	
	// TIMESERIES
  TRY( XLALDestroyREAL4TimeSeries(hderiv_tmp), status);               hderiv_tmp              = NULL;
  TRY( XLALDestroyREAL4TimeSeries(phiStart), status);               phiStart              = NULL;
  TRY( XLALDestroyREAL4TimeSeries(phiEnd), status);        phiEnd       = NULL;
  TRY( XLALDestroyREAL4TimeSeries(AStart), status);        AStart       = NULL;
  TRY( XLALDestroyREAL4TimeSeries(AEnd), status);        AEnd       = NULL;
  
  // PARAMSTRUCS
  TRY( XLALDestroyREAL4Vector( htPPNparams.ppn ), status);
  
  /* COHERENT GW'S */ 
  TRY( XLALDestroyCoherentGW(&htwaveform), status); 
  
  /*********************************************************************
   *
   *  Detach Error Handling
   * 
   ********************************************************************/
  DETATCHSTATUSPTR( status );
  if(ACS->verbose_switch == 1)fprintf(stdout, "Running LALInspiralComputeDerivatives: paramid = %d - done \n", paramid);
  RETURN(status);
  
}

/* <lalVerbatim file="LALInspiralComputeFisherComponentsCP">  */
void LALInspiralComputeFisherComponents (
																				 LALStatus															 *status,
																				 REAL4					                         *Comp,
																				 REAL4TimeSeries                        *hderiv1,
																				 REAL4TimeSeries                        *hderiv2,
																				 REAL8FrequencySeries                   *psd,
																				 UINT4                                  compid,
																				 FisherACS															 *ACS)
/* </lalVerbatim> */
{
	if(ACS->verbose_switch == 1) fprintf(stdout, "Running XLALInspiralComputeFisherComponents: compid = %d \n", compid);
	
  /*********************************************************************
   *
   *  Error handling
   * 
   ********************************************************************/
	
  INITSTATUS( status, "LALInspiralComputeFisherComponents", LALINSPIRALCOMPUTEFISHERC);
  ATTATCHSTATUSPTR(status);
  
  /*********************************************************************
   *
   *  CREATING TEMPORARY VARIABLES
   * 
   ********************************************************************/  
  
  if(ACS->verbose_switch == 1) fprintf(stdout, "... Setting variables and allocating memory \n");
	
  /* COUNTERS */
  INT4 i, k;
  REAL8 invpsd = 0.0;
  
	/* TIME/FREQUENCY SERIES */
  COMPLEX8FrequencySeries *Hfderiv1;
  COMPLEX8FrequencySeries *Hfderiv2;
  COMPLEX8FrequencySeries *integrand; 
  REAL4TimeSeries *temp;
  
  if( 	(Hfderiv1 = (COMPLEX8FrequencySeries *) LALMalloc( sizeof(COMPLEX8FrequencySeries) )) == NULL
		 &&	(Hfderiv2 = (COMPLEX8FrequencySeries *) LALMalloc( sizeof(COMPLEX8FrequencySeries) )) == NULL
		 &&	(integrand = (COMPLEX8FrequencySeries *) LALMalloc( sizeof(COMPLEX8FrequencySeries) )) == NULL
		 &&	(temp = (REAL4TimeSeries *) LALMalloc( sizeof(REAL4TimeSeries) ) ) == NULL )
	{
		TRY( LALFree(Hfderiv1), status);		Hfderiv1	= NULL;
		TRY( LALFree(Hfderiv1), status);		Hfderiv1	= NULL;
		TRY( LALFree(integrand), status);		integrand	= NULL;
		TRY( LALFree(temp), status);		temp	= NULL;
		ABORT( status, LALINSPIRALCOMPUTEFISHERC_EMEM, LALINSPIRALCOMPUTEFISHERC_EMEM_MSG);
	}
	else
	{
		Hfderiv1 = XLALCreateCOMPLEX8FrequencySeries("Hfderiv1", &(ACS->epoch), 0.0, ACS->deltaF_fs, &lalStrainUnit, ACS->N_fs);
		Hfderiv2 = XLALCreateCOMPLEX8FrequencySeries("Hfderiv2", &(ACS->epoch), 0.0, ACS->deltaF_fs, &lalStrainUnit, ACS->N_fs);
		integrand = XLALCreateCOMPLEX8FrequencySeries("integrand", &(ACS->epoch), 0.0, ACS->deltaF_fs, &lalStrainUnit, ACS->N_fs);
		temp = XLALCreateREAL4TimeSeries("temp", &(ACS->epoch), 0.0, ACS->deltaT_ts, &lalStrainUnit, ACS->N_ts);
	}
	
	for(i=0;i<ACS->N_ts;i++)
	{
		if(i<ACS->N_fs)
		{
			Hfderiv1->data->data[i].re *= 0.0;
			Hfderiv1->data->data[i].im *= 0.0;
			Hfderiv2->data->data[i].re *= 0.0;
			Hfderiv2->data->data[i].im *= 0.0;
			integrand->data->data[i].re *= 0.0;
			integrand->data->data[i].im *= 0.0;
		}
		else
		{
			temp->data->data[i] = 0.0;
		}
		
	}
	
  /* FFT PLANS */
  RealFFTPlan *fwdRealPlan1 = NULL;
  RealFFTPlan *fwdRealPlan2 = NULL;
  RealFFTPlan *revRealPlan  = NULL;	
	
  /*********************************************************************
   *
   *  COMPUTING FISHER COMPONENTS USING TWO DERIVATIVES
   * 		- TAKING FFT OF DERIVATIVES
   * 		- PERFORM SCALAR PRODUCT
   * 
   ********************************************************************/
  
  if(ACS->verbose_switch == 1) fprintf(stdout, "... Computing Fisher component \n");
	
  /* Take FFT of derivatives */ 
  
  //Hfderiv1.data = NULL;
  //Hfderiv2.data = NULL;
  //LALCCreateVector( status->statusPtr, &Hfderiv1.data, ACS->N_fs );
  //LALCCreateVector( status->statusPtr, &Hfderiv2.data, ACS->N_fs );
	
  LALCreateForwardRealFFTPlan( status->statusPtr, &fwdRealPlan1, ACS->N_ts, 0 ); 
  LALCreateForwardRealFFTPlan( status->statusPtr, &fwdRealPlan2, ACS->N_ts, 0 ); 
	
  hderiv1->f0 = 0.0;
  hderiv1->deltaT = ACS->deltaT_ts;
	
  //Hfderiv1.deltaF = ACS->deltaF_fs;
	
  hderiv2->f0 = 0.0;
  hderiv2->deltaT = ACS->deltaT_ts;
	
  //Hfderiv2.deltaF = ACS->deltaF_fs;
  
  LALTimeFreqRealFFT( status->statusPtr, Hfderiv1, hderiv1, fwdRealPlan1 );
  LALTimeFreqRealFFT( status->statusPtr, Hfderiv2, hderiv2, fwdRealPlan2 );
  
  /* Rescale by a dynamical range */
  for (i = 0; i < ACS->N_fs; ++i)
  {
    Hfderiv1->data->data[i].re *= ACS->dynRange;
    Hfderiv1->data->data[i].im *= ACS->dynRange;    
		
    Hfderiv2->data->data[i].re *= ACS->dynRange;
    Hfderiv2->data->data[i].im *= ACS->dynRange;    
  }
	
  /* Compute Fisher matrix component */
  //integrand.data = NULL;
  //temp.data = NULL;
  //LALCCreateVector( status->statusPtr, &integrand.data, ACS->N_fs );
  //LALSCreateVector( status->statusPtr, &temp.data, ACS->N_ts );
  
  LALCreateReverseRealFFTPlan( status->statusPtr, &revRealPlan, ACS->N_ts, 0 ); 
	
  for (k = 0; k < ACS->N_fs; k++)
  {
		// check PSD and set integration limits (40.0 - LSO)
    if (psd->data->data[k] == 0.0 || ((REAL8) k)*ACS->deltaF_fs < ACS->fstart || ((REAL8) k)*ACS->deltaF_fs > ACS->fstop	)
      invpsd = 0.0;
    else
      invpsd = 1.0 / ( ACS->dynRange * ACS->dynRange * psd->data->data[k] * 2.0 );
		
    integrand->data->data[k].re = (Hfderiv1->data->data[k].re * Hfderiv2->data->data[k].re + Hfderiv1->data->data[k].im * Hfderiv2->data->data[k].im) * invpsd;
    integrand->data->data[k].im = 0.0;
  }
	
  //integrand.deltaF = ACS->deltaF_fs;
	
  LALFreqTimeRealFFT( status->statusPtr, temp, integrand, revRealPlan );
	
  Comp[compid] = 4.0*temp->data->data[0]; 
	
  /*********************************************************************
   *
   *  CLEANING UP
   * 
   ********************************************************************/
  
  if(ACS->verbose_switch == 1) fprintf(stdout, "... Cleaning up \n");
  
  /* VECTORS */
  XLALDestroyCOMPLEX8FrequencySeries(Hfderiv1);		Hfderiv1 = NULL;
  XLALDestroyCOMPLEX8FrequencySeries(Hfderiv2);		Hfderiv2 = NULL;
  XLALDestroyCOMPLEX8FrequencySeries(integrand);	Hfderiv2 = NULL;
  XLALDestroyREAL4TimeSeries(temp);								temp = NULL;
  
  //LALCDestroyVector(status->statusPtr, &Hfderiv1.data);
  //LALCDestroyVector(status->statusPtr, &Hfderiv2.data);
  //LALCDestroyVector(status->statusPtr, &integrand.data);
  //LALSDestroyVector(status->statusPtr, &temp.data);
  
  /* FFT PLANS */
  LALDestroyRealFFTPlan( status->statusPtr, &fwdRealPlan1); 
  LALDestroyRealFFTPlan( status->statusPtr, &fwdRealPlan2); 
  LALDestroyRealFFTPlan( status->statusPtr, &revRealPlan); 
	
  /*********************************************************************
   *
   *  Detach Error handling
   * 
   ********************************************************************/
  
  DETATCHSTATUSPTR( status );
  RETURN(status);
	
}

/* <lalVerbatim file="LALInspiralInvertMatrixCP">  */
void LALInspiralInvertMatrix (
															LALStatus															*status,
															REAL4                                 *Input,
															REAL4                                 *Inverse,
															FisherACS															*ACS)
/* </lalVerbatim> */
{ 
	if(ACS->verbose_switch == 1) fprintf(stdout, "Running LALInspiralInvertMatrix \n");
	
  /*********************************************************************
   *
   *  Error handling
   * 
   ********************************************************************/
	
  INITSTATUS( status, "LALInspiralInvertMatrix", LALINSPIRALCOMPUTEFISHERC);
  ATTATCHSTATUSPTR(status);
	
  /*********************************************************************
   *
   *  CREATING TEMPORARY VARIABLES
   * 
   ********************************************************************/  
	
  if(ACS->verbose_switch == 1) fprintf(stdout, "... Setting Variables and Allocating Memory \n");
	
	/* COUNTERS */
  INT4 i = 0;
  
  // Initialise matrices A, U, S and V
  gsl_matrix *A     = gsl_matrix_calloc (ACS->N_inputparams, ACS->N_inputparams); // Input matrix
  gsl_matrix *InvS  = gsl_matrix_calloc (ACS->N_inputparams, ACS->N_inputparams); // inverse S
  gsl_matrix *V     = gsl_matrix_calloc (ACS->N_inputparams, ACS->N_inputparams); // V
  gsl_matrix *U     = gsl_matrix_calloc (ACS->N_inputparams, ACS->N_inputparams); // U
  gsl_matrix *InvA  = gsl_matrix_calloc (ACS->N_inputparams, ACS->N_inputparams); // Inverse A
  gsl_matrix *C     = gsl_matrix_calloc (ACS->N_inputparams, ACS->N_inputparams); // temporary storage
  gsl_matrix *I     = gsl_matrix_calloc (ACS->N_inputparams, ACS->N_inputparams); // testing idenity
  gsl_vector *s     = gsl_vector_alloc (ACS->N_inputparams);     // eigenvalues AA^T
  
  // Loading Gamma matrix from input into gsl matrix A
  for (i = 0; i<ACS->N_inputparams * ACS->N_inputparams; i++)
  {
    gsl_matrix_set (A, i/ACS->N_inputparams, i%ACS->N_inputparams, Input[i]);
  }
  
  /*********************************************************************
   *
   *  COMPUTING INVERSE
   * 		- PERFORM SVD
   * 		- CALCULATE INVERSE
   * 
   ********************************************************************/ 
  
  if(ACS->verbose_switch == 1) fprintf(stdout, "... Computing Inverse \n");
	
	// Prepare U for SVD
	gsl_matrix_memcpy(U, A);
	
	// Perform SVD
	if(ACS->verbose_switch == 1) fprintf(stdout, "... Computing Inverse - Performing SVD\n");
	gsl_linalg_SV_decomp_jacobi(U, V, s);  
	
	//printf("Eigenvalues\n");
	//for(i=0;i<ACS->N_inputparams;i++)
	//{
	//printf("%e\t", gsl_vector_get(s, i));
	//}
	//printf("\n");
	
	if(ACS->verbose_switch == 1) fprintf(stdout, "... Computing Inverse - Inverting SVD elements \n");
	// Compute Inverse S
	for (i = 0; i<ACS->N_inputparams; i++)
	{
		gsl_vector_set( s, i, 1./gsl_vector_get( s, i) );
		gsl_matrix_set( InvS, i, i, gsl_vector_get( s, i) );
	}
	
	//fprintf(stdout, "InvS = \n");
	//gsl_matrix_fprintf(stdout, InvS, "%e");
	
	// Tranpose U
	gsl_matrix_transpose(U);
	
	if(ACS->verbose_switch == 1) fprintf(stdout, "... Computing Inverse - Multiplying inverted elements \n");
	
	// Multiply V and InvS
	gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
									1.0, V, InvS,
									0.0, C);
	
	gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
									1.0, C, U,
									0.0, InvA);                             
  
  for(i=0;i<(ACS->N_inputparams*ACS->N_inputparams);i++)
  {
		Inverse[i] = fabs(gsl_matrix_get(InvA, i/ACS->N_inputparams, i%ACS->N_inputparams));
		//fprintf(stdout, "%e",gsl_matrix_get(InvA, i/ACS->N_inputparams, i%ACS->N_inputparams));
		//if(i%ACS->N_inputparams==(ACS->N_inputparams-1)) printf("\n");
		//else printf("\t");
	}
  
  if(ACS->verbose_switch == 1) fprintf(stdout, "... Computing Inverse - done \n");
  
  /*********************************************************************
   *
   *  TESTING ACCURACY
   * 		- A * INVA = 1
   * 
   ********************************************************************/  
  
  if(ACS->verbose_switch == 1) fprintf(stdout, "... Testing accuracy \n");
  
  gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
									1.0, A, InvA,
									0.0, I);      
	
  if(ACS->verbose_switch == 1) fprintf(stdout, "Unit = \n");
  if(ACS->verbose_switch == 1)
  {
		for(i=0;i<(ACS->N_inputparams*ACS->N_inputparams);i++)
		{
			fprintf(stdout, "%e", gsl_matrix_get(I, i/ACS->N_inputparams, i%ACS->N_inputparams));
			if(i%ACS->N_inputparams==(ACS->N_inputparams-1)) fprintf(stdout,"\n");
			else fprintf(stdout,"\t");
		}		
	}
  
  /*********************************************************************
   *
   *  CLEANING UP
   * 
   ********************************************************************/  
  
  if(ACS->verbose_switch == 1) fprintf(stdout, "... Cleaning memory \n");
	
	/* MATRICES */
  gsl_matrix_free(A);
  gsl_matrix_free(U);
  gsl_matrix_free(InvS);
  gsl_matrix_free(InvA);
  gsl_matrix_free(V);
  gsl_matrix_free(C);
  gsl_matrix_free(I);
  gsl_vector_free(s);
  
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
														REAL8FrequencySeries                   *psd,
														PPNConsistencyParamStruc                          *PPNparams,
														REAL8																	 *SNR,
														FisherACS															 *ACS)
/* </lalVerbatim> */
{
	if(ACS->verbose_switch == 1) fprintf(stdout, "Running LALInspiralTestWaveform \n");
	
  /*********************************************************************
   *
   *  Error handling
   * 
   ********************************************************************/ 
  INITSTATUS( status, "LALInspiralComputeSNR", LALINSPIRALCOMPUTEFISHERC);
  ATTATCHSTATUSPTR(status);
  
  /*********************************************************************
   *
   *  Output Files
   * 
   ********************************************************************/   
	
	
  /*********************************************************************
   *
   *  Temporary Variables
   * 
   ********************************************************************/ 
  
	if(ACS->verbose_switch == 1){fprintf(stdout, "... Initialising Variables \n");}
	
	//Counters
	INT4 i,j,k;
	
	// REAL NUMBERS
  REAL8 invpsd;										// Inverse PSD		
  REAL8 SNRtot;										// Total SNR	
  
  // PARAMETERS AND WAVEFORMS
  CoherentGW    	SNRwaveform;		// Holding waveform, loaded onto ht
  PPNConsistencyParamStruc   SNRparams;			// Local copy of input parameters  
	
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
   *  Setting Variables and Allocating Memory
   * 
   ********************************************************************/
  
  if(ACS->verbose_switch==1){fprintf(stdout, "... Setting Variables and Allocating Memory \n");}
	
	// CREATING TIME/FREQUENCY SERIES
	if( (	ht 				= (REAL4TimeSeries *) 				LALMalloc( sizeof(REAL4TimeSeries) ) ) == NULL	&&
		 (	invFT			= (REAL4TimeSeries *) 				LALMalloc( sizeof(REAL4TimeSeries) ) ) == NULL &&
		 ( Hf 				= (COMPLEX8FrequencySeries *) LALMalloc(sizeof(COMPLEX8FrequencySeries))) == NULL &&
		 ( snrpower 	= (COMPLEX8FrequencySeries *) LALMalloc(sizeof(COMPLEX8FrequencySeries))) == NULL			)
	{
		LALFree(ht); 			ht=NULL;
		LALFree(invFT); 	invFT=NULL;
		LALFree(Hf);			Hf = NULL;
		LALFree(snrpower);	snrpower=NULL;
		ABORT(status, LALINSPIRALCOMPUTEFISHERC_EMEM, LALINSPIRALCOMPUTEFISHERC_EMEM_MSG);
	}
	else
	{
		ht 				= XLALCreateREAL4TimeSeries("ht", &(ACS->epoch), 0.0, ACS->deltaT_ts, &lalStrainUnit, ACS->N_ts);
		invFT			=	XLALCreateREAL4TimeSeries("invFT", &(ACS->epoch), 0.0, ACS->deltaT_ts,&lalStrainUnit, ACS->N_ts);
		Hf				= XLALCreateCOMPLEX8FrequencySeries("Hf", &(ACS->epoch), 0.0, ACS->deltaF_fs, &lalStrainUnit, ACS->N_ts/2+1);
		snrpower	= XLALCreateCOMPLEX8FrequencySeries("snrpower", &(ACS->epoch), 0.0, ACS->deltaF_fs, &lalStrainUnit, ACS->N_ts/2+1);
		
		for(j=0;j<ACS->N_ts;j++)
		{
			if(j<ACS->N_ts){ht->data->data[j] 		= 0.0; invFT->data->data[j]	= 0.0;}
			if(j<(ACS->N_ts/2+1)){Hf->data->data[j].re=0.0; Hf->data->data[j].im = 0.0; snrpower->data->data[j].re=0.0, snrpower->data->data[j].im=0.0;}
		}
		
	}
  
  // COPY INPUT PARAMETERS
  SNRparams.ppn = XLALCreateREAL4Vector( PPNparams->ppn->length );
  TRY( XLALCopyPPNConsistencyParamStruc(PPNparams, &SNRparams), status);
  
  /* CLEARING COHERENTGW */
	memset(&SNRwaveform, 0, sizeof(CoherentGW) );  
  
  // FFT
  fwdRealPlan = XLALCreateForwardREAL4FFTPlan(ACS->N_ts,0);
  revRealPlan = XLALCreateReverseREAL4FFTPlan(ACS->N_ts,0);
	
  
  /*********************************************************************
   *
   *  Compute SNR
   * 
   ********************************************************************/
  
  if(ACS->verbose_switch==1){fprintf(stdout, "... Computing SNR - fstart = %e | fstop = %e \n", PPNparams->fStartIn, ACS->fstop);}
  
  // GENERATE WAVEFORM AND COMBINE HP AND HC
  LALGeneratePPNAmpCorConsistency(status->statusPtr, &SNRwaveform, &SNRparams);
  LALInspiralCombinePlusCross(status->statusPtr, SNRwaveform, SNRparams, ht, ACS);
  
  // COPY INTO TIMESERIES
  //for (i=0; i<((INT4)SNRwaveform.f->data->length); i++)
  //{
	//ht->data->data[(ACS->N_ts-SNRwaveform.f->data->length)+i] = SNRwaveform.h->data->data[2*i];
  //}
  
  // PERFORM FFT ON HT
  LALTimeFreqRealFFT( status->statusPtr, Hf, ht, fwdRealPlan );
  
  // Scale using dynamical range
  for (i = 0; i < ACS->N_fs; ++i)
  {
    Hf->data->data[i].re *= ACS->dynRange; Hf->data->data[i].im *= ACS->dynRange; 
    //if(ACS->printall == 1) fprintf(ACS->modsq_out, "%e\t%e\n", i*ACS->deltaF_fs, sqrt(Hf->data->data[i].re*Hf->data->data[i].re+Hf->data->data[i].im*Hf->data->data[i].im)/ACS->dynRange);
  }  
	
	// Calculate Optimal SNR
	for (k = 0; k < ACS->N_ts / 2 + 1; ++k)
  {
    if (psd->data->data[k] == 0.0 || ((REAL8) k)*ACS->deltaF_fs < ACS->fstart || ((REAL8) k)*ACS->deltaF_fs > ACS->fstop )
      invpsd = 0.0;
    else
      invpsd = 1.0 / ( ACS->dynRange * ACS->dynRange * psd->data->data[k]*2.0);
		
    snrpower->data->data[k].re = (Hf->data->data[k].re * Hf->data->data[k].re + Hf->data->data[k].im * Hf->data->data[k].im) * invpsd;
    snrpower->data->data[k].im = 0.0;
  }
  
  // Reserve FT
	LALFreqTimeRealFFT( status->statusPtr, invFT, snrpower, revRealPlan );
	
	// Compute Total SNR
	SNRtot = pow( 4.0*invFT->data->data[0] , 0.5 ); 
  if(ACS->verbose_switch==1){fprintf(stdout, "... Computing SNR - SNRtot = %e\n", SNRtot);}
  *SNR = SNRtot;
  
  /*********************************************************************
   *
   *  Compute SNR
   * 
   ********************************************************************/
	if(ACS->verbose_switch == 1){fprintf(stdout, "... Cleaning up \n");}
  
	XLALDestroyREAL4TimeSeries(ht);							ht = NULL;
	XLALDestroyREAL4TimeSeries(invFT);					invFT = NULL;
	XLALDestroyCOMPLEX8FrequencySeries(Hf);			Hf = NULL;
	XLALDestroyCOMPLEX8FrequencySeries(snrpower); snrpower =NULL;
  XLALDestroyCoherentGW(&SNRwaveform);
  XLALDestroyREAL4Vector(SNRparams.ppn);
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

void LALInspiralComputeDerivatives_linReg(
																					LALStatus																*status,			// Status pointer
																					PPNConsistencyParamStruc														*InputParams,	// Input parameters
																					REAL8																		deltax,				// total width of differentiation
																					INT4																		paramid,			// parameter choice
																					INT4																		icur,				// start waveform element
																					INT4																		iend,					// end waveform element (length)
																					REAL4TimeSeries													*hderiv,			// derivative calculated with linReg algorithm
																					INT4                                    xsamples,
																					FisherACS																*ACS)    
{
	if(ACS->verbose_switch == 1)fprintf(stdout, "Running ComputeDerivatives_linReg - paramid = %d | icur = %d | iend = %d | deltax = %e  | xsamples = %d \n", paramid, icur, iend, deltax, xsamples);
	/* ------------------LAL error handling ----------------------------*/
  INITSTATUS( status, "LALInspiralComputeDerivatives_linReg", LALINSPIRALCOMPUTEFISHERC);
  ATTATCHSTATUSPTR(status);
  
  /*------------------ Check Input parameters ------------------------*/
  ASSERT( InputParams, status, LALINSPIRALCOMPUTEFISHERC_EINPUT, LALINSPIRALCOMPUTEFISHERC_EINPUT_MSG );
  ASSERT( hderiv, status, LALINSPIRALCOMPUTEFISHERC_EINPUT, LALINSPIRALCOMPUTEFISHERC_EINPUT_MSG );	
  
	/*********************************************************************
	 * 
	 * Output Files
	 * 
	 ********************************************************************/  
  
  //FILE *deltaPhi_out;         // phase difference
  //FILE *bandlength_out;       // length of band
  //FILE *phi_out;
	
  //// Create names for output files  
  //char deltaPhi_name[128];
  //char bandlength_name[128];
  //char phi_name[128];
	
  //const char folder[]       = "/data/gravwav/tgfli/Fisher/";  
  //const char deltaPhi_base[]      = "deltaPhi_out_";
  //const char bandlength_base[]  = "bandlength_out_";
  //const char phi_base[] = "phi_out_";
	
  //sprintf(deltaPhi_name, "%s%s%d%s", folder, deltaPhi_base, paramid,".dat");
  //sprintf(bandlength_name,"%s%s%d%s", folder, bandlength_base, paramid, ".dat"); 
  //sprintf(phi_name, "%s%s%d%s", folder, phi_base, paramid, ".dat");
	
  //// Open output files
  //deltaPhi_out = fopen(deltaPhi_name, "w");
  //bandlength_out = fopen(bandlength_name, "w");
  //phi_out = fopen(phi_name, "a");
	
	/*********************************************************************
	 * 
	 * Initiate Variables
	 * 
	 ********************************************************************/
	
	if(ACS->verbose_switch == 1)fprintf(stdout, "... Initialising Variables \n");
	
	/* linReg settings - tune for optimal perfomance */
	
	/* Miscelaneous variables */
	INT4 i = 0;								// Counter	
	INT4 j = 0;								// Counter
	INT4 N = ((INT4)hderiv->data->length);	// Size of TimeSerie
	
	/* Temporary variables for linReg calculations */
	INT4 linReg_count	= 0;			// Counting amount of data points over which linReg takes place
	INT4 hxCur_len		= 0;		// current length of waveform
	INT4 hxPrev_len 	= 0;		// previous length of waveform
	
	
  REAL8 parameter		= 0.0;    // Value of current parameter
	REAL8 curX 				= 0.0;		// current value of x
	REAL8 prevX				= 0.0;		// previous value of x
	REAL8 SumX_min 		= 0.0;		// sum x, minimal point of band
	REAL8 SumXX_min		= 0.0;		// sum x^2, minimal point of band
	REAL8 SumX_max 		= 0.0;		// sum x, maximal point of band
	REAL8 SumXX_max		= 0.0;		// sum x^2, maximal point of band
  REAL8 band_temp   = 0.0;    // calculating width of band
  REAL8 band_cum    = 0.0;
	
	/* Creating TimeSeries */
	if(ACS->verbose_switch == 1)fprintf(stdout, "... Initialising Variables - Creating TimeSeries \n");
	REAL4TimeSeries							*phiStart;			// phase @ xstart
	REAL4TimeSeries							*phiEnd;				// phase @ xend
	REAL4TimeSeries							*hxCur;					// storing current temporary waveform
	REAL4TimeSeries							*hxCur_tmp;			// tmp for shifting waveform in case of time derivative
	REAL4TimeSeries							*hxPrev;				// storing previous temporary waveform
	REAL4TimeSeries							*SumY_min;			// sum y, minimal point of band
	REAL4TimeSeries							*SumXY_min;			// sum x*y, minimal point of band
	REAL4TimeSeries							*SumY_max;			// sum y, maximal point of band
	REAL4TimeSeries							*SumXY_max;			// sum x*y, maximal point of band
	
	/* Create TimeSeries if there is enough memory */
	if( 	(phiStart 	= (REAL4TimeSeries *) LALMalloc( sizeof( REAL4TimeSeries) ) ) == NULL
		 &&	(phiEnd	 		= (REAL4TimeSeries *) LALMalloc( sizeof( REAL4TimeSeries) ) ) == NULL 
		 &&	(hxCur	 		= (REAL4TimeSeries *) LALMalloc( sizeof( REAL4TimeSeries) ) ) == NULL
		 &&	(hxCur_tmp	= (REAL4TimeSeries *) LALMalloc( sizeof( REAL4TimeSeries) ) ) == NULL
		 &&	(hxPrev	 		= (REAL4TimeSeries *) LALMalloc( sizeof( REAL4TimeSeries) ) ) == NULL
		 &&	(SumY_min		= (REAL4TimeSeries *) LALMalloc( sizeof( REAL4TimeSeries) ) ) == NULL	
		 &&	(SumXY_min	= (REAL4TimeSeries *) LALMalloc( sizeof( REAL4TimeSeries) ) ) == NULL
		 &&	(SumY_max		= (REAL4TimeSeries *) LALMalloc( sizeof( REAL4TimeSeries) ) ) == NULL
		 &&	(SumXY_max	= (REAL4TimeSeries *) LALMalloc( sizeof( REAL4TimeSeries) ) ) == NULL)
	{
		TRY( LALFree(phiStart), 	status);			phiStart 	= NULL;
		TRY( LALFree(phiEnd), 		status);			phiEnd 		= NULL;
		TRY( LALFree(hxCur), 			status);			hxCur 		= NULL;
		TRY( LALFree(hxCur_tmp),	status);			hxCur_tmp	= NULL;
		TRY( LALFree(hxPrev), 		status);			hxPrev 		= NULL;
		TRY( LALFree(SumY_min),		status);			SumY_min 	= NULL;
		TRY( LALFree(SumY_min),		status);			SumY_min 	= NULL;
		TRY( LALFree(SumXY_min),	status);			SumXY_min	= NULL;
		TRY( LALFree(SumY_max),		status);			SumY_min 	= NULL;
		TRY( LALFree(SumXY_max),	status);			SumXY_min	= NULL;
		
		ABORT( status, LALINSPIRALCOMPUTEFISHERC_EMEM, LALINSPIRALCOMPUTEFISHERC_EMEM_MSG);
	}
	else
	{    
		// Global TimeSeries variables
    LIGOTimeGPS epoch= LIGOTIMEGPSZERO;
    
    // Creating TimeSeries
    if ( XLALCreateREAL4TimeSeries("test", &epoch, 0.0, InputParams->deltaT, &lalStrainUnit, N) == NULL)
    {
      if(ACS->verbose_switch == 1)fprintf(stdout, "... Initialising Variables - Creating TimeSeries -> Error, destroy TimeSeries and Abort \n");
      ABORT( status, LALINSPIRALCOMPUTEFISHERC_ETSCREATE, LALINSPIRALCOMPUTEFISHERC_ETSCREATE_MSG );
    }
    else
    {
      phiStart  = XLALCreateREAL4TimeSeries("phiStart", &epoch, 0.0, InputParams->deltaT, &lalStrainUnit, N);
      phiEnd   	= XLALCreateREAL4TimeSeries("phiEnd", &epoch, 0.0, InputParams->deltaT, &lalStrainUnit, N);
      hxCur   	= XLALCreateREAL4TimeSeries("hxCur", &epoch, 0.0, InputParams->deltaT, &lalStrainUnit, N);
      hxCur_tmp	= XLALCreateREAL4TimeSeries("hxCur_tmp", &epoch, 0.0, InputParams->deltaT, &lalStrainUnit, N);
      hxPrev   	= XLALCreateREAL4TimeSeries("hxCur", &epoch, 0.0, InputParams->deltaT, &lalStrainUnit, N);
      SumY_min	= XLALCreateREAL4TimeSeries("SumY_min", &epoch, 0.0, InputParams->deltaT, &lalStrainUnit, N);
      SumXY_min	= XLALCreateREAL4TimeSeries("SumXY_min", &epoch, 0.0, InputParams->deltaT, &lalStrainUnit, N);
      SumY_max	= XLALCreateREAL4TimeSeries("SumY_max", &epoch, 0.0, InputParams->deltaT, &lalStrainUnit, N);
      SumXY_max	= XLALCreateREAL4TimeSeries("SumXY_max", &epoch, 0.0, InputParams->deltaT, &lalStrainUnit, N);
		}
		
		if(		phiStart->data 	== NULL
			 &&	phiEnd->data 		== NULL 
			 &&	hxCur->data 		== NULL
			 &&	hxCur_tmp->data	== NULL
			 &&	hxPrev->data 		== NULL
			 &&	SumY_min->data	== NULL
			 &&	SumXY_min->data	== NULL
			 &&	SumY_max->data	== NULL
			 &&	SumXY_max->data	== NULL)
		{
			TRY(XLALDestroyREAL4TimeSeries(phiStart), 	status);
			TRY(XLALDestroyREAL4TimeSeries(phiEnd), 		status);
			TRY(XLALDestroyREAL4TimeSeries(hxCur), 			status);
			TRY(XLALDestroyREAL4TimeSeries(hxCur_tmp),	status);
			TRY(XLALDestroyREAL4TimeSeries(hxPrev),			status);
			TRY(XLALDestroyREAL4TimeSeries(SumY_min),		status);
			TRY(XLALDestroyREAL4TimeSeries(SumXY_min),	status);
			TRY(XLALDestroyREAL4TimeSeries(SumY_max),		status);
			TRY(XLALDestroyREAL4TimeSeries(SumXY_max),	status);
			ABORT( status, LALINSPIRALCOMPUTEFISHERC_ETSCREATE, LALINSPIRALCOMPUTEFISHERC_ETSCREATE_MSG);
		}
		else
		{
			for(i=0; i<N; i++)
			{
				phiStart->data->data[i] 	= 0.0;
				phiEnd->data->data[i]			= 0.0;
				hxCur->data->data[i]			= 0.0;
				hxCur_tmp->data->data[i]	= 0.0;
				hxPrev->data->data[i]			= 0.0;
				SumY_min->data->data[i]		= 0.0;
				SumXY_min->data->data[i]	= 0.0;
				SumY_max->data->data[i]		= 0.0;
				SumXY_max->data->data[i]	= 0.0;
			}
		}
	}
	
	/*********************************************************************
	 * 
	 * linReg Algorithm
   *  - Load in parameters from input
   *  - Increment parameter with deltax (except for dl and tc)
   *    ~ For dl, derivative is -h(t) and hence trivial
   *    ~ For tc, the waveforms are shifted when copied 
   *  - Filter out band edges (if applicable)
   *  - Pass on for condition checking (next block)
	 * 
	 ********************************************************************/		
	
	/* Loading in relevant parameter through paramid */
  switch(paramid)
  {
    case 0:
			parameter = InputParams->mTot_real8;
			if (parameter == 0.0) ABORT( status, LALINSPIRALCOMPUTEFISHERC_EMBAD, LALINSPIRALCOMPUTEFISHERC_EMBAD_MSG );
      break;
    case 1:
			parameter = InputParams->eta_real8; 
			if (parameter == 0.0) ABORT( status, LALINSPIRALCOMPUTEFISHERC_EMBAD, LALINSPIRALCOMPUTEFISHERC_EMBAD_MSG );     
      break;
			//case 2:
      //break;
    case 2:
      parameter = InputParams->phi;
    case 3:
      break;
    case 4:
			parameter = InputParams->cosI;
			break;
		case 5:
			parameter = InputParams->phasePNparams[ACS->testPhaseParam];
			break;
    default:
      ABORT( status, LALINSPIRALCOMPUTEFISHERC_EINPUT, LALINSPIRALCOMPUTEFISHERC_EINPUT_MSG );
      break;
  }
	
	/* Looping over parameter space */
  for(i=0; i<xsamples; i++)
  {
    if( ((100*i/xsamples)%10) == 0){if(ACS->verbose_switch==1){fprintf(stdout, "\r... Creating h(x): %d percent done", 100*i/xsamples);}}
    if(i==(xsamples-1)){if(ACS->verbose_switch==1){fprintf(stdout, "\n");}}
    
    // Creating hxwaveform structure
    CoherentGW hxwaveform;                // h(x) temp waveform
    memset(&hxwaveform, 0, sizeof(CoherentGW) );
    
    // Creating temporary paramsStruc
    PPNConsistencyParamStruc hxPPNparams;
    hxPPNparams.ppn = XLALCreateREAL4Vector( InputParams->ppn->length );
    TRY( XLALCopyPPNConsistencyParamStruc(InputParams, &hxPPNparams), status);
		
    if (paramid == 3)
    {
      prevX = parameter +  (i-1) * deltax;
      curX = parameter +  i * deltax;
    }
    else
    {
      prevX = parameter +  ((i-xsamples/2)-1) * deltax;
      curX = parameter +  (i-xsamples/2) * deltax;
    }
		
    // Generate and copy waveforms
    switch(paramid)
    {
      case 0:
        hxPPNparams.mTot_real8 = curX;
        
				/* REFILL PHI_I*/
				if(ACS->testPhaseParam!=0) hxPPNparams.phasePNparams[0] = -pow(hxPPNparams.eta_real8,-3.0/8.0)*pow(5.0*LAL_MTSUN_SI*hxPPNparams.mTot_real8,-5.0/8.0); 
				if(ACS->testPhaseParam!=1) hxPPNparams.phasePNparams[1] = -(3715.0/8064.0 + 55.0/96.0*hxPPNparams.eta_real8)*pow(hxPPNparams.eta_real8,-5.0/8.0)*pow(5.0*LAL_MTSUN_SI*hxPPNparams.mTot_real8,-3.0/8.0);
				if(ACS->testPhaseParam!=2) hxPPNparams.phasePNparams[2] = 3.0/4.0*LAL_PI*pow(hxPPNparams.eta_real8,-0.75)*pow(5.0*LAL_MTSUN_SI*hxPPNparams.mTot_real8,-0.25); 
				if(ACS->testPhaseParam!=3) hxPPNparams.phasePNparams[3] = -(9275495.0/14450688.0 + 284875.0/258048.0*hxPPNparams.eta_real8 + 1855.0/2048.0*pow(hxPPNparams.eta_real8,2.0))*pow(hxPPNparams.eta_real8,-7.0/8.0)*pow(5.0*LAL_MTSUN_SI*hxPPNparams.mTot_real8,-1.0/8.0);
				if(ACS->testPhaseParam!=4) hxPPNparams.phasePNparams[4] = -1.0/hxPPNparams.eta_real8*(-38645.0/172032.0 + 65.0/2048.0*hxPPNparams.eta_real8)*LAL_PI*log(hxPPNparams.eta_real8/(5.0*LAL_MTSUN_SI*hxPPNparams.mTot_real8));
				if(ACS->testPhaseParam!=5) hxPPNparams.phasePNparams[5] = -1.0/hxPPNparams.eta_real8*(-38645.0/172032.0 + 65.0/2048.0*hxPPNparams.eta_real8)*LAL_PI; 
				if(ACS->testPhaseParam!=6) hxPPNparams.phasePNparams[6] = -(831032450749357.0/57682522275840.0 - 53.0/40.0*LAL_PI*LAL_PI - 107.0/56.0*LAL_GAMMA + 107.0/448.0*log(hxPPNparams.eta_real8/(256*5.0*LAL_MTSUN_SI*hxPPNparams.mTot_real8)) + (-123292747421.0/4161798144.0 + 2255.0/2048.0*LAL_PI*LAL_PI + 385.0/48.0*(-1987.0/3080.0) - 55.0/16.0*(-11831.0/9240.0))*hxPPNparams.eta_real8 + 154565.0/1835008.0*pow(hxPPNparams.eta_real8,2.0) - 1179625.0/1769472.0*pow(hxPPNparams.eta_real8,3.0))*pow(hxPPNparams.eta_real8,-9.0/8.0)*pow(5.0*LAL_MTSUN_SI*hxPPNparams.mTot_real8,1.0/8.0);
				if(ACS->testPhaseParam!=7) hxPPNparams.phasePNparams[7] = -107.0/448.0*pow(hxPPNparams.eta_real8,-9.0/8.0)*pow(5.0*LAL_MTSUN_SI*hxPPNparams.mTot_real8,1.0/8.0);
				if(ACS->testPhaseParam!=8) hxPPNparams.phasePNparams[8] = -(188516689.0/173408256.0 + 488825.0/516096.0*hxPPNparams.eta_real8 - 141769.0/516096.0*pow(hxPPNparams.eta_real8,2.0))*LAL_PI*pow(hxPPNparams.eta_real8,-5.0/4.0)*pow(5.0*LAL_MTSUN_SI*hxPPNparams.mTot_real8,1.0/4.0); 		        
        
        break;
      case 1:
        hxPPNparams.eta_real8 = curX;
        
				/* REFILL PHI_I FOR */
				if(ACS->testPhaseParam!=0) hxPPNparams.phasePNparams[0] = -pow(hxPPNparams.eta_real8,-3.0/8.0)*pow(5.0*LAL_MTSUN_SI*hxPPNparams.mTot_real8,-5.0/8.0); 
				if(ACS->testPhaseParam!=1) hxPPNparams.phasePNparams[1] = -(3715.0/8064.0 + 55.0/96.0*hxPPNparams.eta_real8)*pow(hxPPNparams.eta_real8,-5.0/8.0)*pow(5.0*LAL_MTSUN_SI*hxPPNparams.mTot_real8,-3.0/8.0);
				if(ACS->testPhaseParam!=2) hxPPNparams.phasePNparams[2] = 3.0/4.0*LAL_PI*pow(hxPPNparams.eta_real8,-0.75)*pow(5.0*LAL_MTSUN_SI*hxPPNparams.mTot_real8,-0.25); 
				if(ACS->testPhaseParam!=3) hxPPNparams.phasePNparams[3] = -(9275495.0/14450688.0 + 284875.0/258048.0*hxPPNparams.eta_real8 + 1855.0/2048.0*pow(hxPPNparams.eta_real8,2.0))*pow(hxPPNparams.eta_real8,-7.0/8.0)*pow(5.0*LAL_MTSUN_SI*hxPPNparams.mTot_real8,-1.0/8.0);
				if(ACS->testPhaseParam!=4) hxPPNparams.phasePNparams[4] = -1.0/hxPPNparams.eta_real8*(-38645.0/172032.0 + 65.0/2048.0*hxPPNparams.eta_real8)*LAL_PI*log(hxPPNparams.eta_real8/(5.0*LAL_MTSUN_SI*hxPPNparams.mTot_real8));
				if(ACS->testPhaseParam!=5) hxPPNparams.phasePNparams[5] = -1.0/hxPPNparams.eta_real8*(-38645.0/172032.0 + 65.0/2048.0*hxPPNparams.eta_real8)*LAL_PI; 
				if(ACS->testPhaseParam!=6) hxPPNparams.phasePNparams[6] = -(831032450749357.0/57682522275840.0 - 53.0/40.0*LAL_PI*LAL_PI - 107.0/56.0*LAL_GAMMA + 107.0/448.0*log(hxPPNparams.eta_real8/(256*5.0*LAL_MTSUN_SI*hxPPNparams.mTot_real8)) + (-123292747421.0/4161798144.0 + 2255.0/2048.0*LAL_PI*LAL_PI + 385.0/48.0*(-1987.0/3080.0) - 55.0/16.0*(-11831.0/9240.0))*hxPPNparams.eta_real8 + 154565.0/1835008.0*pow(hxPPNparams.eta_real8,2.0) - 1179625.0/1769472.0*pow(hxPPNparams.eta_real8,3.0))*pow(hxPPNparams.eta_real8,-9.0/8.0)*pow(5.0*LAL_MTSUN_SI*hxPPNparams.mTot_real8,1.0/8.0);
				if(ACS->testPhaseParam!=7) hxPPNparams.phasePNparams[7] = -107.0/448.0*pow(hxPPNparams.eta_real8,-9.0/8.0)*pow(5.0*LAL_MTSUN_SI*hxPPNparams.mTot_real8,1.0/8.0);
				if(ACS->testPhaseParam!=8) hxPPNparams.phasePNparams[8] = -(188516689.0/173408256.0 + 488825.0/516096.0*hxPPNparams.eta_real8 - 141769.0/516096.0*pow(hxPPNparams.eta_real8,2.0))*LAL_PI*pow(hxPPNparams.eta_real8,-5.0/4.0)*pow(5.0*LAL_MTSUN_SI*hxPPNparams.mTot_real8,1.0/4.0); 		                
        
        break;
				//case 2:
        //break;
      case 2:
        hxPPNparams.phi = curX;
      case 3:
        break;
			case 4:
				hxPPNparams.cosI = curX;
				/* MIRROR AROUND -1 IF CROSSING -1 */
				if(hxPPNparams.cosI < -1.0 ){ hxPPNparams.cosI = -2.0 - hxPPNparams.cosI; hxPPNparams.sinI = -hxPPNparams.sinI;}
				/* MIRROR AROUND +1 IF CROSSING +1 */
				if(hxPPNparams.cosI > 1.0 ){ hxPPNparams.cosI = 2.0 - hxPPNparams.cosI; hxPPNparams.sinI = -hxPPNparams.sinI;}
				break;
			case 5:
				hxPPNparams.phasePNparams[ACS->testPhaseParam] = curX;
				break;
      default:
        if(ACS->verbose_switch == 1)fprintf(stdout, "XLALInspiralComputeDerivatives failed: invalid derivative variable \n");
        break;
    }
		
    // Create waveforms
    TRY( LALGeneratePPNAmpCorConsistency( status->statusPtr, &hxwaveform, &hxPPNparams), status);
		if(paramid==3){LALInspiralCombinePlusCross( status->statusPtr, hxwaveform, hxPPNparams, hxCur_tmp, ACS);}
		else{LALInspiralCombinePlusCross( status->statusPtr, hxwaveform, hxPPNparams, hxCur, ACS);}
		
    // Get length and check if vector is long enough to hold waveform
    hxCur_len = hxwaveform.f->data->length;
    
    /* SHIFT WAVEFORM FOR TIME DERIVATIVE */
    if (paramid == 3)
    {
			for(j=0;j+i<N;j++)
			{
				hxCur->data->data[j+i] = hxCur_tmp->data->data[j];
			}
		}
		//LALInspiralWaveTaper(status->statusPtr, hxCur->data, 3);
    
    // if time series length is long enough, copy waveform into time series
    //if( hxCur_len < N )
    //{
		//for (j=0;j<hxCur_len;j++)
		//{
		//if (paramid == 4 && ((N-hxCur_len)+j+i) < N)
		//{
		//hxCur->data->data[(N-hxCur_len)+j+i] =  hxwaveform.h->data->data[2*j];
		////if (j==0) fprintf(stdout, "k = %d, h(t) = %e \n", (N-hxCur_len)+j+i, hxCur->data->data[(N-hxCur_len)+j+i] );
		//}
		//else
		//{
		//hxCur->data->data[(INT4)((N-hxCur_len)+j)] = hxwaveform.h->data->data[2*j];
		//}
		//}
    //}
    //else
    //{
		//ABORT(status, LALINSPIRALCOMPUTEFISHERC_ETSCREATE, LALINSPIRALCOMPUTEFISHERC_ETSCREATE_MSG);
    //}  
    
    /*if(paramid == 0 || paramid == 1){} */
		
		/* PRINT OUTPUT TO FILE */
		/*if(paramid==1)
		 {
		 INT4 i0           = 0;
		 fprintf(ACS->hx_out, "%e\t",curX);
		 for(i0=0;i0<((INT4)hxwaveform.f->data->length);i0+=200)
		 {
		 fprintf(ACS->hx_out, "%e\t", hxCur->data->data[N-1-i0]);
		 }
		 fprintf(ACS->hx_out, "\n");
		 }     */  
    
    if( ((paramid == 2 || paramid == 3 || paramid == 4 || paramid == 5) && ((INT4)i) !=0 ) || ( paramid == 1 && ((INT4)i) != 0 && hxCur_len != hxPrev_len && curX < 0.25) || (paramid == 0 && ((INT4)i) != 0 && hxCur_len != hxPrev_len ) )
    {
      linReg_count++;   
      
      if(linReg_count > 1) band_cum  += (curX - band_temp);
      band_temp = curX;
			
			// Statistics for LinReg
      SumX_min    += curX;
      SumX_max    += prevX;
      SumXX_min  += pow(curX, 2.0);
      SumXX_max  += pow(prevX, 2.0);
      for (j=0; j<((INT4) N); j++)
      {
        SumY_min->data->data[j]  += (*hxCur).data->data[j];
        SumY_max->data->data[j]  += (*hxPrev).data->data[j];
        SumXY_min->data->data[j] += (*hxCur).data->data[j] * curX;
        SumXY_max->data->data[j] += (*hxPrev).data->data[j] * prevX;
        if(linReg_count == 1) 
        {
          if (j<hxCur_len)phiStart->data->data[(N-hxCur_len)+j] = hxwaveform.phi->data->data[j];
        }
        else 
        {
          if (paramid == 3)
          {
            if ( ((N-hxCur_len)+j+i) <N ) phiEnd->data->data[(N-hxCur_len)+j+i] = hxwaveform.phi->data->data[j];
          }
          else
          {
            if (j<hxCur_len) phiEnd->data->data[(N-hxCur_len)+j] = hxwaveform.phi->data->data[j];
          }
        }
      }
    }
    
    // Loop completed, storing cur into old variables, empty cur
    hxPrev_len = hxCur_len;
    
    // move waveform into old waveform and clean new waveform
    for(j=0; j<((INT4) N) ; j++)
    {
      hxPrev->data->data[j] = hxCur->data->data[j];
      hxCur->data->data[j] = 0.0;
      hxCur_tmp->data->data[j] = 0.0;
    }
    
    TRY( XLALDestroyCoherentGW(&hxwaveform), status); 
    TRY( XLALDestroyREAL4Vector( hxPPNparams.ppn ), status);
    
  }
  
  /* Computing derivative */
  if(ACS->verbose_switch == 1)fprintf(stdout, "... Calculating Derivative \n");
  for(i=icur;i<iend;i++)
  {
    if(icur>N) if(ACS->verbose_switch == 1)fprintf(stdout, "... icur too large \n");
    hderiv->data->data[N-1-i] = ((SumXY_min->data->data[N-1-i] - SumX_min * SumY_min->data->data[N-1-i] / linReg_count) / (SumXX_min - pow(SumX_min, 2.0) / linReg_count)
																 + (SumXY_max->data->data[N-1-i] - SumX_max * SumY_max->data->data[N-1-i] / linReg_count) / (SumXX_max - pow(SumX_max, 2.0) / linReg_count))/2.0; 
    
    //fprintf(phi_out, "%e\t%e\t%e\t%e\n",i*InputParams->deltaT, phiStart->data->data[N-1-i], phiEnd->data->data[N-1-i], fabs(phiStart->data->data[N-1-i] - phiEnd->data->data[N-1-i]));
  }
  
  /* Extra checks */
  //fprintf(stdout, "phiStart = %e \t phiEnd = %e \n", phiStart->data->data[N-i0], phiEnd->data->data[N-i0]);
  //fprintf(deltaPhi_out, "%e\t%e \n", parameter, phiStart->data->data[N-i0] - phiEnd->data->data[N-i0]);
  //fprintf(bandlength_out, "%e\t%e \n", parameter, band_cum/((REAL8)(linReg_count-1)));
	
	
	/*********************************************************************
	 * 
	 * Clean Up
	 * 
	 ********************************************************************/	
	if(ACS->verbose_switch == 1)fprintf(stdout, "... Cleaning up \n");
  
  if(ACS->verbose_switch == 1)fprintf(stdout, "... Cleaning up - Closing files\n");
  //fclose(deltaPhi_out);
  //fclose(bandlength_out);
  //fclose(phi_out);
	
	if(ACS->verbose_switch == 1)fprintf(stdout, "... Cleaning up - TimeSeries\n");
	TRY( XLALDestroyREAL4TimeSeries(phiStart),	status);
	TRY( XLALDestroyREAL4TimeSeries(phiEnd),		status);
	TRY( XLALDestroyREAL4TimeSeries(hxCur),			status);
	TRY( XLALDestroyREAL4TimeSeries(hxCur_tmp),	status);
	TRY( XLALDestroyREAL4TimeSeries(SumY_min),	status);
	TRY( XLALDestroyREAL4TimeSeries(SumXY_min),	status);
	TRY( XLALDestroyREAL4TimeSeries(SumY_max),	status);
	TRY( XLALDestroyREAL4TimeSeries(SumXY_max), status);
	
	/*********************************************************************
	 * 
	 * Detach LAL error handling
	 * 
	 ********************************************************************/	
  DETATCHSTATUSPTR( status );
  RETURN(status);
}

void LALInspiralComputeDerivatives_5point(
																					LALStatus																*status,			// Status pointer
																					PPNConsistencyParamStruc														*InputParams,	// Input parameters
																					INT4																		paramid,			// parameter choice
																					REAL8																		deltax,				// differentiation width
																					INT4																		istart,				// starting time
																					INT4																		iend,					// end time
																					REAL4TimeSeries													*hderiv, 
																					FisherACS																*ACS)			// TimeSerie holding derivative
{
	if(ACS->verbose_switch == 1)fprintf(stdout, "Running LALInspiralComputeDerivatives_5point \n");
	
	/*********************************************************************
	 * 
	 * LAL error handling
	 * 
	 *********************************************************************/	
  
  INITSTATUS( status, "LALInpiralComputeDerivatives_5point", LALINSPIRALCOMPUTEFISHERC);
  ATTATCHSTATUSPTR(status);
  
	
	/*********************************************************************
	 * 
	 * Check input
	 * 
	 *********************************************************************/
  
	ASSERT( InputParams, status, LALINSPIRALCOMPUTEFISHERC_EINPUT, LALINSPIRALCOMPUTEFISHERC_EINPUT_MSG );
  ASSERT( hderiv, status, LALINSPIRALCOMPUTEFISHERC_EINPUT, LALINSPIRALCOMPUTEFISHERC_EINPUT_MSG );	
  
	/*********************************************************************
	 * 
	 * Input parameters
	 * 
	 *********************************************************************/
  
	// Global TimeSeries variables
	LIGOTimeGPS epoch= LIGOTIMEGPSZERO;
  
	/*********************************************************************
	 * 
	 * Temporary Variables
	 * 
	 *********************************************************************/	 
	
	if(ACS->verbose_switch == 1) fprintf(stdout, "... Initialising Variables \n");
	
	// COUNTERS	
	INT4 i;
	
	// VARIOUS
	INT4 N		=		hderiv->data->length;		// segment size
	
	// TIME SERIES
	REAL4TimeSeries			*ht;
	REAL4TimeSeries			*m1ht;
	REAL4TimeSeries			*m2ht;
	REAL4TimeSeries			*p1ht;
	REAL4TimeSeries			*p2ht;
	
	REAL4TimeSeries			*ht_tmp;
	REAL4TimeSeries			*m1ht_tmp;
	REAL4TimeSeries			*m2ht_tmp;
	REAL4TimeSeries			*p1ht_tmp;
	REAL4TimeSeries			*p2ht_tmp;
	
	// Coherent GW's
	CoherentGW					htwaveform;
	CoherentGW					m1htwaveform;
	CoherentGW					m2htwaveform;
	CoherentGW					p1htwaveform;
	CoherentGW					p2htwaveform;
	
	// PPN PARAMS
	PPNConsistencyParamStruc htPPNparams;
	PPNConsistencyParamStruc m1htPPNparams;
	PPNConsistencyParamStruc m2htPPNparams;
	PPNConsistencyParamStruc p1htPPNparams;
	PPNConsistencyParamStruc p2htPPNparams;
	
	/*********************************************************************
	 * 
	 * Allocate Memory
	 * 
	 *********************************************************************/	
  if(ACS->verbose_switch == 1) fprintf(stdout, "... Allocating Memory \n");
  
  // TIME SERIES
  if( 		(ht = (REAL4TimeSeries *) LALMalloc( sizeof( REAL4TimeSeries) ) ) == NULL
		 &&	(m1ht = (REAL4TimeSeries *) LALMalloc( sizeof( REAL4TimeSeries) ) ) == NULL	
		 &&	(m2ht = (REAL4TimeSeries *) LALMalloc( sizeof( REAL4TimeSeries) ) ) == NULL
		 &&	(p1ht = (REAL4TimeSeries *) LALMalloc( sizeof( REAL4TimeSeries) ) ) == NULL
		 &&	(p2ht = (REAL4TimeSeries *) LALMalloc( sizeof( REAL4TimeSeries) ) ) == NULL
		 &&	(ht_tmp = (REAL4TimeSeries *) LALMalloc( sizeof( REAL4TimeSeries) ) ) == NULL	
		 &&	(m1ht_tmp = (REAL4TimeSeries *) LALMalloc( sizeof( REAL4TimeSeries) ) ) == NULL	
		 &&	(m2ht_tmp = (REAL4TimeSeries *) LALMalloc( sizeof( REAL4TimeSeries) ) ) == NULL	
		 &&	(p1ht_tmp = (REAL4TimeSeries *) LALMalloc( sizeof( REAL4TimeSeries) ) ) == NULL	
		 &&	(p2ht_tmp = (REAL4TimeSeries *) LALMalloc( sizeof( REAL4TimeSeries) ) ) == NULL	)
	{
		TRY(LALFree(ht), status);			ht = NULL;
		TRY(LALFree(m1ht), status);		m1ht = NULL;
		TRY(LALFree(m2ht), status);		m2ht = NULL;
		TRY(LALFree(p1ht), status);		p1ht = NULL;
		TRY(LALFree(p2ht), status);		p2ht = NULL;
		TRY(LALFree(ht_tmp), status);		ht_tmp = NULL;
		TRY(LALFree(m1ht_tmp), status);		m1ht_tmp = NULL;
		TRY(LALFree(m2ht_tmp), status);		m2ht_tmp = NULL;
		TRY(LALFree(p1ht_tmp), status);		p1ht_tmp = NULL;
		TRY(LALFree(p2ht_tmp), status);		p2ht_tmp = NULL;
		
		ABORT( status, LALINSPIRALCOMPUTEFISHERC_EMEM, LALINSPIRALCOMPUTEFISHERC_EMEM_MSG);
	}
	else
	{
		ht		=		XLALCreateREAL4TimeSeries("ht", &epoch, 0.0, InputParams->deltaT, &lalStrainUnit, N);
		m1ht	=		XLALCreateREAL4TimeSeries("m1", &epoch, 0.0, InputParams->deltaT, &lalStrainUnit, N);
		m2ht	=		XLALCreateREAL4TimeSeries("m2", &epoch, 0.0, InputParams->deltaT, &lalStrainUnit, N);
		p1ht	=		XLALCreateREAL4TimeSeries("p1", &epoch, 0.0, InputParams->deltaT, &lalStrainUnit, N);
		p2ht	=		XLALCreateREAL4TimeSeries("p2", &epoch, 0.0, InputParams->deltaT, &lalStrainUnit, N);
		
		ht_tmp	=		XLALCreateREAL4TimeSeries("ht_tmp", &epoch, 0.0, InputParams->deltaT, &lalStrainUnit, N);
		m1ht_tmp	=		XLALCreateREAL4TimeSeries("m1_tmp", &epoch, 0.0, InputParams->deltaT, &lalStrainUnit, N);
		m2ht_tmp	=		XLALCreateREAL4TimeSeries("m2_tmp", &epoch, 0.0, InputParams->deltaT, &lalStrainUnit, N);
		p1ht_tmp	=		XLALCreateREAL4TimeSeries("p1_tmp", &epoch, 0.0, InputParams->deltaT, &lalStrainUnit, N);
		p2ht_tmp	=		XLALCreateREAL4TimeSeries("p2_tmp", &epoch, 0.0, InputParams->deltaT, &lalStrainUnit, N);
	}
	
	// CLEAR TIMESERIES
	for(i=0;i<N;i++)
	{
		ht->data->data[i] = 0.0;
		m1ht->data->data[i] = 0.0;
		m2ht->data->data[i] = 0.0;
		p1ht->data->data[i] = 0.0;
		p2ht->data->data[i] = 0.0;
		
		ht_tmp->data->data[i] = 0.0;
		m1ht_tmp->data->data[i] = 0.0;
		m2ht_tmp->data->data[i] = 0.0;
		p1ht_tmp->data->data[i] = 0.0;
		p2ht_tmp->data->data[i] = 0.0;
	}
	
	
	// COHERENT GW'S
	memset(&htwaveform, 0, sizeof(CoherentGW) );
	memset(&m1htwaveform, 0, sizeof(CoherentGW) );
	memset(&m2htwaveform, 0, sizeof(CoherentGW) );
	memset(&p1htwaveform, 0, sizeof(CoherentGW) );
	memset(&p2htwaveform, 0, sizeof(CoherentGW) );
	
	// Creating temporary paramsStruc
	htPPNparams.ppn = XLALCreateREAL4Vector( InputParams->ppn->length );
	m1htPPNparams.ppn = XLALCreateREAL4Vector( InputParams->ppn->length );
	m2htPPNparams.ppn = XLALCreateREAL4Vector( InputParams->ppn->length );
	p1htPPNparams.ppn = XLALCreateREAL4Vector( InputParams->ppn->length );
	p2htPPNparams.ppn = XLALCreateREAL4Vector( InputParams->ppn->length );
	
	/*********************************************************************
	 * 
	 * Calculate Derivative using 5-point method
	 * 
	 *********************************************************************/	
	
	if(ACS->verbose_switch == 1) fprintf(stdout, "... Calculating 5-point derivative \n");
	
	// COPY INPUT PARAMETERS
	if(ACS->verbose_switch == 1) fprintf(stdout, "... Calculating 5-point derivative - Copy input parameters \n");	
	TRY( XLALCopyPPNConsistencyParamStruc(InputParams, &htPPNparams), status);
	TRY( XLALCopyPPNConsistencyParamStruc(InputParams, &m1htPPNparams), status);
	TRY( XLALCopyPPNConsistencyParamStruc(InputParams, &m2htPPNparams), status);
	TRY( XLALCopyPPNConsistencyParamStruc(InputParams, &p1htPPNparams), status);
	TRY( XLALCopyPPNConsistencyParamStruc(InputParams, &p2htPPNparams), status);
	
	// SHIFT PARAMETERS ACCORDING TO ABSDELTA/DELTAT
	if(ACS->verbose_switch == 1) fprintf(stdout, "... Calculating 5-point derivative - Shift parameters \n");	
	switch(paramid)
	{
		case 0:
			m1htPPNparams.mTot_real8 = htPPNparams.mTot_real8 - deltax;
			m2htPPNparams.mTot_real8 = htPPNparams.mTot_real8 - 2.0*deltax;
			p1htPPNparams.mTot_real8 = htPPNparams.mTot_real8 + deltax;
			p2htPPNparams.mTot_real8 = htPPNparams.mTot_real8 + 2.0*deltax;
			break;
			
		case 1:
			m1htPPNparams.eta_real8 = htPPNparams.eta_real8 - deltax;
			m2htPPNparams.eta_real8 = htPPNparams.eta_real8 - 2.0*deltax;
			p1htPPNparams.eta_real8 = htPPNparams.eta_real8 + deltax;
			p2htPPNparams.eta_real8 = htPPNparams.eta_real8 + 2.0*deltax;
			break;
			
			//case 2:
			//break;
			
		case 2:
			m1htPPNparams.phi = htPPNparams.phi - deltax;
			m2htPPNparams.phi = htPPNparams.phi - 2.0*deltax;
			p1htPPNparams.phi = htPPNparams.phi + deltax;
			p2htPPNparams.phi = htPPNparams.phi + 2.0*deltax;		
			break;
			
		case 3:
			break;
			
		case 4:
			/* CHECKING BOUNDS OF COSi */
			m1htPPNparams.cosI = htPPNparams.cosI - deltax;
			if(m1htPPNparams.cosI <-1.0){m1htPPNparams.cosI = -2.0 - m1htPPNparams.cosI;} if(m1htPPNparams.cosI > 1.0){m1htPPNparams.cosI = 2.0 - m1htPPNparams.cosI;}
			m2htPPNparams.cosI = htPPNparams.cosI - 2.0*deltax;
			if(m2htPPNparams.cosI <-1.0){m2htPPNparams.cosI = -2.0 - m2htPPNparams.cosI;} if(m2htPPNparams.cosI > 1.0){m2htPPNparams.cosI = 2.0 - m2htPPNparams.cosI;}
			p1htPPNparams.cosI = htPPNparams.cosI + deltax;
			if(p1htPPNparams.cosI <-1.0){p1htPPNparams.cosI = -2.0 - p1htPPNparams.cosI;} if(p1htPPNparams.cosI > 1.0){p1htPPNparams.cosI = 2.0 - p1htPPNparams.cosI;}
			p2htPPNparams.cosI = htPPNparams.cosI + 2.0*deltax;
			if(p2htPPNparams.cosI <-1.0){p2htPPNparams.cosI = -2.0 - p2htPPNparams.cosI;} if(p2htPPNparams.cosI > 1.0){p2htPPNparams.cosI = 2.0 - p2htPPNparams.cosI;}
			break;
			
		case 5:
			m1htPPNparams.phasePNparams[ACS->testPhaseParam] = htPPNparams.phasePNparams[ACS->testPhaseParam] - deltax;
			m2htPPNparams.phasePNparams[ACS->testPhaseParam] = htPPNparams.phasePNparams[ACS->testPhaseParam] - 2.0*deltax;
			p1htPPNparams.phasePNparams[ACS->testPhaseParam] = htPPNparams.phasePNparams[ACS->testPhaseParam] + deltax;
			p2htPPNparams.phasePNparams[ACS->testPhaseParam] = htPPNparams.phasePNparams[ACS->testPhaseParam] + 2.0*deltax;
			break;
			
		default:
			ABORT(status, LALINSPIRALCOMPUTEFISHERC_EINPUT, LALINSPIRALCOMPUTEFISHERC_EINPUT_MSG);
			break;
			
	}
	
	
	// GENERATE WAVEFORMS
	if(ACS->verbose_switch == 1) fprintf(stdout, "... Calculating 5-point derivative - Generate Waveforms \n");	
	LALGeneratePPNAmpCorConsistency( status->statusPtr, &htwaveform, &htPPNparams);
	LALGeneratePPNAmpCorConsistency( status->statusPtr, &m1htwaveform, &m1htPPNparams);
	LALGeneratePPNAmpCorConsistency( status->statusPtr, &m2htwaveform, &m2htPPNparams);
	LALGeneratePPNAmpCorConsistency( status->statusPtr, &p1htwaveform, &p1htPPNparams);
	LALGeneratePPNAmpCorConsistency( status->statusPtr, &p2htwaveform, &p2htPPNparams);
	
	/* COMBINE HP AND HC */
	if(paramid==3)
	{
		LALInspiralCombinePlusCross(status->statusPtr, htwaveform, htPPNparams, ht_tmp, ACS);
		LALInspiralCombinePlusCross(status->statusPtr, m1htwaveform, m1htPPNparams, m1ht_tmp, ACS);
		LALInspiralCombinePlusCross(status->statusPtr, m2htwaveform, m2htPPNparams, m2ht_tmp, ACS);
		LALInspiralCombinePlusCross(status->statusPtr, p1htwaveform, p1htPPNparams, p1ht_tmp, ACS);
		LALInspiralCombinePlusCross(status->statusPtr, p2htwaveform, p2htPPNparams, p2ht_tmp, ACS);		
	}
	else
	{	
		LALInspiralCombinePlusCross(status->statusPtr, htwaveform, htPPNparams, ht, ACS);
		LALInspiralCombinePlusCross(status->statusPtr, m1htwaveform, m1htPPNparams, m1ht, ACS);
		LALInspiralCombinePlusCross(status->statusPtr, m2htwaveform, m2htPPNparams, m2ht, ACS);
		LALInspiralCombinePlusCross(status->statusPtr, p1htwaveform, p1htPPNparams, p1ht, ACS);
		LALInspiralCombinePlusCross(status->statusPtr, p2htwaveform, p2htPPNparams, p2ht, ACS);
	}
	
	/* SHIFT WAVEFORM IN CASE OF PARAMID == 4 */
	if(paramid == 3)
	{
		for(i=0;i<N;i++)
		{
			if(i+2<N) ht->data->data[i+2] 	= ht_tmp->data->data[i];
			if(i<N) m1ht->data->data[i] 		= m1ht_tmp->data->data[i];
			if(i+1<N) m2ht->data->data[i+1] = m2ht_tmp->data->data[i];
			if(i+3<N) p1ht->data->data[i+3] = p1ht_tmp->data->data[i];
			if(i+4<N) p2ht->data->data[i+4] = p2ht_tmp->data->data[i];
		}
	}
	
	
	// COPY WAVEFORMS INTO TIMESERIES (+ SHIFT IN CASE OF PARAMID = 4)
	if(ACS->verbose_switch == 1) fprintf(stdout, "... Calculating 5-point derivative - Copy into TimeSeries \n");	
	//for(i=0;i<N;i++)
	//{
	//if( (paramid == 4) && (N-((INT4)htwaveform.f->data->length)+i < N) )
	//{
	//ht->data->data[(N-htwaveform.f->data->length-2)+i] = htwaveform.h->data->data[2*i];
	//m1ht->data->data[(N-htwaveform.f->data->length-3)+i] = m1htwaveform.h->data->data[2*i];
	//m2ht->data->data[(N-htwaveform.f->data->length-4)+i] = m2htwaveform.h->data->data[2*i];
	//p1ht->data->data[(N-htwaveform.f->data->length-1)+i] = p1htwaveform.h->data->data[2*i];
	//p2ht->data->data[(N-htwaveform.f->data->length)+i] = p2htwaveform.h->data->data[2*i];
	
	//}
	//else
	//{
	//if(i<((INT4)htwaveform.f->data->length)) ht->data->data[(N-htwaveform.f->data->length)+i] = htwaveform.h->data->data[2*i];
	//if(i<((INT4)m1htwaveform.f->data->length)) m1ht->data->data[(N-m1htwaveform.f->data->length)+i] = m1htwaveform.h->data->data[2*i];
	//if(i<((INT4)m2htwaveform.f->data->length)) m2ht->data->data[(N-m2htwaveform.f->data->length)+i] = m2htwaveform.h->data->data[2*i];
	//if(i<((INT4)p1htwaveform.f->data->length)) p1ht->data->data[(N-p1htwaveform.f->data->length)+i] = p1htwaveform.h->data->data[2*i];
	//if(i<((INT4)p2htwaveform.f->data->length)) p1ht->data->data[(N-p2htwaveform.f->data->length)+i] = p2htwaveform.h->data->data[2*i];
	//}
	//}
	
	// CALCULATE DERIVATIVE
	if(ACS->verbose_switch == 1) fprintf(stdout, "... Calculating 5-point derivative - Calculate derivative \n");	
	for(i=istart;i<iend;i++)
	{
		//if(paramid == 2) hderiv->data->data[N-1-i] = -1.0*ht->data->data[N-1-i];
		hderiv->data->data[N-1-i] = (m2ht->data->data[N-1-i]-p2ht->data->data[N-1-i]+8.0*p1ht->data->data[N-1-i]-8.0*m1ht->data->data[N-1-i])/(12.0*deltax);
	}
	
	/*********************************************************************
	 * 
	 * Clean up
	 * 
	 *********************************************************************/	
	
	// TIME SERIES
	TRY( XLALDestroyREAL4TimeSeries(ht),	status);
	TRY( XLALDestroyREAL4TimeSeries(m1ht),	status);
	TRY( XLALDestroyREAL4TimeSeries(m2ht),	status);
	TRY( XLALDestroyREAL4TimeSeries(p1ht),	status);
	TRY( XLALDestroyREAL4TimeSeries(p2ht),	status);
	
	TRY( XLALDestroyREAL4TimeSeries(ht_tmp),	status);
	TRY( XLALDestroyREAL4TimeSeries(m1ht_tmp),	status);
	TRY( XLALDestroyREAL4TimeSeries(m2ht_tmp),	status);
	TRY( XLALDestroyREAL4TimeSeries(p1ht_tmp),	status);
	TRY( XLALDestroyREAL4TimeSeries(p2ht_tmp),	status);
	
	// COHERENT GW'S
	TRY( XLALDestroyCoherentGW(&htwaveform), status); 
	TRY( XLALDestroyCoherentGW(&m1htwaveform), status);
	TRY( XLALDestroyCoherentGW(&m2htwaveform), status);
	TRY( XLALDestroyCoherentGW(&p1htwaveform), status);
	TRY( XLALDestroyCoherentGW(&p2htwaveform), status);
	
	// PPN PARAMS
	TRY( XLALDestroyREAL4Vector( htPPNparams.ppn ), status);
	TRY( XLALDestroyREAL4Vector( m1htPPNparams.ppn ), status);
	TRY( XLALDestroyREAL4Vector( m2htPPNparams.ppn ), status);
	TRY( XLALDestroyREAL4Vector( p1htPPNparams.ppn ), status);
	TRY( XLALDestroyREAL4Vector( p2htPPNparams.ppn ), status);
  
	/*********************************************************************
	 * 
	 * Detach LAL error handling
	 * 
	 *********************************************************************/	  
  
  DETATCHSTATUSPTR(status);
  RETURN(status);
	
}

void LALInspiralComputeDerivatives_3PolyChi2(
																						 LALStatus																*status,			// Status pointer
																						 PPNConsistencyParamStruc								*InputParams,	// Input parameters
																						 REAL8																		deltax,				// total width of differentiation
																						 INT4																		paramid,			// parameter choice
																						 INT4																		icur,				// start waveform element
																						 INT4																		iend,					// end waveform element (length)
																						 REAL4TimeSeries													*hderiv,			// derivative calculated with linReg algorithm
																						 INT4                                    xsamples,
																						 FisherACS																*ACS)    
{
	if(ACS->verbose_switch == 1)fprintf(stdout, "Running ComputeDerivatives_3PolyChi2 - paramid = %d | icur = %d | iend = %d | deltax = %e  | xsamples = %d \n", paramid, icur, iend, deltax, xsamples);
	/* ------------------LAL error handling ----------------------------*/
  INITSTATUS( status, "LALInspiralComputeDerivatives_linReg", LALINSPIRALCOMPUTEFISHERC);
  ATTATCHSTATUSPTR(status);
  
  /*------------------ Check Input parameters ------------------------*/
  ASSERT( InputParams, status, LALINSPIRALCOMPUTEFISHERC_EINPUT, LALINSPIRALCOMPUTEFISHERC_EINPUT_MSG );
  ASSERT( hderiv, status, LALINSPIRALCOMPUTEFISHERC_EINPUT, LALINSPIRALCOMPUTEFISHERC_EINPUT_MSG );	
  
	/*********************************************************************
	 * 
	 * Output Files
	 * 
	 ********************************************************************/  
	
	/*********************************************************************
	 * 
	 * Tweakable parameters
	 * 
	 ********************************************************************/ 
  
	INT4 NPoly	= 3;
	
	/*********************************************************************
	 * 
	 * Initiate Variables
	 * 
	 ********************************************************************/
	
	if(ACS->verbose_switch == 1)fprintf(stdout, "... Initialising Variables \n");
	
	/* linReg settings - tune for optimal perfomance */
	
	/* Miscelaneous variables */
	INT4 i = 0;								// Counter	
	INT4 j = 0;								// Counter
	INT4 k = 0;								// Counter
	INT4 N = ((INT4)hderiv->data->length);	// Size of TimeSerie
	
	INT4 hxCur_len = 0;
	INT4 hxPrev_len = 0;
	
  REAL8 parameter		= 0.0;    // Value of current parameter
	REAL8 curX 				= 0.0;		// current value of x
	REAL8 prevX				= 0.0;		// previous value of x
	REAL8 SumX[7]	= {0.0};
	
	
	/* Creating TimeSeries */
	if(ACS->verbose_switch == 1)fprintf(stdout, "... Initialising Variables - Creating TimeSeries \n");
	REAL4TimeSeries							*hxCur;					// storing current temporary waveform
	REAL4TimeSeries							*hxPrev;				// storing previous temporary waveform
	REAL4TimeSeries							*SumYX[NPoly+1];
	
	hxCur   	= XLALCreateREAL4TimeSeries("hxCur", &ACS->epoch, 0.0, InputParams->deltaT, &lalStrainUnit, N);
	hxPrev   	= XLALCreateREAL4TimeSeries("hxPrev", &ACS->epoch, 0.0, InputParams->deltaT, &lalStrainUnit, N);
	
	for (i=0; i<(NPoly+1); i++) {
		SumYX[i] = XLALCreateREAL4TimeSeries("SumYX", &ACS->epoch, 0.0, InputParams->deltaT, &lalStrainUnit, N);
	}
	
	for(i=0; i<N; i++)	{
		hxCur->data->data[i]			= 0.0;
		hxPrev->data->data[i]			= 0.0;
		for (j=0; j<(NPoly+1); j++) {
			SumYX[j]->data->data[i]	= 0.0;
		}
	}
	
	/*********************************************************************
	 * 
	 * linReg Algorithm
   *  - Load in parameters from input
   *  - Increment parameter with deltax (except for dl and tc)
   *  - Filter out band edges (if applicable)
   *  - Pass on for condition checking (next block)
	 * 
	 ********************************************************************/		
	
	/* Loading in relevant parameter through paramid */
  switch(paramid)
  {
    case 0:
			parameter = InputParams->mTot_real8;
			if (parameter == 0.0) ABORT( status, LALINSPIRALCOMPUTEFISHERC_EMBAD, LALINSPIRALCOMPUTEFISHERC_EMBAD_MSG );
      break;
    case 1:
			parameter = InputParams->eta_real8; 
			if (parameter == 0.0) ABORT( status, LALINSPIRALCOMPUTEFISHERC_EMBAD, LALINSPIRALCOMPUTEFISHERC_EMBAD_MSG );     
      break;
			//case 2:
      //break;
    case 2:
      parameter = InputParams->phi;
    case 3:
      break;
    case 4:
			parameter = InputParams->cosI;
			break;
		case 5:
			parameter = InputParams->phasePNparams[ACS->testPhaseParam];
			break;
    default:
      ABORT( status, LALINSPIRALCOMPUTEFISHERC_EINPUT, LALINSPIRALCOMPUTEFISHERC_EINPUT_MSG );
      break;
  }
	
	/* Looping over parameter space */
  for(i=0; i<xsamples; i++)
  {
    if(ACS->verbose_switch==1){fprintf(stdout, "\r... Creating h(x): %d of %d waveforms computed", i, xsamples);}
    if(i==(xsamples-1)){if(ACS->verbose_switch==1){fprintf(stdout, "\n");}}
    
    // Creating hxwaveform structure
    CoherentGW hxwaveform;                // h(x) temp waveform
    memset(&hxwaveform, 0, sizeof(CoherentGW) );
    
    // Creating temporary paramsStruc
    PPNConsistencyParamStruc hxPPNparams;
    hxPPNparams.ppn = XLALCreateREAL4Vector( InputParams->ppn->length );
    TRY( XLALCopyPPNConsistencyParamStruc(InputParams, &hxPPNparams), status);
		
		prevX = parameter +  ((i-xsamples/2)-1) * deltax;
		curX = parameter +  (i-xsamples/2) * deltax;
		
    // Generate and copy waveforms
    switch(paramid)
    {
      case 0:
        hxPPNparams.mTot_real8 = curX;
        
				/* REFILL PHI_I*/
				if(ACS->testPhaseParam!=0) hxPPNparams.phasePNparams[0] = -pow(hxPPNparams.eta_real8,-3.0/8.0)*pow(5.0*LAL_MTSUN_SI*hxPPNparams.mTot_real8,-5.0/8.0); 
				if(ACS->testPhaseParam!=1) hxPPNparams.phasePNparams[1] = -(3715.0/8064.0 + 55.0/96.0*hxPPNparams.eta_real8)*pow(hxPPNparams.eta_real8,-5.0/8.0)*pow(5.0*LAL_MTSUN_SI*hxPPNparams.mTot_real8,-3.0/8.0);
				if(ACS->testPhaseParam!=2) hxPPNparams.phasePNparams[2] = 3.0/4.0*LAL_PI*pow(hxPPNparams.eta_real8,-0.75)*pow(5.0*LAL_MTSUN_SI*hxPPNparams.mTot_real8,-0.25); 
				if(ACS->testPhaseParam!=3) hxPPNparams.phasePNparams[3] = -(9275495.0/14450688.0 + 284875.0/258048.0*hxPPNparams.eta_real8 + 1855.0/2048.0*pow(hxPPNparams.eta_real8,2.0))*pow(hxPPNparams.eta_real8,-7.0/8.0)*pow(5.0*LAL_MTSUN_SI*hxPPNparams.mTot_real8,-1.0/8.0);
				if(ACS->testPhaseParam!=4) hxPPNparams.phasePNparams[4] = -1.0/hxPPNparams.eta_real8*(-38645.0/172032.0 + 65.0/2048.0*hxPPNparams.eta_real8)*LAL_PI*log(hxPPNparams.eta_real8/(5.0*LAL_MTSUN_SI*hxPPNparams.mTot_real8));
				if(ACS->testPhaseParam!=5) hxPPNparams.phasePNparams[5] = -1.0/hxPPNparams.eta_real8*(-38645.0/172032.0 + 65.0/2048.0*hxPPNparams.eta_real8)*LAL_PI; 
				if(ACS->testPhaseParam!=6) hxPPNparams.phasePNparams[6] = -(831032450749357.0/57682522275840.0 - 53.0/40.0*LAL_PI*LAL_PI - 107.0/56.0*LAL_GAMMA + 107.0/448.0*log(hxPPNparams.eta_real8/(256*5.0*LAL_MTSUN_SI*hxPPNparams.mTot_real8)) + (-123292747421.0/4161798144.0 + 2255.0/2048.0*LAL_PI*LAL_PI + 385.0/48.0*(-1987.0/3080.0) - 55.0/16.0*(-11831.0/9240.0))*hxPPNparams.eta_real8 + 154565.0/1835008.0*pow(hxPPNparams.eta_real8,2.0) - 1179625.0/1769472.0*pow(hxPPNparams.eta_real8,3.0))*pow(hxPPNparams.eta_real8,-9.0/8.0)*pow(5.0*LAL_MTSUN_SI*hxPPNparams.mTot_real8,1.0/8.0);
				if(ACS->testPhaseParam!=7) hxPPNparams.phasePNparams[7] = -107.0/448.0*pow(hxPPNparams.eta_real8,-9.0/8.0)*pow(5.0*LAL_MTSUN_SI*hxPPNparams.mTot_real8,1.0/8.0);
				if(ACS->testPhaseParam!=8) hxPPNparams.phasePNparams[8] = -(188516689.0/173408256.0 + 488825.0/516096.0*hxPPNparams.eta_real8 - 141769.0/516096.0*pow(hxPPNparams.eta_real8,2.0))*LAL_PI*pow(hxPPNparams.eta_real8,-5.0/4.0)*pow(5.0*LAL_MTSUN_SI*hxPPNparams.mTot_real8,1.0/4.0); 		        
        
        break;
				
      case 1:
        hxPPNparams.eta_real8 = curX;
        
				/* REFILL PHI_I */
				if(ACS->testPhaseParam!=0) hxPPNparams.phasePNparams[0] = -pow(hxPPNparams.eta_real8,-3.0/8.0)*pow(5.0*LAL_MTSUN_SI*hxPPNparams.mTot_real8,-5.0/8.0); 
				if(ACS->testPhaseParam!=1) hxPPNparams.phasePNparams[1] = -(3715.0/8064.0 + 55.0/96.0*hxPPNparams.eta_real8)*pow(hxPPNparams.eta_real8,-5.0/8.0)*pow(5.0*LAL_MTSUN_SI*hxPPNparams.mTot_real8,-3.0/8.0);
				if(ACS->testPhaseParam!=2) hxPPNparams.phasePNparams[2] = 3.0/4.0*LAL_PI*pow(hxPPNparams.eta_real8,-0.75)*pow(5.0*LAL_MTSUN_SI*hxPPNparams.mTot_real8,-0.25); 
				if(ACS->testPhaseParam!=3) hxPPNparams.phasePNparams[3] = -(9275495.0/14450688.0 + 284875.0/258048.0*hxPPNparams.eta_real8 + 1855.0/2048.0*pow(hxPPNparams.eta_real8,2.0))*pow(hxPPNparams.eta_real8,-7.0/8.0)*pow(5.0*LAL_MTSUN_SI*hxPPNparams.mTot_real8,-1.0/8.0);
				if(ACS->testPhaseParam!=4) hxPPNparams.phasePNparams[4] = -1.0/hxPPNparams.eta_real8*(-38645.0/172032.0 + 65.0/2048.0*hxPPNparams.eta_real8)*LAL_PI*log(hxPPNparams.eta_real8/(5.0*LAL_MTSUN_SI*hxPPNparams.mTot_real8));
				if(ACS->testPhaseParam!=5) hxPPNparams.phasePNparams[5] = -1.0/hxPPNparams.eta_real8*(-38645.0/172032.0 + 65.0/2048.0*hxPPNparams.eta_real8)*LAL_PI; 
				if(ACS->testPhaseParam!=6) hxPPNparams.phasePNparams[6] = -(831032450749357.0/57682522275840.0 - 53.0/40.0*LAL_PI*LAL_PI - 107.0/56.0*LAL_GAMMA + 107.0/448.0*log(hxPPNparams.eta_real8/(256*5.0*LAL_MTSUN_SI*hxPPNparams.mTot_real8)) + (-123292747421.0/4161798144.0 + 2255.0/2048.0*LAL_PI*LAL_PI + 385.0/48.0*(-1987.0/3080.0) - 55.0/16.0*(-11831.0/9240.0))*hxPPNparams.eta_real8 + 154565.0/1835008.0*pow(hxPPNparams.eta_real8,2.0) - 1179625.0/1769472.0*pow(hxPPNparams.eta_real8,3.0))*pow(hxPPNparams.eta_real8,-9.0/8.0)*pow(5.0*LAL_MTSUN_SI*hxPPNparams.mTot_real8,1.0/8.0);
				if(ACS->testPhaseParam!=7) hxPPNparams.phasePNparams[7] = -107.0/448.0*pow(hxPPNparams.eta_real8,-9.0/8.0)*pow(5.0*LAL_MTSUN_SI*hxPPNparams.mTot_real8,1.0/8.0);
				if(ACS->testPhaseParam!=8) hxPPNparams.phasePNparams[8] = -(188516689.0/173408256.0 + 488825.0/516096.0*hxPPNparams.eta_real8 - 141769.0/516096.0*pow(hxPPNparams.eta_real8,2.0))*LAL_PI*pow(hxPPNparams.eta_real8,-5.0/4.0)*pow(5.0*LAL_MTSUN_SI*hxPPNparams.mTot_real8,1.0/4.0);
        
				break;
				
      case 2:
        hxPPNparams.phi = curX;
				break;
				
      case 3:
        break;
				
			case 4:
				hxPPNparams.cosI = curX;
				/* MIRROR AROUND -1 IF CROSSING -1 */
				if(hxPPNparams.cosI < -1.0 ){ hxPPNparams.cosI = -2.0 - hxPPNparams.cosI; hxPPNparams.sinI = -hxPPNparams.sinI;}
				/* MIRROR AROUND +1 IF CROSSING +1 */
				if(hxPPNparams.cosI > 1.0 ){ hxPPNparams.cosI = 2.0 - hxPPNparams.cosI; hxPPNparams.sinI = -hxPPNparams.sinI;}
				break;
				
			case 5:
				hxPPNparams.phasePNparams[ACS->testPhaseParam] = curX;
				break;
				
      default:
        if(ACS->verbose_switch == 1)fprintf(stdout, "XLALInspiralComputeDerivatives failed: invalid derivative variable \n");
        break;
    }
		
    // Create waveforms
    TRY( LALGeneratePPNAmpCorConsistency( status->statusPtr, &hxwaveform, &hxPPNparams), status);
		LALInspiralCombinePlusCross( status->statusPtr, hxwaveform, hxPPNparams, hxCur, ACS);
		
    // Get length and check if vector is long enough to hold waveform
    hxCur_len = hxwaveform.f->data->length;
		ASSERT(hxCur_len<N, status, LALINSPIRALCOMPUTEFISHERC_ETSCREATE, LALINSPIRALCOMPUTEFISHERC_ETSCREATE_MSG);
    
    if(			(paramid == 2 || paramid == 3 || paramid == 4 || paramid == 5)
			 ||	(paramid == 1 && ((INT4)i) != 0 && hxCur_len != hxPrev_len && curX <= 0.25) 
			 ||	(paramid == 0 && ((INT4)i) != 0 && hxCur_len != hxPrev_len ) )
    {
      for (j=0; j<(2*NPoly+1); j++) {
				SumX[j] += pow(curX-parameter,((double)j));
			}
			
			for (j=0; j<(NPoly+1); j++) {
				for (k=0; k<N; k++) {
					SumYX[j]->data->data[k] += hxCur->data->data[k]*pow(curX-parameter,((double)j));
				}
			}
			
    }
    
    // Loop completed, storing cur into old variables, empty cur
    hxPrev_len = hxCur_len;
    
    // move waveform into old waveform and clean new waveform
    for(k=0; k<((INT4) N) ; k++)
    {
      hxPrev->data->data[k] = hxCur->data->data[k];
      hxCur->data->data[k] = 0.0;
    }
    
    TRY( XLALDestroyCoherentGW(&hxwaveform), status); 
    TRY( XLALDestroyREAL4Vector( hxPPNparams.ppn ), status);
    
  }
  
  /* Computing derivative */
  if(ACS->verbose_switch == 1)fprintf(stdout, "... Calculating Derivative \n");
  for(i=icur;i<iend;i++)
  {		
		hderiv->data->data[i] = (SumX[1]*pow(SumX[5],2.0)*SumYX[0]->data->data[i] - SumX[1]*SumX[4]*SumX[6]*SumYX[0]->data->data[i] - SumX[0]*pow(SumX[5],2.0)*SumYX[1]->data->data[i] + SumX[0]*SumX[4]*SumX[6]*SumYX[1]->data->data[i] + pow(SumX[3],3.0)*SumYX[2]->data->data[i] + SumX[0]*SumX[4]*SumX[5]*SumYX[2]->data->data[i] - SumX[0]*pow(SumX[4],2.0)*SumYX[3]->data->data[i] - pow(SumX[3],2.0)*(SumX[5]*SumYX[0]->data->data[i]+SumX[4]*SumYX[1]->data->data[i]+SumX[2]*SumYX[3]->data->data[i]) + pow(SumX[2],2.0)*(-1.0*SumX[6]*SumYX[1]->data->data[i]+SumX[4]*SumYX[3]->data->data[i]) + SumX[3]*(pow(SumX[4],2.0)*SumYX[0]->data->data[i]+SumX[2]*SumX[6]*SumYX[0]->data->data[i]+2.0*SumX[2]*SumX[5]*SumYX[1]->data->data[i]-SumX[2]*SumX[4]*SumYX[2]->data->data[i]-SumX[1]*SumX[5]*SumYX[2]->data->data[i]-SumX[0]*SumX[6]*SumYX[2]->data->data[i]+SumX[1]*SumX[4]*SumYX[3]->data->data[i]+SumX[0]*SumX[5]*SumYX[3]->data->data[i]) - SumX[2]*(SumX[4]*SumX[5]*SumYX[0]->data->data[i]-SumX[1]*SumX[6]*SumYX[2]->data->data[i]+SumX[1]*SumX[5]*SumYX[3]->data->data[i]))/(pow(SumX[3],4.0) + pow(SumX[2],2.0)*pow(SumX[4],2.0) - SumX[0]*pow(SumX[4],3.0) + pow(SumX[1],2.0)*pow(SumX[5],2.0) - SumX[2]*SumX[5]*(2.0*SumX[1]*SumX[4]+SumX[0]*SumX[5]) - pow(SumX[2],3.0)*SumX[6] - pow(SumX[1],2.0)*SumX[4]*SumX[6] + SumX[0]*SumX[2]*SumX[4]*SumX[6] - pow(SumX[3],2.0)*(3.0*SumX[2]*SumX[4]+2.0*SumX[1]*SumX[5]+SumX[0]*SumX[6]) + 2.0*SumX[3]*(SumX[5]*(pow(SumX[2],2.0)+SumX[0]*SumX[4])+SumX[1]*(pow(SumX[4],2.0)+SumX[2]*SumX[6])) );
  }
	
	
	/*********************************************************************
	 * 
	 * Clean Up
	 * 
	 ********************************************************************/	
	if(ACS->verbose_switch == 1)fprintf(stdout, "... Cleaning up \n");
  
  if(ACS->verbose_switch == 1)fprintf(stdout, "... Cleaning up - Closing files\n");
	
	if(ACS->verbose_switch == 1)fprintf(stdout, "... Cleaning up - TimeSeries\n");
	TRY( XLALDestroyREAL4TimeSeries(hxCur),			status);
	TRY( XLALDestroyREAL4TimeSeries(hxPrev),		status);
	for (i=0; i<(NPoly+1); i++) {
		XLALDestroyREAL4TimeSeries(SumYX[i]);
	}
	
	/*********************************************************************
	 * 
	 * Detach LAL error handling
	 * 
	 ********************************************************************/	
  DETATCHSTATUSPTR( status );
  RETURN(status);
}


void LALInspiralCombinePlusCross(
																 LALStatus															*status,			// Status pointer
																 CoherentGW														GravWav,			// Input waveform, data[2i] = hp, data[2i+1] = hc
																 PPNConsistencyParamStruc													params,				// 
																 REAL4TimeSeries												*ht,
																 FisherACS															*ACS)					// combined h(t) from hp and hc
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
  
  INITSTATUS( status, "LALInpiralCombinePlusCross", LALINSPIRALCOMPUTEFISHERC);
  ATTATCHSTATUSPTR(status);	
  
	/*********************************************************************
	 * 
	 * Check input parameters
	 * 
	 *********************************************************************/	  
  
	ASSERT( ht, status, LALINSPIRALCOMPUTEFISHERC_EINPUT, LALINSPIRALCOMPUTEFISHERC_EINPUT_MSG );
	ASSERT( GravWav.h, status, LALINSPIRALCOMPUTEFISHERC_EINPUT, LALINSPIRALCOMPUTEFISHERC_EINPUT_MSG );
	
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
	
	if(		(hp = (REAL4TimeSeries *) LALMalloc( sizeof(REAL4TimeSeries) ) ) == NULL
		 &&	(hc = (REAL4TimeSeries *) LALMalloc( sizeof(REAL4TimeSeries) ) ) == NULL)
	{
		LALFree(hp);			hp = NULL;
		LALFree(hc);			hc = NULL;
		ABORT(status, LALINSPIRALCOMPUTEFISHERC_EMEM, LALINSPIRALCOMPUTEFISHERC_EMEM_MSG);
	}
	else
	{
		hp = XLALCreateREAL4TimeSeries("hp", &(ht->epoch), 0.0, ht->deltaT, &lalStrainUnit, ht->data->length);
		hc = XLALCreateREAL4TimeSeries("hc", &(ht->epoch), 0.0, ht->deltaT, &lalStrainUnit, ht->data->length);
	}
	
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
  
  if(ht->data->length > GravWav.f->data->length)
  {
		for(i=0;i<GravWav.f->data->length;i++)
		{
			hp->data->data[ht->data->length-GravWav.f->data->length+i] = GravWav.h->data->data[2*i];
			hc->data->data[ht->data->length-GravWav.f->data->length+i] = GravWav.h->data->data[2*i+1];
		}
	}
	else
	{
		for(i=0;i<ht->data->length;i++)
		{
			hp->data->data[i] = GravWav.h->data->data[2*i];
			hc->data->data[i] = GravWav.h->data->data[2*i+1];
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
  
	LALDetector detector = lalCachedDetectors[LAL_LHO_4K_DETECTOR]; /* H1 */
	LALSource source;
	LIGOTimeGPS epoch;
	LALDetAndSource det_source; 
	LALDetAMResponse det_resp; 									
  
  switch(ACS->coordinate)
  {
		case 0:
			theta = params.position.latitude;
			phi		= params.position.longitude;
			psi		= params.psi;
			
			Fp = 1.0/2.0 * ( 1 + pow(cos(theta),2.0) ) * cos(2.0*phi) * cos(2.0*psi)
			- cos(theta) * sin(2.0*phi) * sin(2.0*psi);
			
			
			//fprintf(stdout, "Fp = %e \n", Fp);
			
			Fc =  1.0/2.0 * (1 + pow(cos(theta),2.0) ) * cos(2.0*phi) * sin(2.0*psi)
			+ cos(theta) * sin(2.0*phi) * cos(2.0*psi);
			
			//fprintf(stdout, "Fc = %e \n", Fc);		
			break;
			
		case 1:
			/* LOAD DETECTOR INFORMATION */
			
			
			/* The position and polarisation of the binary */
			
			source.equatorialCoords.longitude = params.position.longitude;
			source.equatorialCoords.latitude = params.position.latitude;
			
			/* Specify EQUATORIAL coords so that long = R.A. and lat=dec */
			source.equatorialCoords.system = COORDINATESYSTEM_EQUATORIAL;
			source.orientation = params.psi;
			
			/* Time of the event */
			
			epoch = params.epoch; 
			
			/* Create a detector and source structure */
			
			det_source.pSource		= &source; /* Set source */
			det_source.pDetector	= &detector; /* detector */	
			
			/* LALDETECTOR RESPONSE */
			
			
			/* Compute det_resp */
			LALComputeDetAMResponse(status->statusPtr, &det_resp, &det_source, &epoch);
			
			Fp = det_resp.plus;
			Fc = det_resp.cross;  		
			break;
			
		default:
			printf("ERROR - INVALID COORDINATE SYSTEM \n");
			break;
			
	}
  
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

void InvertMatrixSVD (
											gsl_matrix        *A,
											gsl_matrix				*InvA,
											int								N)
{ 
	
  /*********************************************************************
   *
   *  CREATING TEMPORARY VARIABLES
   * 
   ********************************************************************/  
	
	/* COUNTERS */
  int i = 0;
  
  // Initialise matrices A, U, S and V
  gsl_matrix *InvS  = gsl_matrix_calloc (N, N); // inverse S
  gsl_matrix *V     = gsl_matrix_calloc (N, N); // V
  gsl_matrix *U     = gsl_matrix_calloc (N, N); // U
  gsl_matrix *C     = gsl_matrix_calloc (N, N); // temporary storage
  gsl_matrix *I     = gsl_matrix_calloc (N, N); // testing idenity
  gsl_vector *s     = gsl_vector_alloc (N);     // eigenvalues AA^T
	
  //printf("INPUT \n");
	//for(i=0;i<(N*N);i++)
	//{
	//printf("%e", gsl_matrix_get(A, i/N, i%N));
	//if(i%N==(N-1)) printf("\n");
	//else printf("\t");
	//}  
	//printf("\n");
  
  /*********************************************************************
   *
   *  COMPUTING INVERSE
   * 		- PERFORM SVD
   * 		- CALCULATE INVERSE
   * 
   ********************************************************************/ 
	
	// Prepare U for SVD
	gsl_matrix_memcpy(U, A);
	
	// Perform SVD
	gsl_linalg_SV_decomp_jacobi(U, V, s);  
	
	// Compute Inverse S
	for (i = 0; i<N; i++)
	{
		gsl_vector_set( s, i, 1./gsl_vector_get( s, i) );
		gsl_matrix_set( InvS, i, i, gsl_vector_get( s, i) );
	}
	
	//printf("EIGENVECTORS \n");
	//for(i=0;i<N;i++)
	//{
	//printf("%e", gsl_vector_get(s,i));
	//if(i==(N-1)) printf("\n");
	//else printf("\t");
	//}
	//printf("\n");
	
	// Tranpose U
	gsl_matrix_transpose(U);
	
	// Multiply V and InvS
	gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
									1.0, V, InvS,
									0.0, C);
	
	gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
									1.0, C, U,
									0.0, InvA);                             
  
  //printf("INVERSE \n");
	//for(i=0;i<(N*N);i++)
	//{
	//printf("%e", gsl_matrix_get(InvA, i/N, i%N));
	//if(i%N==(N-1)) printf("\n");
	//else printf("\t");
	//}  
	//printf("\n");  
  
  /*********************************************************************
   *
   *  TESTING ACCURACY
   * 		- A * INVA = 1
   * 
   ********************************************************************/  
  
  gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
									1.0, A, InvA,
									0.0, I);      
	
	//printf("UNIT\n");
	//for(i=0;i<(N*N);i++)
	//{
	//printf("%e", gsl_matrix_get(I, i/N, i%N));
	//if(i%N==(N-1)) printf("\n");
	//else printf("\t");
	//}
	//printf("\n");
  
  /*********************************************************************
   *
   *  CLEANING UP
   * 
   ********************************************************************/  
	
	/* MATRICES */
  //gsl_matrix_free(A);
  gsl_matrix_free(U);
  gsl_matrix_free(InvS);
  //gsl_matrix_free(InvA);
  gsl_matrix_free(V);
  gsl_matrix_free(C);
  gsl_matrix_free(I);
  gsl_vector_free(s);
  
  /*********************************************************************
   *
   *  Detach Error handling
   * 
   ********************************************************************/
  
  return;
}

void LALSavitskyGolayFilter(
														LALStatus						*status,
														INT4								R,		/* Polynomial fitter order */
														INT4								M,		/* time steps forward/backward sampled */
														REAL4TimeSeries			*y,		/* Time Series to be filtered */
														REAL4TimeSeries			*Output, /* Output of Filter */
														INT4 								D			/* Order of derivative outputted */
														)
{
	
	
  /*********************************************************************
   *
   * Savitsky-Golay Filter
   * - Based on least square fitting of a polynomial of Rth order
   * - Smoothens function by extrapolating from m neighbouring points
   * - See Abraham Savitsky and Marcel J. E. Golay 
   * 		"Smoothing and differentiation of data by simplified least squares procedures"
   * 
   ********************************************************************/ 
	
	/*********************************************************************
	 * 
	 * LAL error handling
	 * 
	 *********************************************************************/	
  
  INITSTATUS( status, "LALSavitskyGolayFilter", LALINSPIRALCOMPUTEFISHERC);
  ATTATCHSTATUSPTR(status);	 
	
	/*********************************************************************
	 * 
	 * Read input
	 * 
	 *********************************************************************/	   
	
  /* LOAD IN TIMESERIES PARAMETERS */
  REAL8 dt = y->deltaT;
  INT4 N = y->data->length;
  //printf("Input parameters: R = %d | M = %d | N = %d | dt = %e \n", R, M, N, dt);	
	
  /*********************************************************************
   *
   * Create Temporary Variables
   * 
   ********************************************************************/  
  
  //printf("Initialising Variables \n");
  
  /* COUNTERS */
  int i,j,k;
  
  /* factorial of D (used for derivatives) */
  INT4 factorial = 1;
  
  /* MATRICES AND VECTORS */
	gsl_matrix *m     	= gsl_matrix_calloc (R+1, 2*M+1);   /* m_ij = j^i */
	gsl_matrix *U				= gsl_matrix_calloc (R+1, R+1);		/* U_ii = deltaT^i */
	gsl_vector *a				= gsl_vector_calloc (R+1);			/* a_j, for y(t_i) = Sum_{j=0}^{R} a_j t_i^j */
	gsl_matrix *c				= gsl_matrix_calloc (R+1, 2*M+1);		/* c_ij = U_-1 (m m^T)^-1 m, in a_j = c_ji y_i */ 
	gsl_vector *ym			= gsl_vector_calloc (2*M+1);	/* y_m = [y_-M, ... , y_M] */
	gsl_matrix *tmr			= gsl_matrix_calloc (R+1, 2*M+1);		/* t_m^r = U*m */
	
	/* COMBINED MATRICES AND VECTORS */
	gsl_matrix *mT			= gsl_matrix_calloc (2*M+1, R+1);		/* m^T */
	gsl_matrix *mmT			= gsl_matrix_calloc (R+1, R+1);		/* mm^T */
	gsl_matrix *InvmmT	= gsl_matrix_calloc (R+1, R+1);		/* (mm^T)^-1 */
	gsl_matrix *InvmmTm	= gsl_matrix_calloc (R+1, 2*M+1);		/* (mm^T)^-1 m */
	gsl_matrix *InvU		= gsl_matrix_calloc (R+1, R+1);		/* U^-1 */
	
  /*********************************************************************
   *
   * Filling matrices
   * 
   ********************************************************************/ 
	//printf("Filling parameters \n");
	
  /* m_ij = j^i */
  //printf("Filling parameters -  m %dx%d\n", m->size1, m->size2);
  for(i=0;i<(R+1);i++)
  {
		for(j=0;j<(2*M+1);j++)
		{
			gsl_matrix_set(m, i, j, pow((j-M), i));
		}
	}
	
  //printf("m %dx%d\n", m->size1, m->size2);
	//for(i=0;i<((R+1)*(2*M+1));i++)
	//{
	//printf("%e", gsl_matrix_get(m, i/(2*M+1), i%(2*M+1)));
	//if(i%(2*M+1)==(2*M)) printf("\n");
	//else printf("\t");
	//}  
	//printf("\n");	 
	
  /* U_ii = deltaT^i */
  //printf("Filling parameters -  U %dx%d \n", U->size1, U->size2);
  for(i=0;i<(R+1);i++)
  {
		gsl_matrix_set(U, i, i, pow(dt,i));
	} 
	
  //printf("U %dx%d\n", U->size1, U->size2);
	//for(i=0;i<((R+1)*(R+1));i++)
	//{
	//printf("%e", gsl_matrix_get(U, i/(R+1), i%(R+1)));
	//if(i%(R+1)==(R)) printf("\n");
	//else printf("\t");
	//}  
	//printf("\n");	 	
	
	/* m^T */
	//printf("Filling parameters -  mT %dx%d\n", mT->size1, mT->size2);
	for(i=0;i<(R+1); i++)
	{
		for(j=0;j<(2*M+1);j++)
		{
			gsl_matrix_set(mT, j, i, gsl_matrix_get(m, i, j));
		}
	}
	
  //printf("mT %dx%d\n", mT->size1, mT->size2);
	//for(i=0;i<((2*M+1)*(R+1));i++)
	//{
	//printf("%e", gsl_matrix_get(mT, i/(R+1), i%(R+1)));
	//if(i%(R+1)==(R)) printf("\n");
	//else printf("\t");
	//}  
	//printf("\n");	 	
	
	/* mm^T */
	//printf("Filling parameters -  mmT %dx%d\n", mmT->size1, mmT->size2);
	gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
									1.0, m, mT,
									0.0, mmT);	
	
  //printf("mmT %dx%d\n", mmT->size1, mmT->size2);
	//for(i=0;i<((R+1)*(R+1));i++)
	//{
	//printf("%e", gsl_matrix_get(mmT, i/(R+1), i%(R+1)));
	//if(i%(R+1)==(R)) printf("\n");
	//else printf("\t");
	//}  
	//printf("\n");							
	
	/* (mm^T)^-1 */
	//printf("Filling parameters -  InvmmT %dx%d\n", InvmmT->size1, InvmmT->size2);
	InvertMatrixSVD(mmT, InvmmT, R+1);
	
	/* U^-1*/
	//printf("Filling parameters -  InvU %dx%d\n", InvU->size1, InvU->size2);
	//InvertMatrixSVD(U, InvU, R+1);
	
	for(i=0;i<(R+1);i++)
	{
		//printf("%e | %e \n", 1.0/gsl_matrix_get(U, i, i), gsl_matrix_get(InvU, i, i));
		gsl_matrix_set(InvU, i, i, 1.0/gsl_matrix_get(U, i, i));
	}
	
	/* (mm^T)^-1 m */
	//printf("Filling parameters -  InvmmTm %dx%d \n", InvmmTm->size1, InvmmTm->size2 );
	gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
									1.0, InvmmT, m,
									0.0, InvmmTm);	
	
  //printf("InvmmTm \n");
	//for(i=0;i<((R+1)*(2*M+1));i++)
	//{
	//printf("%e", gsl_matrix_get(InvmmTm, i/(2*M+1), i%(2*M+1)));
	//if(i%(2*M+1)==(2*M)) printf("\n");
	//else printf("\t");
	//}  
	//printf("\n");	
	
	/* c_ij = U_-1 (m m^T)^-1 m */
	//printf("Filling parameters -  c \n");
	gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
									1.0, InvU, InvmmTm,
									0.0, c);
	
  //printf("c \n");
	//for(i=0;i<((c->size1)*(c->size2));i++)
	//{
	//printf("%e", gsl_matrix_get(c, i/(c->size2), i%(c->size2)));
	//if(i%(c->size2)==(c->size2-1)) printf("\n");
	//else printf("\t");
	//}  
	//printf("\n");	
	
	/* t_m^r = U*m */
	//printf("%dx%d -> (%dx%d)x(%dx%d)\n", tmr->size1, tmr->size2, U->size1, U->size2, m->size1, m->size2);
	//printf("Filling parameters -  tmr \n");
	gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
									1.0, U, m,
									0.0, tmr);	
	
	
	/*********************************************************************
   *
   * Set polynomial prefactors and smooth function
   * 
   ********************************************************************/ 
  
  for(i=D;i>0;i--){factorial *= i;}
  
  //printf("Smoothing \n");
  for(k=M;k<(N-M);k++)
  {
		/* READ DATA POINTS INTO VECTOR */
		for(j=0;j<2*M+1;j++)
		{
			gsl_vector_set(ym, j, y->data->data[k-M+j]);		
		}
		
		/* a = c*y */
		gsl_blas_dgemv( CblasNoTrans, 
									 1.0, c, ym, 
									 0.0, a );
		
		Output->data->data[k] = factorial*gsl_vector_get(a, D);
		
		/* CLEAN UP IN THE LOOP */
		gsl_vector_set_zero (ym);
		gsl_vector_set_zero (a);
	}
	
  /*********************************************************************
   *
   * Output to file
   * 
   ********************************************************************/  
  //printf("Write to file \n");
  //FILE *smoothOut;
  //smoothOut = fopen("smoothOut.dat", "w");
  
  //for(k=0;k<N;k++)
  //{
	//fprintf(smoothOut, "%e\t%e\n", k*dt, Output[k]);
	//}
	//fclose(smoothOut);
	
  /*********************************************************************
   *
   * Clean Up
   *  
   ********************************************************************/  
	
	//printf("Cleaning up \n");
	gsl_matrix_free(m);
	gsl_matrix_free(U);
	gsl_vector_free(a);
	gsl_matrix_free(c);  
	gsl_vector_free(ym);
	gsl_matrix_free(tmr);
	gsl_matrix_free(mT);
	gsl_matrix_free(mmT);
	gsl_matrix_free(InvmmT);
	gsl_matrix_free(InvmmTm);
	gsl_matrix_free(InvU);
	
	DETATCHSTATUSPTR(status);
	RETURN(status);
}


void LALStrainInGeocentricTime(
															 LALStatus				*status,
															 LALDetAndSource	*det_source, 
															 CoherentGW				waveform,
															 REAL4TimeSeries	*output)
{
	/*********************************************************************
	 * 
	 * LAL error handling
	 * 
	 *********************************************************************/	
  
  INITSTATUS( status, "LALStrainInGeocentricTime", LALINSPIRALCOMPUTEFISHERC);
  ATTATCHSTATUSPTR(status);	
	
	/*********************************************************************
	 * 
	 * Check Input
	 * 
	 *********************************************************************/

	
		
	/*********************************************************************
	 * 
	 * Define Variables
	 * 
	 *********************************************************************/
	
	INT4 i;
	
	REAL8 TimeDiffFromGC = 0.0;	// Compute time difference between detector (to be specified) and the geocentre
	INT4 TimeDiffFromGCInt = 0;
	LALDetAMResponse det_resp;	// Structure holding the results of LALComputeDetAMResponse; det_resp.plus and det_resp.cross are the beam pattern functions
	LIGOTimeGPS *det_time = NULL;
	memcpy(det_time, &(waveform.h->epoch), sizeof(LIGOTimeGPS));
	
	REAL4TimeSeries *hp, *hc;
	hp = XLALCreateREAL4TimeSeries("hp",&(output->epoch), output->f0, output->deltaT, &lalStrainUnit, output->data->length);
	hc = XLALCreateREAL4TimeSeries("hc",&(output->epoch), output->f0, output->deltaT, &lalStrainUnit, output->data->length);
	
	//hp = XLALCreateREAL4TimeSeries("hp", &(ht->epoch), 0.0, ht->deltaT, &lalStrainUnit, ht->data->length);
	
	/*********************************************************************
	 * 
	 * Compute Time Difference and correct epoch in coherentGW
	 * 
	 *********************************************************************/
	
	TimeDiffFromGC = XLALTimeDelayFromEarthCenter( det_source->pDetector->location, det_source->pSource->equatorialCoords.longitude, det_source->pSource->equatorialCoords.latitude, &(waveform.h->epoch));
	
	det_time = XLALGPSAdd(&(waveform.h->epoch), TimeDiffFromGC);
	
	TimeDiffFromGCInt = TimeDiffFromGC / output->deltaT;
	
	for (i=0; i<((INT4)output->data->length)-1; i++) {
		
		if( ( (i+TimeDiffFromGCInt) < (((INT4)output->data->length)-1) ) && ( i+TimeDiffFromGCInt > 0 ) )
		{
			hp->data->data[((INT4)(0.9*output->data->length-1))-i+TimeDiffFromGCInt] = waveform.h->data->data[((INT4)waveform.h->data->length-1)-2*i];
			hc->data->data[((INT4)(0.9*output->data->length-1))-i+TimeDiffFromGCInt] = waveform.h->data->data[((INT4)waveform.h->data->length-1)-2*i+1];
		}
		else exit(-1);
	}
	
	
	/*********************************************************************
	 * 
	 * Compute beam pattern function & combine hplus and hcross
	 * 
	 *********************************************************************/
	
	LALComputeDetAMResponse(status->statusPtr, &det_resp, det_source, det_time );
	
	for (i=0; i<((INT4) output->data->length)-1 ; i++) {
		output->data->data[i] = det_resp.plus*hp->data->data[i] + det_resp.cross*hc->data->data[i];
	}
	
	
	/*********************************************************************
	 * 
	 * Output to file (test)
	 * 
	 *********************************************************************/
	
	/*********************************************************************
	 * 
	 * Detach function and returning status
	 * 
	 *********************************************************************/	
	DETATCHSTATUSPTR(status);
	RETURN(status);
}


/*  version tracker
 v1: Base version with both Interband and AdaptiveDelta Method
 v2: Write LinReg based on Interband
 remove AdaptiveDelta as it was shown to perform worse
 v3: Include adaptive deltax according to a phase check
 put algorithm in different function for looping ability
 v4: Reinstate 5 point derivative, mainly for the calculation of the 
 time derivative.
 v5: Include Algorithm Control Structure (for versatile code tuning)
 v6: Beam pattern function for h+ and hx
 RNG for inclination angle	
 v7: -
 v8: Include LSO dependent length of waveform
 Include Savitsky-Golay smoothing filter
 v9: Include derivatives for PN phase parameters
 v10: Include higher order chi2 derivatives
 
 */
