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

#if 0
<lalVerbatim file="FisherTestCV">
Author: Tjonnie Li, Chris Van Den Broeck
$Id$
</lalVerbatim>

<lalLaTeX>

\subsection{Program \texttt{FisherTest.c}}
\label{ss:FisherTest.c}

Calculates the Fisher matrix at a single point in parameter space.

\subsubsection*{Usage}

\subsubsection*{Description}

\subsubsection*{Exit codes}

\subsubsection*{Algorithm}

\subsubsection*{Uses}
\begin{verbatim}
lalDebugLevel
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{FisherTestCV}}

</lalLaTeX>
#endif


#include <lal/LALInspiralComputeFisher.h>

NRCSID(FISHERTESTC,"$Id$");

int lalDebugLevel = 1;

int main( INT4 argc, char **argv )
{
  fprintf(stdout, "Running FisherTest \n"); 
  /*********************************************************************
   *
   *  Error Handling
   * 
   ********************************************************************/
  
	static LALStatus status;
  
  /*********************************************************************
   *
   * Get seed from argument of executable and set up Random Number
   * Generator (RNG)
   * 
   ********************************************************************/	

	/* GET SEED FROM ARGUMENT */
	
	if (argc!=2)
  {
    fprintf(stdout, "Usage: %s random seed\n", argv[0]);
    exit(0);
  }
  INT4 seed = atoi(argv[1]);
  
  /* INITIATE GSL RANDOM NUMBER GENERATOR */
  
	const gsl_rng_type *type;           // RNG type
  gsl_rng *p;                         // Generator
  gsl_rng_env_setup();                // Setup environment
  gsl_rng_default_seed = seed;		    // vary generation sequence
  type = gsl_rng_default;             // set RNG type to default
  p = gsl_rng_alloc (type);           // Set RNG type  
  
  /*********************************************************************
   *
   * Input parameters for waveform generation
   * 
   ********************************************************************/	  
  
	// WAVEFORM PARAMETERS
  REAL8 m1            = 2.0;           // Component Mass 1 
  REAL8 m2            = 3.0; 					// Component Mass 2
		//m1*((1.0-2.0*eta)+sqrt(1.0-4*eta))/(2*eta);        
  REAL8 inc           = LAL_PI/7.0;     // Inclination angle
  REAL8 tc            = 0.0;            // End time of waveform
  REAL8 phic          = 0.0;            // Phase of coalescenes
  REAL8 dist          = 1.0e8*LAL_PC_SI;// Distance to source
  REAL8 fStartIn      = 40.0;           // Starting frequency
  REAL8 fStopIn       = 0.0;            // No early termination
	INT4 phaseOrder     = 6;              // Expansion order PPN
  INT4 ampOrder       = 3;              // Amplitude Order    
  
  // LOADING PARAMETER INTO STRUCTURE
  PPNParamStruc params;
  params.position.longitude = acos(gsl_ran_flat (p, -1.0, 1.0));
  params.position.latitude = gsl_ran_flat (p, 0, LAL_TWOPI);
  params.position.system = COORDINATESYSTEM_GEOGRAPHIC;
  params.psi			= acos(gsl_ran_flat(p,-1.0,1.0));
  params.deltaT = 1.0/16384.0;
  params.mTot_real8 = m1 + m2;
  params.mTot = m1 + m2;
  params.eta_real8  = m1 * m2 / pow(m1 + m2, 2.);
  params.inc = inc;
  params.cosI = cos(inc);
  params.sinI = sin(inc);
  params.tc  = tc;
  params.phi = phic;
  params.d   = dist;
  params.fStartIn = fStartIn;
  params.fStopIn  = fStopIn;
  params.ampOrder = ampOrder;
  
  INT4 i;
  params.ppn = XLALCreateREAL4Vector( phaseOrder+1 );
  for (i=0; i<=(phaseOrder); i++)
  {
    if(phaseOrder > 0 && i==1) {params.ppn->data[i] = 0.0;}
    else {params.ppn->data[i] = 1.0;}
    
	}
	
  // PSD MODEL
  void (*noisemodel)(LALStatus*,REAL8*,REAL8) = LALLIGOIPsd;	
	
  /*********************************************************************
   *
   * Algorithm Control Structure
   * 
   ********************************************************************/	
	
	/* ALGORITHM CONTROL STRUCTURE */
	FisherACS testACS;
	
	/* GENERAL CONTROL PARAMETERS */
	testACS.verbose_switch		=	1;
	testACS.printall					= 1;
	//testACS.epoch							= LIGOTIMEGPSZERO;
	testACS.N_inputparams			=	6;
	testACS.dynRange					= 1E+27;
	
	/* TIME/FREQUENCY SERIES */
	testACS.N_ts							= 1E+5;
	testACS.N_fs							= testACS.N_ts/2+1;
	testACS.deltaT_ts					= 1.0/16384.0;
	testACS.deltaF_fs					= 1.0/((REAL8)testACS.N_ts * testACS.deltaT_ts);
	
	/* DIFFERENTIATION PARAMETERS */
	testACS.xsamples_tot						= 400;
	testACS.linReg_countMax				= testACS.xsamples_tot/2;
	testACS.linReg_countMin				= 20;
	testACS.linReg_loopMax					= 100;
	testACS.deltaPhi_max					= LAL_PI/25;
	testACS.reldelta							= 5E-6;
	testACS.absdelta_adjustUp			= 1.5;
	testACS.absdelta_adjustDown		= 0.5;
	testACS.xsamples_adjustUp			= 2;
	testACS.xsamples_adjustDown		= 5;
	
	/* DETECTOR FREQUENCY WINDOW */
	testACS.fstart								= 40.0;
	testACS.fstop									= (pow(LAL_C_SI, 3.0))/(pow(6.0,1.5)*LAL_PI*LAL_G_SI*(m1+m2)*LAL_MSUN_SI);		
	
  /*********************************************************************
   *
   *  Output Files
   * 
   ********************************************************************/  
  
  if(testACS.verbose_switch==1){fprintf(stdout, "... Creating Output Files \n");}
  
  /* FOLDER FOR OUTPUT FILES */
	char folder[128] = "/home/spxcv/src/lalsuite/lalinspiral/test/";
  
  /* NAME OUTPUT FILES */
	const char powerspectral[]		= "psd.dat";
	const char modsq[]						= "modsq.dat";	
  const char derivatives[] 			= "derivatives.dat";
  const char fourierderivs[] 		= "fourierderivatives.dat";
  const char fisher[]						= "fisher";
  const char covariance[]				= "covariance";
  const char hx[]								= "hx";
  
  /* CREATE OUTPUT FILES */
  FILE *psd_out;
  FILE *modsq_out=NULL;
  FILE *derivatives_out=NULL; 
  FILE *fourierderivs_out=NULL;  
  FILE *fisher_out=NULL;
  FILE *cov_out=NULL;  
  FILE *hx_out=NULL;
  
	char psd_name[128];
  char modsq_name[128];
	char derivatives_name[128];
	char fourierderivs_name[128];
	char fisher_name[128];
  char cov_name[128];
  char hx_name[128];
  
  sprintf(psd_name, "%s%s", folder, powerspectral);
	sprintf(modsq_name, "%s%s", folder, modsq);
	sprintf(derivatives_name, "%s%s", folder, derivatives);
	sprintf(fourierderivs_name, "%s%s", folder, fourierderivs);
	sprintf(fisher_name, "%s%s%d%s", folder, fisher, seed, ".dat");
	sprintf(cov_name, "%s%s%d%s", folder, covariance, seed, ".dat");
	sprintf(hx_name, "%s%s%d%s", folder, hx, seed, ".dat");
	
	if(testACS.printall == 1)
	{
		psd_out						= fopen(psd_name, "w");
		modsq_out					= fopen(modsq_name, "w");
		derivatives_out		= fopen(derivatives_name, "w");
		fourierderivs_out	= fopen(fourierderivs_name, "w");  
		hx_out						= fopen(hx_name, "w");
	}
  fisher_out				= fopen(fisher_name, "a");
  cov_out						= fopen(cov_name, "a");
	
	/* INCLUDE OUTPUT FILES INTO ACS */
  testACS.psd_out						= psd_out;
	testACS.modsq_out					= modsq_out;
  testACS.derivatives_out		= derivatives_out;
  testACS.fourierderivs_out	= fourierderivs_out;  
  testACS.fisher_out				= fisher_out;
  testACS.cov_out						= cov_out;
  testACS.hx_out						= hx_out;
  
  /* WRITE PARAMETER INTO OUTPUT FILES */
  fprintf(testACS.fisher_out, "%e\t%e\t%e\t%e", m1, m2, m1 + m2, m1 * m2 / pow(m1 + m2, 2.));
  fprintf(testACS.cov_out, "%e\t%e\t%e\t%e", m1, m2, m1 + m2, m1 * m2 / pow(m1 + m2, 2.));  	
  
  /*********************************************************************
   *
   *  Determine waveform length to adjust for 
   *  optimal sampling frequency
   * 
   ********************************************************************/   
  
  if(testACS.verbose_switch==1){fprintf(stdout, "... Determine optimal sampling frequency \n");}  
  
  // CREATE WAVEFORM TO DETERMINE OPTIMAL SAMPLING FREQUENCY
  CoherentGW tmpwaveform;
  CoherentGW tmp2waveform;
  PPNParamStruc tmpparams;
  
  ///* TEST COMBINING HP AND HC */
  //REAL4TimeSeries *ht_test;
  //ht_test = XLALCreateREAL4TimeSeries("ht", &(testACS.epoch), 0.0, testACS.deltaT_ts, &lalStrainUnit, testACS.N_ts);
  
  ///* CLEAN WAVEFORM */
	//for(k=0;k<ht_test->data->length;k++)
	//{
		//ht_test->data->data[k] = 0.0;
	//}
  
  memset(&tmpwaveform, 0, sizeof(CoherentGW));
  memset(&tmp2waveform, 0, sizeof(CoherentGW));
  tmpparams.ppn = XLALCreateREAL4Vector( params.ppn->length );
  XLALCopyPPNParamStruc(&params, &tmpparams);
  
  fprintf(stdout, "... Ready for first function call\n");
  LALGeneratePPNAmpCorInspiral( &status, &tmpwaveform, &tmpparams);
  fprintf(stdout, "... Past first function call\n");    
  LALGeneratePPNAmpCorInspiral( &status, &tmp2waveform, &params);
  
  for(i=0;i<((INT4)tmpwaveform.f->data->length);i++)
  {
		if(tmpwaveform.h->data->data[2*i] != tmp2waveform.h->data->data[2*i] && tmpwaveform.h->data->data[2*i+1] != tmp2waveform.h->data->data[2*i+1] && tmpwaveform.phi->data->data[i] != tmp2waveform.phi->data->data[i]) printf("copy failed for %d \n", i);
		
	}    
	
  /*********************************************************************
   *
   *  Create PSD
   * 
   ********************************************************************/    
  
  if(testACS.verbose_switch==1){fprintf(stdout, "... Creating PSD \n");}
  
  // CREATE FREQUENCY SERIES FOR PSD COMPUTING
  REAL8FrequencySeries *psd;
  if ( ( psd = (REAL8FrequencySeries *) LALMalloc( sizeof(REAL8FrequencySeries) ) ) == NULL)
  {
      LALFree(psd); psd = NULL;
      if(testACS.verbose_switch==1){fprintf(stdout, "... Creating PSD - error creating Frequency Series \n");}
  }
  else
  {
    psd = XLALCreateREAL8FrequencySeries("PSD_LIGO",  &(testACS.epoch), 0.0, testACS.deltaF_fs, &lalStrainUnit, testACS.N_fs);
  }
  
  // Computing LIGO PSD
  LALNoiseSpectralDensity( &status, psd->data, noisemodel, psd->deltaF );

	INT4 j;
	if(testACS.printall == 1)
	{
		for (j=0;j<testACS.N_fs;j++){fprintf(testACS.psd_out, "%e\t%e\n", j*testACS.deltaF_fs, pow(psd->data->data[j],0.5));}
	}	
	
  /*********************************************************************
   *
   *  Clean Up
   * 
   ********************************************************************/ 
  
  if(testACS.verbose_switch==1){fprintf(stdout, "... Cleaning Up \n");}
  
  /* GSL RNG */
  //if(testACS.verbose_switch==1){fprintf(stdout, "... Cleaning Up - RNG \n");}
  gsl_rng_free (p);
  
  /* TEMPORARY STRUCTURES */
  //if(testACS.verbose_switch==1){fprintf(stdout, "... Cleaning Up - InputParams \n");}
  XLALDestroyREAL4Vector(params.ppn); params.ppn = NULL;
  //if(testACS.verbose_switch==1){fprintf(stdout, "... Cleaning Up - PSD \n");}
  XLALDestroyREAL8FrequencySeries(psd); psd = NULL;

	/* CLOSING OUTPUT FILES */
	//if(testACS.verbose_switch==1){fprintf(stdout, "... Cleaning Up - Output Files \n");}
  if(testACS.printall == 1)
  {
		fclose(psd_out);
		fclose(modsq_out);
		fclose(derivatives_out);
		fclose(fourierderivs_out);
		fclose(hx_out);
	}
  fclose(fisher_out);
  fclose(cov_out);
  //fclose(length_out);	
  
  /*********************************************************************
   *
   * Closing program
   * 
   ********************************************************************/	    

  /*if(testACS.verbose_switch==1)*/fprintf(stdout, "Finished FisherTest \n");
  REPORTSTATUS( &status );
	return status.statusCode;
}
