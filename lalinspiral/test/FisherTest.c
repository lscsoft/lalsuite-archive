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

/* PLEASE NOTE THAT THIS VERSION IS STILL BEING DEVELOPED. DOCUMENTATION
 * AND CLEARER CODE WILL BE UPDATED AS SOON AS POSSIBLE
 * IF YOU HAVE ANY QUESTION, PLEASE EMAIL TGFLI@NIKHEF.NL */ 

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

double ET_B(double f)
{
  /* Calculate detector sensitivity using ET_B */
	
  double x;       // scaled frequency
  double Output;  // output power spectrum
	
  x = f/100.;
  Output = pow(2.39*1.e-27*pow(x,-15.64) + 0.349*pow(x,-2.145) + 
							 1.76*pow(x,-0.12) + 0.409*pow(x,1.10),2.);
  
  /* Sh *= pow(10.,-50.); */ // ? scaling factor?
  return( Output );  
	
}

double ET_C(double f)
{
  /* Calculate PSD for ET_C (times 10^50)*/
	
  double Output;
  double fb[23] = {10000.,4600.,1900.,210.,110.,80.,51.,38.,29.8,22.,18.5,9.,7.5,5.7,4.,3.,2.4,2.,1.65,1.4,1.26,1.12,1.};
	
  if(f>=fb[1])
	{ Output = 6.405e-3*pow(f,0.9969)-5.092e-2; }
  else if(fb[1]>f && f>=fb[2])
	{ Output = 7.458e-3*pow(f,0.9794)-1.793e-1; }
  else if(fb[2]>f && f>=fb[3])
	{ Output = -7.646e-30*pow(f,9)+1.54e-25*pow(f,8)-1.331e-21*pow(f,7)+6.458e-18*pow(f,6)-1.936e-14*pow(f,5)+3.719e-11*pow(f,4)-4.606e-8*pow(f,3)+3.614e-5*pow(f,2)-1.095e-2*f+4.397e-0; }
  else if(fb[3]>f && f>=fb[4])
	{ Output = -3.317e-7*pow(f,3)+2.027e-4*pow(f,2)-4.354e-2*f+6.601e-0; }
  else if(fb[4]>f && f>=fb[5])
	{ Output = -2.935e-21*pow(f,9)+6.757e-18*pow(f,8)-6.794e-15*pow(f,7)+3.916e-12*pow(f,6)-1.426e-9*pow(f,5)+3.414e-7*pow(f,4)-5.399e-5*pow(f,3)+5.52e-3*pow(f,2)-3.433e-1*f+1.389e1; }
  else if(fb[5]>f && f>=fb[6])
	{ Output = -2.998e-5*pow(f,3)+7.391e-3*pow(f,2)-6.303e-1*f +2.279e1; }
  else if(fb[6]>f && f>=fb[7])
	{ Output = 1.27e-4*pow(f,3)-1.118e-2*pow(f,2)+2.716e-2*f +1.674e1; }
  else if(fb[7]>f && f>=fb[8])
	{ Output = 3.151e-3*pow(f,3)-3.676e-1*pow(f,2)+1.4e1*f -1.655e2; }
  else if(fb[8]>f && f>=fb[9])
	{ Output = -3.173e-3*pow(f,3)+3.028e-1*pow(f,2)-9.138e0*f + 9.6e1; }
  else if(fb[9]>f && f>=fb[10])
	{ Output = -5.873e-4*pow(f,4)+5.781e-2*pow(f,3)-2.06e-0*pow(f,2)+3.134e1*f-1.627e2; }
  else if(fb[10]>f && f>=fb[11])
	{ Output = 5.104e-3*pow(f,3)-2.817e-1*pow(f,2)+4.761e-0*f-1.466e1; }
  else if(fb[11]>f && f>=fb[12])
	{ Output = -1.368e-1*pow(f,3)+3.581e-0*pow(f,2)-3.031e1*f+9.151e1; }
  else if(fb[12]>f && f>=fb[13])
	{ Output = -6.333e-1*pow(f,3)+1.509e1*pow(f,2)-1.193e2*f+3.210e2; }
  else if(fb[13]>f && f>=fb[14])
	{ Output = 1.362e4*pow(f,(-3.952))-3.004e-3; }
  else if(fb[14]>f && f>=fb[15])
	{ Output = 1.135e4*pow(f,(-3.806))-1.16e-0; }
  else if(fb[15]>f && f>=fb[16])
	{ Output = 8.52e3*pow(f,(-3.53))-4.264e-0; }
  else if(fb[16]>f && f>=fb[17])
	{ Output = 6.899e3*pow(f,(-3.276))-9.181e-0; }
  else if(fb[17]>f && f>=fb[18])
	{ Output = 7.056e3*pow(f,(-3.521))+9.042e1; }
  else if(fb[18]>f && f>=fb[19])
	{ Output = 2.717e4*pow(f,4)-1.842e5*pow(f,3)+4.709e5*pow(f,2)-5.402e5*f+2.3667e5; }
  else if(fb[19]>f && f>=fb[20])
	{ Output = -1.98e5*pow(f,3)+8.486e5*pow(f,2)-1.219e6*f+5.8894e5; }
  else if(fb[20]>f && f>=fb[21])
	{ Output = -1.042e6*pow(f,3)+4.045e6*pow(f,2)-5.257e6*f+2.29051e6; }
  else if(fb[21]>f)
	{ Output = 4.741e4*pow(f,(-11.72))+1.906e2; }
  else
	{ printf("Frequency %f out of range.\n",f);
	  Output = 42; }
  
  Output = pow(Output,2);
	
  return(Output);
}

int main( INT4 argc, char **argv )
{
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
  fprintf(stdout, "Running FisherTest - seed %d \n", atoi(argv[1])); 
  
  /* INITIATE GSL RANDOM NUMBER GENERATOR */
  
	const gsl_rng_type *type;           // RNG type
  gsl_rng *p;                         // Generator
  gsl_rng_env_setup();                // Setup environment
  gsl_rng_default_seed = seed;		    // vary generation sequence
  type = gsl_rng_default;             // set RNG type to default
  p = gsl_rng_alloc (type);           // Set RNG type
	
  /*********************************************************************
   *
   *  Input parameters, could also be read in from files
   * 
   ********************************************************************/
  //
  // WAVEFORM PARAMETERS
  REAL4 m1            = 9.0;           // Component Mass 1 
  REAL4 m2            = 90.0; 					// Component Mass 2
  REAL4 dist          = 3.0e8*LAL_PC_SI;// Distance to source
  REAL4 latitude			= LAL_PI/6.0; //acos(gsl_ran_flat (p, -1.0, 1.0));						// latitude
	REAL4 longitude			= LAL_PI/6.0; //gsl_ran_flat (p, 0, LAL_TWOPI);						// longitude
  REAL4 inc           = LAL_PI/3.0; //acos(gsl_ran_flat (p, -1.0, 1.0));     // Inclination angle
  REAL4 psi						= LAL_PI/4.0; //gsl_ran_flat (p, 0, LAL_TWOPI);						// Polarisation angle
  REAL4 tc            = 0.0;            // End time of waveform
  REAL4 phic          = 0.0;            // Phase of coalescenes
  REAL4 fStartIn      = 20.0;           // Starting frequency
  REAL8 fStopIn       = 0.0;            // No early termination
	INT4 phaseOrder     = 6;              // Expansion order PPN
  INT4 ampOrder       = 5;              // Amplitude Order  
  INT4 gpsSeconds			= 0;							// GPS time, seconds
  INT4 gpsNanoSeconds = 0;							// GPS time, nanoseconds
  
  /*********************************************************************
   *
   *  Read parameters from file
   * 
   ********************************************************************/  
  
  //INT4 q = 0;
  //INT4 scanitems = 1;
  INT4 testPhaseParam = 2;
  INT4 psd_switch = 0;
	
  
  //FILE *inputParams;
  //inputParams = fopen("/home/tgfli/lscsoft/ampCorPPN_filter/src/lalsuite/lal/packages/bank/src/ET_params.dat", "r"); 
  
  //while(scanitems > 0)
  //{
	//scanitems=fscanf(inputParams, "%d\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%d\t%d\t%e\n", &gpsSeconds, &m1, &m2, &dist, &longitude, &latitude, &inc, &phic, &psi, &testPhaseParam, &psd_switch, &fStartIn);
	//if(scanitems <=0) {printf("Nothing to read \n"); exit(0);}
	//if(q==(seed))
	//{
	//dist *= 1.0E6*LAL_PC_SI; 
	//break;
	//}
	//q++;
	//}
  //fclose(inputParams);
  
	//if(testACS.verbose_switch==1)printf("input: m1 = %1.6e | m2 = %1.6e | dist = %1.6e | theta = %1.6e | phi = %1.6e | iota = %1.6e | psi = %1.6e \n", m1, m2, dist, latitude, longitude, inc, psi);  
  
  /*********************************************************************
   *
   *  Algorithm Control Structure
   * 
   ********************************************************************/      
	
	/* ALGORITHM CONTROL STRUCTURE */
	FisherACS testACS;
	
	/* GENERAL CONTROL PARAMETERS */
	testACS.seed							= seed;
	testACS.verbose_switch		=	1;
	testACS.printall					= 1;
	testACS.SavitskyGolay_switch = 1;
	testACS.N_inputparams			=	6;
	testACS.dynRange					= 1E+27;
	testACS.coordinate				= 0; // 0: detector frame, 1 equitorial frame
	
	/* TIME/FREQUENCY SERIES (TO BE ADJUSTED BELOW) */
	testACS.N_ts							= 1E+5;
	testACS.N_fs							= testACS.N_ts/2+1;
	testACS.deltaT_ts					= 1.0/16384.0;
	testACS.deltaF_fs					= 1.0/((REAL8)testACS.N_ts * testACS.deltaT_ts);
	
	/* DIFFERENTIATION PARAMETERS */
	testACS.xsamples_tot						= 400;
	testACS.linReg_countMax				= testACS.xsamples_tot/2;
	testACS.linReg_countMin				= 20;
	testACS.linReg_loopMax				= 100;
	testACS.deltaPhi_max					= LAL_PI/25;
	testACS.reldelta							= 5E-6;
	testACS.absdelta_adjustUp			= 1.5;
	testACS.absdelta_adjustDown		= 0.5;
	testACS.xsamples_adjustUp			= 2;
	testACS.xsamples_adjustDown		= 5;
	
	/* DETECTOR FREQUENCY WINDOW */
	REAL4 adjustLSO	= (ampOrder+2.0)/2.0;
	
	testACS.fstart								= fStartIn;
	testACS.fstop									= adjustLSO*(pow(LAL_C_SI, 3.0))/(pow(6.0,1.5)*LAL_PI*LAL_G_SI*(m1+m2)*LAL_MSUN_SI);	
	
	/* SAVITSKY GOLAY FILTER */
	testACS.SG_R									= 7; // Rth order polynomial
	testACS.SG_M									= 4; // M SURROUNDING points (2M+1 points)
	
	testACS.testPhaseParam				= testPhaseParam;
	testACS.mass2phi_switch				= 0;
	
  if(testACS.verbose_switch==1)printf("... LSO = %e\n",testACS.fstop);	
	
  // PSD MODEL
  void (*noisemodel)(LALStatus*,REAL8*,REAL8) = LALAdvLIGOPsd;
  
  /*********************************************************************
   *
   *  Declare and set temporary variables
   * 
   ********************************************************************/
  
  if(testACS.verbose_switch==1){fprintf(stdout, "... Declaring variables \n");}
  
  // COUNTERS
  INT4 i,j;
  
  // LOADING PARAMETER INTO STRUCTURE
  PPNConsistencyParamStruc params;
  params.position.longitude = longitude; /*acos(gsl_ran_flat (p, -1.0, 1.0));*/
  params.position.latitude 	= latitude;/*gsl_ran_flat (p, 0, LAL_TWOPI);*/
  params.position.system 		= COORDINATESYSTEM_GEOGRAPHIC;
  params.epoch.gpsSeconds		= gpsSeconds;
  params.epoch.gpsNanoSeconds = gpsNanoSeconds;
  params.psi								= psi; /*acos(gsl_ran_flat(p,-1.0,1.0));*/
  params.deltaT 						= testACS.deltaT_ts;
  params.mTot_real8 				= m1 + m2;
  params.eta_real8  				= m1 * m2 / pow(m1 + m2, 2.);
  params.inc 								= inc;
  params.cosI 							= cos(inc);
  params.sinI 							= sin(inc);
  params.tc  								= tc;
  params.phi 								= phic;
  params.d   								= dist;
  params.fStartIn 					= fStartIn;
  params.fStopIn  					= fStopIn;
  params.ampOrder 					= ampOrder;
  
	params.phasePNparams[0] = -pow(params.eta_real8,-3.0/8.0)*pow(5.0*LAL_MTSUN_SI*params.mTot_real8,-5.0/8.0); 
	params.phasePNparams[1] = -(3715.0/8064.0 + 55.0/96.0*params.eta_real8)*pow(params.eta_real8,-5.0/8.0)*pow(5.0*LAL_MTSUN_SI*params.mTot_real8,-3.0/8.0);
	params.phasePNparams[2] = 3.0/4.0*LAL_PI*pow(params.eta_real8,-0.75)*pow(5.0*LAL_MTSUN_SI*params.mTot_real8,-0.25); 
	params.phasePNparams[3] = -(9275495.0/14450688.0 + 284875.0/258048.0*params.eta_real8 + 1855.0/2048.0*pow(params.eta_real8,2.0))*pow(params.eta_real8,-7.0/8.0)*pow(5.0*LAL_MTSUN_SI*params.mTot_real8,-1.0/8.0);
	params.phasePNparams[4] = -1.0/params.eta_real8*(-38645.0/172032.0 + 65.0/2048.0*params.eta_real8)*LAL_PI*log(params.eta_real8/(5.0*LAL_MTSUN_SI*params.mTot_real8));
	params.phasePNparams[5] = -1.0/params.eta_real8*(-38645.0/172032.0 + 65.0/2048.0*params.eta_real8)*LAL_PI; 
	params.phasePNparams[6] = -(831032450749357.0/57682522275840.0 - 53.0/40.0*LAL_PI*LAL_PI - 107.0/56.0*LAL_GAMMA + 107.0/448.0*log(params.eta_real8/(256*5.0*LAL_MTSUN_SI*params.mTot_real8)) + (-123292747421.0/4161798144.0 + 2255.0/2048.0*LAL_PI*LAL_PI + 385.0/48.0*(-1987.0/3080.0) - 55.0/16.0*(-11831.0/9240.0))*params.eta_real8 + 154565.0/1835008.0*pow(params.eta_real8,2.0) - 1179625.0/1769472.0*pow(params.eta_real8,3.0))*pow(params.eta_real8,-9.0/8.0)*pow(5.0*LAL_MTSUN_SI*params.mTot_real8,1.0/8.0);
	params.phasePNparams[7] = -107.0/448.0*pow(params.eta_real8,-9.0/8.0)*pow(5.0*LAL_MTSUN_SI*params.mTot_real8,1.0/8.0);
	params.phasePNparams[8] = -(188516689.0/173408256.0 + 488825.0/516096.0*params.eta_real8 - 141769.0/516096.0*pow(params.eta_real8,2.0))*LAL_PI*pow(params.eta_real8,-5.0/4.0)*pow(5.0*LAL_MTSUN_SI*params.mTot_real8,1.0/4.0);
	/*
	 params.phi0 = -pow(params.eta_real8,-3.0/8.0)*pow(5.0*LAL_MTSUN_SI*params.mTot_real8,-5.0/8.0); 
	 params.phi2 = -(3715.0/8064.0 + 55.0/96.0*params.eta_real8)*pow(params.eta_real8,-5.0/8.0)*pow(5.0*LAL_MTSUN_SI*params.mTot_real8,-3.0/8.0);
	 params.phi3 = 3.0/4.0*LAL_PI*pow(params.eta_real8,-0.75)*pow(5.0*LAL_MTSUN_SI*params.mTot_real8,-0.25); 
	 params.phi4 = -(9275495.0/14450688.0 + 284875.0/258048.0*params.eta_real8 + 1855.0/2048.0*pow(params.eta_real8,2.0))*pow(params.eta_real8,-7.0/8.0)*pow(5.0*LAL_MTSUN_SI*params.mTot_real8,-1.0/8.0);
	 params.phi5 = -1.0/params.eta_real8*(-38645.0/172032.0 + 65.0/2048.0*params.eta_real8)*LAL_PI*log(params.eta_real8/(5.0*LAL_MTSUN_SI*params.mTot_real8));
	 params.phi5l = -1.0/params.eta_real8*(-38645.0/172032.0 + 65.0/2048.0*params.eta_real8)*LAL_PI; 
	 params.phi6 = -(831032450749357.0/57682522275840.0 - 53.0/40.0*LAL_PI*LAL_PI - 107.0/56.0*LAL_GAMMA + 107.0/448.0*log(params.eta_real8/(256*5.0*LAL_MTSUN_SI*params.mTot_real8)) + (-123292747421.0/4161798144.0 + 2255.0/2048.0*LAL_PI*LAL_PI + 385.0/48.0*(-1987.0/3080.0) - 55.0/16.0*(-11831.0/9240.0))*params.eta_real8 + 154565.0/1835008.0*pow(params.eta_real8,2.0) - 1179625.0/1769472.0*pow(params.eta_real8,3.0))*pow(params.eta_real8,-9.0/8.0)*pow(5.0*LAL_MTSUN_SI*params.mTot_real8,1.0/8.0);
	 params.phi6l = -107.0/448.0*pow(params.eta_real8,-9.0/8.0)*pow(5.0*LAL_MTSUN_SI*params.mTot_real8,1.0/8.0);
	 params.phi7 = -(188516689.0/173408256.0 + 488825.0/516096.0*params.eta_real8 - 141769.0/516096.0*pow(params.eta_real8,2.0))*LAL_PI*pow(params.eta_real8,-5.0/4.0)*pow(5.0*LAL_MTSUN_SI*params.mTot_real8,1.0/4.0);
	 */
  params.ppn = XLALCreateREAL4Vector( phaseOrder+1 );
  for (i=0; i<=(phaseOrder); i++)
  {
    if(phaseOrder > 0 && i==1) {params.ppn->data[i] = 0.0;}
    else {params.ppn->data[i] = 1.0;}
  }
  
  if(testACS.verbose_switch==1){fprintf(stdout, "... m1 = %e | m2 = %e | mtot = %e | eta = %e | phiOrder = %d | ampOrder = %d \n", m1, m2, params.mTot_real8, params.eta_real8, phaseOrder, ampOrder);}
  
  /*********************************************************************
   *
   *  Determine waveform length to adjust for 
   *  optimal sampling frequency
   * 
   ********************************************************************/    
  
  if(testACS.verbose_switch==1){fprintf(stdout, "... Determine optimal sampling frequency \n");}
  
  /* COUNTER */
  INT4 k;
  
  /* TOTAL LENGTH WAVEFORM (CONSTANT) */
  REAL4 deltaT_total = 0.0;
  
  // CREATE WAVEFORM TO DETERMINE OPTIMAL SAMPLING FREQUENCY
  CoherentGW tmpwaveform;
  PPNConsistencyParamStruc tmpparams;
  
  /* TEST COMBINING HP AND HC */
  REAL4TimeSeries *ht_test;
  ht_test = XLALCreateREAL4TimeSeries("ht", &(testACS.epoch), 0.0, testACS.deltaT_ts, &lalStrainUnit, testACS.N_ts);
  
  /* CLEAN WAVEFORM */
	for(k=0;k<((INT4)ht_test->data->length);k++)
	{
		ht_test->data->data[k] = 0.0;
	}
  
  memset(&tmpwaveform, 0, sizeof(CoherentGW));
  tmpparams.ppn = XLALCreateREAL4Vector( params.ppn->length );
  XLALCopyPPNConsistencyParamStruc(&params, &tmpparams);
  
  LALGeneratePPNAmpCorConsistency( &status, &tmpwaveform, &params);
  LALInspiralCombinePlusCross( &status, tmpwaveform, params, ht_test, &testACS);
	
  if(testACS.verbose_switch==1){fprintf(stdout, "... Determine optimal sampling frequency - wavelength %d, deltaT = %e \n", tmpwaveform.f->data->length, testACS.deltaT_ts);}
  deltaT_total = tmpwaveform.f->data->length*testACS.deltaT_ts;
  
  /* SET MINIMUM LENGTH OF WAVEFORM TO ACCOMODATE LENGTH SHIFTS */
  if( (2*(deltaT_total*testACS.fstop-1))>1E3 ) testACS.N_ts = 2*(deltaT_total*testACS.fstop-1);
  else testACS.N_ts = 1E3;
  
  testACS.deltaT_ts = deltaT_total/((REAL4)testACS.N_ts);
  
  /* ADD 5% EXTRA SPACE FOR SHIFTS */
  testACS.N_ts *= 1.05;
  
  /* SETTING RELATED PARAMETERS */
  testACS.N_fs = testACS.N_ts/2+1;
  testACS.deltaF_fs = 1.0/((REAL8)testACS.N_ts * testACS.deltaT_ts);
  params.deltaT = testACS.deltaT_ts;
  
  XLALDestroyREAL4Vector(tmpparams.ppn);
	XLALDestroyCoherentGW(&tmpwaveform); 
	
  if(testACS.verbose_switch==1){fprintf(stdout, "... Determine optimal sampling frequency - UPDATED: segmentlength %d, deltaT = %e \n", testACS.N_ts, testACS.deltaT_ts);}
  
	
  /*********************************************************************
   *
   *  Output Files
   * 
   ********************************************************************/  
  
  if(testACS.verbose_switch==1){fprintf(stdout, "... Creating Output Files \n");}
  
  /* FOLDER FOR OUTPUT FILES */
	char folder[128] = "/home/tjonnie/Data/Fisher/test/";
	testACS.folder = folder;
  
  /* NAME OUTPUT FILES */
  const char fisher[]						= "fisher";
  const char covariance[]				= "covariance";
  
  /* CREATE OUTPUT FILES */
  FILE *fisher_out;
  FILE *cov_out;  
	
	char fisher_name[4096];
  char cov_name[4096];
  
	sprintf(fisher_name, "%s%s%d%s", folder, fisher, seed, ".dat");
	sprintf(cov_name, "%s%s%d%s", folder, covariance, seed, ".dat");
	
	/* GENERAL OUTPUT */
  fisher_out				= fopen(fisher_name, "w");
  cov_out						= fopen(cov_name, "w");
	
  
  /*********************************************************************
   *
   *  Create PSD
   * 
   ********************************************************************/    
  
  if(testACS.verbose_switch==1){fprintf(stdout, "... Creating PSD \n");}
  
  /* PSD */
  REAL8FrequencySeries *psd;
  
  // CREATE FREQUENCY SERIES FOR PSD COMPUTING
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
  if(psd_switch==0)LALNoiseSpectralDensity( &status, psd->data, noisemodel, psd->deltaF );
  
  /* ET SENSITIVITY CURVE */
	if(psd_switch==1) for(j=0;j<testACS.N_fs;j++){psd->data->data[j] = ET_B(j*testACS.deltaF_fs)*1E-50;}
	if(psd_switch==2) for(j=0;j<testACS.N_fs;j++){psd->data->data[j] = ET_C(j*testACS.deltaF_fs)*1E-50;}
	
  
  /*********************************************************************
   *
   *  Print input parameters into file
   * 
   ********************************************************************/ 
	
	if(testACS.verbose_switch==1)fprintf(stdout,"...Printing input to file \n");
	
	FILE *print_input;
  char print_name[128];
	const char print_base[]		 = "input_params_";
	sprintf(print_name,"%s%s%d%s", folder, print_base, testACS.seed, ".dat");
	print_input = fopen(print_name, "w");
	
	/* CREATE TIME STAMP */
	time_t     now;
	struct tm  *ts;
	
	now = time(NULL);
	
	ts = localtime(&now);
	
	fprintf(print_input, "Input Parameters - %d \n", testACS.seed);
	fprintf(print_input, "date: %d-%d-%d \n", ts->tm_mday, ts->tm_mon, ts->tm_year+1900);
	fprintf(print_input, "time: %d:%d:%d \n", ts->tm_hour, ts->tm_min, ts->tm_sec);
	fprintf(print_input, "\n");
	fprintf(print_input, "PPN parameters \n");
	fprintf(print_input, "\n");
	fprintf(print_input, "latitude\t = %e \n", params.position.latitude);
	fprintf(print_input, "longitude\t = %e \n", params.position.longitude);
	fprintf(print_input, "polarization\t = %e \n", params.psi);
	fprintf(print_input, "epoch\t\t = %d.%d \n", params.epoch.gpsSeconds, params.epoch.gpsNanoSeconds);
	fprintf(print_input, "m1\t\t = %e \n", m1);
	fprintf(print_input, "m2\t\t = %e \n", m2);
	fprintf(print_input, "total mass\t = %e \n", params.mTot_real8);
	fprintf(print_input, "symmetric mass\t = %e \n", params.eta_real8);
	fprintf(print_input, "distance\t = %e \n", params.d);
	fprintf(print_input, "inclination\t = %e \n", params.inc);
	fprintf(print_input, "cosi\t\t = %e \n", params.cosI);
	fprintf(print_input, "sini\t\t = %e \n", params.sinI);
	fprintf(print_input, "coal. phase\t = %e \n", params.phi);
	fprintf(print_input, "deltaT\t\t = %e \n", params.deltaT);
	fprintf(print_input, "fStart\t\t = %e \n", params.fStartIn);
	fprintf(print_input, "fStop\t\t = %e \n", params.fStopIn);
	//fprintf(print_input, "maxLength \t\t = \n");
	fprintf(print_input, "PPN\t\t = %d \n", phaseOrder);
	fprintf(print_input, "APN\t\t = %d \n", params.ampOrder);
	fprintf(print_input, "phi0\t\t = %e \n", params.phasePNparams[0]);
	fprintf(print_input, "phi2\t\t = %e \n", params.phasePNparams[1]);
	fprintf(print_input, "phi3\t\t = %e \n", params.phasePNparams[2]);
	fprintf(print_input, "phi4\t\t = %e \n", params.phasePNparams[3]);
	fprintf(print_input, "phi5\t\t = %e \n", params.phasePNparams[4]);
	fprintf(print_input, "phi5l\t\t = %e \n", params.phasePNparams[5]);
	fprintf(print_input, "phi6\t\t = %e \n", params.phasePNparams[6]);
	fprintf(print_input, "phi6l\t\t = %e \n", params.phasePNparams[7]);
	fprintf(print_input, "phi7\t\t = %e \n", params.phasePNparams[8]);
	fprintf(print_input, "\n");
	//	fprintf(print_input, "\n");
	//	fprintf(print_input, "Output parameters \n");
	//	fprintf(print_input, "\n");
	//	fprintf(print_input, "tc \t\t = \n");
	//	fprintf(print_input, "dfdt \t\t = \n");
	//	fprintf(print_input, "fStartActual \t\t = \n");
	//	fprintf(print_input, "fStopActual \t\t =\n");
	//	fprintf(print_input, "waveform length \t\t = \n");
	//	fprintf(print_input, "termination \t\t = \n");
	
	fprintf(print_input, "Algorithm Control Structure \n");
	fprintf(print_input, "\n");
	fprintf(print_input, "verbose switch\t = %d\n", testACS.verbose_switch);
	fprintf(print_input, "print switch\t = %d\n", testACS.printall);
	fprintf(print_input, "coordinates\t = %d \n", testACS.coordinate);
	fprintf(print_input, "SG switch\t = %d \n", testACS.SavitskyGolay_switch);
	fprintf(print_input, "folder\t\t = %s \n", testACS.folder);
	fprintf(print_input, "N parameters\t = %d \n", testACS.N_inputparams);
	fprintf(print_input, "N time series\t = %d \n", testACS.N_ts);
	fprintf(print_input, "N freq series\t = %d \n", testACS.N_fs);
	fprintf(print_input, "deltaT\t\t = %e \n", testACS.deltaT_ts);
	fprintf(print_input, "deltaF\t\t = %e \n", testACS.deltaF_fs);
	fprintf(print_input, "xsamples\t = %d \n", testACS.xsamples_tot);
	fprintf(print_input, "linReg_countMax\t= %d \n", testACS.linReg_countMax);
	fprintf(print_input, "linReg_countMin\t= %d \n", testACS.linReg_countMin);
	fprintf(print_input, "linReg_loopMax\t= %d \n", testACS.linReg_loopMax);
	fprintf(print_input, "deltaPhi_max\t= %e \n", testACS.deltaPhi_max);
	fprintf(print_input, "reldelta\t= %e \n",testACS.reldelta);
	fprintf(print_input, "absdelta_Up\t= %e \n", testACS.absdelta_adjustUp);
	fprintf(print_input, "absdelta_Down\t= %e \n", testACS.absdelta_adjustDown);
	fprintf(print_input, "xsamples_Up\t= %d \n", testACS.xsamples_adjustUp);
	fprintf(print_input, "xsamples_Down\t= %d \n", testACS.xsamples_adjustDown);
	fprintf(print_input, "test param\t= %d \n", testACS.testPhaseParam);
	fprintf(print_input, "SG R\t\t= %d \n", testACS.SG_R);
	fprintf(print_input, "SG M\t\t= %d \n", testACS.SG_M);
	fprintf(print_input, "psd switch\t= %d \n", psd_switch);	
	
	fclose(print_input);	 
	
	
  /*********************************************************************
   *
   *  Test Savitzky-Golay Filter
   * 
   ********************************************************************/   
	
//  CoherentGW SGTestWaveform;
//  CoherentGW SGTestWaveform2;
//  PPNConsistencyParamStruc SGTestParams;
//  
//  XLALCopyPPNConsistencyParamStruc(&params, &SGTestParams);
//  
//  REAL4TimeSeries *SGTest_ht;
//  REAL4TimeSeries *SGTest_ht2;
//	REAL4TimeSeries *phiStart;
//	REAL4TimeSeries *phiEnd;
//  REAL4TimeSeries *SGTest_ht_SG;
//  REAL4TimeSeries *SGTest_dhdt_SG;
//  REAL4TimeSeries *SGTest_dhdt_IB;
//  
//  SGTest_ht = XLALCreateREAL4TimeSeries("ht", &(testACS.epoch), 0.0, testACS.deltaT_ts, &lalStrainUnit, testACS.N_ts);
//  SGTest_ht2 = XLALCreateREAL4TimeSeries("ht", &(testACS.epoch), 0.0, testACS.deltaT_ts, &lalStrainUnit, testACS.N_ts);
//	phiStart = XLALCreateREAL4TimeSeries("phiStart", &(testACS.epoch), 0.0, testACS.deltaT_ts, &lalStrainUnit, testACS.N_ts);
//  phiEnd = XLALCreateREAL4TimeSeries("phiEnd", &(testACS.epoch), 0.0, testACS.deltaT_ts, &lalStrainUnit, testACS.N_ts);
//  SGTest_ht_SG = XLALCreateREAL4TimeSeries("ht SG", &(testACS.epoch), 0.0, testACS.deltaT_ts, &lalStrainUnit, testACS.N_ts);
//  SGTest_dhdt_SG = XLALCreateREAL4TimeSeries("dhdt SG", &(testACS.epoch), 0.0, testACS.deltaT_ts, &lalStrainUnit, testACS.N_ts);
//  SGTest_dhdt_IB = XLALCreateREAL4TimeSeries("dhdt IB", &(testACS.epoch), 0.0, testACS.deltaT_ts, &lalStrainUnit, testACS.N_ts);
//  
//  for(i=0;i<testACS.N_ts;i++)
//  {
//		SGTest_ht->data->data[i] = 0.0;
//		SGTest_ht2->data->data[i] = 0.0;
//		phiStart->data->data[i] = 0.0;
//		phiEnd->data->data[i] = 0.0;
//		SGTest_dhdt_SG->data->data[i] = 0.0;
//		SGTest_dhdt_IB->data->data[i] = 0.0;
//		
//	}
//  
//  LALInspiralComputeDerivatives(&status, SGTest_dhdt_IB, &(SGTestParams), 0, &(testACS));	
//  
//  LALGeneratePPNAmpCorConsistency( &status, &(SGTestWaveform), &(SGTestParams));
//  SGTestParams.mTot_real8 *= 1.01;
//  SGTestParams.phasePNparams[0] = -pow(SGTestParams.eta_real8,-3.0/8.0)*pow(5.0*LAL_MTSUN_SI*SGTestParams.mTot_real8,-5.0/8.0); 
//	SGTestParams.phasePNparams[1] = -(3715.0/8064.0 + 55.0/96.0*SGTestParams.eta_real8)*pow(SGTestParams.eta_real8,-5.0/8.0)*pow(5.0*LAL_MTSUN_SI*SGTestParams.mTot_real8,-3.0/8.0);
//	SGTestParams.phasePNparams[2] = 3.0/4.0*LAL_PI*pow(SGTestParams.eta_real8,-0.75)*pow(5.0*LAL_MTSUN_SI*SGTestParams.mTot_real8,-0.25); 
//	SGTestParams.phasePNparams[3] = -(9275495.0/14450688.0 + 284875.0/258048.0*SGTestParams.eta_real8 + 1855.0/2048.0*pow(SGTestParams.eta_real8,2.0))*pow(SGTestParams.eta_real8,-7.0/8.0)*pow(5.0*LAL_MTSUN_SI*SGTestParams.mTot_real8,-1.0/8.0);
//	SGTestParams.phasePNparams[4] = -1.0/SGTestParams.eta_real8*(-38645.0/172032.0 + 65.0/2048.0*SGTestParams.eta_real8)*LAL_PI*log(SGTestParams.eta_real8/(5.0*LAL_MTSUN_SI*SGTestParams.mTot_real8));
//	SGTestParams.phasePNparams[5] = -1.0/SGTestParams.eta_real8*(-38645.0/172032.0 + 65.0/2048.0*SGTestParams.eta_real8)*LAL_PI; 
//	SGTestParams.phasePNparams[6] = -(831032450749357.0/57682522275840.0 - 53.0/40.0*LAL_PI*LAL_PI - 107.0/56.0*LAL_GAMMA + 107.0/448.0*log(SGTestParams.eta_real8/(256*5.0*LAL_MTSUN_SI*SGTestParams.mTot_real8)) + (-123292747421.0/4161798144.0 + 2255.0/2048.0*LAL_PI*LAL_PI + 385.0/48.0*(-1987.0/3080.0) - 55.0/16.0*(-11831.0/9240.0))*SGTestParams.eta_real8 + 154565.0/1835008.0*pow(SGTestParams.eta_real8,2.0) - 1179625.0/1769472.0*pow(SGTestParams.eta_real8,3.0))*pow(SGTestParams.eta_real8,-9.0/8.0)*pow(5.0*LAL_MTSUN_SI*SGTestParams.mTot_real8,1.0/8.0);
//	SGTestParams.phasePNparams[7] = -107.0/448.0*pow(SGTestParams.eta_real8,-9.0/8.0)*pow(5.0*LAL_MTSUN_SI*SGTestParams.mTot_real8,1.0/8.0);
//	SGTestParams.phasePNparams[8] = -(188516689.0/173408256.0 + 488825.0/516096.0*SGTestParams.eta_real8 - 141769.0/516096.0*pow(SGTestParams.eta_real8,2.0))*LAL_PI*pow(SGTestParams.eta_real8,-5.0/4.0)*pow(5.0*LAL_MTSUN_SI*SGTestParams.mTot_real8,1.0/4.0);
//  LALGeneratePPNAmpCorConsistency( &status, &(SGTestWaveform2), &(SGTestParams));
//  LALInspiralCombinePlusCross(&status, SGTestWaveform, SGTestParams, SGTest_ht, &(testACS));
//  LALInspiralCombinePlusCross(&status, SGTestWaveform2, SGTestParams, SGTest_ht2, &(testACS));
//  
//  LALSavitskyGolayFilter(&status, testACS.SG_R, testACS.SG_M, SGTest_ht, SGTest_ht_SG, 0);
//  LALSavitskyGolayFilter(&status, 8, 4, SGTest_ht, SGTest_dhdt_SG, 1);
//	
//	for (i=0; i<testACS.N_ts; i++) {
//		if(i<((INT4)SGTestWaveform.f->data->length)) phiStart->data->data[testACS.N_ts-1-i] = SGTestWaveform.phi->data->data[SGTestWaveform.f->data->length-1-i];
//		if(i<((INT4)SGTestWaveform2.f->data->length)) phiEnd->data->data[testACS.N_ts-1-i] = SGTestWaveform2.phi->data->data[SGTestWaveform2.f->data->length-1-i];
//	}
//  
//  /* PRINT OUTPUT TO FILE */
//  FILE *SGTest_out;
//  FILE *phi_out;
//  SGTest_out = fopen("/home/tjonnie/Programs/lscsoft/ampCorPPN_filter/src/lalsuite/lalinspiral/test/SGTest.dat", "w");
//  phi_out = fopen("/home/tjonnie/Programs/lscsoft/ampCorPPN_filter/src/lalsuite/lalinspiral/test/phi_all.dat", "w");
//  
//  for(i=0;i<testACS.N_ts;i++)
//  {
//		fprintf(SGTest_out, "%e\t%e\t%e\t%e\t%e\n", i*testACS.deltaT_ts, SGTest_ht->data->data[i],  SGTest_ht_SG->data->data[i], SGTest_dhdt_SG->data->data[i], SGTest_dhdt_IB->data->data[i]);		
//		fprintf(phi_out, "%e\t%e\t%e\t%e\n", i*testACS.deltaT_ts, phiStart->data->data[i], phiEnd->data->data[i], phiEnd->data->data[i]-phiStart->data->data[i]);
//	}
//  
//  fclose(SGTest_out);
//  fclose(phi_out);
  
	
  /*********************************************************************
   *
   *  Test Computation of time derivative
   * 
   ********************************************************************/ 		
	
	//REAL4TimeSeries *hderiv_test;
	//hderiv_test = XLALCreateREAL4TimeSeries("hderiv_test", &(testACS.epoch), 0.0, testACS.deltaT_ts, &lalStrainUnit, testACS.N_ts);
	
	//for(j=0;j<testACS.N_ts;j++)
	//{
		//hderiv_test->data->data[j] =0.0;
	//}
	
	//LALInspiralComputeDerivatives(&status, hderiv_test, &params, 3, &testACS);
	
	//FILE *hderiv_out;
	//hderiv_out = fopen("/home/tjonnie/Programs/lscsoft/ampCorPPN_filter/src/lalsuite/lalinspiral/test/testderiv.dat","w");
	
	//for(j=0;j<testACS.N_ts;j++)
	//{
		//fprintf(hderiv_out, "%e\t%e\n", j*testACS.deltaT_ts, hderiv_test->data->data[j]);
	//}
	
	//fclose(hderiv_out);
	
	//XLALDestroyREAL4TimeSeries(hderiv_test);
	
  /*********************************************************************
   *
   *  Compute Fisher Matrix
   * 
   ********************************************************************/ 	
	
	PPNConsistencyParamStruc FisherParams1;
	PPNConsistencyParamStruc FisherParams2;
	
	FisherParams1.ppn = XLALCreateREAL4Vector( params.ppn->length );
	FisherParams2.ppn = XLALCreateREAL4Vector( params.ppn->length );
  XLALCopyPPNConsistencyParamStruc(&params, &FisherParams1);
  XLALCopyPPNConsistencyParamStruc(&params, &FisherParams2);
  
  FisherParams2.position.longitude = FisherParams1.position.longitude+LAL_PI/4.0;
	
  if(testACS.verbose_switch==1){fprintf(stdout, "... Computing Fisher Matrix \n");}
  REAL4 Fisher1[testACS.N_inputparams*testACS.N_inputparams];
  //REAL4 Fisher2[testACS.N_inputparams*testACS.N_inputparams];
  REAL4 Fishertot[testACS.N_inputparams*testACS.N_inputparams];
  REAL4 Covariance[testACS.N_inputparams*testACS.N_inputparams];
  LALInspiralComputeFisherMatrix(&status, psd, &FisherParams1, Fisher1, &testACS);
  //LALInspiralComputeFisherMatrix(&status, psd, &FisherParams2, Fisher2, &testACS);
  
  /* WRITE PARAMETER INTO OUTPUT FILES */
  fprintf(fisher_out, "%d\t%d\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t", testACS.seed, params.epoch.gpsSeconds, params.mTot_real8, params.eta_real8, params.phasePNparams[0], params.phasePNparams[1], params.phasePNparams[testACS.testPhaseParam], testACS.SNR, params.d, params.position.longitude, params.position.latitude, params.inc, params.phi, params.psi);
 	fprintf(cov_out, "%d\t%d\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t", testACS.seed, params.epoch.gpsSeconds, params.mTot_real8, params.eta_real8, params.phasePNparams[0], params.phasePNparams[1], params.phasePNparams[testACS.testPhaseParam],testACS.SNR, params.d, params.position.longitude, params.position.latitude, params.inc, params.phi, params.psi);
  
  /* OUTPUT FISHER TO FILE AND SCREEN */
  if(testACS.verbose_switch==1)fprintf(stdout, "FISHER MATRIX \n");
  for (i = 0; i < testACS.N_inputparams*testACS.N_inputparams; ++i)
  {
    Fishertot[i] =  Fisher1[i];// + Fisher2[i];
    fprintf(fisher_out, "%e", Fishertot[i]);
    if(testACS.verbose_switch==1)fprintf(stdout, "%e", Fishertot[i]);
    
    if(i==(testACS.N_inputparams*testACS.N_inputparams-1))
    {
			fprintf(fisher_out, "\n");
			if(testACS.verbose_switch==1)fprintf(stdout, "\n");
		}
    else
    {
			fprintf(fisher_out, "\t");
			if(testACS.verbose_switch==1)fprintf(stdout, "\t");
		}
	}
	
	/* INVERTING FISHER TO GET COVARIANCE */
  LALInspiralInvertMatrix(&status, Fishertot, Covariance, &testACS); 
	
	if(testACS.verbose_switch==1)fprintf(stdout, "INVERSE MATRIX \n");
  for (i = 0; i < testACS.N_inputparams*testACS.N_inputparams; ++i)
  {
    fprintf(cov_out, "%e", Covariance[i]);
    if(testACS.verbose_switch==1)fprintf(stdout, "%e", Covariance[i]);
    
    if(i==(testACS.N_inputparams*testACS.N_inputparams-1))
    {
			fprintf(cov_out, "\n");
			if(testACS.verbose_switch==1)fprintf(stdout, "\n");
		}
    else
    {
			fprintf(cov_out, "\t");
			if(testACS.verbose_switch==1)fprintf(stdout, "\t");
		}
	}
	
  /*********************************************************************
   *
   *  Run the Matrix inversion individually
   * 
   ********************************************************************/ 
  
  //REAL4 testArray[49];
  //REAL4 tempArray[49];
	//testArray[	0	]=	2.309848E-01	;
	//testArray[	1	]=	-3.736297E-04	;
	//testArray[	2	]=	-1.598149E+02	;
	//testArray[	3	]=	6.731071E-03	;
	//testArray[	4	]=	-1.991984E-21	;
	//testArray[	5	]=	-1.846048E-22	;
	//testArray[	6	]=	-3.370836E-23	;
	//testArray[	7	]=	-3.736297E-04	;
	//testArray[	8	]=	8.228739E-07	;
	//testArray[	9	]=	6.149972E-01	;
	//testArray[	10	]=	-1.099523E-05	;
	//testArray[	11	]=	4.735530E-24	;
	//testArray[	12	]=	-3.642607E-25	;
	//testArray[	13	]=	-3.273877E-25	;
	//testArray[	14	]=	-1.598149E+02	;
	//testArray[	15	]=	6.149972E-01	;
	//testArray[	16	]=	9.578669E+05	;
	//testArray[	17	]=	1.140922E+00	;
	//testArray[	18	]=	1.915651E-02	;
	//testArray[	19	]=	6.631099E-03	;
	//testArray[	20	]=	-5.556514E-03	;
	//testArray[	21	]=	6.731071E-03	;
	//testArray[	22	]=	-1.099523E-05	;
	//testArray[	23	]=	1.140922E+00	;
	//testArray[	24	]=	3.314605E-04	;
	//testArray[	25	]=	-7.405499E-23	;
	//testArray[	26	]=	-1.779479E-23	;
	//testArray[	27	]=	-1.218393E-23	;
	//testArray[	28	]=	-1.991984E-21	;
	//testArray[	29	]=	4.735530E-24	;
	//testArray[	30	]=	1.915651E-02	;
	//testArray[	31	]=	-7.405499E-23	;
	//testArray[	32	]=	1.323647E-06	;
	//testArray[	33	]=	4.581856E-07	;
	//testArray[	34	]=	-3.839355E-07	;
	//testArray[	35	]=	-1.846048E-22	;
	//testArray[	36	]=	-3.642607E-25	;
	//testArray[	37	]=	6.631099E-03	;
	//testArray[	38	]=	-1.779479E-23	;
	//testArray[	39	]=	4.581856E-07	;
	//testArray[	40	]=	1.586027E-07	;
	//testArray[	41	]=	-1.329008E-07	;
	//testArray[	42	]=	-3.370836E-23	;
	//testArray[	43	]=	-3.273877E-25	;
	//testArray[	44	]=	-5.556514E-03	;
	//testArray[	45	]=	-1.218393E-23	;
	//testArray[	46	]=	-3.839355E-07	;
	//testArray[	47	]=	-1.329008E-07	;
	//testArray[	48	]=	1.113639E-07	;
	
	//LALInspiralInvertMatrix(&status, testArray, tempArray, testACS);
	
  /*********************************************************************
   *
   *  PLOTTING LENGTH CHANGE WRT TO MTOT AND ETA
   * 
   ********************************************************************/ 	
  
  //FILE *length_out;
  //length_out = fopen("/data/gravwav/tgfli/Fisher/length_eta0p15.dat", "w");
  //INT4 k;
  
  //for(k=1;k<101;k++)
  //{
	//CoherentGW length_waveform;
	//PPNConsistencyParamStruc length_params;
	
	//memset(&length_waveform, 0, sizeof(CoherentGW));
	//length_params.ppn = XLALCreateREAL4Vector( params.ppn->length );
	//XLALCopyPPNConsistencyParamStruc(&params, &length_params);
	
	//length_params.mTot_real8 = 1.0*k;
	//length_params.eta_real8 = 0.15;
	
	//LALGeneratePPNAmpCorInspiral( &status, &length_waveform, &length_params);
	
	//fprintf(stdout, "\r%e\t%e", length_params.mTot_real8, length_waveform.f->data->length*length_params.deltaT);
	//fprintf(length_out, "%e\t%e\n", length_params.mTot_real8, length_waveform.f->data->length*length_params.deltaT);
	
	//XLALDestroyCoherentGW(&length_waveform); 
	//XLALDestroyREAL4Vector( length_params.ppn );
	//}
	//fclose(length_out);
	
	
  /*********************************************************************
   *
   *  Clean Up
   * 
   ********************************************************************/ 
  
  if(testACS.verbose_switch==1){fprintf(stdout, "... Cleaning Up \n");}
  
  /* GSL RNG */
  gsl_rng_free (p);
  
  /* TEMPORARY STRUCTURES */
  XLALDestroyREAL4Vector(params.ppn);
  XLALDestroyREAL8FrequencySeries(psd); psd = NULL;
	
	/* CLOSING OUTPUT FILES */
	fclose(fisher_out);
  fclose(cov_out);  
	
  if(testACS.verbose_switch==1)fprintf(stdout, "Finished FisherTest \n");
  REPORTSTATUS( &status );
	return status.statusCode;
}
