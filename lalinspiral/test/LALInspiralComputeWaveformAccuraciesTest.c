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
#include <stdlib.h>
#include <lal/LALStdlib.h>
#include <lal/LALInspiralComputeWaveformAccuracies.h>

#define EPOCH (315187200000000000LL) /* about Jan. 1, 1990 */


/* A function to convert INT8 nanoseconds to LIGOTimeGPS. */
void
I8ToLIGOTimeGPS( LIGOTimeGPS *output, INT8 input );


NRCSID( LALINSPIRALCOMPUTEWAVEFORMACCURACIESTESTC, "$Id$" );

extern int lalDebugLevel;

int main( int argc, char **argv )
{
  // LAL ERROR REPORTING SYSTEM
  static LALStatus stat;
  lalDebugLevel = 0;

	argc--;
	printf("Running %s\n", argv[0]);
	
  /*********************************************************************
   *
   *  CREATE PARAMETER STRUCTURE
   * 
   ********************************************************************/
   
   INT4 i;
   INT4 j;
   
   REAL8 longitude = 0.0;
   REAL8 latitude = 0.0;
   INT4 gpsSeconds = 800000000;
   INT4 gpsNanoSeconds = 0;
   REAL8 psi = 0.0;
   REAL8 deltaT = 6.1035e-5; //1.0/16384.0;
   REAL8 m1 = 1.35;
   REAL8 m2 = 1.35;
   REAL8 inc = 0.0;
   REAL8 tc = 0.0;
   REAL8 phic = 0.0;
   REAL8 dist = 100.0*LAL_PC_SI*1.0e6;
   REAL8 fStartIn = 40.0;
   //REAL8 fStopIn = 2000.0;
   INT4 ampOrder = 0;
   INT4 phaseOrder = 6;
   
  /*******************************************************************
   * INPUT SETUP                                                     *
   *******************************************************************/ 
   
  PPNConsistencyParamStruc seedParams;
  seedParams.position.longitude = longitude; /*acos(gsl_ran_flat (p, -1.0, 1.0));*/
  seedParams.position.latitude 	= latitude;/*gsl_ran_flat (p, 0, LAL_TWOPI);*/
  seedParams.position.system 		= COORDINATESYSTEM_EQUATORIAL;
  seedParams.epoch.gpsSeconds		= gpsSeconds;
  seedParams.epoch.gpsNanoSeconds = gpsNanoSeconds;
  seedParams.lengthIn = 0;
  seedParams.psi								= psi; /*acos(gsl_ran_flat(p,-1.0,1.0));*/
  seedParams.deltaT 						= deltaT;
  seedParams.mTot_real8 				= m1 + m2;
  seedParams.mTot 							= m1 + m2;
  seedParams.eta_real8  				= m1*m2/( seedParams.mTot_real8*seedParams.mTot_real8 );
  seedParams.eta  				= m1*m2/( seedParams.mTot_real8*seedParams.mTot_real8 );;
  seedParams.delta							= pow( 1.0 - 4.0*seedParams.eta_real8, 0.5 );
  seedParams.inc 								= inc;
  seedParams.tc  								= tc;
  seedParams.phi 								= phic;
  seedParams.d   								= dist;
  seedParams.fStartIn 					= fStartIn;
  seedParams.fStopIn  					= pow(LAL_C_SI, 3.0)/(pow(6.0,1.5)*LAL_PI*LAL_G_SI*(seedParams.mTot_real8)*LAL_MSUN_SI);
  seedParams.ampOrder 					= ampOrder;
  REAL8 seedDeltas[10] = {0.0};
  LALPopulatePhasePNparams(&seedParams, seedDeltas);
  
  seedParams.ppn = XLALCreateREAL4Vector( phaseOrder+1 );
  for (i=0; i<=(phaseOrder); i++)
  {
    if(phaseOrder > 0 && i==1) {seedParams.ppn->data[i] = 0.0;}
    else {seedParams.ppn->data[i] = 1.0;}
  }
  
  /*********************************************************************
   *
   *  SETUP PARAMETERS FOR TAYLORF2
   * 
   ********************************************************************/  
  
  InspiralTemplate seedParamsF2;
  
  seedParamsF2.approximant 	= /*TaylorF2Test;*/AmpCorPPNTest;
  //seedParamsF2.approximant 	= TaylorF2Test;
  seedParamsF2.order				= LAL_PNORDER_THREE_POINT_FIVE;
  seedParamsF2.ampOrder			= LAL_PNORDER_NEWTONIAN;
  seedParamsF2.mass1				= 5.0;			// SHOULD BE IN SOLAR MASSES, SEE 127 AND ONWARDS IN LALINSPIRALSETUP.C
  seedParamsF2.mass2				= 5.0;
  seedParamsF2.fCutoff			= pow(LAL_C_SI, 3.0)/(pow(6.0,1.5)*LAL_PI*LAL_G_SI*(seedParamsF2.mass1 + seedParamsF2.mass2)*LAL_MSUN_SI);
  seedParamsF2.fLower				= 40.0;
  seedParamsF2.tSampling		= 4096.0;	// SAMPLING RATE
  seedParamsF2.distance			= 100.0;		// IN MPC
  seedParamsF2.signalAmplitude = 1.0; // TO BE FILLED BY SEPERATE FUNCTION? CHECK IN LALNEST.C
  seedParamsF2.startPhase		= 0.0;
  seedParamsF2.startTime		= 0.0;
  seedParamsF2.ieta					= 1;
  seedParamsF2.massChoice		= m1Andm2;
  LALInspiralParameterCalc(&stat,&seedParamsF2);
  LALInspiralRestrictedAmplitude(&stat,&seedParamsF2);
  
  expnCoeffs ak;
  expnFunc expnFunction;
  
  memset(&ak,0,sizeof(expnCoeffs));
  
  LALInspiralSetup(&stat, &ak, &seedParamsF2);
  LALInspiralChooseModel(&stat, &expnFunction, &ak, &seedParamsF2);
  
  printf("tc = %e\n",seedParamsF2.tC );
  
  /* FROM LALINSPIRALSETUP.C
		00120    ieta determines the nature of the waveforms:
		00121    ieta=0 testmass waveforms
		00122    ieta=1 comparable mass waveforms.
		00123 
	*/
  /*
  LALInspiralParameterCalc(&stat,&seedParamsF2);
  REPORTSTATUS( &stat );
	LALInspiralRestrictedAmplitude(&stat,&seedParamsF2);
	REPORTSTATUS( &stat );
	
	printf("%e\n", seedParamsF2.signalAmplitude);
	printf("%e\n", seedParamsF2.fCutoff);
  
  REAL4FrequencySeries *signalvec;
  INT4 NF2 = 1000;
  REAL8 df = seedParamsF2.tSampling/(2*NF2);
  signalvec 				= XLALCreateREAL4FrequencySeries("exact_ht", &(seedParams.epoch), 0.0, df, &lalStrainUnit, 2*NF2);
  REAL8 dphis[9] = {0.0};
  
  // TEST F2 WAVEFORM OUTPUT
  LALInspiralStationaryPhaseApprox2Test (&stat,signalvec->data,&seedParamsF2,dphis);
  REPORTSTATUS( &stat );
  
  FILE *F2GW_file;
  F2GW_file = fopen("f2gw.dat", "w");
  
  REAL8 sigRe = 0.0;
  REAL8 sigIm = 0.0;
  REAL8 dynRange = 1.0E27;
  
  for(i=1;i<NF2/2-1;i++){
  	sigRe = signalvec->data->data[i] * dynRange;
  	sigIm = signalvec->data->data[2*NF2-i] * dynRange;
  	fprintf(F2GW_file, "%e\t%e\t%e\t%e\n", i*df, sigRe, sigIm,  sqrt(sigRe*sigRe+sigIm*sigIm)/dynRange);
  }
  
  fclose(F2GW_file);
  */
  //exit(-1);
  
  /*********************************************************************
   *
   *  COMPUTE SNR
   * 
   ********************************************************************/                               
	PPNConsistencyParamStruc SNRParams;
	//FILE *SNR_file;
	//SNR_file = fopen("SNR_mtot.dat","w");
	
	SNRParams.ppn = XLALCreateREAL4Vector(seedParams.ppn->length);
	XLALCopyPPNConsistencyParamStruc (&seedParams,&SNRParams); 
	//SNRParams.mTot_real8 = 3.0+i;
	//SNRParams.eta_real8 = 0.1+i/100.0;
	//SNRParams.d = (100.0+50*i)*LAL_PC_SI*1.0e6;
	//SNRParams.fStopIn = pow(LAL_C_SI, 3.0)/(pow(6.0,1.5)*LAL_PI*LAL_G_SI*(SNRParams.mTot_real8)*LAL_MSUN_SI);
	
	REAL8 SNRDeltas[10] = {0.0};
																			 
	LALPopulatePhasePNparams(&SNRParams, SNRDeltas); 
	
	printf("computing SNR \n");
	REAL8 SNR = 0.0;
	LALInspiralComputeSNR (&stat,&seedParamsF2,SNRDeltas,&SNR);  
	
	printf("SNR = %e\n", SNR);
	
		// CLOSING PROGRAM
  REPORTSTATUS( &stat );
  return stat.statusCode;
														
	//fprintf(SNR_file, "%e\t%e\n", SNRParams.mTot_real8, SNR);														
	//memset(&SNRParams, 0, sizeof(PPNConsistencyParamStruc) );
	
	//fclose(SNR_file);
  
  /*********************************************************************
   *
   *  COMPUTE DHINDH
   * 
   ********************************************************************/
	
	PPNConsistencyParamStruc modelParams;
	PPNConsistencyParamStruc exactParams;
	
	FILE *accuracies_file;
	accuracies_file = fopen("accuracies_phi0phi2.dat","w");
  
  for(i=0;i<21;i++){
  	for(j=0;j<21;j++){
			modelParams.ppn = XLALCreateREAL4Vector(seedParams.ppn->length);
			XLALCopyPPNConsistencyParamStruc (&seedParams,&modelParams); 
			exactParams.ppn = XLALCreateREAL4Vector(seedParams.ppn->length);
			XLALCopyPPNConsistencyParamStruc (&seedParams,&exactParams); 
			
			REAL8 exactDeltas[10] = {0.0};
			
			REAL8 modelDeltas[10] = {0.0};
			modelDeltas[0] = -0.0001+0.00001*i;
			modelDeltas[2] = -0.001+0.0001*j;
																					 
			LALPopulatePhasePNparams(&exactParams, exactDeltas);
			LALPopulatePhasePNparams(&modelParams, modelDeltas);
			
			REAL8 dhindh = 0.0;
			LALInspiralComputeDhInDh (&stat,&exactParams,&modelParams,fStartIn,seedParams.fStopIn,&dhindh);
			fprintf(accuracies_file, "%e", dhindh*dhindh/(SNR*SNR));
			fprintf(accuracies_file, "\t");
			printf("%e\t%e\t%e\t%e\t%e\n", modelDeltas[0], modelDeltas[2], dhindh*dhindh, SNR*SNR, dhindh*dhindh/(SNR*SNR));
			
			memset(&modelParams, 0, sizeof(PPNConsistencyParamStruc) );
			memset(&exactParams, 0, sizeof(PPNConsistencyParamStruc) );
		}
		fprintf(accuracies_file, "\n");

	}													
	// CLOSING PROGRAM
  REPORTSTATUS( &stat );
  return stat.statusCode;
}

/* A function to convert INT8 nanoseconds to LIGOTimeGPS. */
void
I8ToLIGOTimeGPS( LIGOTimeGPS *output, INT8 input )
{
  INT8 s = input / 1000000000LL;
  output->gpsSeconds = (INT4)( s );
  output->gpsNanoSeconds = (INT4)( input - 1000000000LL*s );
  return;
}

/* </lalVerbatim> */
