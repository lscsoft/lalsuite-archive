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
#include <stdlib.h>
#include <lal/LALStdlib.h>
#include "LALCalibrationErrors.h"
#include <lal/LALDatatypes.h>

NRCSID( LALCALIBRATIONERRORSC, "$Id$" );
extern int enable_calamp;
extern int enable_calfreq;
extern int enable_gradual_cal;
extern REAL8 calibration_percent;
extern REAL8 calibration_out_max;
//static LALStatus stat;
extern double *CalAmpFacs;

void SampleCalibrationErrorsAmplitude(REAL8 *logF, INT4 length, INT4 IFO, INT4 seed, REAL8 *errors){
    INT4 i;
    //UINT4 length;
    
    // SETTING UP GSL RANDOM GENERATOR
    const gsl_rng_type *type;           // RNG type
  	gsl_rng *p;                         // Generator
  	gsl_rng_env_setup();                // Setup environment
  	gsl_rng_default_seed = seed;          // vary generation sequence
  	type = gsl_rng_default;             // set RNG type to default
  	p = gsl_rng_alloc (type);           // Set RNG type
    
    REAL8 stddev[3]={0.0};
    switch (IFO) {
        case 1:
            stddev[0]=0.104;
            stddev[1]=0.154;
            stddev[2]=0.242;
            break;
        case 2:
            stddev[0]=0.144;
            stddev[1]=0.139;
            stddev[2]=0.138;
            break;
        case 3:
            stddev[0]=0.10;
            stddev[1]=0.10;
            stddev[2]=0.20;
            break;
        default:
            fprintf(stderr,"Unknown IFO! Valid codes are H1, L1, V1. Aborting\n");
            exit(-1);
            break;
    }
    
    //length = sizeof(logF)/sizeof(*logF);
    //REAL8 errors[length];
    for (i=0; i<length; i++) {
        if (logF[i]>=log10(1.0) && logF[i]<log10(2000.0)) {
            errors[i]=gsl_ran_gaussian(p, stddev[0]);
        } 
        else if (logF[i]>=log10(2000.0) && logF[i]<log10(4000.0)){
            errors[i]=gsl_ran_gaussian(p, stddev[1]);
        } 
        else if (logF[i]>=log10(4000.0) && logF[i]<=log10(6000.0)){
            errors[i]=gsl_ran_gaussian(p, stddev[2]);
				}
        errors[i]=1.05; // The errors represent a ratio between amplitudes, and then must be centered around 1.0 
        }		
    return;    
}

/* function to return the random phase calibration errors in the logfrequency array */

void SampleCalibrationErrorsPhase(REAL8 *logF, INT4 length, INT4 IFO, INT4 seed, REAL8 *errors){
    INT4 i;
    //UINT4 length;
    
    // SETTING UP GSL RANDOM GENERATOR
		const gsl_rng_type *type;           // RNG type
  	gsl_rng *p;                         // Generator
  	gsl_rng_env_setup();                // Setup environment
  	gsl_rng_default_seed = seed;    // vary generation sequence
  	type = gsl_rng_default;             // set RNG type to default
  	p = gsl_rng_alloc (type);           // Set RNG type
    

    REAL8 stddev[6]={0.0};
    switch (IFO) {
        case 1:  //H1
            stddev[0]=4.5;  // 1-500
            stddev[1]=4.5;  // 500-1000
            stddev[2]=4.5;   // 1k-2k
            stddev[3]=4.9;   // 2k -2.8k
            stddev[4]=4.9;   //2.8k-4k
            stddev[5]=5.8;   // >4k 
            break;
        case 2: //L1
            stddev[0]=4.2;  // 1-500
            stddev[1]=4.2;  // 500-1000
            stddev[2]=4.2;   // 1k-2k
            stddev[3]=3.6;   // 2k -2.8k
            stddev[4]=3.6;   //2.8k-4k
            stddev[5]=3.3;   // >4k 
            break;
        case 3: //V1
            stddev[0]=2.2918;  // 1-500
            stddev[1]=0.5729;  // 500-1000
            stddev[2]=6.87;   // 1k-2k
            stddev[3]=6.87;   // 2k -2.8k
            stddev[4]=360.0*7e-06;   //2.8k-4k
            stddev[5]=360.0*7e-06;   // >4k 
            break;
        default:
            fprintf(stderr,"Unknown IFO! Valid codes are H1, L1, V1. Aborting\n");
            exit(-1);
            break;
    }
    
    for (i=0; i<length; i++) {
        if (logF[i]>=log10(1.0) && logF[i]<log10(500.0) && IFO==3) {
            errors[i]=gsl_ran_gaussian(p, (stddev[0]+0.00286479*pow(10.0,logF[i])));
        }
        else if (logF[i]>=log10(1.0) && logF[i]<log10(500.0)) {
            errors[i]=gsl_ran_gaussian(p, stddev[0]);
        }
        else if (logF[i]>=log10(500.0) && logF[i]<log10(1000.0) && IFO==3) {
            errors[i]=gsl_ran_gaussian(p, (stddev[1]+0.0063*pow(10.0,logF[i])));
        } 
        else if (logF[i]>=log10(500.0) && logF[i]<log10(1000.0)) {
            errors[i]=gsl_ran_gaussian(p, stddev[1]);
        }
        else if (logF[i]>=log10(1000.0) && logF[i]<log10(2000.0)) {
            errors[i]=gsl_ran_gaussian(p, stddev[2]);
        } 
        else if (logF[i]>=log10(2000.0) && logF[i]<log10(2800.0)) {
            errors[i]=gsl_ran_gaussian(p, stddev[3]);
        }
        else if (logF[i]>=log10(2800.0) && logF[i]<log10(4000.0) && IFO==3){
            errors[i]=gsl_ran_gaussian(p, (stddev[4]*pow(10.0,logF[i])));
        } 
        else if (logF[i]>=log10(2800.0) && logF[i]<log10(4000.0)){
            errors[i]=gsl_ran_gaussian(p, stddev[4]);
        } 
        else if (logF[i]>=log10(4000.0) && logF[i]<=log10(6000.0) && IFO==3){
            errors[i]=gsl_ran_gaussian(p, (stddev[5]*pow(10.0,logF[i])));
				}
         else if (logF[i]>=log10(4000.0) && logF[i]<=log10(6000.0)){
            errors[i]=gsl_ran_gaussian(p, stddev[5]);
				}
        //errors[i]*=LAL_PI/180.0;			
         errors[i]=0.0;
        //printf("error[%i] | logf = %e | error = %e\n", i, logF[i], errors[i]);
    }
    return; /* this is in radians ! */   
}

void FitNoiseRealisation(LALStatus *status,INT4	R,	INT4 N,	REAL8    *y,REAL8 dlogf,REAL8	*D)
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
  
  INITSTATUS( status, "LALSavitskyGolayFilter", LALCALIBRATIONERRORSC);
  ATTATCHSTATUSPTR(status);	 
	
	/*********************************************************************
	 * 
	 * Read input
	 * 
	 *********************************************************************/	   
	
  /* LOAD IN TIMESERIES PARAMETERS */
  INT4 M = (N-1)/2;
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
  INT4 Factorial = 1;
  
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
		gsl_matrix_set(U, i, i, pow(dlogf,i));
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
  
  //printf("Smoothing \n");
	/* READ DATA POINTS INTO VECTOR */
	for(j=0;j<2*M+1;j++)
	{
		gsl_vector_set(ym, j, y[j]);		
	}
	
	/* a = c*y */
	gsl_blas_dgemv( CblasNoTrans, 
								 1.0, c, ym, 
								 0.0, a );
	
	for(k=0; k<R; k++)
	{
		D[k] = Factorial*gsl_vector_get(a, k);
	}
	
	gsl_vector_set_zero (ym);
	gsl_vector_set_zero (a);
	
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

void InvertMatrixSVD (gsl_matrix *A,gsl_matrix	*InvA,	int	N)
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

REAL8 ConvertCoefficientsToFunction(REAL8 *coeff, REAL8 f, REAL8 cen)
{
	REAL8 logF = log10(f)-cen; // FIT USED CEN AS CENTRAL POINT!
	REAL8 output = 0.0;
	
	INT4 i;
	
	for(i=0;i<7;i++)
	{
		output += coeff[i]*pow(logF, (REAL8) i);
	}
	
	return output;

}

/* function to return a frequency array logarithmic spaced */
REAL8 GenerateFrequencySamples(REAL8 f_min, REAL8 f_max, UINT4 Length){
    REAL8 logFreq[Length];
    UINT4 i;
    REAL8 step=(log(f_max)-log(f_min))/Length;
    for (i=0; i<Length; i++) {
        logFreq[i]=log(f_min)-step*i;
    }
    return *logFreq;
}
/* function to return the random amplitude calibration errors in the logfrequency array */


void CreateErrorStreams(LALMCMCInput *inputMCMC,CHAR *IFOname, int i,int seed){
/* This only fills the MCMC structures that hold the amplitude and phase errors. This function DOES NOT modify the actual noise and WF */
//injTime = injTable->geocent_end_time.gpsSeconds + 1.0E-9 * injTable->geocent_end_time.gpsNanoSeconds;
    int IFOnum=0;
    UINT4 j;
    REAL8 fstart = 1.0;
    REAL8 fend = 4000.0;
    INT4 length= 13;
    REAL8 logfstart = log10(fstart);
    REAL8 logfend = log10(fend);
    REAL8 deltalogf = (logfend-logfstart)/(length-1);
    REAL8 logfcur[length];
    INT4 I;
    for(I = 0; I<length; I++){
        logfcur[I] = logfstart + deltalogf*I;
    }

    REAL8 phaseErr[length];
    REAL8 ampErr[length];
    INT4 FitOrder = 7;
    REAL8 phaseFit[FitOrder];
    REAL8 ampFit[FitOrder];
    REAL8 cen = logfcur[(length-1)/2]; 
    static LALStatus stat;
    if(!strcmp(IFOname,"H1")){IFOnum =1;}
    if(!strcmp(IFOname,"L1")){IFOnum =2;}
    if(!strcmp(IFOname,"V1")){IFOnum =3;}
    switch(IFOnum) {
		case 1:
            SampleCalibrationErrorsPhase(logfcur, length, 1, seed, phaseErr);
            SampleCalibrationErrorsAmplitude(logfcur, length, 1, seed+4101, ampErr);
			break;
        case 2:
            SampleCalibrationErrorsPhase(logfcur, length, 2, seed, phaseErr);
            SampleCalibrationErrorsAmplitude(logfcur, length, 2, seed+4101, ampErr);
			break;
		case 3:
            SampleCalibrationErrorsPhase(logfcur, length, 3, seed, phaseErr);
            SampleCalibrationErrorsAmplitude(logfcur, length, 3, seed+4101, ampErr);
            break;
        default:
			fprintf(stderr,"Unknown interferometer %s. Valid codes: H1 L1 V1\n",IFOname); exit(-1);
    }
    FitNoiseRealisation(&stat,FitOrder,length,phaseErr,deltalogf,phaseFit);
    FitNoiseRealisation(&stat,FitOrder,length,ampErr,deltalogf,ampFit);
    for(j=0;j<inputMCMC->invspec[i]->data->length;j++){
    inputMCMC->calibAmplitude[i]->data->data[j]= ConvertCoefficientsToFunction(ampFit,j*inputMCMC->deltaF,cen);
    inputMCMC->calibPhase[i]->data->data[j]=ConvertCoefficientsToFunction(phaseFit,j*inputMCMC->deltaF,cen);
    }
}

void ApplyCalibrationErrorsToData(LALMCMCInput inputMCMC, COMPLEX16FrequencySeries *CalibNoise,CHAR *IFOname, int i,int seed){
/* Modify the noise datastream  */
//injTime = injTable->geocent_end_time.gpsSeconds + 1.0E-9 * injTable->geocent_end_time.gpsNanoSeconds;
    int IFOnum=0;
    UINT4 j;
    REAL8 fstart = 1.0;
    REAL8 fend = 4000.0;
    INT4 length= 13;
    REAL8 logfstart = log10(fstart);
    REAL8 logfend = log10(fend);
    REAL8 deltalogf = (logfend-logfstart)/(length-1);
    REAL8 logfcur[length];
    INT4 I;
    for(I = 0; I<length; I++){
        logfcur[I] = logfstart + deltalogf*I;
    }

    REAL8 phaseErr[length];
    REAL8 ampErr[length];
    INT4 FitOrder = 7;
    REAL8 phaseFit[FitOrder];
    REAL8 ampFit[FitOrder];
    REAL8 cen = logfcur[(length-1)/2]; 
    static LALStatus stat;
    if(!strcmp(IFOname,"H1")){IFOnum =1;}
    if(!strcmp(IFOname,"L1")){IFOnum =2;}
    if(!strcmp(IFOname,"V1")){IFOnum =3;}
    switch(IFOnum) {
		case 1:
            SampleCalibrationErrorsPhase(logfcur, length, 1, seed, phaseErr);
            SampleCalibrationErrorsAmplitude(logfcur, length, 1, seed+4101, ampErr);
			break;
        case 2:
            SampleCalibrationErrorsPhase(logfcur, length, 2, seed, phaseErr);
            SampleCalibrationErrorsAmplitude(logfcur, length, 2, seed+4101, ampErr);
			break;
		case 3:
            SampleCalibrationErrorsPhase(logfcur, length, 3, seed, phaseErr);
            SampleCalibrationErrorsAmplitude(logfcur, length, 3, seed+4101, ampErr);
            break;
        default:
			fprintf(stderr,"Unknown interferometer %s. Valid codes: H1 L1 V1\n",IFOname); exit(-1);
    }
    FitNoiseRealisation(&stat,FitOrder,length,phaseErr,deltalogf,phaseFit);
    FitNoiseRealisation(&stat,FitOrder,length,ampErr,deltalogf,ampFit);
    /* Modify the fake noise PSD */
    for(j=0;j<inputMCMC.invspec[i]->data->length;j++){
        if(enable_calamp){  inputMCMC.invspec[i]->data->data[j]/=((REAL8)CalAmpFacs[i]*(REAL8)CalAmpFacs[i]);}
        else if(enable_calfreq){ inputMCMC.invspec[i]->data->data[j]/=((ConvertCoefficientsToFunction(ampFit,j*inputMCMC.deltaF,cen)*ConvertCoefficientsToFunction(ampFit,j*inputMCMC.deltaF,cen)));}
    }

    /* Modifiy the noise stream in the frequency domain */
    if(enable_calamp){
        for(j=0;j<inputMCMC.invspec[i]->data->length;j++) {
            inputMCMC.stilde[i]->data->data[j].re*=(REAL8)CalAmpFacs[i];
            inputMCMC.stilde[i]->data->data[j].im*=(REAL8)CalAmpFacs[i];
        }
    }
    if(enable_calfreq){
        fprintf(stderr,"Calling calibpolar from ApplyCalibrationErrorsToData \n");
        CalibPolar(inputMCMC.stilde[i],CalibNoise,IFOname, seed);	
        for(j=0;j<inputMCMC.invspec[i]->data->length;j++){
            inputMCMC.stilde[i]->data->data[j].re = CalibNoise->data->data[j].re;
            inputMCMC.stilde[i]->data->data[j].im = CalibNoise->data->data[j].im;
        }
    }
}


void ApplyCalibrationErrorsToWaveform(COMPLEX16FrequencySeries *injF,COMPLEX16FrequencySeries *CalibInj,CHAR *IFOname, int i,int seed){
    UINT4 j;
    REAL8 fstart = 1.0;
    REAL8 fend = 4000.0;
    INT4 length= 13;
    REAL8 logfstart = log10(fstart);
    REAL8 logfend = log10(fend);
    REAL8 deltalogf = (logfend-logfstart)/(length-1);
    REAL8 logfcur[length];
    INT4 I;
    for(I = 0; I<length; I++){
        logfcur[I] = logfstart + deltalogf*I;
    }
        
    if(enable_calamp){
        fprintf(stderr,"Calling calamp \n");
        for(j=0;j<injF->data->length;j++) {
            
            injF->data->data[j].re*=(REAL8)CalAmpFacs[i];
            injF->data->data[j].im*=(REAL8)CalAmpFacs[i];
        }
    }
    if(enable_calfreq){
        fprintf(stderr,"Calling calibpolar from ApplyCalibrationErrorsToWaveform \n");
        CalibPolar(injF,CalibInj,IFOname,seed);    
        for(j=0;j<injF->data->length;j++){
            injF->data->data[j].re = CalibInj->data->data[j].re;
            injF->data->data[j].im = CalibInj->data->data[j].im;
        }
    }
}

void CalibPolar(COMPLEX16FrequencySeries *injF, COMPLEX16FrequencySeries *calibInjF, CHAR *IFOname,int seed){

    fprintf(stderr, "seed in calib polar %i \n",seed);
    REAL8 amplitude=0.0;
    REAL8 phase=0.0;
    REAL8 deltaf=0.0;
    UINT4 j;
    deltaf=injF->deltaF;
    REAL8 fstart = 1.0;
    REAL8 fend = 4000.0;
    INT4 length= 13;
    REAL8 logfstart = log10(fstart);
    REAL8 logfend = log10(fend);
    REAL8 deltalogf = (logfend-logfstart)/(length-1);
    REAL8 logfcur[length];
    INT4 I;
    static LALStatus stat;
    for(I = 0; I<length; I++){
        logfcur[I] = logfstart + deltalogf*I;
        }
    REAL8 phaseErr[length];
    REAL8 ampErr[length];
    INT4 FitOrder = 7;
    REAL8 phaseFit[FitOrder];
    REAL8 ampFit[FitOrder];
    REAL8 cen = logfcur[(length-1)/2]; 


    int IFO;
    if(!strcmp(IFOname,"H1")){IFO =1;}
	if(!strcmp(IFOname,"L1")){IFO =2;}
	if(!strcmp(IFOname,"V1")){IFO =3;}
	switch(IFO) {
        case 1:
            SampleCalibrationErrorsPhase(logfcur, length, 1, seed, phaseErr);
            SampleCalibrationErrorsAmplitude(logfcur, length, 1, seed+4101, ampErr);
			break;
		case 2:
            SampleCalibrationErrorsPhase(logfcur, length, 2, seed, phaseErr);
            SampleCalibrationErrorsAmplitude(logfcur, length, 2, seed+4101, ampErr);
            break;
        case 3:
            SampleCalibrationErrorsPhase(logfcur, length, 3, seed, phaseErr);
            SampleCalibrationErrorsAmplitude(logfcur, length, 3, seed+4101, ampErr);
			break;
        default:
            fprintf(stderr,"Unknown interferometer %s. Valid codes: H1 L1 V1\n",IFOname); exit(-1);
    }
    FitNoiseRealisation(&stat,FitOrder,length,phaseErr,deltalogf,phaseFit);
    FitNoiseRealisation(&stat,FitOrder,length,ampErr,deltalogf,ampFit);
    //FILE *calibout;
    //char caliboutname[100];
    //sprintf(caliboutname,"phase_errors_%s.dat",IFOname);
    //calibout=fopen(caliboutname,"w");
	for(j=0;j<injF->data->length;j++){
        if(!enable_calamp){
            amplitude=ConvertCoefficientsToFunction(ampFit,j*deltaf,cen)*sqrt(pow(injF->data->data[j].re,2.0)+pow(injF->data->data[j].im,2.0));
        }
        else {
            amplitude=sqrt(pow(injF->data->data[j].re,2.0)+pow(injF->data->data[j].im,2.0));
        }
        //fprintf(calibout,"%g\t%g\n",j*deltaf,ConvertCoefficientsToFunction(phaseFit,j*deltaf,cen));
        phase=ConvertCoefficientsToFunction(phaseFit,j*deltaf,cen)+atan2(injF->data->data[j].im,injF->data->data[j].re);
        calibInjF->data->data[j].re=amplitude*cos(phase);
        calibInjF->data->data[j].im=amplitude*sin(phase);
		//fprintf(calibout,"%g\t%g\t%g\n",j*deltaf,amplitude,phase);
        
    }
    
	//fclose(calibout);
}

/*
void OLD_CalibPolar(COMPLEX16FrequencySeries *injF, COMPLEX16FrequencySeries *calibInjF, CHAR *IFOname){
    
	REAL8 amplitude=0.0;
        REAL8 phase=0.0;
        REAL8 deltaf=0.0;
        UINT4 j;
        //FILE *calibout;
        char caliboutname[100];
        if(isWavesDir == 1){
            fprintf(stderr,"waves directory is present\n");
            fprintf(stderr,"Writing calibrated waves \n");
            sprintf(caliboutname,"./waves/calibwave_%s_%9.0f.dat",IFOname,InjTime);}
        else {
			fprintf(stderr,"waves directory is not present\n");
            fprintf(stderr,"Writing calibrated waves on the run  directory.\n");
            sprintf(caliboutname,"calibwave_%s_%9.0f.dat",IFOname,InjTime);}
 
        calibout=fopen(caliboutname,"w");
 /                 
      deltaf=injF->deltaF;
		int IFO;
		if(!strcmp(IFOname,"H1")){IFO =1;}
		if(!strcmp(IFOname,"L1")){IFO =2;}
		if(!strcmp(IFOname,"V1")){IFO =3;}
		switch(IFO) {
			case 1:
				R_A=&Amp_H1;
				R_PH=&Ph_H1;
				break;
			case 2:
				R_A=&Amp_L1;
				R_PH=&Ph_L1;
				break;
			case 3:
				R_A=&Amp_V1;
				R_PH=&Ph_V1;
				break;
			default:
				fprintf(stderr,"Unknown interferometer %s. Valid codes: H1 L1 V1\n",IFOname); exit(-1);
		}
		for(j=0;j<injF->data->length;j++){
            	if(!enable_calamp){
                amplitude=R_A(j*deltaf)*sqrt(pow(injF->data->data[j].re,2.0)+pow(injF->data->data[j].im,2.0));
                }
                else {
                amplitude=sqrt(pow(injF->data->data[j].re,2.0)+pow(injF->data->data[j].im,2.0));
                }
              	phase=R_PH(j*deltaf)+atan2(injF->data->data[j].im,injF->data->data[j].re);
                calibInjF->data->data[j].re=amplitude*cos(phase);
               	calibInjF->data->data[j].im=amplitude*sin(phase);
		//fprintf(calibout,"%g\t%g\t%g\n",j*deltaf,amplitude,phase);
       		}
	//fclose(calibout);
       	}
*/
// Until we have something for V1, we use L1 quantities.

/*
REAL8 Amp_H1(REAL8 f){
		double output = 1.0;

		if(f>60.0 && f<=100.0)
			output = 0.000144921*f+0.953962+2.19779*pow(f,-1.0);

		if(f>100.0 && f<=150.0)
			output = -1.65116e-05*f+0.991484+0.07191*pow(f,-1.0);

		if(f>150.0 && f<=318.0)
			output = 1.42451e-05*f+0.98561+0.271245*pow(f,-1.0);

		if(f>318.0 && f<=500.0)
			output = 3.04006e-05*f+0.977116+1.35158*pow(f,-1.0);

		// return a constant 
		//output= output-1.0;
                // return 1.0+(calibration_out_max*output-1.0)*calibration_percent;  // this with gradual_cal
		//return 1+(sqrt(output*output + calibration_systematic_H1_AM*calibration_systematic_H1_AM))*calibration_percent;
                return 1.035;
}
REAL8 Amp_L1(REAL8 f){
		double output = 1.0;
                //new fit updated by Salvatore V. 01/05/11
		if(f>60.0 && f<=150.0)
output = -1.95484+0.170278*f-0.388292e-02*pow(f,2)+0.448703e-04*pow(f,3)-2.71089e-07*pow(f,4)+7.35001e-10*pow(f,5)-1.11474e-13*pow(f,6)-2.309e-15*pow(f,7);
	        if(f>150.0 && f<=500.0)
output = 0.808333+0.453077e-02*f-0.390407e-04*pow(f,2)+1.77562e-07*pow(f,3)-4.72241e-10*pow(f,4)+7.34368e-13*pow(f,5)-6.15246e-16*pow(f,6)+2.12688e-19*pow(f,7);
		// return a constant
//		output = output-1.0;
//                return 1.0+(output*calibration_out_max-1.0)*calibration_percent;  // this with gradual_cal
//		return 1+(sqrt(output*output + calibration_systematic_L1_AM*calibration_systematic_L1_AM))*calibration_percent;
		return 1.15;
}
REAL8 Amp_V1(REAL8 f){
		double output = 1.0;
                //new fit updated by Salvatore V. 01/05/11
                if(log10(f)>1 && log10(f)<=1.48443)
                { 
		output = -0.876892 + 6.65445*pow(log10(f),1.0) - 8.81424*pow(log10(f), 2.0) + 5.12016*pow(log10(f),3.0) - 1.09537*pow(log10(f),4.0);
                }
	 	else if(log10(f)>1.48443 && log10(f)<=2.02130)
	        {
		output = -16.7486 + 41.0096*pow(log10(f),1.0) - 35.1488*pow(log10(f), 2.0) + 13.2676*pow(log10(f),3.0) - 1.86709*pow(log10(f),4.0);
	        }
		else if(log10(f)>2.02130 && log10(f)<=2.45144)
	        {
		output = 172.466866 - 304.746*pow(log10(f),1.0) + 202.191*pow(log10(f), 2.0) - 59.3845*pow(log10(f),3.0) + 6.51706*pow(log10(f),4.0);
               	} 
	 	else if(log10(f)>2.45144 && log10(f)<= 3.08153)
	        {
		output = 43.05292 - 60.222*pow(log10(f),1.0) + 32.2545*pow(log10(f), 2.0) - 7.66187*pow(log10(f),3.0) + 0.681219*pow(log10(f),4.0);
	        }
		
                //return a constant
		output = output-1.0;
//                return 1.0+(output*calibration_out_max-1.0)*calibration_percent;  // this with gradual_cal
		return 1.0 +(sqrt(output*output + calibration_systematic_V1_AM*calibration_systematic_V1_AM))*calibration_percent;;

}
REAL8 Ph_H1(REAL8 f){
		double output = 0.0;

		if(f>60.0 && f<=80.0)
			output = 114.005-6.23854*f+0.127996*pow(f,2)-0.00116878*pow(f,3)+4.00732e-06*pow(f,4);
		if(f>80.0 && f<=500.0)
			output = -0.0701154 +0.0170887*log(0.914066*f)-15.5936*pow(f,-1);
                // convert in rads /
//                return (LAL_PI*output/180.0)*calibration_percent*calibration_out_max;
		//return (LAL_PI*sqrt(output*output + calibration_systematic_H1_PH*calibration_systematic_H1_PH)/180.0)*calibration_percent;
                return (LAL_PI/180.0)*2.6;
}
REAL8 Ph_L1(REAL8 f){
		double output = 0.0;

		if(f>60.0 && f<=110.0)
			output = -69.493+4.77314*f-0.123966*pow(f,2)+0.0015403*pow(f,3)-9.24682e-06*pow(f,4)+2.16121e-08*pow(f,5);

		if(f>110.0 && f<=500.0)
			output = -0.040558+0.315112+0.012216*f-0.00022649*pow(f,2)+9.75241e-07*pow(f,3)-1.72514e-09*pow(f,4)+1.11536e-12*pow(f,5);

                // convert in rads /
//		return  (LAL_PI*sqrt(output*output+calibration_systematic_L1_PH*calibration_systematic_L1_PH)/180.0)*calibration_percent;
//                return (LAL_PI*output/180.0)*calibration_percent*calibration_out_max;
 		return (LAL_PI/180.0)*3.8;
}
REAL8 Ph_V1(REAL8 f){
		double output = 0.0;

                if(log10(f)>1 && log10(f)<=1.64080)
            	{
		output = 5.1621 - 16.1787*pow(log10(f),1.0) + 18.7131*pow(log10(f), 2.0) - 9.45203*pow(log10(f),3.0) + 1.75489*pow(log10(f),4.0);
	        }
		else if(log10(f)>1.64080 && log10(f)<=2.15234)
	        {
		output = -76.8328 + 164.194*pow(log10(f),1.0) - 131.131*pow(log10(f), 2.0) + 46.3381*pow(log10(f),3.0) - 6.10831*pow(log10(f),4.0);
	        }
		else if(log10(f)>2.15234 && log10(f)<=3.2)
	        {
		output = -0.392164e-02 + 2.66324 - 3.48059*pow(log10(f),1.0) + 1.76342*pow(log10(f), 2.0) - 0.415683*pow(log10(f),3.0) + 0.0393829*pow(log10(f),4.0);
	        }
               // Virgo data are already in rads /
//                return output*calibration_percent*calibration_out_max;
		return output*calibration_percent;
}
*/

/* </lalVerbatim> */
