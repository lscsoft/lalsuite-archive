#include <stdlib.h>
#include <gsl/gsl_vector.h>



/*#include <lal/LALDatatypes.h>
#include <lal/LALSimInspiral.h>
#include <lal/TimeSeries.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/TimeSeries.h>
#include <lal/Units.h>

#include <lal/LALSimIMR.h>
#include <lal/Date.h>

#include <lal/SeqFactories.h>*/

#include "LALSimIMREOBNRv2.h"
#include "LALSimIMRSpinEOB.h"
#include "LALSimFindAttachTime.h"
#include "LALSimIMRSpinEOBHamiltonian.c"





double  XLALSimLocateOmegaTime(
    REAL8Array *dynamicsHi,
    unsigned int numdynvars,
    unsigned int retLenHi,
    SpinEOBParams   seobParams,
    SpinEOBHCoeffs  seobCoeffs,
    REAL8 m1,
    REAL8 m2,
    int *found
        )
{        
    /* 
    * Locate merger point (max omega), 
    * WaveStep 1.1: locate merger point */
    int debugPK = 0;
    int debugRD = 1;
    FILE *out = NULL; 
    gsl_spline    *spline = NULL;
    gsl_interp_accel *acc = NULL;
    
    if (debugPK) {debugRD = 0;}
  
    unsigned int peakIdx, i, j;
    REAL8Vector *values = NULL;
    REAL8Vector *dvalues = NULL;
    REAL8Vector *omegaHi = NULL;
    
    if ( !(values = XLALCreateREAL8Vector( numdynvars )) )
    {
        XLAL_ERROR(  XLAL_ENOMEM );
    }
    if ( !(dvalues = XLALCreateREAL8Vector( numdynvars )) )
    {
        XLAL_ERROR(  XLAL_ENOMEM );
    }
    if ( !(omegaHi = XLALCreateREAL8Vector( retLenHi )) )
    {
        XLAL_ERROR(  XLAL_ENOMEM );
    }
    REAL8 rdotvec[3] = {0,0,0};
    REAL8 rvec[3] = {0,0,0};
    REAL8 rcrossrdot[3] = {0,0,0};
    REAL8Vector timeHi;
  
    timeHi.length = retLenHi;
    timeHi.data = dynamicsHi->data;
     
    double omega  = 0.0;  
    double magR;
    double time1, time2, omegaDeriv, omegaDerivMid, tPeakOmega;
  
    if(debugPK) {
        out = fopen( "omegaHi.dat", "w" );
        printf("length of values = %d, retLenHi = %d\n", values->length, retLenHi);
        fflush(NULL);
    }
    if(debugRD) {
        out = fopen( "omegaHi.dat", "w" );
    }
  
    for ( i = 0; i < retLenHi; i++ )
    {
        for ( j = 0; j < values->length; j++ )
            { values->data[j] = *(dynamicsHi->data+(j+1)*retLenHi+i); }
    
        /* Calculate dr/dt */
        memset( dvalues->data, 0, numdynvars*sizeof(dvalues->data[0]));
        if( XLALSpinHcapRvecDerivative( 0, values->data, dvalues->data, 
            &seobParams) != XLAL_SUCCESS )
        {
                printf(
                    " Calculation of dr/dt failed while computing omegaHi time series\n");
                XLAL_ERROR( XLAL_EFUNC );
        }
    
        /* Calculare r x dr/dt */
        for (j=0; j<3; j++){
            rvec[j] = values->data[j];
            rdotvec[j] = dvalues->data[j];
        }
    
        //memcpy(rdotvec, dvalues->data, 3*sizeof(REAL8));
        //rvec[0] = posVecxHi.data[i]; rvec[1] = posVecyHi.data[i]; 
        //rvec[2] = posVeczHi.data[i];
        cross_product( rvec, rdotvec, rcrossrdot );        
   
        /* Calculate omega = |r x dr/dt| / r*r */
        magR = sqrt(inner_product(rvec, rvec));
        omegaHi->data[i] = sqrt(inner_product(rcrossrdot, rcrossrdot)) / (magR*magR); 
   
        if(debugPK || debugRD){
            fprintf( out, "%.16e\t%.16e\n", timeHi.data[i], omegaHi->data[i]);
        }
    }
  
    // Searching for crude omega_max (extremum)
    peakIdx = 0;
    *found = 0;
    for ( i = 1, peakIdx = 0; i < retLenHi-1; i++ ){
        omega = omegaHi->data[i];
        if (omega >= omegaHi->data[i-1] && omega > omegaHi->data[i+1]){
            peakIdx = i;
            *found = 1;
            if (debugPK){
                printf("PK: Crude peak of Omega is at idx = %d. t = %f,  OmegaPeak = %.16e\n", 
                    peakIdx, timeHi.data[peakIdx], omega);
                fflush(NULL);
            }
        }  
    }
    
    if(debugPK) {
        fclose(out);
        if (peakIdx ==0){
            printf("Stas: peak of orbital frequency was not found. peakIdx = %d, retLenHi = %d, i at exit = %d\n", peakIdx, retLenHi, i);
            fflush(NULL);
        }
    }
    if(debugRD) {
        fclose(out);
    }
  
    // refining the omega_max search (if it is found)
    tPeakOmega = 0.0;
    if(peakIdx != 0){
        spline = gsl_spline_alloc( gsl_interp_cspline, retLenHi );
        acc    = gsl_interp_accel_alloc();
        time1 = timeHi.data[peakIdx-2];
        gsl_spline_init( spline, timeHi.data, omegaHi->data, retLenHi );
        omegaDeriv = gsl_spline_eval_deriv( spline, time1, acc );
   
        if ( omegaDeriv > 0. ) { time2 = timeHi.data[peakIdx+2]; }
        else{
            time2 = time1;
            peakIdx = peakIdx-2;
	        time1 = timeHi.data[peakIdx-2];	      
	        omegaDeriv = gsl_spline_eval_deriv( spline, time1, acc );
        }
   
        do
        {
            tPeakOmega = ( time1 + time2 ) / 2.;
	        omegaDerivMid = gsl_spline_eval_deriv( spline, tPeakOmega, acc );
	   
	        if ( omegaDerivMid * omegaDeriv < 0.0 ) { time2 = tPeakOmega; }
	        else
	        {
		        omegaDeriv = omegaDerivMid;
		        time1 = tPeakOmega;
		    }
            if (debugPK){
                printf("Stas: searching for orbital max: %f, %f, %f, %f \n", time1, time2, omegaDeriv, omegaDerivMid);
            }
        } while ( time2 - time1 > 1.0e-5 );
        if(debugPK) {
          printf( "Estimation of the orbital peak is now at time %.16e \n", tPeakOmega);
          fflush(NULL);
        }
    }
  
    if(*found == 0 || debugRD){
        if(debugPK){
            printf("Stas: We couldn't find the maximum of orbital frequency, search for maximum of A(r)/r^2 \n");
        }   
        REAL8 rad, rad2, m1PlusetaKK, bulk, logTerms, deltaU, u, u2, u3, u4, u5;
        REAL8 listAOverr2[retLenHi];
        REAL8 Aoverr2;
        REAL8Vector *sigmaStar = NULL;
        REAL8Vector *sigmaKerr = NULL;
        if ( !(sigmaStar = XLALCreateREAL8Vector( 3 )) )
        {
          XLALDestroyREAL8Vector( sigmaStar );
          XLAL_ERROR( XLAL_ENOMEM );
        }
        if ( !(sigmaKerr = XLALCreateREAL8Vector( 3 )) )
        {
          XLALDestroyREAL8Vector( sigmaStar );
          XLAL_ERROR( XLAL_ENOMEM );
        }
        REAL8Vector s1Vec, s2Vec;
        s1Vec.length = s2Vec.length = 3;
        REAL8 s1Data[3], s2Data[3];
        REAL8 mTotal = m1 + m2;
        REAL8 a;
        REAL8 eta = m1*m2/(mTotal*mTotal);
        
        if(debugPK || debugRD){ 
            out = fopen( "OutAofR.dat", "w" );
        }
        for ( i = 0; i < retLenHi; i++ )
        {
            for ( j = 0; j < values->length; j++ )
            {
                values->data[j] = *(dynamicsHi->data+(j+1)*retLenHi+i);
            }
            for( j = 0; j < 3; j++ )
            {
                //s1DataNorm[k] = values->data[k+6];
                //s2DataNorm[k] = values->data[k+9];
                s1Data[j] = values->data[j+6] * mTotal * mTotal;
                s2Data[j] = values->data[j+9] * mTotal * mTotal;
            }
            s1Vec.data = s1Data;
            s2Vec.data = s2Data;
            XLALSimIMRSpinEOBCalculateSigmaStar( sigmaStar, m1, m2, &s1Vec, &s2Vec );
            XLALSimIMRSpinEOBCalculateSigmaKerr( sigmaKerr, m1, m2, &s1Vec, &s2Vec );
            
            seobParams.a = a = sqrt(inner_product(sigmaKerr->data, sigmaKerr->data));
            m1PlusetaKK = -1. + eta * seobCoeffs.KK;
            rad2 =  values->data[0]*values->data[0] + values->data[1]*values->data[1] + values->data[2]*values->data[2];
            rad = sqrt(rad2);
            u = 1./rad;
            u2 = u*u;
            u3 = u2*u;
            u4 = u2*u2;
            u5 = u4*u;
            bulk = 1./(m1PlusetaKK*m1PlusetaKK) + (2.*u)/m1PlusetaKK + a*a*u2;
            logTerms = 1. + eta*seobCoeffs.k0 + eta*log(1. + seobCoeffs.k1*u + seobCoeffs.k2*u2 + seobCoeffs.k3*u3 + seobCoeffs.k4*u4 + seobCoeffs.k5*u5 + seobCoeffs.k5l*u5*log(u));
            deltaU = bulk*logTerms;
            listAOverr2[i] = deltaU / rad2;
            if(debugPK || debugRD){
                fprintf(out, "%3.10f %3.10f\n", timeHi.data[i], listAOverr2[i]);
            }
            
        }
        if(debugPK || debugRD) fclose(out);
        if (*found == 0){
            // searching formaximum of A(r)/r^2
            peakIdx = 0;
            *found = 0;
            for ( i = 1, peakIdx = 0; i < retLenHi-1; i++ ){
                Aoverr2 = listAOverr2[i];
                if (Aoverr2 >= listAOverr2[i-1] && Aoverr2 > listAOverr2[i+1]){
                    peakIdx = i;
                    tPeakOmega = timeHi.data[i];
                    *found = 1;
                    if (debugPK){
                        printf("PK: Peak of A(r)/r^2 is at idx = %d. t = %f, Peak ampl. = %.16e\n", 
                            peakIdx, timeHi.data[peakIdx], Aoverr2);
                        fflush(NULL);
                    }
                }  
            }
        }
    
        if(debugPK) {
            if (peakIdx ==0){
                printf("Stas: peak of A(r)/r^2 was not found. \
                    peakIdx = %d, retLenHi = %d, i at exit = %d\n", peakIdx, retLenHi, i);
                fflush(NULL);
            }
        }
        XLALDestroyREAL8Vector(sigmaStar);
        XLALDestroyREAL8Vector(sigmaKerr);
    }
    if (spline != NULL)
        gsl_spline_free(spline);
    if (acc != NULL)
        gsl_interp_accel_free(acc);
    XLALDestroyREAL8Vector( values );
    XLALDestroyREAL8Vector( dvalues );
    XLALDestroyREAL8Vector( omegaHi );
    if (*found == 0){
        return(timeHi.data[retLenHi-1]);
    }
    else{
        return(tPeakOmega);
    }
}
  
double XLALSimLocateAmplTime(
    REAL8Vector *timeHi, 
    COMPLEX16Vector *hP22,
    int *found)
{
    int debugPK = 0;
    int debugRD = 1;
    FILE *out = NULL; 
    gsl_spline    *spline = NULL;
    gsl_interp_accel *acc = NULL;
    if (debugPK) {debugRD = 0;}
    
    // First we search for the maximum (extremum) of amplitude
    unsigned int i, peakIdx; 
    double minoff = 0.0;
    double maxoff = 30.0;
    unsigned int Nps = timeHi->length; 
    // this definesthe search interval for maximum (we might use min0ff= 0.051 instead)
    double tMax = timeHi->data[Nps-1] - minoff;
    double tMin = timeHi->data[Nps-1] - maxoff;
    double AmplN, AmplO;
    double tAmpMax, AmpMax, tAmp;
    tAmpMax = 0.;
    REAL8 Ampl[Nps];
    
    if(debugPK || debugRD) {
            out = fopen( "AmpPHi.dat", "w" );
    }
    AmplO = sqrt(creal(hP22->data[0])*creal(hP22->data[0]) + cimag(hP22->data[0])*cimag(hP22->data[0]));
    Ampl[0] = AmplO;
    peakIdx = 0;
    for (i=0; i<Nps-1; i++){
        AmplN = sqrt(creal(hP22->data[i+1])*creal(hP22->data[i+1]) + cimag(hP22->data[i+1])*cimag(hP22->data[i+1]));
        //Ampl = sqrt(hreP22->data[i]*hreP22->data[i] + himP22->data[i]*himP22->data[i]);
        if(debugPK || debugRD){
            fprintf(out, "%3.10f %3.10f\n", timeHi->data[i], Ampl[i]);
        }
        if (Ampl[i] >= AmplO && Ampl[i] >AmplN){
            tAmp = timeHi->data[i]; 
            if (tAmp >=tMin && tAmp <= tMax ){
                *found = 1;
                tAmpMax = tAmp;
                AmpMax = Ampl[i];
                peakIdx = i;
            }else{
                if (debugPK){
                    printf("Stas dismissing time %3.10f outside limits %3.10f, %3.10f \n", 
                        tAmp, tMin, tMax);
                }
            }        
        }
        AmplO = Ampl[i];
        Ampl[i+1] = AmplN;                
    }
    
    if (debugPK) 
    {
        fclose(out);
        if (*found ==0){
            printf("Stas: peak of 2,2 mode in P-frame was not found. peakIdx = %d, retLenHi = %d, i at exit = %d\n", peakIdx, Nps, i);
            fflush(NULL);
        }else{
            printf("Stas: we have found maximum of amplitude %3.10f at t = %3.10f \n", AmpMax, tAmpMax);
        }
    }
    if (debugRD) 
    {
        fclose(out);
    }

    if (*found ==0 || debugRD){
        // we haven't found the maximum of amplitude -> search for minimum of derivative (extremum)
        spline = gsl_spline_alloc( gsl_interp_cspline, Nps );
        acc    = gsl_interp_accel_alloc();
        gsl_spline_init( spline, timeHi->data, Ampl, Nps );
        
        double AmplDerivO = gsl_spline_eval_deriv(spline, timeHi->data[1], acc);
        double AmplDeriv = AmplDerivO;
        double AmplDerivN;
        if(debugPK || debugRD) {
                out = fopen( "DotAmpPHi.dat", "w" );
        }
        printf("stas wrinting the file\n");
        for (i=1; i<Nps-2; i++){
            if(debugPK || debugRD){
                fprintf(out, "%3.10f %3.10f\n", timeHi->data[i], AmplDeriv);
            }
            AmplDerivN = gsl_spline_eval_deriv( spline, timeHi->data[i+1], acc );
            if (*found == 0){
                if (AmplDeriv <= AmplDerivO && AmplDeriv < AmplDerivN){
                    tAmp = timeHi->data[i];
                    if (tAmp >=tMin && tAmp <= tMax  && *found==0){
                        *found = 1;
                        tAmpMax = tAmp;
                        AmpMax = AmplDeriv;
                        peakIdx = i;
                        //break;
                    }else{
                        if (debugPK){
                            printf("Stas dismissing time %3.10f outside limits %3.10f, %3.10f \n", 
                                tAmp, tMin, tMax);
                        }
                    }                              
                }
            }
            AmplDerivO = AmplDeriv;
            AmplDeriv = AmplDerivN;
        }
        
        if (debugPK) 
        {
            fclose(out);
            if (*found ==0){
                printf("Stas: peak of 2,2 mode in P-frame was not found. peakIdx = %d, retLenHi = %d, i at exit = %d\n", peakIdx, Nps, i);
                fflush(NULL);
            }else{
                printf("Stas: we have found maximum of amplitude %3.10f at t = %3.10f \n", AmpMax, tAmpMax);
            }
        }
         if (debugRD) 
        {
            fclose(out);
        }

        
    }
    if (spline != NULL)
        gsl_spline_free(spline);
    if (acc != NULL)
        gsl_interp_accel_free(acc);
    if (*found == 0){
        return(timeHi->data[Nps-1]);
    }
    else{
        return(tAmpMax);
    }
    
    
    
   
    
}    
