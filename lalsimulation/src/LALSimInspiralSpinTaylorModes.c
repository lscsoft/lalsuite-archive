
/*
 * Copyright (C) 2013 R. O'Shaughnessy
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

#include <complex.h>
#include <math.h>

#include <gsl/gsl_const.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_odeiv.h>

#include <lal/LALSimInspiral.h>
#include <lal/LALSimInspiralPrecess.h>
#include <lal/LALConstants.h>
#include <lal/LALStdlib.h>
#include <lal/TimeSeries.h>
#include <lal/FrequencySeries.h>
#include <lal/Units.h>
#include <lal/SphericalHarmonics.h>

/* Appends the start and end time series together, skipping the redundant first
 * sample of end.  Frees end before returning a pointer to the result, which is
 * the resized start series.  */
static REAL8TimeSeries *appendTSandFree(REAL8TimeSeries *start, 
        REAL8TimeSeries *end) {
    unsigned int origlen = start->data->length;
    start = XLALResizeREAL8TimeSeries(start, 0, 
            start->data->length + end->data->length - 1);
    
    memcpy(start->data->data + origlen, end->data->data+1, 
            (end->data->length-1)*sizeof(REAL8));

    XLALGPSAdd(&(start->epoch), -end->deltaT*(end->data->length - 1));

    XLALDestroyREAL8TimeSeries(end);

    return start;        
}

SphHarmTimeSeries * XLALSimInspiralSpinTaylorCorotatingModesDriver(
	REAL8 phiRef,                   /**< orbital phase at reference pt. */
	REAL8 v0,                       /**< tail gauge term (default = 1) */
	REAL8 deltaT,                   /**< sampling interval (s) */
	REAL8 m1,                       /**< mass of companion 1 (kg) */
	REAL8 m2,                       /**< mass of companion 2 (kg) */
	REAL8 fStart,                   /**< start GW frequency (Hz) */
	REAL8 fRef,                     /**< reference GW frequency (Hz) */
	REAL8 r,                        /**< distance of source (m) */
	REAL8 s1x,                      /**< initial value of S1x */
	REAL8 s1y,                      /**< initial value of S1y */
	REAL8 s1z,                      /**< initial value of S1z */
	REAL8 s2x,                      /**< initial value of S2x */
	REAL8 s2y,                      /**< initial value of S2y */
	REAL8 s2z,                      /**< initial value of S2z */
	REAL8 lnhatx,                   /**< initial value of LNhatx */
	REAL8 lnhaty,                   /**< initial value of LNhaty */
	REAL8 lnhatz,                   /**< initial value of LNhatz */
	REAL8 e1x,                      /**< initial value of E1x */
	REAL8 e1y,                      /**< initial value of E1y */
	REAL8 e1z,                      /**< initial value of E1z */
	REAL8 lambda1,                  /**< (tidal deformability of mass 1) / (mass of body 1)^5 (dimensionless) */
	REAL8 lambda2,                  /**< (tidal deformability of mass 2) / (mass of body 2)^5 (dimensionless) */
	REAL8 quadparam1,               /**< phenom. parameter describing induced quad. moment of body 1 (=1 for BHs, ~2-12 for NSs) */
	REAL8 quadparam2,               /**< phenom. parameter describing induced quad. moment of body 2 (=1 for BHs, ~2-12 for NSs) */
	LALSimInspiralSpinOrder spinO,  /**< twice PN order of spin effects */
	LALSimInspiralTidalOrder tideO, /**< twice PN order of tidal effects */
	int phaseO,                     /**< twice PN phase order */
	int amplitudeO,                 /**< twice PN amplitude order */
                                                          Approximant approx,              /**< PN approximant (SpinTaylorT1/T2/T4) */
                                                       int lmax
	)
{
	SphHarmTimeSeries *hlm = NULL;
    REAL8TimeSeries *V, *Phi, *S1x, *S1y, *S1z, *S2x, *S2y, *S2z;
    REAL8TimeSeries *LNhatx, *LNhaty, *LNhatz, *E1x, *E1y, *E1z;
    int status, n;
    unsigned int i;
    REAL8 fS, fE, phiShift;
    /* The Schwarzschild ISCO frequency - for sanity checking fRef */
    REAL8 fISCO = pow(LAL_C_SI,3) / (pow(6.,3./2.)*LAL_PI*(m1+m2)*LAL_G_SI);

    /* Sanity check fRef value */
    if( fRef < 0. )
    {
        XLALPrintError("XLAL Error - %s: fRef = %f must be >= 0\n", 
                __func__, fRef);
        XLAL_ERROR(XLAL_EINVAL);
    }
    if( fRef != 0. && fRef < fStart )
    {
        XLALPrintError("XLAL Error - %s: fRef = %f must be > fStart = %f\n", 
                __func__, fRef, fStart);
        XLAL_ERROR(XLAL_EINVAL);
    }
    if( fRef >= fISCO )
    {
        XLALPrintError("XLAL Error - %s: fRef = %f must be < Schwar. ISCO=%f\n",
                __func__, fRef, fISCO);
        XLAL_ERROR(XLAL_EINVAL);
    }

    /* if fRef=0, just integrate from start to end. Let phiRef=phiC */
    if( fRef < LAL_REAL4_EPS )
    {
        fS = fStart;
        fE = 0.;
        /* Evolve the dynamical variables */
        n = XLALSimInspiralSpinTaylorPNEvolveOrbit(&V, &Phi,
                &S1x, &S1y, &S1z, &S2x, &S2y, &S2z, 
                &LNhatx, &LNhaty, &LNhatz, &E1x, &E1y, &E1z, 
                deltaT, m1, m2, fS, fE, s1x, s1y, s1z, s2x, s2y, s2z, 
                lnhatx, lnhaty, lnhatz, e1x, e1y, e1z, lambda1, lambda2,
                quadparam1, quadparam2, spinO, tideO, phaseO, approx);
        if( n < 0 )
            XLAL_ERROR(XLAL_EFUNC);

        /* Apply phase shift so orbital phase ends with desired value */
        phiShift = phiRef - Phi->data->data[Phi->data->length-1];
        for( i=0; i < Phi->data->length; i++)
        {
            Phi->data->data[i] += phiShift;
        }
    }
    /* if fRef=fStart, just integrate from start to end. Let phiRef=phiStart */
    else if( abs(fRef - fStart) < LAL_REAL4_EPS )
    {
        fS = fStart;
        fE = 0.;
        /* Evolve the dynamical variables */
        n = XLALSimInspiralSpinTaylorPNEvolveOrbit(&V, &Phi,
                &S1x, &S1y, &S1z, &S2x, &S2y, &S2z, 
                &LNhatx, &LNhaty, &LNhatz, &E1x, &E1y, &E1z, 
                deltaT, m1, m2, fS, fE, s1x, s1y, s1z, s2x, s2y, s2z, 
                lnhatx, lnhaty, lnhatz, e1x, e1y, e1z, lambda1, lambda2,
                quadparam1, quadparam2, spinO, tideO, phaseO, approx);
        if( n < 0 )
            XLAL_ERROR(XLAL_EFUNC);

        /* Apply phase shift so orbital phase starts with desired value */
        phiShift = phiRef - Phi->data->data[0];
        for( i=0; i < Phi->data->length; i++)
        {
            Phi->data->data[i] += phiShift;
        }
    }
    else /* Start in middle, integrate backward and forward, stitch together */
    {
        REAL8TimeSeries *V1=NULL, *Phi1=NULL, *S1x1=NULL, *S1y1=NULL, *S1z1=NULL, *S2x1=NULL, *S2y1=NULL, *S2z1=NULL;
        REAL8TimeSeries *LNhatx1=NULL, *LNhaty1=NULL, *LNhatz1=NULL, *E1x1=NULL, *E1y1=NULL, *E1z1=NULL;
        REAL8TimeSeries *V2=NULL, *Phi2=NULL, *S1x2=NULL, *S1y2=NULL, *S1z2=NULL, *S2x2=NULL, *S2y2=NULL, *S2z2=NULL;
        REAL8TimeSeries *LNhatx2=NULL, *LNhaty2=NULL, *LNhatz2=NULL, *E1x2=NULL, *E1y2=NULL, *E1z2=NULL;

        /* Integrate backward to fStart */
        fS = fRef;
        fE = fStart;
        n = XLALSimInspiralSpinTaylorPNEvolveOrbit(&V1, &Phi1,
                &S1x1, &S1y1, &S1z1, &S2x1, &S2y1, &S2z1, 
                &LNhatx1, &LNhaty1, &LNhatz1, &E1x1, &E1y1, &E1z1,
                deltaT, m1, m2, fS, fE, s1x, s1y, s1z, s2x, s2y,
                s2z, lnhatx, lnhaty, lnhatz, e1x, e1y, e1z, lambda1, lambda2,
                quadparam1, quadparam2, spinO, tideO, phaseO, approx);
        
        /* Apply phase shift so orbital phase has desired value at fRef */
        phiShift = phiRef - Phi1->data->data[Phi1->data->length-1];
        for( i=0; i < Phi1->data->length; i++)
        {
            Phi1->data->data[i] += phiShift;
        }

        /* Integrate forward to end of waveform */
        fS = fRef;
        fE = 0.;
        n = XLALSimInspiralSpinTaylorPNEvolveOrbit(&V2, &Phi2,
                &S1x2, &S1y2, &S1z2, &S2x2, &S2y2, &S2z2, 
                &LNhatx2, &LNhaty2, &LNhatz2, &E1x2, &E1y2, &E1z2,
                deltaT, m1, m2, fS, fE, s1x, s1y, s1z, s2x, s2y,
                s2z, lnhatx, lnhaty, lnhatz, e1x, e1y, e1z, lambda1, lambda2,
                quadparam1, quadparam2, spinO, tideO, phaseO, approx);
        
        /* Apply phase shift so orbital phase has desired value at fRef */
        phiShift = phiRef - Phi2->data->data[0];
        for( i=0; i < Phi2->data->length; i++)
        {
            Phi2->data->data[i] += phiShift;
        }

        /* Stitch 2nd set of vectors onto 1st set. Free 2nd set. */
        V = appendTSandFree(V1, V2); 
        Phi = appendTSandFree(Phi1, Phi2);
        S1x = appendTSandFree(S1x1, S1x2);
        S1y = appendTSandFree(S1y1, S1y2);
        S1z = appendTSandFree(S1z1, S1z2);
        S2x = appendTSandFree(S2x1, S2x2);
        S2y = appendTSandFree(S2y1, S2y2);
        S2z = appendTSandFree(S2z1, S2z2);
        LNhatx = appendTSandFree(LNhatx1, LNhatx2);
        LNhaty = appendTSandFree(LNhaty1, LNhaty2);
        LNhatz = appendTSandFree(LNhatz1, LNhatz2);
        E1x = appendTSandFree(E1x1, E1x2);
        E1y = appendTSandFree(E1y1, E1y2);
        E1z = appendTSandFree(E1z1, E1z2);
    }


    /* Use the dynamical variables to build the polarizations */
    /* Currently uses NO SPIN AMPLITUDE CORRECTIONS*/
    int m, l;
    COMPLEX16TimeSeries *hxx;
    for(l=2; l<=lmax; l++){
        for(m=-l; m<=l; m++){
            hxx = XLALCreateSimInspiralPNModeCOMPLEX16TimeSeries(V, Phi,
                v0, m1, m2, r, amplitudeO, l, m);
            if ( !hxx ){
                XLAL_ERROR_NULL(XLAL_EFUNC);
            }
            XLALSphHarmTimeSeriesAddMode(hxx, l, m, hlm);
            XLALDestroyCOMPLEX16TimeSeries(hxx);
        }
    }


    /* Destroy vectors of dynamical variables, check for errors then exit */
    XLALDestroyREAL8TimeSeries(V);
    XLALDestroyREAL8TimeSeries(Phi);
    XLALDestroyREAL8TimeSeries(S1x);
    XLALDestroyREAL8TimeSeries(S1y);
    XLALDestroyREAL8TimeSeries(S1z);
    XLALDestroyREAL8TimeSeries(S2x);
    XLALDestroyREAL8TimeSeries(S2y);
    XLALDestroyREAL8TimeSeries(S2z);
    XLALDestroyREAL8TimeSeries(LNhatx);
    XLALDestroyREAL8TimeSeries(LNhaty);
    XLALDestroyREAL8TimeSeries(LNhatz);
    XLALDestroyREAL8TimeSeries(E1x);
    XLALDestroyREAL8TimeSeries(E1y);
    XLALDestroyREAL8TimeSeries(E1z);
    if( status < 0 )
        XLAL_ERROR(XLAL_EFUNC);

    return hlm;
}
