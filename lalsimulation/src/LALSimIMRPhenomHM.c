/*
 * Copyright (C) 2017 Sebastian Khan, Francesco Pannarale
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

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif


#include <math.h>
#include <stdio.h>

#include <lal/LALSimIMR.h>
#include <lal/SphericalHarmonics.h>

/* This allows us to reuse internal IMRPhenomD functions without making those functions XLAL */
#include "LALSimIMRPhenomD_internals.c"
#include "LALSimRingdownCW.c"


#include "LALSimIMRPhenomHM.h"

#ifndef _OPENMP
#define omp ignore
#endif

//FP: make PhenomDQuantities, etc. EXTERNAL?

/* List of modelled modes */
/* NOTE: l=2, m=2 should be the first mode in the list called ModeArray. Or the code will break */
static int NMODES = NMODES_MAX;
// static const int ModeArray[NMODES_MAX][2] = { {2,2}, {2,1}, {3,3}, {4,4}, {5,5} };
static const int ModeArray[NMODES_MAX][2] = { {2,2}, {2,1}, {3,3}, {3,2}, {4,4}, {4,3} };
// static const int ModeArray[NMODES_MAX][2] = { {2,2}, {2,1}, {3,3}, {3,2}, {4,4} };
// static const int ModeArray[NMODES_MAX][2] = { {2,2}, {2,1}, {3,3}, {4,4}, {4,3} };
// static const int ModeArray[NMODES_MAX][2] = { {2,2}, {2,1} };
// static const int ModeArray[NMODES_MAX][2] = { {2,2} };

/* List of phase shifts: the index is the azimuthal number m */
static const double cShift[7] = {0.0,
                                 LAL_PI_2 /* i shift */,
                                 0.0,
                                 -LAL_PI_2 /* -i shift */,
                                 LAL_PI /* 1 shift */,
                                 LAL_PI_2 /* -1 shift */,
                                 0.0};

/* Dimensionless frequency of last data point in waveform */
#define Mf_CUT_HM 0.5
// #define Mf_CUT_HM 0.285
// #define Mf_CUT_HM 0.2
/* Activates amplitude part of the model */
#define AmpFlagTrue 1
#define AmpFlagFalse 0
#define MfAtScale_wf_amp 0.0001
/* Scale factor multiplying the ringdown frequency for mode (l,m).
* Used to estimate the end frequency for each mode.
*/
#define F_FACTOR 5. /*NOTE: this can cause little jumps in the strain because of sudden termination of some modes */


/* START: newer functions  */

int init_useful_mf_powers(UsefulMfPowers *p, REAL8 number)
{
    XLAL_CHECK(0 != p, XLAL_EFAULT, "p is NULL");
    XLAL_CHECK(!(number < 0) , XLAL_EDOM, "number must be non-negative");

    // consider changing pow(x,1/6.0) to cbrt(x) and sqrt(x) - might be faster
    p->itself = number;
    p->sqrt = sqrt(number);
    //p->sixth = pow(number, 1/6.0); //FP
    p->sixth = cbrt(p->sqrt);
    p->m_sixth = 1.0/p->sixth;
    p->third = p->sixth * p->sixth;
    p->two_thirds = p->third * p->third;
    p->four_thirds = number * p->third;
    p->five_thirds = number * p->two_thirds;
    p->two = number * number;
    p->seven_thirds = p->third * p->two;
    p->eight_thirds = p->two_thirds * p->two;
    p->m_seven_sixths = p->m_sixth/number;
    p->m_five_sixths = p->m_seven_sixths*p->third;
    p->m_sqrt = 1./p->sqrt;

    return XLAL_SUCCESS;
}

/* TODO: at this point it is probably unnecessary to have UsefulPowers and
 * UsefulMfPowers: one struct should suffice */

/* Given a full UsefulMfPowers variable, extract a UsefulPowers one from it */
int downsize_useful_mf_powers(UsefulPowers *out, UsefulMfPowers *in);
int downsize_useful_mf_powers(UsefulPowers *out, UsefulMfPowers *in)
{
	XLAL_CHECK(0 != in, XLAL_EFAULT, "in is NULL");

	out->third          = in->third;
	out->two_thirds     = in->two_thirds;
	out->four_thirds    = in->four_thirds;
	out->five_thirds    = in->five_thirds;
	out->two            = in->two;
	out->seven_thirds   = in->seven_thirds;
	out->eight_thirds   = in->eight_thirds;
    out->inv            = 1./in->itself;
    out->m_seven_sixths = in->m_seven_sixths;
    out->m_third        = 1./in->third;
    out->m_two_thirds   = out->m_third*out->m_third;
    out->m_five_thirds  = out->inv*out->m_two_thirds;

	return XLAL_SUCCESS;
}

/* Precompute a bunch of PhenomD related quantities and store them filling in a PhenomDStorage variable */
int init_PhenomD_Storage(PhenomDStorage* p, const REAL8 m1, const REAL8 m2, const REAL8 chi1z, const REAL8 chi2z)
{
    XLAL_CHECK(0 != p, XLAL_EFAULT, "p is NULL");

    p->m1 = m1; /* Inherits units from m1 in function arguments */
    p->m2 = m2; /* Inherits units from m2 in function arguments */
    p->Mtot = m1+m2; /* Inherits units from m1 and m2 in function arguments */
    p->eta = m1*m2/(p->Mtot*p->Mtot);
    p->Inv1MinusEradRational0815 = 1.0/(1.0-EradRational0815(p->eta, chi1z, chi2z));
    p->finspin = XLALSimIMRPhenomDFinalSpin(m1, m2, chi1z, chi2z); /* dimensionless final spin */
    if (p->finspin > 1.0) XLAL_ERROR(XLAL_EDOM, "PhenomD fring function: final spin > 1.0 not supported\n");

    /* Precompute Lionel's QNM higher mode fits (both real and imaginary parts of
    the complex ringdown frequency).
    Annoying, but faster than looking up ell's and mm's.
    WARNING: bug prone, as one may access PhenomHMfring and PhenomHMfdamp
    for ell's and mm's not in ModeArray. */
    complex double ZZ;
    const REAL8 inv2Pi = 0.5/LAL_PI;
    for( int j=0; j<NMODES; j++ ){
        int ell = ModeArray[j][0];
        int mm = ModeArray[j][1];
        ZZ = CW07102016( KAPPA(p->finspin, ell, mm), ell, mm, 0 );
        /* lm mode ringdown frequency (real part of ringdown), geometric units */
        const REAL8 Mf_RD = inv2Pi * creal(ZZ); /* GW ringdown frequency, converted from angular frequency */
        p->PhenomHMfring[ell][mm] = Mf_RD * p->Inv1MinusEradRational0815; /* scale by predicted final mass */
        /* lm mode ringdown damping time (imaginary part of ringdown), geometric units */
        const REAL8 f_damp = inv2Pi * cimag(ZZ); /* this is the 1./tau in the complex QNM */
        p->PhenomHMfdamp[ell][mm] = f_damp * p->Inv1MinusEradRational0815; /* scale by predicted final mass */
    }
    p->Mf_RD_22 = p->PhenomHMfring[2][2];
    p->Mf_DM_22 = p->PhenomHMfdamp[2][2];

    for( int j=0; j<NMODES; j++ ){
        int ell = ModeArray[j][0];
        int mm = ModeArray[j][1];
        p->Rholm[ell][mm] = p->Mf_RD_22/p->PhenomHMfring[ell][mm];
        p->Taulm[ell][mm] = p->PhenomHMfdamp[ell][mm]/p->Mf_DM_22;
        // printf("p->Rholm[%i][%i] = %f, p->Taulm[%i][%i] = %f\n",ell, mm, p->Rholm[ell][mm],ell, mm, p->Taulm[ell][mm]);
    }

    /* A bunch of useful powers used in XLALSimIMRPhenomHMPNAmplitudeLeadingOrder */
    /* pow_Mf_wf_prefactor are the coefficients from the leading order amplitude terms for each l,m mode we consider.  */
    /* If A_lm(f) = alpha_lm * f^klm then pow_Mf_wf_prefactor = alpha_lm */
    /* note The factors of pi's normally multiplying with the frequency are factored into these terms. */
    /* note delta == sqrt(1. - 4.*eta) */
    /* note delta2 == 1. - 4.*eta */
    REAL8 sqrteta = sqrt(p->eta);
    REAL8 Seta = sqrt( 1.0 - 4.0 * p->eta );
    REAL8 delta = Seta;
    REAL8 delta2 = 1.0 - 4.0 * p->eta ;
    REAL8 Blm_ans = 0.;
    for( int j=0; j<NMODES; j++ ){
        int ell = ModeArray[j][0];
        int mm = ModeArray[j][1];

        if ( ell==2 ) {
            if (mm==2 ) {
                Blm_ans = 1.0;
            } else { // mm==1
                Blm_ans = delta * pow(LAL_PI, 1.0/3.0) /  3.0;
            }
        } else if ( ell==3 ) {
            if ( mm==3 ) {
                Blm_ans = (3.0/4.0) * sqrt(15.0/14.0) * pow(LAL_PI, 1.0/3.0) * delta;
            }
            else { // mm==2
                Blm_ans = sqrt(5.0/63.0) * pow(LAL_PI, 2.0/3.0) * (delta2 + p->eta);
            }
        } else if ( ell==4 ) {
            if ( mm==4 ) {
                Blm_ans = sqrt(320.0/567.0) * pow(LAL_PI, 2.0/3.0) * (delta2 + p->eta);
            } else { // mm==3
                Blm_ans = sqrt(81.0/1120.0) * LAL_PI * (delta2 + 2*p->eta) * delta;
            }
        } else if ( ell==5 ) {
            if ( mm==5 ) {
                //CHECK ME
                Blm_ans = (625.0/96.0) * (1.0 / sqrt(66.0)) * LAL_PI * (delta2 + 2*p->eta) * delta;
            } else { // mm==4
            //NOT IMPLEMENTED
            Blm_ans = 0.0;
            }
        } else if ( ell==6 ) {
            if ( mm==6 ) {
                //CHECK ME
                Blm_ans = (9.0/5.0) * sqrt(6.0) * pow(LAL_PI, 1.0/6.0) / sqrteta;
            }
            else { // mm==5
                //NOT IMPLEMENTED
                Blm_ans = 0.0;
            }
        } else {
            //NOT IMPLEMENTED
            Blm_ans = 0.0;
        }
        p->Blm_prefactor[ell][mm] = Blm_ans;
    }

    return XLAL_SUCCESS;
}

//mathematica function postfRDflm
//double XLALIMRPhenomHMpostfRDflm(REAL8 Mf, REAL8 Mf_RD_22, REAL8 Mf_RD_lm, const INT4 AmpFlag, const INT4 ell, const INT4 mm, PhenomDStorage* PhenomDQuantities){
double XLALIMRPhenomHMTrd(REAL8 Mf, REAL8 Mf_RD_22, REAL8 Mf_RD_lm, const INT4 AmpFlag, const INT4 ell, const INT4 mm, PhenomDStorage* PhenomDQuantities);
double XLALIMRPhenomHMTrd(REAL8 Mf, REAL8 Mf_RD_22, REAL8 Mf_RD_lm, const INT4 AmpFlag, const INT4 ell, const INT4 mm, PhenomDStorage* PhenomDQuantities){
    double ans = 0.0;
    if ( AmpFlag==1 ) {
        /* For amplitude */
        ans = Mf-Mf_RD_lm+Mf_RD_22; /*Used for the Amplitude as an approx fix for post merger powerlaw slope */
    } else {
        /* For phase */
        REAL8 Rholm = PhenomDQuantities->Rholm[ell][mm];
        ans = Rholm * Mf;          /* Used for the Phase */
    }

    return ans;
}
// mathematica function Trd
// domain mapping function - ringdown
//UNUSED double XLALIMRPhenomHMTrd(REAL8 Mf, REAL8 Mf_RD_22, REAL8 Mf_RD_lm, const INT4 AmpFlag, const INT4 ell, const INT4 mm, PhenomDStorage* PhenomDQuantities){
//    return XLALIMRPhenomHMpostfRDflm(Mf, Mf_RD_22, Mf_RD_lm, AmpFlag, ell, mm, PhenomDQuantities);
//}

// mathematica function Ti
// domain mapping function - inspiral
double XLALIMRPhenomHMTi(REAL8 Mf, const INT4 mm);
double XLALIMRPhenomHMTi(REAL8 Mf, const INT4 mm){
    return 2.0 * Mf / mm;
}

void XLALIMRPhenomHMSlopeAmAndBm(double *Am, double *Bm, const INT4 mm, REAL8 fi, REAL8 fr, REAL8 Mf_RD_22, REAL8 Mf_RD_lm, const INT4 AmpFlag, const INT4 ell, PhenomDStorage* PhenomDQuantities);
void XLALIMRPhenomHMSlopeAmAndBm(double *Am, double *Bm, const INT4 mm, REAL8 fi, REAL8 fr, REAL8 Mf_RD_22, REAL8 Mf_RD_lm, const INT4 AmpFlag, const INT4 ell, PhenomDStorage* PhenomDQuantities){
    REAL8 Trd = XLALIMRPhenomHMTrd(fr, Mf_RD_22, Mf_RD_lm, AmpFlag, ell, mm, PhenomDQuantities);
    REAL8 Ti = XLALIMRPhenomHMTi(fi, mm);

    //Am = ( Trd[fr]-Ti[fi] )/( fr - fi );
    *Am = ( Trd-Ti )/( fr - fi );

    //Bm = Ti[fi] - fi*Am;
    *Bm = Ti - fi*(*Am);
}

int XLALIMRPhenomHMMapParams(REAL8 *a, REAL8 *b, REAL8 flm, REAL8 fi, REAL8 fr, REAL8 Ai, REAL8 Bi, REAL8 Am, REAL8 Bm, REAL8 Ar, REAL8 Br);
int XLALIMRPhenomHMMapParams(REAL8 *a, REAL8 *b, REAL8 flm, REAL8 fi, REAL8 fr, REAL8 Ai, REAL8 Bi, REAL8 Am, REAL8 Bm, REAL8 Ar, REAL8 Br){
    // Define function to output map params used depending on
    if ( flm > fi ){
        if ( flm > fr ){
            *a = Ar;
            *b = Br;
        } else {
            *a = Am;
            *b = Bm;
        }
    }
    else {
        *a = Ai;
        *b = Bi;
    };
    return XLAL_SUCCESS;
}

int XLALIMRPhenomHMFreqDomainMapParams(REAL8 *a, REAL8 *b, REAL8 *fi, REAL8 *fr, REAL8 *f1, const REAL8 flm, const INT4 ell, const INT4 mm, PhenomDStorage *PhenomDQuantities, const int AmpFlag);
int XLALIMRPhenomHMFreqDomainMapParams( REAL8 *a,/**< [Out]  */
                                        REAL8 *b,/**< [Out]  */
                                        REAL8 *fi,/**< [Out]  */
                                        REAL8 *fr,/**< [Out]  */
                                        REAL8 *f1,/**< [Out]  */
                                        const REAL8 flm, /**< input waveform frequency */
                                        const INT4 ell, /**< spherical harmonics ell mode */
                                        const INT4 mm, /**< spherical harmonics m mode */
                                        PhenomDStorage *PhenomDQuantities, /**< Stores quantities in order to calculate them only once */
                                        const int AmpFlag /**< is ==1 then computes for amplitude, if ==0 then computes for phase */)
{

    /*check output points are NULL*/
    XLAL_CHECK(a != NULL, XLAL_EFAULT);
    XLAL_CHECK(b != NULL, XLAL_EFAULT);
    XLAL_CHECK(fi != NULL, XLAL_EFAULT);
    XLAL_CHECK(fr != NULL, XLAL_EFAULT);
    XLAL_CHECK(f1 != NULL, XLAL_EFAULT);

    /* Account for different f1 definition between PhenomD Amplitude and Phase derivative models */
    REAL8 Mf_1_22  = 0.; /* initalise variable */
    if ( AmpFlag==1 ) {
        /* For amplitude */
        Mf_1_22  = AMP_fJoin_INS; /* inspiral joining frequency from PhenomD [amplitude model], for the 22 mode */
    } else {
        /* For phase */
        Mf_1_22  = PHI_fJoin_INS; /* inspiral joining frequency from PhenomD [phase model], for the 22 mode */
    }

    *f1 = Mf_1_22;

    REAL8 Mf_RD_22 = PhenomDQuantities->Mf_RD_22;
    REAL8 Mf_RD_lm = PhenomDQuantities->PhenomHMfring[ell][mm];

    // Define a ratio of QNM frequencies to be used for scaling various quantities
    REAL8 Rholm = PhenomDQuantities->Rholm[ell][mm];

    // Given experiments with the l!=m modes, it appears that the QNM scaling rather than the PN scaling may be optimal for mapping f1
    REAL8 Mf_1_lm = Mf_1_22 / Rholm;

    /* Define transition frequencies */
    *fi = Mf_1_lm;
    *fr = Mf_RD_lm;

    /*Define the slope and intercepts of the linear transformation used*/
    REAL8 Ai = 2.0/mm;
    REAL8 Bi = 0.0;
    REAL8 Am;
    REAL8 Bm;
    XLALIMRPhenomHMSlopeAmAndBm(&Am, &Bm, mm, *fi, *fr, Mf_RD_22, Mf_RD_lm, AmpFlag, ell, PhenomDQuantities);

    REAL8 Ar = 1.0;
    REAL8 Br = 0.0;
    if ( AmpFlag==1 ) {
        /* For amplitude */
        Br = -Mf_RD_lm+Mf_RD_22;
    } else {
        /* For phase */
        Ar = Rholm;
    }

    /* Define function to output map params used depending on */
    int ret = XLALIMRPhenomHMMapParams(a, b, flm, *fi, *fr, Ai, Bi, Am, Bm, Ar, Br);
    if (ret != XLAL_SUCCESS){
        XLALPrintError("XLAL Error - XLALIMRPhenomHMMapParams failed in XLALIMRPhenomHMMapParams (1)\n");
        XLAL_ERROR(XLAL_EDOM);
    }

    // *f2lm = 0.5 * ( *fr + *fi );

    return XLAL_SUCCESS;
}

/**
 * XLALSimIMRPhenomHMFreqDomainMap
 * Input waveform frequency in Geometric units (Mflm)
 * and computes what frequency this corresponds
 * to scaled to the 22 mode.
 */
double XLALSimIMRPhenomHMFreqDomainMap(REAL8 Mflm, const INT4 ell, const INT4 mm, PhenomDStorage* PhenomDQuantities, const int AmpFlag);
double XLALSimIMRPhenomHMFreqDomainMap(REAL8 Mflm,
                                        const INT4 ell,
                                        const INT4 mm,
                                        PhenomDStorage* PhenomDQuantities,
                                        const int AmpFlag)
{

    /* Mflm here has the same meaning as Mf_wf in XLALSimIMRPhenomHMFreqDomainMapHM (old deleted function). */
    REAL8 a = 0.;
    REAL8 b = 0.;
    /* Following variables not used in this funciton but are returned in XLALIMRPhenomHMFreqDomainMapParams */
    REAL8 fi = 0.;
    REAL8 fr = 0.;
    REAL8 f1 = 0.;
    int ret = XLALIMRPhenomHMFreqDomainMapParams(&a, &b, &fi, &fr, &f1, Mflm, ell, mm, PhenomDQuantities, AmpFlag);
    if (ret != XLAL_SUCCESS){
        XLALPrintError("XLAL Error - XLALIMRPhenomHMFreqDomainMapParams failed in XLALSimIMRPhenomHMFreqDomainMapParams\n");
        XLAL_ERROR(XLAL_EDOM);
    }
    REAL8 Mf22 = a * Mflm + b;
    return Mf22;
}

/* END: newer functions  */




/*
 * For a given frequency and ell and m spherical harmonic mode
 * return the frequency scaled to give the leading order PN
 * amplitude for the given ell and m modes.
 * Calcuated from mathematica function: FrequencyPower[f, {ell, m}] / FrequencyPower[f, {2, 2}]
 * FrequencyPower function just returns the leading order PN term in the amplitude.
 */
 /**
   * If A_lm(f) = alpha_lm * f^klm then this function
   * returns f^klm / f^k22
   */
double XLALSimIMRPhenomHMPNFrequencyScale( UsefulPowers *p, REAL8 Mf, INT4 ell, INT4 mm);
double XLALSimIMRPhenomHMPNFrequencyScale( UsefulPowers *p,
                                           REAL8 Mf,
                                           INT4 ell,
                                           INT4 mm ) {

    /* Initialise answer */
    REAL8 ans = 0.0;

    //FP: just compute the power in here
    if ( ell==2 ) {
        if ( mm==2 ) {
            ans = 1.0;
        }
        else { //mm==1
            ans = p->third;
        }
    } else if ( ell==3 ) {
        if ( mm==3 ) {
          ans = p->third;
        }
        else { //mm==2
          ans = p->two_thirds;
        }
    } else if ( ell==4 ) {
        if ( mm==4 ) {
          ans = p->two_thirds;
        }
        else { //mm==3
          ans = Mf;
        }
    } else if ( ell==5 ) {
        if ( mm==5 ) {
          ans = Mf;
        }
        else { //mm==4
          ans = p->four_thirds;
        }
    } else if ( ell==6 ) {
        if ( mm==6 ) {
          ans = p->four_thirds;
        }
        else { //mm==5
          ans = p->five_thirds;
        }
    } else {
        XLALPrintError("XLAL Error - requested ell = %i and m = %i mode not available, check documentation for available modes\n", ell, mm);
        XLAL_ERROR(XLAL_EDOM);
    }

    return ans;

}

/**
  * If A_lm(f) = alpha_lm * f^klm then this function
  * returns f^(klm) where klm is the appropriate exponent for th l,m mode
  */
double XLALSimIMRPhenomHMPNAmplitudeLeadingOrderFpow(INT4 ell, INT4 mm, REAL8 Mf);
double XLALSimIMRPhenomHMPNAmplitudeLeadingOrderFpow(INT4 ell, INT4 mm, REAL8 Mf) {
    /* Initialise answer */
    REAL8 ans = 0.0;

    UsefulMfPowers powers_of_Mf;
    int errcode = XLAL_SUCCESS;
    errcode = init_useful_mf_powers(&powers_of_Mf, Mf);
    XLAL_CHECK(errcode == XLAL_SUCCESS, errcode, "init_useful_mf_powers failed for Mf");

    //FP: some of these can be computed directly here rather than for each mm and ll
    if ( ell==2 ) {
        if ( mm==2 ) {
          ans = powers_of_Mf.m_seven_sixths;
        } else { //mm==1
          ans = powers_of_Mf.m_five_sixths;
        }
    } else if ( ell==3 ) {
        if ( mm==3 ) {
          ans = powers_of_Mf.m_five_sixths;
        }
        else { //mm==2
          ans = powers_of_Mf.m_sqrt;
        }
    } else if ( ell==4 ) {
        if ( mm==4 ) {
          ans = powers_of_Mf.m_sqrt;
        }
        else { //mm==3
          ans = powers_of_Mf.m_sixth;
        }
    } else if ( ell==5 ) {
        if ( mm==5 ) {
          ans = powers_of_Mf.m_sixth;
        }
        else { //mm==4
          ans = powers_of_Mf.sixth;
        }
    } else if ( ell==6 ) {
        if ( mm==6 ) {
          ans = powers_of_Mf.sixth;
        }
        else { //mm==5
          ans = powers_of_Mf.sqrt;
        }
    } else {
        XLALPrintError("XLAL Error - requested ell = %i and m = %i mode not available, check documentation for available modes\n", ell, mm);
        XLAL_ERROR(XLAL_EDOM);
    }

    return ans;
}

double XLALSimIMRPhenomHMAmplitude(double Mf_wf, int ell, int mm, IMRPhenomDAmplitudeCoefficients *pAmp, AmpInsPrefactors * amp_prefactors, PhenomDStorage * PhenomDQuantities);
double XLALSimIMRPhenomHMAmplitude( double Mf_wf,
                                    int ell,
                                    int mm,
                                    IMRPhenomDAmplitudeCoefficients *pAmp,
                                    AmpInsPrefactors * amp_prefactors,
                                    PhenomDStorage * PhenomDQuantities
                                  )
{
    double Mf_22 =  XLALSimIMRPhenomHMFreqDomainMap(Mf_wf, ell, mm, PhenomDQuantities, AmpFlagTrue);//FP PhenomHMfring[ell][mm], Rholm[ell][mm], Mf_RD_22

    UsefulPowers powers_of_Mf_22;
    int errcode = XLAL_SUCCESS;
    errcode = init_useful_powers(&powers_of_Mf_22, Mf_22);
    XLAL_CHECK(errcode == XLAL_SUCCESS, errcode, "init_useful_powers failed for Mf_22");

    double PhenDamp = IMRPhenDAmplitude(Mf_22, pAmp, &powers_of_Mf_22, amp_prefactors);

    /* coefficients of leading order PN amplitude */
    double Blm_prefac = PhenomDQuantities->Blm_prefactor[ell][mm];

    /* ratio of frequency term in leadering order PN */
    /* at the scaled frequency */
    double f_frac = XLALSimIMRPhenomHMPNFrequencyScale( &powers_of_Mf_22, Mf_22, ell, mm);

    double Blm = Blm_prefac * f_frac;

    /* (m/2)^klm NOTE in the paper this is (2/m)^(-klm) i.e. inverted. */
    double m_over_2_pow_klm = XLALSimIMRPhenomHMPNAmplitudeLeadingOrderFpow(ell, mm, mm/2.0);

    double betalm = Blm * m_over_2_pow_klm;

    double HMamp = PhenDamp * betalm;
    return HMamp;

}


/*TODO Should also probably add LALDict as a argument here.*/
int XLALSimIMRPhenomHMPhasePreComp(HMPhasePreComp *q, const INT4 ell, const INT4 mm, const REAL8 eta, const REAL8 chi1z, const REAL8 chi2z, PhenomDStorage *PhenomDQuantities, LALDict *extraParams);
int XLALSimIMRPhenomHMPhasePreComp(HMPhasePreComp *q, const INT4 ell, const INT4 mm, const REAL8 eta, const REAL8 chi1z, const REAL8 chi2z, PhenomDStorage *PhenomDQuantities, LALDict *extraParams)
{
    /*
    { ai,bi,fi,fr,f1,f2lm,fRD22,fRDlm,Ti,Trd,Tm } = FreqDomainMapParams[ 0.0001, mode, eta,chi1,chi2,False ];
    { a2lm,b2lm,fi,fr,f1,f2lm,fRD22,fRDlm,Ti,Trd,Tm } = FreqDomainMapParams[ f2lm+0.0001, mode, eta,chi1,chi2,False ];
    { ar,br,fi,fr,f1,f2lm,fRD22,fRDlm,Ti,Trd,Tm } = FreqDomainMapParams[ fr+0.0001, mode, eta,chi1,chi2,False ];

    PhDBconst = IMRPhenDPhasev2[ a2lm*fi   + b2lm, eta,chi1,chi2, mode ]/a2lm;
    PhDCconst = IMRPhenDPhasev2[ a2lm*f2lm + b2lm, eta,chi1,chi2, mode ]/a2lm;

    PhDBAterm = IMRPhenDPhasev2[ ai*fi + bi, eta,chi1,chi2, mode ]/ai;

    Return[{ai,bi,a2lm,b2lm,ar,br,fi,f2lm,fr,PhDBconst,PhDCconst,PhDBAterm}];
    */

    REAL8 ai = 0.0;
    REAL8 bi = 0.0;
    REAL8 am = 0.0;
    REAL8 bm = 0.0;
    REAL8 ar = 0.0;
    REAL8 br = 0.0;
    REAL8 fi = 0.0;
    REAL8 f1 = 0.0;
    REAL8 fr = 0.0;

    const INT4 AmpFlag = 0;

    /* NOTE: As long as Mfshit + f2lm isn't >= fr then the value of the shift is arbitrary. */
    const REAL8 Mfshift = 0.0001;

    int ret = XLALIMRPhenomHMFreqDomainMapParams(&ai, &bi, &fi, &fr, &f1, Mfshift, ell, mm, PhenomDQuantities, AmpFlag);
    if (ret != XLAL_SUCCESS){
        XLALPrintError("XLAL Error - XLALIMRPhenomHMFreqDomainMapParams failed in XLALIMRPhenomHMFreqDomainMapParams - inspiral\n");
        XLAL_ERROR(XLAL_EDOM);
    }

    q->ai = ai;
    q->bi = bi;

    ret = XLALIMRPhenomHMFreqDomainMapParams(&am, &bm, &fi, &fr, &f1, fi+Mfshift, ell, mm, PhenomDQuantities, AmpFlag);
    if (ret != XLAL_SUCCESS){
        XLALPrintError("XLAL Error - XLALIMRPhenomHMFreqDomainMapParams failed in XLALIMRPhenomHMFreqDomainMapParams - intermediate\n");
        XLAL_ERROR(XLAL_EDOM);
    }

    q->am = am;
    q->bm = bm;

    ret = XLALIMRPhenomHMFreqDomainMapParams(&ar, &br, &fi, &fr, &f1, fr+Mfshift, ell, mm, PhenomDQuantities, AmpFlag);
    if (ret != XLAL_SUCCESS){
        XLALPrintError("XLAL Error - XLALIMRPhenomHMFreqDomainMapParams failed in XLALIMRPhenomHMFreqDomainMapParams - merger-ringdown\n");
        XLAL_ERROR(XLAL_EDOM);
    }

    q->ar = ar;
    q->br = br;

    q->fi = fi;
    q->fr = fr;

    // printf("q->ai = %f\n", q->ai);
    // printf("q->bi = %f\n", q->bi);
    // printf("q->am = %f\n", q->am);
    // printf("q->bm = %f\n", q->bm);
    // printf("q->ar = %f\n", q->ar);
    // printf("q->br = %f\n", q->br);


    IMRPhenomDPhaseCoefficients *pPhi;
    pPhi = XLALMalloc(sizeof(IMRPhenomDPhaseCoefficients));
    ComputeIMRPhenomDPhaseCoefficients(pPhi, eta, chi1z, chi2z, PhenomDQuantities->finspin, extraParams);
    if (!pPhi) XLAL_ERROR(XLAL_EFUNC);

    PNPhasingSeries *pn = NULL;
    XLALSimInspiralTaylorF2AlignedPhasing(&pn, PhenomDQuantities->m1, PhenomDQuantities->m2, chi1z, chi2z, extraParams);
    if (!pn) XLAL_ERROR(XLAL_EFUNC);

    // Subtract 3PN spin-spin term below as this is in LAL's TaylorF2 implementation
    // (LALSimInspiralPNCoefficients.c -> XLALSimInspiralPNPhasing_F2), but
    // was not available when PhenomD was tuned.
    REAL8 testGRcor=1.0;
    testGRcor += XLALSimInspiralWaveformParamsLookupNonGRDChi6(extraParams);

    pn->v[6] -= (Subtract3PNSS(PhenomDQuantities->m1, PhenomDQuantities->m2, PhenomDQuantities->Mtot, eta, chi1z, chi2z) * pn->v[0])* testGRcor;


    PhiInsPrefactors phi_prefactors;
    int status = init_phi_ins_prefactors(&phi_prefactors, pPhi, pn);
    XLAL_CHECK(XLAL_SUCCESS == status, status, "init_phi_ins_prefactors failed");

    double Rholm = PhenomDQuantities->Rholm[ell][mm];
    double Taulm = PhenomDQuantities->Taulm[ell][mm];

    /* Compute coefficients to make phase C^1 continuous (phase and first derivative) */
    ComputeIMRPhenDPhaseConnectionCoefficients(pPhi, pn, &phi_prefactors, Rholm, Taulm);

    REAL8 PhDBMf = am*fi + bm;
    UsefulPowers powers_of_PhDBMf;
    status = init_useful_powers(&powers_of_PhDBMf, PhDBMf);
    XLAL_CHECK(XLAL_SUCCESS == status, status, "init_useful_powers for powers_of_PhDBMf failed");
    q->PhDBconst = IMRPhenDPhase(PhDBMf, pPhi, pn, &powers_of_PhDBMf, &phi_prefactors, Rholm, Taulm)/am;

    REAL8 PhDCMf = ar*fr + br;
    UsefulPowers powers_of_PhDCMf;
    status = init_useful_powers(&powers_of_PhDCMf, PhDCMf);
    XLAL_CHECK(XLAL_SUCCESS == status, status, "init_useful_powers for powers_of_PhDCMf failed");
    q->PhDCconst = IMRPhenDPhase(PhDCMf, pPhi, pn, &powers_of_PhDCMf, &phi_prefactors, Rholm, Taulm)/ar;

    REAL8 PhDBAMf = ai*fi + bi;
    UsefulPowers powers_of_PhDBAMf;
    status = init_useful_powers(&powers_of_PhDBAMf, PhDBAMf);
    XLAL_CHECK(XLAL_SUCCESS == status, status, "init_useful_powers for powers_of_PhDBAMf failed");
    q->PhDBAterm = IMRPhenDPhase(PhDBAMf, pPhi, pn, &powers_of_PhDBAMf, &phi_prefactors, Rholm, Taulm)/ai;
    LALFree(pPhi);
    LALFree(pn);

    return XLAL_SUCCESS;

}

double XLALSimIMRPhenomHMPhase( double Mf_wf, int mm, HMPhasePreComp *q, PNPhasingSeries *pn, IMRPhenomDPhaseCoefficients *pPhi, PhiInsPrefactors *phi_prefactors, double Rholm, double Taulm );
double XLALSimIMRPhenomHMPhase( double Mf_wf, /**< input frequency in geometric units*/
                                int mm,
                                HMPhasePreComp *q,
                                PNPhasingSeries *pn,
                                IMRPhenomDPhaseCoefficients *pPhi,
                                PhiInsPrefactors *phi_prefactors,
                                double Rholm,
                                double Taulm
                              )
{

    // const INT4 AmpFlagFalse = 0; /* FIXME: Could make this a global variable too */
    // double Mf_22 = XLALSimIMRPhenomHMFreqDomainMap( Mf_wf, ell, mm, eta, chi1z, chi2z, AmpFlagFalse );

    // UsefulPowers powers_of_Mf_22;
    // status = init_useful_powers(&powers_of_Mf_22, Mf_22);
    // XLAL_CHECK(XLAL_SUCCESS == status, status, "init_useful_powers for powers_of_Mf_22 failed");

    /* phi_lm(f) = m * phi_22(f_22) / 2.0 */
    // double PhenDphase = mm * IMRPhenDPhase(Mf_22, pPhi, pn, &powers_of_Mf_22, &phi_prefactors, Rholm, Taulm) / 2.0;

    REAL8 Mf = 0.0;
    REAL8 Mfr = 0.0;
    REAL8 retphase = 0.0;
    REAL8 tmpphaseC = 0.0;
    int status = 0;
    // printf("Mf_wf = %f, IMRPhenDPhaseA(Mf = q->ai * Mf_wf + q->bi;) = %f, IMRPhenDPhaseB(Mf = q->am*Mf_wf + q->bm) = %f, IMRPhenDPhaseC(Mf = q->ar*Mf_wf + q->br) = %f \n",Mf_wf, q->ai * Mf_wf + q->bi, q->am*Mf_wf + q->bm, q->ar*Mf_wf + q->br);
    // This if ladder is in the mathematica function HMPhase. PhenomHMDev.nb
    if ( !(Mf_wf > q->fi) ){ // in mathematica -> IMRPhenDPhaseA
        // printf("IMRPhenDPhaseA\n");
        Mf = q->ai * Mf_wf + q->bi;
        UsefulPowers powers_of_Mf;
        status = init_useful_powers(&powers_of_Mf, Mf);
        XLAL_CHECK(XLAL_SUCCESS == status, status, "init_useful_powers for powers_of_Mf failed");
        retphase = IMRPhenDPhase(Mf, pPhi, pn, &powers_of_Mf, phi_prefactors, Rholm, Taulm) / q->ai;
    } else if ( !(Mf_wf > q->fr) ){ // in mathematica -> IMRPhenDPhaseB
        // printf("IMRPhenDPhaseB\n");
        Mf = q->am*Mf_wf + q->bm;
        UsefulPowers powers_of_Mf;
        status = init_useful_powers(&powers_of_Mf, Mf);
        XLAL_CHECK(XLAL_SUCCESS == status, status, "init_useful_powers for powers_of_Mf failed");
        retphase = IMRPhenDPhase(Mf, pPhi, pn, &powers_of_Mf, phi_prefactors, Rholm, Taulm) / q->am - q->PhDBconst + q->PhDBAterm;
    } else if ( Mf_wf > q->fr ) { // in mathematica -> IMRPhenDPhaseC
        // printf("IMRPhenDPhaseC\n");
        Mfr = q->am*q->fr + q->bm;
        UsefulPowers powers_of_Mfr;
        status = init_useful_powers(&powers_of_Mfr, Mfr);
        XLAL_CHECK(XLAL_SUCCESS == status, status, "init_useful_powers for powers_of_Mfr failed");
        tmpphaseC = IMRPhenDPhase(Mfr, pPhi, pn, &powers_of_Mfr, phi_prefactors, Rholm, Taulm) / q->am - q->PhDBconst + q->PhDBAterm;
        /* tmpphaseC = IMRPhenDPhaseB[fr] in mathematica */

        Mf = q->ar*Mf_wf + q->br;
        UsefulPowers powers_of_Mf;
        status = init_useful_powers(&powers_of_Mf, Mf);
        XLAL_CHECK(XLAL_SUCCESS == status, status, "init_useful_powers for powers_of_Mf failed");
        retphase = IMRPhenDPhase(Mf, pPhi, pn, &powers_of_Mf, phi_prefactors, Rholm, Taulm) / q->ar - q->PhDCconst + tmpphaseC;
    } else {
        XLALPrintError("XLAL_ERROR - should not get here - in function XLALSimIMRPhenomHMPhase");
        XLAL_ERROR(XLAL_EDOM);
    }

    /*
     * Phase shift due to leading order complex amplitude
     * [L.Blancet, arXiv:1310.1528 (Sec. 9.5)]
     * "Spherical hrmonic modes for numerical relativity"
     */
    /*TODO: Need to have a list at the begining and a function to check the input
    lm mode to see if it is one that is included in the model.*/
    retphase += cShift[mm];

    // LALFree(pPhi);
    // LALFree(pn);

    return retphase;
}

/********************* Function to add modes for frequency-domain structures ********************/

/*
 * Helper function to add a mode to hplus, hcross in Fourier domain
 * - copies the function XLALSimAddMode, which was done only for TD structure
 * This function was lifted from the EOBNRv2HM_ROM code
 */
static INT4 FDAddMode(COMPLEX16FrequencySeries *hptilde, COMPLEX16FrequencySeries *hctilde, COMPLEX16FrequencySeries *hlmtilde, REAL8 theta, REAL8 phi, INT4 l, INT4 m, INT4 sym);
static INT4 FDAddMode(COMPLEX16FrequencySeries *hptilde, COMPLEX16FrequencySeries *hctilde, COMPLEX16FrequencySeries *hlmtilde, REAL8 theta, REAL8 phi, INT4 l, INT4 m, INT4 sym) {
  /* Deleted the definition of the string 'func': usage ? */
  COMPLEX16 Y;
  UINT4 j;
  COMPLEX16 hlmtildevalue;

  // printf("ell = %i, m = %i", l, m);
  // printf("length hlmtilde = %i  \n\n\n\n\n", hlmtilde->data->length);
  // printf("length hptilde = %i  \n\n\n\n\n", hptilde->data->length);
  // printf("length hctilde = %i  \n\n\n\n\n", hctilde->data->length);

  /* Checks LAL_CHECK_VALID_SERIES and LAL_CHECK_CONSISTENT_TIME_SERIES removed
   * - they do not seem available for frequency series ? */
  INT4 minus1l; /* (-1)^l */
  if ( l%2 ) minus1l = -1;
  else minus1l = 1;
  if ( sym ) { /* Equatorial symmetry: add in -m mode */
    Y = XLALSpinWeightedSphericalHarmonic(theta, phi, -2, l, m);
    COMPLEX16 Ymstar = conj(XLALSpinWeightedSphericalHarmonic(theta, phi, -2, l, -m));
    COMPLEX16 factorp = 0.5*(Y + minus1l*Ymstar);
    COMPLEX16 factorc = I*0.5*(Y - minus1l*Ymstar);
    //COMPLEX16* datap = hptilde->data->data;
    //COMPLEX16* datac = hctilde->data->data;
    for ( j = 0; j < hlmtilde->data->length; ++j ) {
      hlmtildevalue = (hlmtilde->data->data[j]);
      //datap[j] += factorp*hlmtildevalue;
      //datac[j] += factorc*hlmtildevalue;
      hptilde->data->data[j] += factorp*hlmtildevalue;
      hctilde->data->data[j] += factorc*hlmtildevalue;
    }
  }
  else { /* not adding in the -m mode */
    Y = XLALSpinWeightedSphericalHarmonic(theta, phi, -2, l, m);
    COMPLEX16 factorp = 0.5*Y;
    COMPLEX16 factorc = I*factorp;
    //COMPLEX16* datap = hptilde->data->data;
    //COMPLEX16* datac = hctilde->data->data;
    for ( j = 0; j < hlmtilde->data->length; ++j ) {
      hlmtildevalue = (hlmtilde->data->data[j]);
      //datap[j] += factorp*hlmtildevalue;
      //datac[j] += factorc*hlmtildevalue;
      hptilde->data->data[j] += factorp*hlmtildevalue;
      hctilde->data->data[j] += factorc*hlmtildevalue;
    }
  }

  return 0;
}

/*
 *
 * New functions in for the 'clean' code
 *
 *
 *
 */

/*
 * This returns the hlm of a single (l,m) mode at a single frequency Mf
 */
static COMPLEX16 IMRPhenomHMSingleModehlm(
        int ell,
        int mm,
        double Mf,
        HMPhasePreComp *z,
        IMRPhenomDAmplitudeCoefficients *pAmp,
        AmpInsPrefactors *amp_prefactors,
        PNPhasingSeries *pn,
        IMRPhenomDPhaseCoefficients *pPhi,
        PhiInsPrefactors *phi_prefactors,
        double Rholm,
        double Taulm,
        double phi_precalc, /**< 0.5*phi22(fref) - phi0*/
        PhenomDStorage *PhenomDQuantities
) {

    /*
     * In this function we should pass the phenom model parameters
     * and the (l,m) mode to return a given hlm, not summed with
     * spherical harmonics.
     * Can be evaluated at a single geometric frequency (Mf).
     */

    //FP: is all of PhenomDQuantities necessary?
    //FP: PhenomHMfring[ell][mm], Rholm[ell][mm], Mf_RD_22,
    //FP: pow_Mf_wf_prefactor[ell][mm]
    double HMamp = XLALSimIMRPhenomHMAmplitude( Mf, ell, mm, pAmp, amp_prefactors, PhenomDQuantities );
    double HMphase = XLALSimIMRPhenomHMPhase( Mf, mm, z, pn, pPhi, phi_prefactors, Rholm, Taulm );

    // printf("f(Hz) = %f, Mf = %f, HMamp = %.10e, HMphase = %f\n", Mf / (PhenomDQuantities->Mtot * LAL_MTSUN_SI), Mf, HMamp, HMphase);

    /* Compute reference phase at reference frequency */

    /* Factor of m spherical harmonic mode b/c phi0 is orbital phase */
    /* NOTE: Does HMphaseRef already have the mm scaling? as it's m*(phi0 + phiref) */
    HMphase -= mm * phi_precalc;
    // HMphase -= mm * phi_precalc;

    // Debug lines
    // double Seta = sqrt(1.0 - 4.0*pPhi->eta);
    // double m1 = 0.5 * (1.0 + Seta);
    // double m2 = 0.5 * (1.0 - Seta);
    // const REAL8 M_sec = (m1+m2) * LAL_MTSUN_SI; // Conversion factor Hz -> dimensionless frequency
    // if ( Mf > 0.05 && Mf < 0.051 ){
    //     printf("Mf = %.9f Hz = %f mode l = %i, m = %i       HMphase -= phi_precalc = %.9f\n", Mf, Mf/M_sec, ell, mm, HMphase);
    // }

    COMPLEX16 hlm = HMamp * cexp(-I * HMphase); //FP
    // printf("Mf = %f     hlm  =  %f + i %f\nx",Mf, creal(hlm), cimag(hlm));

    return hlm;
}

/* Given the final frequency in Mf and the total mass, calculate the final frequency in Hz */
/* ComputeIMRPhenomHMfmax:
 * given frequency parameters determine ending frequency
 * if f_max = 0 then return Mf (in Hz)
 * if f_max > Mf(in Hz) then return Mf(in Hz)
 * if f_max < Mf(in Hz) then return f_max
 */
static REAL8 ComputeIMRPhenomHMfmax(REAL8 Mf, REAL8 f_min, REAL8 f_max, REAL8 M);
static REAL8 ComputeIMRPhenomHMfmax(REAL8 Mf    /**< geometric frequency */,
                                    REAL8 f_min /**< low frequency in Hz */,
                                    REAL8 f_max /**< end frequency in Hz */,
                                    REAL8 M     /**< total mass (Msun) */
                                   ){

      const REAL8 M_sec = M * LAL_MTSUN_SI; // Conversion factor Hz -> dimensionless frequency
      const REAL8 fCut = Mf/M_sec; // convert Mf -> Hz

      /*
       * Somewhat arbitrary end point for the waveform.
       * Chosen so that the end of the waveform is well after the ringdown.
       */
      if (!(fCut > f_min))
          XLAL_ERROR(XLAL_EDOM, "(fCut = %g Hz) <= f_min = %g\n", fCut, f_min);

      /* Default f_max to Cut */
      REAL8 f_max_prime = f_max; // initialise f_max_prime
      f_max_prime = f_max ? f_max : fCut; // if f_max is zero then default to fCut otherwise use f_max
      f_max_prime = (f_max_prime > fCut) ? fCut : f_max_prime; // if user has asked for a frequency > fCut then return fCut, otherwise use f_max_prime.
      if (!(f_max_prime > f_min))
          XLAL_ERROR(XLAL_EDOM, "f_max <= f_min\n");

      return f_max_prime;
}

/* Compute t0 as the time of the peak of the 22 mode */
static REAL8 Computet0(REAL8 eta, REAL8 chi1z, REAL8 chi2z, REAL8 finspin, LALDict *extraParams);
static REAL8 Computet0(REAL8 eta, REAL8 chi1z, REAL8 chi2z, REAL8 finspin, LALDict *extraParams){

    if (extraParams==NULL)
      extraParams=XLALCreateDict();


    IMRPhenomDPhaseCoefficients *pPhi;
    pPhi = XLALMalloc(sizeof(IMRPhenomDPhaseCoefficients));
    ComputeIMRPhenomDPhaseCoefficients(pPhi, eta, chi1z, chi2z, finspin, extraParams);
    if (!pPhi) XLAL_ERROR(XLAL_EFUNC);

    IMRPhenomDAmplitudeCoefficients *pAmp;
    pAmp = XLALMalloc(sizeof(IMRPhenomDAmplitudeCoefficients));
    ComputeIMRPhenomDAmplitudeCoefficients(pAmp, eta, chi1z, chi2z, finspin);
    if (!pAmp) XLAL_ERROR(XLAL_EFUNC);

    // double Rholm = XLALSimIMRPhenomHMRholm(eta, chi1z, chi2z, ell, mm);
    // double Taulm = XLALSimIMRPhenomHMTaulm(eta, chi1z, chi2z, ell, mm);

    //time shift so that peak amplitude is approximately at t=0
    //For details see https://www.lsc-group.phys.uwm.edu/ligovirgo/cbcnote/WaveformsReview/IMRPhenomDCodeReview/timedomain
    //NOTE: All modes will have the same time offset. So we use the 22 mode.
    //If we just use the 22 mode then we pass 1.0, 1.0 into DPhiMRD.
    const REAL8 t0 = DPhiMRD(pAmp->fmaxCalc, pPhi, 1.0, 1.0);

    LALFree(pPhi);
    LALFree(pAmp);

    return t0;
}

int EnforcePrimaryIsm1(REAL8 *m1, REAL8 *m2, REAL8 *chi1z, REAL8 *chi2z);
int EnforcePrimaryIsm1(REAL8 *m1, REAL8 *m2, REAL8 *chi1z, REAL8 *chi2z){
    REAL8 chi1z_tmp, chi2z_tmp, m1_tmp, m2_tmp;
    if (*m1>*m2) {
       chi1z_tmp = *chi1z;
       chi2z_tmp = *chi2z;
       m1_tmp   = *m1;
       m2_tmp   = *m2;
   } else { /* swap spins and masses */
       chi1z_tmp = *chi2z;
       chi2z_tmp = *chi1z;
       m1_tmp   = *m2;
       m2_tmp   = *m1;
    }
    *m1 = m1_tmp;
    *m2 = m2_tmp;
    *chi1z = chi1z_tmp;
    *chi2z = chi2z_tmp;

    if (*m1 < *m2)
        XLAL_ERROR(XLAL_EDOM, "XLAL_ERROR in EnforcePrimaryIsm1. When trying\
 to enfore that m1 should be the larger mass.\
 After trying to enforce this m1 = %f and m2 = %f\n", *m1, *m2);

    return XLAL_SUCCESS;
}


/*
 * The goal of this function is to compute hlm for a list of modes
 * and output them in the data type SphHarmFrequencySeries.
 */

/* NOTE: name changed from hlmsphharmfreqseries to hlms */
/* NOTE: Some code duplication in order to keep this function XLAL */
int XLALIMRPhenomHMMultiModehlm(SphHarmFrequencySeries **hlms, REAL8 m1Msun, REAL8 m2Msun, REAL8 chi1z, REAL8 chi2z, REAL8 deltaF, REAL8 f_min, REAL8 f_max, REAL8 fRef_in, REAL8 phi0, REAL8 distance, LALDict *extraParams);
int XLALIMRPhenomHMMultiModehlm(
    SphHarmFrequencySeries **hlms, /**< [out] can access multiple modes with units of Mf */
    REAL8 m1Msun,                  /**< primary mass in Msun */
    REAL8 m2Msun,                  /**< secondary mass in Msun */
    REAL8 chi1z,                   /**< primary spin parameter */
    REAL8 chi2z,                   /**< secondary spin parameter */
    REAL8 deltaF,                  /**< Hz */
    REAL8 f_min,                   /**< Hz */
    REAL8 f_max,                   /**< Hz */
    REAL8 fRef_in,                 /**< reference frequency in Hz */
    REAL8 phi0,                    /**< reference orbital phase */
    REAL8 distance,                /**< distance to source in SI */
    LALDict *extraParams           /**< LALDict structure */
) {

    /* Powers of pi */
    int errcode = XLAL_SUCCESS;

    if (extraParams==NULL)
      extraParams=XLALCreateDict();

    //FP: these are known, so turn them in to #define's
    errcode = init_useful_powers(&powers_of_pi, LAL_PI);
    XLAL_CHECK(XLAL_SUCCESS == errcode, errcode, "init_useful_powers() failed: failed to initiate useful powers of pi.");

    /* Here masses are in Msun */
    errcode = EnforcePrimaryIsm1(&m1Msun, &m2Msun, &chi1z, &chi2z);
    XLAL_CHECK(XLAL_SUCCESS == errcode, errcode, "EnforcePrimaryIsm1 failed");

    /* Compute quantities/parameters related to PhenomD only once and store them */
    PhenomDStorage PhenomDQuantities;
    errcode = init_PhenomD_Storage(&PhenomDQuantities, m1Msun, m2Msun, chi1z, chi2z);
    XLAL_CHECK(XLAL_SUCCESS == errcode, errcode, "init_PhenomD_Storage failed");

    const REAL8 M = PhenomDQuantities.Mtot;
    const REAL8 eta = PhenomDQuantities.eta;
    const REAL8 M_sec = M * LAL_MTSUN_SI;

    if (eta > 0.25 || eta < 0.0)
        XLAL_ERROR(XLAL_EDOM, "Unphysical eta. Must be between 0. and 0.25\n");

    if (fabs(chi1z) > 1.0 || fabs(chi2z) > 1.0)
        XLAL_ERROR(XLAL_EDOM, "Spins outside the range [-1,1] are not supported\n");

    /* If no reference frequency given, set it to the starting GW frequency */
    REAL8 fRef = (fRef_in == 0.0) ? f_min : fRef_in;

    /* Compute phenomD amp coefficients */
    IMRPhenomDAmplitudeCoefficients *pAmp;
    pAmp = (IMRPhenomDAmplitudeCoefficients *) XLALMalloc(sizeof(IMRPhenomDAmplitudeCoefficients));
    ComputeIMRPhenomDAmplitudeCoefficients(pAmp, eta, chi1z, chi2z, PhenomDQuantities.finspin);
    if (!pAmp) XLAL_ERROR(XLAL_EFUNC);
    AmpInsPrefactors amp_prefactors;
    errcode = init_amp_ins_prefactors(&amp_prefactors, pAmp);
    XLAL_CHECK(XLAL_SUCCESS == errcode, errcode, "init_amp_ins_prefactors() failed.");

    /* Compute phenomD phase coefficients */
    IMRPhenomDPhaseCoefficients *pPhi;
    pPhi = (IMRPhenomDPhaseCoefficients *) XLALMalloc(sizeof(IMRPhenomDPhaseCoefficients));
    ComputeIMRPhenomDPhaseCoefficients(pPhi, eta, chi1z, chi2z, PhenomDQuantities.finspin, extraParams);
    if (!pPhi) XLAL_ERROR(XLAL_EFUNC);

    XLALSimInspiralWaveformParamsInsertPNSpinOrder(extraParams,LAL_SIM_INSPIRAL_SPIN_ORDER_35PN);
    PNPhasingSeries *pn = NULL;
    // FP: Malloc?
    XLALSimInspiralTaylorF2AlignedPhasing(&pn, m1Msun, m2Msun, chi1z, chi2z, extraParams);
    if (!pn) XLAL_ERROR(XLAL_EFUNC);

    // Subtract 3PN spin-spin term below as this is in LAL's TaylorF2 implementation
    // (LALSimInspiralPNCoefficients.c -> XLALSimInspiralPNPhasing_F2), but
    REAL8 testGRcor=1.0;
    testGRcor += XLALSimInspiralWaveformParamsLookupNonGRDChi6(extraParams);

    // Was not available when PhenomD was tuned.
    pn->v[6] -= (Subtract3PNSS(m1Msun, m2Msun, M, eta, chi1z, chi2z) * pn->v[0])* testGRcor;

    //FP: Malloc and free this too?
    PhiInsPrefactors phi_prefactors;
    //phi_prefactors = (PhiInsPrefactors *) XLALMalloc(sizeof(PhiInsPrefactors));
    int status = init_phi_ins_prefactors(&phi_prefactors, pPhi, pn);
    XLAL_CHECK(XLAL_SUCCESS == status, status, "init_phi_ins_prefactors failed");

    /* Compute the amplitude pre-factor */
    const REAL8 amp0 = M * LAL_MRSUN_SI * M_sec / distance;

    LIGOTimeGPS ligotimegps_zero = LIGOTIMEGPSZERO; // = {0, 0}

    /* Coalesce at t=0 */
    REAL8 InvDeltaF = 1./deltaF;
    /* Shift by overall length in time */
    XLAL_CHECK ( XLALGPSAdd(&ligotimegps_zero, - InvDeltaF), XLAL_EFUNC, "Failed to shift coalescence time to t=0, tried to apply shift of -1.0/deltaF with deltaF=%g.", deltaF);

    const REAL8 MfRef = M_sec * fRef;

    /*
     * In this function we should setup all the variables that we can
     * precompute to generate the PhenomD amplitude and phase functions
     * We then loop over the static function 'IMRPhenomHMSingleModehlm'
     * to generate the modes. We then sum them up at the end.
     */

    /* Now we have all the PhenomD model parameters, which actually correspond
     * to the (l,m)=(2,2) parameters, we can feed these into the function
     * IMRPhenomHMSingleModehlm to generate any mode we want. */

    /*
     * NOTE: I'm not sure what Mf should be used for the reference time...
     * I think it should be the scaled one. And it should come from the amplitude
     */


     /* Size of the final output hptilde and hctilde.
      All hlms have the same size as the final hptilde and hctilde
      However only they are non-zero in different ranges
      which depends on the ringdown frequency of each mode
      for optimisation purposes */
     size_t n = NextPow2(f_max * InvDeltaF) + 1;

     const REAL8 t0 = Computet0(eta, chi1z, chi2z, PhenomDQuantities.finspin, extraParams);

     double phi_22_at_MfRef = 0.0;
     double phi_const = 0.0;

     for( int j=0; j<NMODES; j++ ){

         int ell = ModeArray[j][0];
         int mm = ModeArray[j][1];

        //  printf("ell = %i\n", ell);
        //  printf("mm = %i\n", mm);


         /* we use a variable ending frequency for each mode so that we don't have to compute all the modes to the same high frequency unnecessarily */
        //  REAL8 f_max_prime = ComputeHigherModeEndingDimensionlessFrequency(ell, mm, &PhenomDQuantities) / M_sec; // Mf -> Hz
        //  REAL8 f_max_prime = F_FACTOR * PhenomDQuantities.PhenomHMfring[ell][mm] / M_sec;
        /* ComputeIMRPhenomHMfmax:
         * given frequency parameters determine ending frequency
         * if f_max = 0 then return the default proposed_f_end_lm_mf (in Hz)
         * if f_max > proposed_f_end_lm_mf(in Hz) then return proposed_f_end_lm_mf(in Hz)
         * if f_max < proposed_f_end_lm_mf(in Hz) then return f_max
         */
         REAL8 proposed_f_end_lm_mf = F_FACTOR * PhenomDQuantities.PhenomHMfring[ell][mm]; // (in Mf)

         REAL8 f_max_prime = ComputeIMRPhenomHMfmax(proposed_f_end_lm_mf, f_min, f_max, M);
        // REAL8 f_max_prime = 0.;
        //  REAL8 proposed_f_end_lm_hz = proposed_f_end_lm_mf / M_sec; // Mf -> Hz
        //  if (proposed_f_end_lm_hz < f_max) {
        //      f_max_prime = proposed_f_end_lm_hz;
        //  } else {
        //      f_max_prime = f_max;
        //  }
        //  printf("proposed_f_end_lm_hz = %f \n", proposed_f_end_lm_hz);
        //  printf("f_max_prime = %f \n", f_max_prime);


         /* Range that will have actual non-zero waveform values generated */
         /* These ranges depend on the l,m mode for optimisation purposes */
         size_t ind_min_lm = (size_t) (f_min * InvDeltaF);
         size_t ind_max_lm = (size_t) (f_max_prime * InvDeltaF);
         XLAL_CHECK ( !(ind_max_lm>n) && !(ind_min_lm>ind_max_lm), XLAL_EDOM, "minimum freq index %zu and maximum freq index %zu do not fulfill 0<=ind_min_lm<=ind_max_lm<=hptilde->data>length=%zu.", ind_min_lm, ind_max_lm, n);

         double Rholm = PhenomDQuantities.Rholm[ell][mm];
         double Taulm = PhenomDQuantities.Taulm[ell][mm];

         /* Compute coefficients to make phase C^1 continuous (phase and first derivative) */
         ComputeIMRPhenDPhaseConnectionCoefficients(pPhi, pn, &phi_prefactors, Rholm, Taulm);

         /* PhenomHM pre-computations */
         /* NOTE: Need to make this an input and NOT part of the frequency loop! */
         HMPhasePreComp z;
         int ret = XLALSimIMRPhenomHMPhasePreComp(&z, ell, mm, eta, chi1z, chi2z, &PhenomDQuantities, extraParams);
         if (ret != XLAL_SUCCESS){
             XLALPrintError("XLAL Error - XLALSimIMRPhenomHMPhasePreComp failed\n");
             XLAL_ERROR(XLAL_EDOM);
         }

         //The size of this COMPLEX16FrequencySeries *hlm needs to be the same as the final hptilde, hctilde
         //The for loop over frequencies needs to be over the indices for the particular (l,m) mode

         /* We loop over (l,m) and use a temporary hlm frequency series to store the results of a single mode */
         COMPLEX16FrequencySeries *hlm = NULL;
         hlm = XLALCreateCOMPLEX16FrequencySeries("hlm: single mode", &ligotimegps_zero, 0.0, deltaF, &lalStrainUnit, n);
         memset(hlm->data->data, 0, n * sizeof(COMPLEX16));

         /* NOTE: Do I need this bit? */
         /* XLALUnitMultiply(hlm->sampleUnits, hlm->sampleUnits, &lalSecondUnit); */

         /* begin computing reference phase shift*/
         /* NOTE: Only works if l=2 and m=2 is first in the list */
         if (ell == 2 && mm == 2 && j==0){
             /* compute the reference phase constant to add to each mode */
             /* NOTE: This only needs to be done once but I'm too lazy to optimise. */
             /* We enforce that phi_22(fref) = 2.0 * phiRef
              * Where phiRef is the inpute reference orbital phase. */
             /* Evaluating XLALSimIMRPhenomHMPhase for the ell=2, mm=2 mode */
             phi_22_at_MfRef = XLALSimIMRPhenomHMPhase( MfRef, 2, &z, pn, pPhi, &phi_prefactors, 1.0, 1.0 );
             phi_const = 0.5 * phi_22_at_MfRef - phi0;
             /* phase_lm (f) -= m*phi_const */
         } else if (ell == 2 && mm == 2 && j!=0){
             XLALPrintError("l=2, m=2 should be the first mode in the list called ModeArray.");
         }


         /* end compute reference phase shift */



         /* correct the time shift at the correctly scaled reference frequency */
        //  double Mf_22_REF =  XLALSimIMRPhenomHMFreqDomainMap(MfRef, ell, mm, &PhenomDQuantities, AmpFlagFalse);



         /* Now generate the waveform for a single (l,m) mode, i.e. compute hlm*/
         /* Loop over frequency */
         REAL8 M_sec_dF = M_sec * deltaF;
         COMPLEX16 It0 = I*t0;
         #pragma omp parallel for
         for (size_t i = ind_min_lm; i < ind_max_lm; i++)
         {
            REAL8 Mf = i * M_sec_dF; /* geometric frequency */
            /* TODO: fix phase and time offsets */
            /* TODO: inclusion of reference frequency */
            // phi -= t0*(Mf-MfRef) + phi_precalc;
            // ((*htilde)->data->data)[i] = amp0 * amp * cexp(-I * phi);
            /* construct hlm at single frequency point and return */
            // (hlm->data->data)[i] = amp0 * IMRPhenomHMSingleModehlm(eta, chi1z, chi2z, ell, mm, Mf, MfRef, phi0, &z);

            //FP: is all of PhenomDQuantities necessary?
            //FP: PhenomHMfring[ell][mm], Rholm[ell][mm], Mf_RD_22,
            //FP: pow_Mf_wf_prefactor[ell][mm]
            (hlm->data->data)[i] = amp0 * IMRPhenomHMSingleModehlm(ell, mm, Mf, &z, pAmp, &amp_prefactors, pn, pPhi, &phi_prefactors, Rholm, Taulm, phi_const, &PhenomDQuantities);
            /* NOTE: The frequency used in the time shift term is the fourier variable of the gravitational wave frequency. i.e., Not rescaled. */
            /* NOTE: normally the t0 term is multiplied by 2pi but the 2pi has been absorbed into the t0. */
            // (hlm->data->data)[i] *= cexp(-I * LAL_PI*t0*(Mf-MfRef)*(2.0-mm) );
            // (hlm->data->data)[i] *= cexp(-I * t0*(Mf-MfRef)*(2.0-mm) );
            (hlm->data->data)[i] *= cexp(-It0*(Mf-MfRef)); //FIXME: 2.0-mm is gone, is this intended? If yes, delete this comment.

            // double Mf_22 =  XLALSimIMRPhenomHMFreqDomainMap(Mf, ell, mm, &PhenomDQuantities, AmpFlagFalse);
            // (hlm->data->data)[i] *= cexp(-It0*(Mf_22-Mf_22_REF));

            // (hlm->data->data)[i] *= cexp(-It0*(Mf-MfRef)*mm/2.0);
            // (hlm->data->data)[i] *= cexp(-It0 * ( (Mf*mm/2.0) - MfRef ));


            /* From phenomD for referene*/
            // REAL8 amp = IMRPhenDAmplitude(Mf, pAmp, &powers_of_f, &amp_prefactors);
            // REAL8 phi = IMRPhenDPhase(Mf, pPhi, pn, &powers_of_f, &phi_prefactors);
            // phi -= t0*(Mf-MfRef) + phi_precalc;
            //  ((*htilde)->data->data)[i] = amp0 * amp * cexp(-I * phi);
         }
         *hlms = XLALSphHarmFrequencySeriesAddMode(*hlms, hlm, ell, mm);
         /* Destroy hlm in (l,m) loop */
         XLALDestroyCOMPLEX16FrequencySeries(hlm);
     }

     LALFree(pAmp);
     LALFree(pPhi);
     LALFree(pn);

    return XLAL_SUCCESS;
}

/* This function is a wrapper of XLALIMRPhenomHMMultiModehlm.
 * It takes the SphHarmFrequencySeries **hlms as inputs and
 * compbines them with the Spherical harmonics at a given
 * inclination angle and azimuthal angle. It returns the
 * hptilde and hctilde.
 */
 /*TODO: Add LALDict as an argument to XLALIMRPhenomHMMultiModeStrain */
 /*TODO: Change name from XLALIMRPhenomHMMultiModeStrain to XLALIMRPhenomHM */
int XLALIMRPhenomHMMultiModeStrain(
    COMPLEX16FrequencySeries **hptilde, /**< [out] */
    COMPLEX16FrequencySeries **hctilde, /**< [out] */
    REAL8 m1,                           /**< primary mass in SI */
    REAL8 m2,                           /**< secondary mass in SI */
    REAL8 chi1z,                        /**< primary spin parameter */
    REAL8 chi2z,                        /**< secondary spin parameter */
    REAL8 deltaF,                       /**< Hz */
    REAL8 f_min,                        /**< Hz */
    REAL8 f_max,                        /**< Hz */
    REAL8 fRef_in,                      /**< reference frequency in Hz */
    REAL8 phi0,                         /**< reference orbital phase */
    REAL8 inclination,                  /**< inclination... */
    REAL8 distance,                     /**< distance to source in SI */
    LALDict *extraParams                /**< LALDict structure */
) {

    /* Sanity checks on input parameters */
    if (m1 < 0) XLAL_ERROR(XLAL_EDOM, "m1 must be positive\n");
    if (m2 < 0) XLAL_ERROR(XLAL_EDOM, "m2 must be positive\n");
    if (fabs(chi1z) > 1.0 || fabs(chi2z) > 1.0 )
        XLAL_ERROR(XLAL_EDOM, "Spins outside the range [-1,1] are not supported\n");
    if (deltaF < 0) XLAL_ERROR(XLAL_EDOM, "deltaF must be positive\n");
    if (f_min < 0) XLAL_ERROR(XLAL_EDOM, "f_min must be positive\n");
    if (f_max < 0) XLAL_ERROR(XLAL_EDOM, "f_max must be greater than 0\n");

    /* External: SI; internal: solar masses */
    m1 /= LAL_MSUN_SI;
    m2 /= LAL_MSUN_SI;

    int ret = EnforcePrimaryIsm1(&m1, &m2, &chi1z, &chi2z);
    XLAL_CHECK(XLAL_SUCCESS == ret, ret, "EnforcePrimaryIsm1 failed");

    /* Mass ratio >= 1 convention */
    const REAL8 q = (m1 > m2) ? (m1 / m2) : (m2 / m1);

    if (q > MAX_ALLOWED_MASS_RATIO)
      XLAL_PRINT_WARNING("Warning: The model is not supported for high mass ratio, see MAX_ALLOWED_MASS_RATIO\n");

    const REAL8 M = m1 + m2; /* total mass (Msun) */
    REAL8 eta = m1 * m2 / (M * M);

    if (eta > 0.25)
        nudge(&eta, 0.25, 1e-6);
    if (eta > 0.25 || eta < 0.0)
        XLAL_ERROR(XLAL_EDOM, "Unphysical eta. Must be between 0. and 0.25\n");

    /* If no reference frequency given, set it to the starting GW frequency */
    REAL8 fRef = (fRef_in == 0.0) ? f_min : fRef_in;

    /* Given the final frequency in Mf and the total mass, calculate the final frequency in Hz */
    /* ComputeIMRPhenomHMfmax:
     * given frequency parameters determine ending frequency
     * if f_max = 0 then return the default Mf_CUT_HM (in Hz)
     * if f_max > Mf_CUT_HM(in Hz) then return Mf_CUT_HM(in Hz)
     * if f_max < Mf_CUT_HM(in Hz) then return f_max
     */
    //f_max_strain : this is the final frequency to be used in hptilde and hctilde
    REAL8 f_max_strain = ComputeIMRPhenomHMfmax(Mf_CUT_HM, f_min, f_max, M);
    /* The aim of this line is to determine the ending frequency for waveform
     * generation.
     * We only want to compute each mode up to approx 1.2*fRDlm
     * For frequencies higher than this we fill with zeros.
     */

    /* Evaluate XLALIMRPhenomHMMultiModehlm */

    // SphHarmFrequencySeries *hlms=NULL;
    SphHarmFrequencySeries **hlms=XLALMalloc(sizeof(SphHarmFrequencySeries));
    *hlms=NULL;

    /*TODO: Add LALDict as an argument to XLALIMRPhenomHMMultiModehlm */
    ret = XLALIMRPhenomHMMultiModehlm(hlms, m1, m2, chi1z, chi2z, deltaF, f_min, f_max_strain, fRef, phi0, distance, extraParams);
    XLAL_CHECK(XLAL_SUCCESS == ret, ret, "XLALIMRPhenomHMMultiModehlm(&hlms) failed");

    LIGOTimeGPS ligotimegps_zero = LIGOTIMEGPSZERO; // = {0, 0}

    /* Coalesce at t=0 */
    /* Shift by overall length in time */
    REAL8 InvDeltaF = 1./deltaF;
    XLAL_CHECK ( XLALGPSAdd(&ligotimegps_zero, - InvDeltaF), XLAL_EFUNC, "Failed to shift coalescence time to t=0, tried to apply shift of -1.0/deltaF with deltaF=%g.", deltaF);

    /* Compute array sizes */
    size_t n = NextPow2(f_max_strain * InvDeltaF) + 1;

    /* Allocate hptilde and hctilde */
    *hptilde = XLALCreateCOMPLEX16FrequencySeries("hptilde: FD waveform", &ligotimegps_zero, 0.0, deltaF, &lalStrainUnit, n);
    if (!(hptilde) ) XLAL_ERROR(XLAL_EFUNC);
    memset((*hptilde)->data->data, 0, n * sizeof(COMPLEX16));
    XLALUnitDivide(&(*hptilde)->sampleUnits, &(*hptilde)->sampleUnits, &lalSecondUnit);

    *hctilde = XLALCreateCOMPLEX16FrequencySeries("hctilde: FD waveform", &ligotimegps_zero, 0.0, deltaF, &lalStrainUnit, n);
    if (!(hctilde) ) XLAL_ERROR(XLAL_EFUNC);
    memset((*hctilde)->data->data, 0, n * sizeof(COMPLEX16));
    XLALUnitDivide(&(*hctilde)->sampleUnits, &(*hctilde)->sampleUnits, &lalSecondUnit);

    /* Adding the modes to form hplus, hcross
     * - use of a function that copies XLALSimAddMode but for Fourier domain structures */
    INT4 sym; /* sym will decide whether to add the -m mode (when equatorial symmetry is present) */
    for( int i=0; i<NMODES; i++){
      INT4 ell = ModeArray[i][0];
      INT4 mm = ModeArray[i][1];

    //   printf("ell = %i\n", ell);
    //   printf("mm = %i\n", mm);

      COMPLEX16FrequencySeries* hlm = XLALSphHarmFrequencySeriesGetMode(*hlms, ell, mm);
      if (!(hlm)) XLAL_ERROR(XLAL_EFUNC);

      /* We test for hypothetical m=0 modes */
      if ( mm==0 ) {
          sym = 0;
      } else {
          sym = 1;
      }
    FDAddMode( *hptilde, *hctilde, hlm, inclination, 0., ell, mm, sym); /* The phase \Phi is set to 0 - assumes phiRef is defined as half the phase of the 22 mode h22 (or the first mode in the list), not for h = hplus-I hcross */
    // FDAddMode( *hptilde, *hctilde, hlm, inclination, phi0, ell, mm, sym); /* Added phi0 here as a quick fix for the reference phase. not sure if it should be m * phi0 or m/2*phi0 . */
    }



    XLALDestroySphHarmFrequencySeries(*hlms);
    XLALFree(hlms);

    // The following lines zero-pad the data to the desired f_max
    // Very important for when performing the inverse FFT
    const REAL8 M_sec = M * LAL_MTSUN_SI; // Add M_sec to PhenomDQuantities?
    REAL8 f_CUT_Hz = Mf_CUT_HM/M_sec;

    //TODO: Add something like this from phenomD
    if (f_max_strain < f_max) {
      // The user has requested a higher f_max than Mf=fCut.
      // Resize the frequency series to fill with zeros beyond the cutoff frequency.
      size_t n_len = (*hptilde)->data->length;
      size_t n_full = NextPow2(f_max / deltaF) + 1; // we actually want to have the length be a power of 2 + 1
      *hptilde = XLALResizeCOMPLEX16FrequencySeries(*hptilde, 0, n_full);
      XLAL_CHECK ( *hptilde, XLAL_ENOMEM, "Failed to resize waveform COMPLEX16FrequencySeries of length %zu (for internal f_CUT_Hz=%f) to new length %zu (for user-requested f_max=%f).", n_len, f_CUT_Hz, n_full, f_max );
      *hctilde = XLALResizeCOMPLEX16FrequencySeries(*hctilde, 0, n_full);
      XLAL_CHECK ( *hctilde, XLAL_ENOMEM, "Failed to resize waveform COMPLEX16FrequencySeries of length %zu (for internal f_CUT_Hz=%f) to new length %zu (for user-requested f_max=%f).", n_len, f_CUT_Hz, n_full, f_max );
    }

    // NOTE: SK: HERE I SWAP hplus with hcross to conform with LAL phase convension
    #pragma omp parallel for
    for (size_t i = 0; i < (*hptilde)->data->length; i++)
    {
       ((*hptilde)->data->data)[i] = I*((*hptilde)->data->data)[i];
       ((*hctilde)->data->data)[i] = -I*((*hctilde)->data->data)[i];
    }



    return XLAL_SUCCESS;
}


/* Convenience function to compute hlm for a single mode - from which
 * the Alm and philm (amplitude and phase) of a particular mode can be
 * obtained. Useful to compare to NR */
int XLALSimIMRPhenomHMSingleModehlm(COMPLEX16FrequencySeries **hlmtilde, /**< [out] */
                                    REAL8 m1Msun, /**< mass1 in units of solar masses */
                                    REAL8 m2Msun, /**< mass2 in units of solar masses */
                                    REAL8 chi1z,
                                    REAL8 chi2z,
                                    REAL8 deltaF,
                                    REAL8 f_min,
                                    REAL8 f_max,
                                    REAL8 fRef_in,
                                    REAL8 phi0,
                                    REAL8 distance,
                                    INT4 ell,
                                    INT4 mm,
                                    LALDict *extraParams){

    int errcode = XLAL_SUCCESS;
    errcode = init_useful_powers(&powers_of_pi, LAL_PI);
    XLAL_CHECK(XLAL_SUCCESS == errcode, errcode, "init_useful_powers() failed: failed to initiate useful powers of pi.");

    if (extraParams==NULL)
      extraParams=XLALCreateDict();

    /* Here masses are in Msun */
    errcode = EnforcePrimaryIsm1(&m1Msun, &m2Msun, &chi1z, &chi2z);
    XLAL_CHECK(XLAL_SUCCESS == errcode, errcode, "EnforcePrimaryIsm1 failed");

    // printf(" m1Msun= %f\n",m1Msun);
    // printf(" m2Msun= %f\n",m2Msun);
    // printf(" chi1z= %f\n",chi1z);
    // printf(" chi2z= %f\n",chi2z);

    PhenomDStorage PhenomDQuantities;
    errcode = init_PhenomD_Storage(&PhenomDQuantities, m1Msun, m2Msun, chi1z, chi2z);
    XLAL_CHECK(XLAL_SUCCESS == errcode, errcode, "init_PhenomD_Storage failed");
    const REAL8 M = PhenomDQuantities.Mtot;//m1Msun + m2Msun;
    REAL8 eta = PhenomDQuantities.eta;

    if (eta > 0.25)
        nudge(&eta, 0.25, 1e-6);
    if (eta > 0.25 || eta < 0.0)
        XLAL_ERROR(XLAL_EDOM, "Unphysical eta. Must be between 0. and 0.25\n");


    const REAL8 M_sec = M * LAL_MTSUN_SI; // Add M_sec to PhenomDQuantities?

    if (eta > 0.25 || eta < 0.0)
        XLAL_ERROR(XLAL_EDOM, "Unphysical eta. Must be between 0. and 0.25\n");

    if (fabs(chi1z) > 1.0 || fabs(chi2z) > 1.0)
        XLAL_ERROR(XLAL_EDOM, "Spins outside the range [-1,1] are not supported\n");

    // if no reference frequency given, set it to the starting GW frequency
    REAL8 fRef = (fRef_in == 0.0) ? f_min : fRef_in;

    //f_max_strain : this is the final frequency to be used in hptilde and hctilde
    REAL8 f_max_strain = ComputeIMRPhenomHMfmax(Mf_CUT_HM, f_min, f_max, M);


    /* Compute PhenomD amp and phase coefficients*/

    IMRPhenomDAmplitudeCoefficients *pAmp;
    pAmp = XLALMalloc(sizeof(IMRPhenomDAmplitudeCoefficients));
    ComputeIMRPhenomDAmplitudeCoefficients(pAmp, eta, chi1z, chi2z, PhenomDQuantities.finspin);

    // printf("pAmp->chi1 = %f\n", pAmp->chi1);

    if (!pAmp) XLAL_ERROR(XLAL_EFUNC);
    AmpInsPrefactors amp_prefactors;
    errcode = init_amp_ins_prefactors(&amp_prefactors, pAmp);
    XLAL_CHECK(XLAL_SUCCESS == errcode, errcode, "init_amp_ins_prefactors() failed.");
    IMRPhenomDPhaseCoefficients *pPhi;
    pPhi = XLALMalloc(sizeof(IMRPhenomDPhaseCoefficients));
    ComputeIMRPhenomDPhaseCoefficients(pPhi, eta, chi1z, chi2z, PhenomDQuantities.finspin, extraParams);
    if (!pPhi) XLAL_ERROR(XLAL_EFUNC);
    XLALSimInspiralWaveformParamsInsertPNSpinOrder(extraParams,LAL_SIM_INSPIRAL_SPIN_ORDER_35PN);
    PNPhasingSeries *pn = NULL;
    XLALSimInspiralTaylorF2AlignedPhasing(&pn, m1Msun, m2Msun, chi1z, chi2z, extraParams);
    if (!pn) XLAL_ERROR(XLAL_EFUNC);


    // printf("pPhi->beta1 = %f\n", pPhi->beta1);

    // Subtract 3PN spin-spin term below as this is in LAL's TaylorF2 implementation
    // (LALSimInspiralPNCoefficients.c -> XLALSimInspiralPNPhasing_F2), but
    REAL8 testGRcor=1.0;
    testGRcor += XLALSimInspiralWaveformParamsLookupNonGRDChi6(extraParams);

    // was not available when PhenomD was tuned.
    pn->v[6] -= (Subtract3PNSS(m1Msun, m2Msun, M, eta, chi1z, chi2z) * pn->v[0])* testGRcor;

    PhiInsPrefactors phi_prefactors;
    errcode = init_phi_ins_prefactors(&phi_prefactors, pPhi, pn);
    XLAL_CHECK(XLAL_SUCCESS == errcode, errcode, "init_phi_ins_prefactors failed");
    /* NOTE: There seems to be a problem here, with the phi_prefactors as in IMRPhenomHMMultiModehlm they work fine here but in this
    function they do not. */
    /*FIXME */
    /*FIXME */
    /*FIXME */
    /*FIXME */
    /*FIXME */
    /*FIXME */
    /*FIXME */

    /* Compute the amplitude pre-factor */
    const REAL8 amp0 = M * LAL_MRSUN_SI * M * LAL_MTSUN_SI / distance;

    LIGOTimeGPS ligotimegps_zero = LIGOTIMEGPSZERO; // = {0, 0}

    /* Coalesce at t=0 */
    REAL8 InvDeltaF = 1./deltaF;
    // shift by overall length in time
    XLAL_CHECK ( XLALGPSAdd(&ligotimegps_zero, -InvDeltaF), XLAL_EFUNC, "Failed to shift coalescence time to t=0, tried to apply shift of -1.0/deltaF with deltaF=%g.", deltaF);

    /* compute array sizes */
    size_t n = NextPow2(f_max_strain * InvDeltaF) + 1;
    /* range that will have actual non-zero waveform values generated */
    size_t ind_min = (size_t) (f_min * InvDeltaF);
    size_t ind_max = (size_t) (f_max_strain * InvDeltaF);
    XLAL_CHECK ( !(ind_max>n) && !(ind_min>ind_max), XLAL_EDOM, "minimum freq index %zu and maximum freq index %zu do not fulfill 0<=ind_min<=ind_max<=hptilde->data>length=%zu.", ind_min, ind_max, n);


    /* Allocate hptilde and hctilde */
    *hlmtilde = XLALCreateCOMPLEX16FrequencySeries("hlmtilde: FD waveform", &ligotimegps_zero, 0.0, deltaF, &lalStrainUnit, n);
    if (!(hlmtilde) ) XLAL_ERROR(XLAL_EFUNC);
    memset((*hlmtilde)->data->data, 0, n * sizeof(COMPLEX16));
    XLALUnitDivide(&(*hlmtilde)->sampleUnits, &(*hlmtilde)->sampleUnits, &lalSecondUnit);

    const REAL8 MfRef = M_sec * fRef;

    /*
     * In this function we should setup
     * all the variables that we can precompute
     * to generate the PhenomD amplitude and phase functions
     * We then loop over the static function 'IMRPhenomHMSingleModehlm'
     * to generate the modes.
     * We then sum them up at the end.
     */

     /* Now we have all the PhenomD model parameters, which actually correspond
      * to the (l,m)=(2,2) parameters, we can feed these into the function
      * IMRPhenomHMSingleModehlm to generate any mode we want. */

     //TODO: Turn t0 computation into a function
     /* NOTE: We could compute the t0 here as the time of the peak of the
        22 mode and make that the assumtion.*/


     /*
      * NOTE: I'm not sure what Mf should be used for the reference time... I think it should be the scaled one. And it should come from the amplitude
      */

     const REAL8 t0 = Computet0(eta, chi1z, chi2z, PhenomDQuantities.finspin, extraParams);

     double Rholm = PhenomDQuantities.Rholm[ell][mm];
     double Taulm = PhenomDQuantities.Taulm[ell][mm];

     // Compute coefficients to make phase C^1 continuous (phase and first derivative)
     ComputeIMRPhenDPhaseConnectionCoefficients(pPhi, pn, &phi_prefactors, Rholm, Taulm);
     // printf("pPhi->C1Int = %f\n", pPhi->C1Int);
     /* compute phenomHM pre computations */
     /* NOTE: Need to make this an input and NOT part of the frequency loop! */
     HMPhasePreComp z;
     errcode = XLALSimIMRPhenomHMPhasePreComp(&z, ell, mm, eta, chi1z, chi2z, &PhenomDQuantities, extraParams);
     if (errcode != XLAL_SUCCESS){
         XLALPrintError("XLAL Error - XLALSimIMRPhenomHMPhasePreComp failed\n");
         XLAL_ERROR(XLAL_EDOM);
     }

      /* compute the reference phase constant to add to each mode */
      /* NOTE: This only needs to be done once but I'm too lazy to optimise. */
      /* We enforce that phi_22(fref) = 2.0 * phiRef
       * Where phiRef is the inpute reference orbital phase. */
      /* Evaluating XLALSimIMRPhenomHMPhase for the ell=2, mm=2 mode */
      double phi_22_at_MfRef = XLALSimIMRPhenomHMPhase( MfRef, 2, &z, pn, pPhi, &phi_prefactors, 1.0, 1.0 );
      double phi_const = 0.5 * phi_22_at_MfRef - phi0;
      /* phase_lm (f) -= m*phi_const */
      // printf("phi_const = %f\n",phi_const);

     /* NOTE: Do I need this bit? */
     /* XLALUnitMultiply(hlm->sampleUnits, hlm->sampleUnits, &lalSecondUnit); */

     /* LOOP OVER FREQUENCY */
     /* Now generate the waveform for a single (l,m) mode */
     REAL8 M_sec_dF = M_sec * deltaF;
     COMPLEX16 It0 = I*t0;
     #pragma omp parallel for
     for (size_t i = ind_min; i < ind_max; i++)
     {
        REAL8 Mf = i * M_sec_dF; /* geometric frequency */
        /* now we can compute the hlm */
        /* TODO: fix phase and time offsets */
        /* TODO: inclusion of reference frequency */
        // phi -= t0*(Mf-MfRef) + phi_precalc;
        // ((*htilde)->data->data)[i] = amp0 * amp * cexp(-I * phi);
        /* construct hlm at single frequency point and return */
        // (hlm->data->data)[i] = amp0 * IMRPhenomHMSingleModehlm(eta, chi1z, chi2z, ell, mm, Mf, MfRef, phi0, &z);

        ((*hlmtilde)->data->data)[i] = amp0 * IMRPhenomHMSingleModehlm(ell, mm, Mf, &z, pAmp, &amp_prefactors, pn, pPhi, &phi_prefactors, Rholm, Taulm, phi_const, &PhenomDQuantities);

        /* NOTE: The frequency used in the time shift term is the fourier variable of the gravitational wave frequency. i.e., Not rescaled. */
        // ((*hlmtilde)->data->data)[i] *= cexp(-It0*(Mf-MfRef)*(2.0-mm) );
        ((*hlmtilde)->data->data)[i] *= cexp(-It0*(Mf-MfRef)); //FIXME: 2.0-mm is gone, is this intended? If yes, delete this comment.
        /*from phenomD for referene*/
        //  REAL8 amp = IMRPhenDAmplitude(Mf, pAmp, &powers_of_f, &amp_prefactors);
        //  REAL8 phi = IMRPhenDPhase(Mf, pPhi, pn, &powers_of_f, &phi_prefactors);
        //  phi -= t0*(Mf-MfRef) + phi_precalc;
        //  ((*htilde)->data->data)[i] = amp0 * amp * cexp(-I * phi);
     }

     // TODO: CHECK THE SIGN OF THE PHASE IN THIS FUNCTION because at the of XLALIMRPhenomHMMultiModeStrain I swap hplus with hcross to get the right phase convention
     // TODO: CHECK THE SIGN OF THE PHASE IN THIS FUNCTION because at the of XLALIMRPhenomHMMultiModeStrain I swap hplus with hcross to get the right phase convention
     // TODO: CHECK THE SIGN OF THE PHASE IN THIS FUNCTION because at the of XLALIMRPhenomHMMultiModeStrain I swap hplus with hcross to get the right phase convention

     LALFree(pAmp);
     LALFree(pPhi);
     LALFree(pn);

    return XLAL_SUCCESS;
}




// Taken from LALSimIMRPhenomP.c
// This function determines whether x and y are approximately equal to a relative accuracy epsilon.
// Note that x and y are compared to relative accuracy, so this function is not suitable for testing whether a value is approximately zero.
static bool approximately_equal(REAL8 x, REAL8 y, REAL8 epsilon) {
  return !gsl_fcmp(x, y, epsilon);
}

// If x and X are approximately equal to relative accuracy epsilon then set x = X.
// If X = 0 then use an absolute comparison.
// Taken from LALSimIMRPhenomP.c
static void nudge(REAL8 *x, REAL8 X, REAL8 epsilon) {
  if (X != 0.0) {
    if (approximately_equal(*x, X, epsilon)) {
      XLAL_PRINT_INFO("Nudging value %.15g to %.15g\n", *x, X);
      *x = X;
    }
  }
  else {
    if (fabs(*x - X) < epsilon)
      *x = X;
  }
}
