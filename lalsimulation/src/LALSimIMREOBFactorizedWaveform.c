/*
*  Copyright (C) 2010 Craig Robinson, Yi Pan
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


/**
 * \author Craig Robinson, Yi Pan
 *
 * \brief The functions contained within this file pre-compute the various
 * coefficients which are required for calculating the factorized waveform
 * in EOBNRv2. Note that for some of the higher modes, the coefficients
 * are changed for the generation of the waveform compared to the generation
 * of the flux. Thus we have a function which adds these additional 
 * contributions to the already computed coefficients.
 */

#include <math.h>
#include <complex.h>
#include "LALSimIMREOBNRv2.h"

/* Include static functions */
#include "LALSimInspiraldEnergyFlux.c"
#include "LALSimIMREOBNewtonianMultipole.c" 
#include "LALSimIMREOBNQCCorrection.c"

#ifndef _LALSIMIMRFACTORIZEDWAVEFORM_C
#define _LALSIMIMRFACTORIZEDWAVEFORM_C

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

UNUSED static int GetGammaFuncPolyFitCoeffs(INT4 l, REAL8 *rc, REAL8 *ic, INT4 *n);
UNUSED static int GammaFunctionPolyFit(INT4 l, REAL8 hathatk, REAL8 *lnr, REAL8 *arg);

/**
 * Constant which comes up in some of the EOB models. Its value is
 * (94/3 -41/32*pi*pi)
 */
#define ninty4by3etc 18.687902694437592603


static inline REAL8 XLALCalculateA5( REAL8 eta );

static inline REAL8 XLALCalculateA6( REAL8 eta );

static
REAL8 XLALCalculateEOBD( REAL8    r,
                         REAL8	eta) UNUSED;


/**
 * Calculates the a5 parameter in the A potential function in EOBNRv2
 */
static inline
REAL8 XLALCalculateA5( const REAL8 eta /**<< Symmetric mass ratio */
                     )
{
  return - 5.82827 - 143.486 * eta + 447.045 * eta * eta;
}

/**
 * Calculates the a6 parameter in the A potential function in EOBNRv2
 */
static inline
REAL8 XLALCalculateA6( const REAL8 UNUSED eta /**<< Symmetric mass ratio */
                     )
{
  return 184.0;
}


/**
 * Function to pre-compute the coefficients in the EOB A potential function
 */
UNUSED static
int XLALCalculateEOBACoefficients(
          EOBACoefficients * const coeffs, /**<< A coefficients (populated in function) */
          const REAL8              eta     /**<< Symmetric mass ratio */
          )
{
  REAL8 eta2, eta3;
  REAL8 a4, a5, a6;

  eta2 = eta*eta;
  eta3 = eta2 * eta;

  /* Note that the definitions of a5 and a6 DO NOT correspond to those in the paper */
  /* Therefore we have to multiply the results of our a5 and a6 finctions by eta. */

  a4 = ninty4by3etc * eta;
  a5 = XLALCalculateA5( eta ) * eta;
  a6 = XLALCalculateA6( eta ) * eta;

  coeffs->n4 =  -64. + 12.*a4 + 4.*a5 + a6 + 64.*eta - 4.*eta2;
  coeffs->n5 = 32. -4.*a4 - a5 - 24.*eta;
  coeffs->d0 = 4.*a4*a4 + 4.*a4*a5 + a5*a5 - a4*a6 + 16.*a6
             + (32.*a4 + 16.*a5 - 8.*a6) * eta + 4.*a4*eta2 + 32.*eta3;
  coeffs->d1 = 4.*a4*a4 + a4*a5 + 16.*a5 + 8.*a6 + (32.*a4 - 2.*a6)*eta + 32.*eta2 + 8.*eta3;
  coeffs->d2 = 16.*a4 + 8.*a5 + 4.*a6 + (8.*a4 + 2.*a5)*eta + 32.*eta2;
  coeffs->d3 = 8.*a4 + 4.*a5 + 2.*a6 + 32.*eta - 8.*eta2;
  coeffs->d4 = 4.*a4 + 2.*a5 + a6 + 16.*eta - 4.*eta2;
  coeffs->d5 = 32. - 4.*a4 - a5 - 24. * eta;

  return XLAL_SUCCESS;
}

/**
 * This function calculates the EOB A function which using the pre-computed
 * coefficients which should already have been calculated.
 */
static
REAL8 XLALCalculateEOBA( const REAL8 r,                     /**<< Orbital separation (in units of total mass M) */
                         EOBACoefficients * restrict coeffs /**<< Pre-computed coefficients for the A function */
                       )
{

  REAL8 r2, r3, r4, r5;
  REAL8 NA, DA;

  /* Note that this function uses pre-computed coefficients,
   * and assumes they have been calculated. Since this is a static function,
   * so only used here, I assume it is okay to neglect error checking
   */

  r2 = r*r;
  r3 = r2 * r;
  r4 = r2*r2;
  r5 = r4*r;


  NA = r4 * coeffs->n4
     + r5 * coeffs->n5;

  DA = coeffs->d0
     + r  * coeffs->d1
     + r2 * coeffs->d2
     + r3 * coeffs->d3
     + r4 * coeffs->d4
     + r5 * coeffs->d5;

  return NA/DA;
}

/**
 * Calculated the derivative of the EOB A function with respect to 
 * r, using the pre-computed A coefficients
 */
static
REAL8 XLALCalculateEOBdAdr( const REAL8 r,                     /**<< Orbital separation (in units of total mass M) */
                            EOBACoefficients * restrict coeffs /**<< Pre-computed coefficients for the A function */
                          )
{
  REAL8 r2, r3, r4, r5;

  REAL8 NA, DA, dNA, dDA, dA;

  r2 = r*r;
  r3 = r2 * r;
  r4 = r2*r2;
  r5 = r4*r;

  NA = r4 * coeffs->n4
     + r5 * coeffs->n5;

  DA = coeffs->d0
     + r  * coeffs->d1
     + r2 * coeffs->d2
     + r3 * coeffs->d3
     + r4 * coeffs->d4
     + r5 * coeffs->d5;

  dNA = 4. * coeffs->n4 * r3
      + 5. * coeffs->n5 * r4;

  dDA = coeffs->d1
      + 2. * coeffs->d2 * r
      + 3. * coeffs->d3 * r2
      + 4. * coeffs->d4 * r3
      + 5. * coeffs->d5 * r4;

  dA = dNA * DA - dDA * NA;

  return dA / (DA*DA);
}

/**
 * Calculate the EOB D function.
 */
static REAL8 XLALCalculateEOBD( REAL8   r, /**<< Orbital separation (in units of total mass M) */
                         REAL8 eta  /**<< Symmetric mass ratio */
                       )
{
	REAL8  u, u2, u3;

	u = 1./r;
	u2 = u*u;
	u3 = u2*u;

	return 1./(1.+6.*eta*u2+2.*eta*(26.-3.*eta)*u3);
}


/**
 * Function to calculate the EOB effective Hamiltonian for the
 * given values of the dynamical variables. The coefficients in the
 * A potential function should already have been computed.
 * Note that the pr used here is the tortoise co-ordinate.
 */
static
REAL8 XLALEffectiveHamiltonian( const REAL8 eta,          /**<< Symmetric mass ratio */
                                const REAL8 r,            /**<< Orbital separation */
                                const REAL8 pr,           /**<< Tortoise co-ordinate */
                                const REAL8 pp,           /**<< Momentum pphi */
                                EOBACoefficients *aCoeffs /**<< Pre-computed coefficients in A function */
                              )
{

        /* The pr used in here is the tortoise co-ordinate */
        REAL8 r2, pr2, pp2, z3, eoba;

        r2   = r * r;
        pr2  = pr * pr;
        pp2  = pp * pp;

        eoba = XLALCalculateEOBA( r, aCoeffs );
        z3   = 2. * ( 4. - 3. * eta ) * eta;
        return sqrt( pr2 + eoba * ( 1.  + pp2/r2 + z3*pr2*pr2/r2 ) );
}


/**
 * Function which calculates the various coefficients used in the generation
 * of the factorized waveform. These coefficients depend only on the symmetric
 * mass ratio eta. It should be noted that this function calculates the 
 * coefficients used in calculating the flux. For generating the waveforms
 * themselves, the coefficients have additional terms added which are calculated
 * using XLALModifyFacWaveformCoefficients(). THe non-spinning parts of these
 * coefficients can be found in Pan et al, arXiv:1106.1021v1 [gr-qc].
 */
UNUSED static int XLALSimIMREOBCalcFacWaveformCoefficients(
          FacWaveformCoeffs * const coeffs, /**<< Structure containing coefficients (populated in function) */
          const REAL8               eta     /**<< Symmetric mass ratio */
          )
{

  REAL8 eta2 = eta*eta;
  REAL8 eta3 = eta2 * eta;

  REAL8 dM, dM2; //dM3;

  REAL8 a = 0;
  REAL8 a2 = 0;
  REAL8 a3 = 0;
  REAL8 chiS = 0;
  REAL8 chiA = 0;

  /* Combination which appears a lot */
  REAL8 m1Plus3eta, m1Plus3eta2, m1Plus3eta3;

  dM2 = 1. - 4.*eta;
  
  /* Check that deltaM has a reasonable value */
  if ( dM2 < 0 )
  {
    XLALPrintError( "eta seems to be < 0.25 - this isn't allowed!\n" );
    XLAL_ERROR( XLAL_EINVAL );
  }

  dM  = sqrt( dM2 );
  //dM3 = dM2 * dM;

  m1Plus3eta  = - 1. + 3.*eta;
  m1Plus3eta2 = m1Plus3eta * m1Plus3eta;
  m1Plus3eta3 = m1Plus3eta * m1Plus3eta2;

  /* Initialize all coefficients to zero */
  /* This is important, as we will not set some if dM is zero */
  memset( coeffs, 0, sizeof( FacWaveformCoeffs ) );


  /* l = 2 */

  coeffs->delta22vh3 = 7./3.;
  coeffs->delta22vh6 = (-4.*a)/3. + (428.*LAL_PI)/105.;
  coeffs->delta22v8 = (20.*a)/63.;
  coeffs->delta22vh9 = -2203./81. + (1712.*LAL_PI*LAL_PI)/315.;
  coeffs->delta22v5  = - 24.*eta;

  coeffs->rho22v2   = -43./42. + (55.*eta)/84.;
  coeffs->rho22v3   = (-2.*(chiS + chiA*dM - chiS*eta))/3.;
  coeffs->rho22v4   = -20555./10584. + (chiS*chiS + 2.*chiA*chiS*dM + chiA*chiA*dM2)/2.
       - (33025.*eta)/21168. + (19583.*eta2)/42336.;
  coeffs->rho22v5   = (-34.*a)/21.;
  coeffs->rho22v6   = 1556919113./122245200. + (89.*a2)/252. - (48993925.*eta)/9779616. 
       - (6292061.*eta2)/3259872. + (10620745.*eta3)/39118464.
       + (41.*eta*LAL_PI*LAL_PI)/192.;
  coeffs->rho22v6l  = - 428./105.;
  coeffs->rho22v7   = (18733.*a)/15876. + a*a2/3.;
  coeffs->rho22v8   = -387216563023./160190110080. + (18353.*a2)/21168. - a2*a2/8.;
  coeffs->rho22v8l  =  9202./2205.;
  coeffs->rho22v10  = -16094530514677./533967033600.;
  coeffs->rho22v10l =  439877./55566.;


  if ( dM2 )
  {
    coeffs->delta21vh3 = 2./3.;
    coeffs->delta21vh6 = (-17.*a)/35. + (107.*LAL_PI)/105.;
    coeffs->delta21vh7 = (3.*a2)/140.;
    coeffs->delta21vh9 = -272./81. + (214.*LAL_PI*LAL_PI)/315.;
    coeffs->delta21v5  = - 493. * eta /42.;

    coeffs->rho21v1   = (-3.*(chiS+chiA/dM))/(4.);
    //coeffs->rho21v2   = -59./56 - (9.*chiAPlusChiSdM*chiAPlusChiSdM)/(32.*dM2) + (23.*eta)/84.;
    /*coeffs->rho21v3   = (-567.*chiA*chiA*chiA - 1701.*chiA*chiA*chiS*dM
                        + chiA*(-4708. + 1701.*chiS*chiS - 2648.*eta)*(-1. + 4.*eta)
                        + chiS* dM3 *(4708. - 567.*chiS*chiS
                        + 1816.*eta))/(2688.*dM3);*/
    coeffs->rho21v2   = -59./56. + (23.*eta)/84. - 9./32.*a2;
    coeffs->rho21v3   = 1177./672.*a - 27./128.*a3;
    coeffs->rho21v4   = -47009./56448.- (865.*a2)/1792. - (405.*a2*a2)/2048. - (10993.*eta)/14112.
                        + (617.*eta2)/4704.;
    coeffs->rho21v5   = (-98635.*a)/75264. + (2031.*a*a2)/7168. - (1701.*a2*a3)/8192.;
    coeffs->rho21v6   = 7613184941./2607897600.+ (9032393.*a2)/1806336. + (3897.*a2*a2)/16384.
                        - (15309.*a3*a3)/65536.; 
    coeffs->rho21v6l  = - 107./105.;
    coeffs->rho21v7   = (-3859374457.*a)/1159065600. - (55169.*a3)/16384.
                        + (18603.*a2*a3)/65536. - (72171.*a2*a2*a3)/262144.;
    coeffs->rho21v7l  =  107.*a/140.;
    coeffs->rho21v8   = -1168617463883./911303737344.;
    coeffs->rho21v8l  = 6313./5880.;
    coeffs->rho21v10  = -63735873771463./16569158860800.; 
    coeffs->rho21v10l = 5029963./5927040.;
  }

  /* l = 3 */
  if ( dM2 )
  {
    coeffs->delta33vh3 = 13./10.;
    coeffs->delta33vh6 = (-81.*a)/20. + (39.*LAL_PI)/7.;
    coeffs->delta33vh9 = -227827./3000. + (78.*LAL_PI*LAL_PI)/7.;
    coeffs->delta33v5  = - 80897.*eta / 2430.;

    coeffs->rho33v2 = -7./6. + (2.*eta)/3.;
    coeffs->rho33v3 = (chiS*dM*(-4. + 5.*eta) + chiA*(-4. + 19.*eta))/(6.*dM);
    coeffs->rho33v4 = -6719./3960. + a2/2. - (1861.*eta)/990. + (149.*eta2)/330.;
    coeffs->rho33v5 = (-4.*a)/3.;
    coeffs->rho33v6 = 3203101567./227026800. + (5.*a2)/36.;
    coeffs->rho33v6l = - 26./7.;
    coeffs->rho33v7 = (5297.*a)/2970. + a*a2/3.;
    coeffs->rho33v8 = -57566572157./8562153600.;
    coeffs->rho33v8l = 13./3.;
  }

  coeffs->delta32vh3 = (10. + 33.*eta)/(-15.*m1Plus3eta);
  coeffs->delta32vh4 = 4.*a;
  coeffs->delta32vh6 = (-136.*a)/45. + (52.*LAL_PI)/21.;
  coeffs->delta32vh9 = -9112./405. + (208.*LAL_PI*LAL_PI)/63.;

  coeffs->rho32v   = (4.*chiS*eta)/(-3.*m1Plus3eta);
  coeffs->rho32v2  = (-4.*a2*eta2)/(9.*m1Plus3eta2) + (328. - 1115.*eta
                        + 320.*eta2)/(270.*m1Plus3eta);
  coeffs->rho32v3  = (2.*(45.*a*m1Plus3eta3
                        - a*eta*(328. - 2099.*eta + 5.*(733. + 20.*a2)*eta2
                        - 960.*eta3)))/(405.*m1Plus3eta3);
  coeffs->rho32v4  = a2/3. + (-1444528.
                        + 8050045.*eta - 4725605.*eta2 - 20338960.*eta3
                        + 3085640.*eta2*eta2)/(1603800.*m1Plus3eta2);
  coeffs->rho32v5  = (-2788.*a)/1215.;
  coeffs->rho32v6  = 5849948554./940355325. + (488.*a2)/405.;
  coeffs->rho32v6l =  - 104./63.;
  coeffs->rho32v8  = -10607269449358./3072140846775.;
  coeffs->rho32v8l = 17056./8505.;

  if ( dM2 )
  {
    coeffs->delta31vh3 = 13./30.;
    coeffs->delta31vh6 = (61.*a)/20. + (13.*LAL_PI)/21.;
    coeffs->delta31vh7 = (-24.*a2)/5.;
    coeffs->delta31vh9 = -227827./81000. + (26.*LAL_PI*LAL_PI)/63.;
    coeffs->delta31v5  = - 17.*eta/10.; 
 
    coeffs->rho31v2  = -13./18. - (2.*eta)/9.;
    coeffs->rho31v3  = (chiA*(-4. + 11.*eta) + chiS*dM*(-4. + 13.*eta))/(6.*dM);
    coeffs->rho31v4  = 101./7128.
                        - (5.*a2)/6. - (1685.*eta)/1782. - (829.*eta2)/1782.;
    coeffs->rho31v5  = (4.*a)/9.;
    coeffs->rho31v6  = 11706720301./6129723600. - (49.*a2)/108.;
    coeffs->rho31v6l =  - 26./63.;
    coeffs->rho31v7  = (-2579.*a)/5346. + a*a2/9.;
    coeffs->rho31v8  = 2606097992581./4854741091200.;
    coeffs->rho31v8l = 169./567.;
  }

  /* l = 4 */
  
  coeffs->delta44vh3 = (112. + 219.*eta)/(-120.*m1Plus3eta);
  coeffs->delta44vh6 = (-464.*a)/75. + (25136.*LAL_PI)/3465.;

  coeffs->rho44v2 = (1614. - 5870.*eta + 2625.*eta2)/(1320.*m1Plus3eta);
  coeffs->rho44v3 = (chiA*(10. - 39.*eta)*dM + chiS*(10. - 41.*eta
                        + 42.*eta2))/(15.*m1Plus3eta);
  coeffs->rho44v4 = a2/2. + (-511573572.
                        + 2338945704.*eta - 313857376.*eta2 - 6733146000.*eta3
                        + 1252563795.*eta2*eta2)/(317116800.*m1Plus3eta2);
  coeffs->rho44v5 = (-69.*a)/55.;
  coeffs->rho44v6 = 16600939332793./1098809712000. + (217.*a2)/3960.;
  coeffs->rho44v6l = - 12568./3465.;

  if ( dM2 )
  {
    coeffs->delta43vh3 = (486. + 4961.*eta)/(810.*(1. - 2.*eta));
    coeffs->delta43vh4 = (11.*a)/4.;
    coeffs->delta43vh6 = 1571.*LAL_PI/385.;

    coeffs->rho43v   = (5.*(chiA - chiS*dM)*eta)/(8.*dM*(-1. + 2.*eta));
    coeffs->rho43v2  = (222. - 547.*eta + 160.*eta2)/(176.*(-1. + 2.*eta));
    coeffs->rho43v4  = -6894273./7047040. + (3.*a2)/8.;
    coeffs->rho43v5  = (-12113.*a)/6160.;
    coeffs->rho43v6  = 1664224207351./195343948800.;
    coeffs->rho43v6l = - 1571./770.;
  }

  coeffs->delta42vh3 = (7.*(1. + 6.*eta))/(-15.*m1Plus3eta);
  coeffs->delta42vh6 = (212.*a)/75. + (6284.*LAL_PI)/3465.;

  coeffs->rho42v2  = (1146. - 3530.*eta + 285.*eta2)/(1320.*m1Plus3eta);
  coeffs->rho42v3  = (chiA*(10. - 21.*eta)*dM + chiS*(10. - 59.*eta
                        + 78.*eta2))/(15.*m1Plus3eta);
  coeffs->rho42v4  = a2/2. + (-114859044. + 295834536.*eta + 1204388696.*eta2 - 3047981160.*eta3
                        - 379526805.*eta2*eta2)/(317116800.*m1Plus3eta2);
  coeffs->rho42v5  = (-7.*a)/110.;
  coeffs->rho42v6  = 848238724511./219761942400. + (2323.*a2)/3960.;
  coeffs->rho42v6l = - 3142./3465.;

  if ( dM2 )
  {
    coeffs->delta41vh3 = (2. + 507.*eta)/(10.*(1. - 2.*eta));
    coeffs->delta41vh4 = (11.*a)/12.;
    coeffs->delta41vh6 = 1571.*LAL_PI/3465.;

    coeffs->rho41v   = (5.*(chiA - chiS*dM)*eta)/(8.*dM*(-1. + 2.*eta));
    coeffs->rho41v2  = (602. - 1385.*eta + 288.*eta2)/(528.*(-1. + 2.*eta));
    coeffs->rho41v4  = -7775491./21141120. + (3.*a2)/8.;
    coeffs->rho41v5  = (-20033.*a)/55440. - (5*a*a2)/6.;
    coeffs->rho41v6  = 1227423222031./1758095539200.;
    coeffs->rho41v6l = - 1571./6930.;
  }

  /* l = 5 */
  if ( dM2 )
  {
    coeffs->delta55vh3 = (96875. + 857528.*eta)/(131250.*(1 - 2*eta));

    coeffs->rho55v2 = (487. - 1298.*eta + 512.*eta2)/(390.*(-1. + 2.*eta));
    coeffs->rho55v3 = (-2.*a)/3.;
    coeffs->rho55v4 = -3353747./2129400. + a2/2.;
    coeffs->rho55v5 = - 241. * a / 195.;
  }

  coeffs->delta54vh3 = 8./15.;
  coeffs->delta54vh4 = 12.*a/5.;

  coeffs->rho54v2 = (-17448. + 96019.*eta - 127610.*eta2
                        + 33320.*eta3)/(13650.*(1. - 5.*eta + 5.*eta2));
  coeffs->rho54v3 = (-2.*a)/15.;
  coeffs->rho54v4 = -16213384./15526875. + (2.*a2)/5.;

  if ( dM2 )
  {
    coeffs->delta53vh3 = 31./70.;

    coeffs->rho53v2 = (375. - 850.*eta + 176.*eta2)/(390.*(-1. + 2.*eta));
    coeffs->rho53v3 = (-2.*a)/3.;
    coeffs->rho53v4 = -410833./709800. + a2/2.;
    coeffs->rho53v5 = - 103.*a/325.;
  }

  coeffs->delta52vh3 = 4./15.;
  coeffs->delta52vh4 = 6.*a/5.;

  coeffs->rho52v2 = (-15828. + 84679.*eta - 104930.*eta2
                        + 21980.*eta3)/(13650.*(1. - 5.*eta + 5.*eta2));
  coeffs->rho52v3 = (-2.*a)/15.;
  coeffs->rho52v4 = -7187914./15526875. + (2.*a2)/5.;

  if ( dM2 )
  {
    coeffs->delta51vh3 = 31./210.;

    coeffs->rho51v2 = (319. - 626.*eta + 8.*eta2)/(390.*(-1. + 2.*eta));
    coeffs->rho51v3 = (-2.*a)/3.;
    coeffs->rho51v4 = -31877./304200. + a2/2.;
    coeffs->rho51v5 = 139.*a/975.;
  }

  /* l = 6 */

  coeffs->delta66vh3 = 43./70.;
  
  coeffs->rho66v2 = (-106. + 602.*eta - 861.*eta2
                        + 273.*eta3)/(84.*(1. - 5.*eta + 5.*eta2));
  coeffs->rho66v3 = (-2.*a)/3.;
  coeffs->rho66v4 = -1025435./659736. + a2/2.;

  if ( dM2 )
  {
    coeffs->delta65vh3 = 10./21.;
    
    coeffs->rho65v2 = (-185. + 838.*eta - 910.*eta2
                        + 220.*eta3)/(144.*(dM2 + 3.*eta2));
    coeffs->rho65v3 = - 2.*a/9.;
  }

  coeffs->delta64vh3 = 43./105.;
  
  coeffs->rho64v2 = (-86. + 462.*eta - 581.*eta2
                        + 133.*eta3)/(84.*(1. - 5.*eta + 5.*eta2));
  coeffs->rho64v3 = (-2.*a)/3.;
  coeffs->rho64v4 = -476887./659736. + a2/2.;

  if ( dM2 )
  {
    coeffs->delta63vh3 = 2./7.;

    coeffs->rho63v2 = (-169. + 742.*eta - 750.*eta2
                        + 156.*eta3)/(144.*(dM2 + 3.*eta2));
    coeffs->rho63v3 = - 2.*a/9.;
  }

  coeffs->delta62vh3 = 43./210.;

  coeffs->rho62v2 = (-74. + 378.*eta - 413.*eta2
                        + 49.*eta3)/(84.*(1. - 5.*eta + 5.*eta2));
  coeffs->rho62v3 = (-2.*a)/3.;
  coeffs->rho62v4 = -817991./3298680. + a2/2.;

  if ( dM2 )
  {
    coeffs->delta61vh3 = 2./21.;

    coeffs->rho61v2 = (-161. + 694.*eta - 670.*eta2
                        + 124.*eta3)/(144.*(dM2 + 3.*eta2));
    coeffs->rho61v3 = - 2. * a / 9.;
  }

  /* l = 7 */
  if ( dM2 )
  {
    coeffs->delta77vh3 = 19./36.;

    coeffs->rho77v2 = (-906. + 4246.*eta - 4963.*eta2
                        + 1380.*eta3)/(714.*(dM2 + 3.*eta2));
    coeffs->rho77v3 = - 2.*a/3.;
  }

  coeffs->rho76v2 = (2144. - 16185.*eta + 37828.*eta2 - 29351.*eta3
                        + 6104.*eta2*eta2) / (1666.*(-1 + 7*eta - 14*eta2
                        + 7*eta3));

  if ( dM2 )
  {
    coeffs->delta75vh3 = 95./252.;

    coeffs->rho75v2 = (-762. + 3382.*eta - 3523.*eta2
                        + 804.*eta3)/(714.*(dM2 + 3.*eta2));
    coeffs->rho75v3 = - 2.*a/3.;
  }

  coeffs->rho74v2 = (17756. - 131805.*eta + 298872.*eta2 - 217959.*eta3
                        + 41076.*eta2*eta2) / (14994.*(-1. + 7.*eta - 14.*eta2
                        + 7.*eta3));

  if ( dM2 )
  {
    coeffs->delta73vh3 = 19./84.;

    coeffs->rho73v2 = (-666. + 2806.*eta - 2563.*eta2
                        + 420.*eta3)/(714.*(dM2 + 3.*eta2));
    coeffs->rho73v3 = - 2.*a/3.;
  }

  coeffs->rho72v2 = (16832. - 123489.*eta + 273924.*eta2 - 190239.*eta3
                        + 32760.*eta2*eta2) /(14994.*(-1. + 7.*eta - 14.*eta2
                        + 7.*eta3));

  if ( dM2 )
  {
    coeffs->delta71vh3 = 19./252.;

    coeffs->rho71v2 = (-618. + 2518.*eta - 2083.*eta2
                        + 228.*eta3)/(714.*(dM2 + 3.*eta2));
    coeffs->rho71v3 = - 2.*a/3.;
  }

  /* l = 8 */
  
  coeffs->rho88v2 = (3482. - 26778.*eta + 64659.*eta2 - 53445.*eta3
                        + 12243.*eta2*eta2) / (2736.*(-1. + 7.*eta - 14.*eta2
                        + 7.*eta3));

  if ( dM2 )
  {
    coeffs->rho87v2 = (23478. - 154099.*eta + 309498.*eta2 - 207550.*eta3
                        + 38920*eta2*eta2) / (18240.*(-1 + 6*eta - 10*eta2
                        + 4*eta3));
  }

  coeffs->rho86v2 = (1002. - 7498.*eta + 17269.*eta2 - 13055.*eta3
                        + 2653.*eta2*eta2) / (912.*(-1. + 7.*eta - 14.*eta2
                        + 7.*eta3));

  if ( dM2 )
  {
    coeffs->rho85v2 = (4350. - 28055.*eta + 54642.*eta2 - 34598.*eta3
                        + 6056.*eta2*eta2) / (3648.*(-1. + 6.*eta - 10.*eta2
                        + 4.*eta3));
  }

  coeffs->rho84v2 = (2666. - 19434.*eta + 42627.*eta2 - 28965.*eta3
                        + 4899.*eta2*eta2) / (2736.*(-1. + 7.*eta - 14.*eta2
                        + 7.*eta3));

  if ( dM2 )
  {
    coeffs->rho83v2 = (20598. - 131059.*eta + 249018.*eta2 - 149950.*eta3
                        + 24520.*eta2*eta2) / (18240.*(-1. + 6.*eta - 10.*eta2
                        + 4.*eta3));
  }

  coeffs->rho82v2 = (2462. - 17598.*eta + 37119.*eta2 - 22845.*eta3
                        + 3063.*eta2*eta2) / (2736.*(-1. + 7.*eta - 14.*eta2
                        + 7.*eta3));

  if ( dM2 )
  {
    coeffs->rho81v2 = (20022. - 126451.*eta + 236922.*eta2 - 138430.*eta3
                        + 21640.*eta2*eta2) / (18240.*(-1. + 6.*eta - 10.*eta2
                        + 4.*eta3));
  }

  /* All relevant coefficients should be set, so we return */

  return XLAL_SUCCESS;
}


/**
 * Function which adds the additional terms required for waveform generation
 * to the factorized waveform coefficients. Note that this function only calculates
 * additional terms not present in the flux, so the factorized waveform coefficients
 * SHOULD ALREADY HAVE BEEN CALCULATED using XLALCalcFacWaveformCoefficients() prior
 * to calling this function.
 */
UNUSED static int XLALSimIMREOBModifyFacWaveformCoefficients( 
                                       FacWaveformCoeffs * const coeffs, /**<< Structure containing coefficients */
                                       const REAL8 eta                   /**<< Symmetric mass ratio */
                                     )
{

  if ( !coeffs )
  {
    XLAL_ERROR( XLAL_EINVAL );
  }

  /* Tweak the relevant coefficients for the generation of the waveform */
  coeffs->rho21v6 += -5. * eta;
  coeffs->rho33v6 += -20. * eta;
  coeffs->rho44v6 += -15. * eta;
  coeffs->rho55v6 += 4. * eta;

  coeffs->delta21v7 += 30. * eta;
  coeffs->delta33v7 += -10. * eta;
  coeffs->delta44v5 += -70. * eta;
  coeffs->delta55v5 += 40. * eta;

  return XLAL_SUCCESS;
}

/**
 * Computes the non-Keplerian correction to the velocity as determined from the
 * frequency obtained assuming a circular orbit. In the early stages of the evolution,
 * this should be a number close to 1.
 */
static REAL8
nonKeplerianCoefficient(
                   REAL8Vector * restrict values, /**<< Dynamics r, phi, pr, pphi */
                   const REAL8       eta,         /**<< Symmetric mass ratio */
                   EOBACoefficients *coeffs       /**<< Pre-computed A coefficients */
                   )
{

  REAL8 r    = values->data[0];
  REAL8 pphi = values->data[3];

  REAL8 A  = XLALCalculateEOBA( r, coeffs );
  REAL8 dA = XLALCalculateEOBdAdr( r, coeffs );

  return 2. * (1. + 2. * eta * ( -1. + sqrt( (1. + pphi*pphi/(r*r)) * A ) ) )
          / ( r*r * dA );
}

/**
 * Computes the factorized waveform according to the prescription
 * given in Pan et al, arXiv:1106.1021v1 [gr-qc], for a given
 * mode l,m, for the given values of the dynamics at that point.
 * The function returns XLAL_SUCCESS if everything works out properly,
 * otherwise XLAL_FAILURE will be returned.
 */
UNUSED static int  XLALSimIMREOBGetFactorizedWaveform( 
                                COMPLEX16   * restrict hlm,    /**<< The value of hlm (populated by the function) */
                                REAL8Vector * restrict values, /**<< Vector containing dynamics r, phi, pr, pphi for a given point */
                                const REAL8 v,                 /**<< Velocity (in geometric units) */
                                const INT4  l,                 /**<< Mode l */
                                const INT4  m,                 /**<< Mode m */
                                EOBParams   * restrict params  /**<< Structure containing pre-computed coefficients, etc. */
                                )
{

  /* Status of function calls */
  INT4 status;
  INT4 i;

  REAL8 eta;
  REAL8 r, pr, pp, Omega, v2, vh, vh3, k, hathatk, eulerlogxabs;
  REAL8 Hreal, Heff, Slm, deltalm, rholm, rholmPwrl;
  COMPLEX16 Tlm;
  COMPLEX16 hNewton;
  gsl_sf_result lnr1, arg1, z2;

  /* Non-Keplerian velocity */
  REAL8 vPhi;

  /* Pre-computed coefficients */
  FacWaveformCoeffs *hCoeffs = params->hCoeffs;

  if ( abs(m) > (INT4) l )
  {
    XLAL_ERROR( XLAL_EINVAL );
  }


  eta = params->eta;

  /* Check our eta was sensible */
  if ( eta > 0.25 )
  {
    XLALPrintError("Eta seems to be > 0.25 - this isn't allowed!\n" );
    XLAL_ERROR( XLAL_EINVAL );
  }
  else if ( eta == 0.25 && m % 2 )
  {
    /* If m is odd and dM = 0, hLM will be zero */
    memset( hlm, 0, sizeof( COMPLEX16 ) );
    return XLAL_SUCCESS;
  }

  r  = values->data[0];
  pr = values->data[2];
  pp = values->data[3];

  Heff  = XLALEffectiveHamiltonian( eta, r, pr, pp, params->aCoeffs );
  Hreal = sqrt( 1.0 + 2.0 * eta * ( Heff - 1.0) );
  v2    = v * v;
  Omega = v2 * v;
  vh3   = Hreal * Omega;
  vh    = cbrt(vh3);
  eulerlogxabs = LAL_GAMMA + log( 2.0 * (REAL8)m * v );


  /* Calculate the non-Keplerian velocity */
  /* given by Eq. (18) of Pan et al, PRD84, 124052(2011) */
  /* psi given by Eq. (19) of Pan et al, PRD84, 124052(2011) */
  /* Assign temporarily to vPhi */
  vPhi = nonKeplerianCoefficient( values, eta, params->aCoeffs );
  /* Assign rOmega value temporarily to vPhi */
  vPhi  = r * cbrt(vPhi);
  /* Assign rOmega * Omega to vPhi */
  vPhi *= Omega;

  /* Calculate the newtonian multipole */
  status = XLALSimIMREOBCalculateNewtonianMultipole( &hNewton, vPhi * vPhi, vPhi/Omega,
            values->data[1], (UINT4)l, m, params );
  if ( status == XLAL_FAILURE )
  {
    XLAL_ERROR( XLAL_EFUNC );
  }

  /* Calculate the source term */
  if ( ( (l+m)%2 ) == 0)
  {
    Slm = Heff;
  }
  else
  {
    Slm = v * pp;
  }

  /* Calculate the Tail term */
  k  = m * Omega;
  hathatk = Hreal * k;
  if ( fabs(hathatk) <= 0.18*l ) {
    /* Use the polynomial fit to the Gamma function along with the identity for adding a positive integer within the argument */
    GammaFunctionPolyFit(l,hathatk,&(lnr1.val),&(arg1.val));
  } else {
    XLAL_CALLGSL( status = gsl_sf_lngamma_complex_e( l+1.0, -2.0*hathatk, &lnr1, &arg1 ) );
    if (status != GSL_SUCCESS)
    {
      XLALPrintError("Error in GSL function\n" );
      XLAL_ERROR( XLAL_EFUNC );
    }
  }
  XLAL_CALLGSL( status = gsl_sf_fact_e( l, &z2 ) );
  if ( status != GSL_SUCCESS)
  {
    XLALPrintError("Error in GSL function\n" );
    XLAL_ERROR( XLAL_EFUNC );
  }
  Tlm = cexp( ( lnr1.val + LAL_PI * hathatk ) + I * (
        arg1.val + 2.0 * hathatk * log(4.0*k/sqrt(LAL_E)) ) );
  Tlm /= z2.val;

  /* Calculate the residue phase and amplitude terms */
  switch( l )
  {
    case 2:
      switch( abs(m) )
      {
        case 2:
          deltalm = vh3*(hCoeffs->delta22vh3 + vh3*(hCoeffs->delta22vh6
            + vh*vh*(hCoeffs->delta22vh9*vh)))
            + hCoeffs->delta22v5 *v*v2*v2 + hCoeffs->delta22v8 *v2*v2*v2*v2;
          rholm  = 1. + v2*(hCoeffs->rho22v2 + v*(hCoeffs->rho22v3
            + v*(hCoeffs->rho22v4
            + v*(hCoeffs->rho22v5 + v*(hCoeffs->rho22v6
            + hCoeffs->rho22v6l*eulerlogxabs + v*(hCoeffs->rho22v7
            + v*(hCoeffs->rho22v8 + hCoeffs->rho22v8l*eulerlogxabs
            + (hCoeffs->rho22v10 + hCoeffs->rho22v10l * eulerlogxabs)*v2)))))));
          break;
        case 1:
          deltalm = vh3*(hCoeffs->delta21vh3 + vh3*(hCoeffs->delta21vh6
            + vh*(hCoeffs->delta21vh7 + (hCoeffs->delta21vh9)*vh*vh)))
            + hCoeffs->delta21v5*v*v2*v2 + hCoeffs->delta21v7*v2*v2*v2*v;
          rholm  = 1. + v*(hCoeffs->rho21v1
            + v*( hCoeffs->rho21v2 + v*(hCoeffs->rho21v3 + v*(hCoeffs->rho21v4
            + v*(hCoeffs->rho21v5 + v*(hCoeffs->rho21v6 + hCoeffs->rho21v6l*eulerlogxabs
            + v*(hCoeffs->rho21v7 + hCoeffs->rho21v7l * eulerlogxabs
            + v*(hCoeffs->rho21v8 + hCoeffs->rho21v8l * eulerlogxabs
            + (hCoeffs->rho21v10 + hCoeffs->rho21v10l * eulerlogxabs)*v2))))))));
          break;
        default:
          XLAL_ERROR( XLAL_EINVAL );
          break;
      }
      break;
    case 3:
      switch (m)
      {
        case 3:
          deltalm = vh3*(hCoeffs->delta33vh3 + vh3*(hCoeffs->delta33vh6 + hCoeffs->delta33vh9*vh3))
            + hCoeffs->delta33v5*v*v2*v2 + hCoeffs->delta33v7*v2*v2*v2*v;
          rholm  = 1. + v2*(hCoeffs->rho33v2 + v*(hCoeffs->rho33v3 + v*(hCoeffs->rho33v4
            + v*(hCoeffs->rho33v5 + v*(hCoeffs->rho33v6 + hCoeffs->rho33v6l*eulerlogxabs
            + v*(hCoeffs->rho33v7 + (hCoeffs->rho33v8 + hCoeffs->rho33v8l*eulerlogxabs)*v))))));
          break;
        case 2:
          deltalm = vh3*(hCoeffs->delta32vh3 + vh*(hCoeffs->delta32vh4 + vh*vh*(hCoeffs->delta32vh6
            + hCoeffs->delta32vh9*vh3)));
          rholm  = 1. + v*(hCoeffs->rho32v
            + v*(hCoeffs->rho32v2 + v*(hCoeffs->rho32v3 + v*(hCoeffs->rho32v4 + v*(hCoeffs->rho32v5
            + v*(hCoeffs->rho32v6 + hCoeffs->rho32v6l*eulerlogxabs
            + (hCoeffs->rho32v8 + hCoeffs->rho32v8l*eulerlogxabs)*v2))))));
          break;
        case 1:
          deltalm = vh3*(hCoeffs->delta31vh3 + vh3*(hCoeffs->delta31vh6
            + vh*(hCoeffs->delta31vh7 + hCoeffs->delta31vh9*vh*vh)))
            + hCoeffs->delta31v5*v*v2*v2;
          rholm  = 1. + v2*(hCoeffs->rho31v2 + v*(hCoeffs->rho31v3 + v*(hCoeffs->rho31v4
            + v*(hCoeffs->rho31v5 + v*(hCoeffs->rho31v6 + hCoeffs->rho31v6l*eulerlogxabs
            + v*(hCoeffs->rho31v7 + (hCoeffs->rho31v8 + hCoeffs->rho31v8l*eulerlogxabs)*v))))));
          break;
        default:
          XLAL_ERROR( XLAL_EINVAL );
          break;
      }
      break;
    case 4:
      switch (m)
      {
        case 4:
          deltalm = vh3*(hCoeffs->delta44vh3 + hCoeffs->delta44vh6 *vh3)
            + hCoeffs->delta44v5*v2*v2*v;
          rholm  = 1. + v2*(hCoeffs->rho44v2
            + v*( hCoeffs->rho44v3 + v*(hCoeffs->rho44v4
            + v*(hCoeffs->rho44v5 + (hCoeffs->rho44v6
            + hCoeffs->rho44v6l*eulerlogxabs)*v))));
          break;
        case 3:
          deltalm = vh3*(hCoeffs->delta43vh3 + vh*(hCoeffs->delta43vh4
            + hCoeffs->delta43vh6*vh*vh));
          rholm  = 1. + v*(hCoeffs->rho43v
            + v*(hCoeffs->rho43v2
            + v2*(hCoeffs->rho43v4 + v*(hCoeffs->rho43v5
            + (hCoeffs->rho43v6 + hCoeffs->rho43v6l*eulerlogxabs)*v))));
          break;
        case 2:
          deltalm = vh3*(hCoeffs->delta42vh3 + hCoeffs->delta42vh6*vh3);
          rholm  = 1. + v2*(hCoeffs->rho42v2
            + v*(hCoeffs->rho42v3 + v*(hCoeffs->rho42v4 + v*(hCoeffs->rho42v5
            + (hCoeffs->rho42v6 + hCoeffs->rho42v6l*eulerlogxabs)*v))));
          break;
        case 1:
          deltalm = vh3*(hCoeffs->delta41vh3 + vh*(hCoeffs->delta41vh4
            + hCoeffs->delta41vh6*vh*vh));
          rholm  = 1. + v*(hCoeffs->rho41v
            + v*(hCoeffs->rho41v2
            + v2*(hCoeffs->rho41v4 + v*(hCoeffs->rho41v5
            + (hCoeffs->rho41v6 +  hCoeffs->rho41v6l*eulerlogxabs)*v))));
          break;
        default:
          XLAL_ERROR( XLAL_EINVAL );
          break;
      }
      break;
    case 5:
      switch (m)
      {
        case 5:
          deltalm = hCoeffs->delta55vh3*vh3 + hCoeffs->delta55v5*v2*v2*v;
          rholm  = 1. + v2*( hCoeffs->rho55v2
            + v*(hCoeffs->rho55v3 + v*(hCoeffs->rho55v4
            + v*(hCoeffs->rho55v5 + hCoeffs->rho55v6*v))));
          break;
        case 4:
          deltalm = vh3*(hCoeffs->delta54vh3 + hCoeffs->delta54vh4*vh);
          rholm  = 1. + v2*(hCoeffs->rho54v2 + v*(hCoeffs->rho54v3
            + hCoeffs->rho54v4*v));
          break;
        case 3:
          deltalm = hCoeffs->delta53vh3 * vh3;
          rholm  = 1. + v2*(hCoeffs->rho53v2
            + v*(hCoeffs->rho53v3 + v*(hCoeffs->rho53v4 + hCoeffs->rho53v5*v)));
          break;
        case 2:
          deltalm = vh3*(hCoeffs->delta52vh3 + hCoeffs->delta52vh4*vh);
          rholm  = 1. + v2*(hCoeffs->rho52v2 + v*(hCoeffs->rho52v3
            + hCoeffs->rho52v4*v));
          break;
        case 1:
          deltalm = hCoeffs->delta51vh3 * vh3;
          rholm  = 1. + v2*(hCoeffs->rho51v2
            + v*(hCoeffs->rho51v3 + v*(hCoeffs->rho51v4 + hCoeffs->rho51v5*v)));
          break;
        default:
          XLAL_ERROR( XLAL_EINVAL );
          break;
      }
      break;
    case 6:
      switch (m)
      {
        case 6:
          deltalm = hCoeffs->delta66vh3*vh3;
          rholm  = 1. + v2*(hCoeffs->rho66v2 + v*(hCoeffs->rho66v3
            + hCoeffs->rho66v4*v));
          break;
        case 5:
          deltalm = hCoeffs->delta65vh3*vh3;
          rholm  = 1. + v2*(hCoeffs->rho65v2 + hCoeffs->rho65v3*v);
          break;
        case 4:
          deltalm = hCoeffs->delta64vh3 * vh3;
          rholm  = 1. + v2*(hCoeffs->rho64v2 + v*(hCoeffs->rho64v3
            + hCoeffs->rho64v4*v));
          break;
        case 3:
          deltalm = hCoeffs->delta63vh3 * vh3;
          rholm  = 1. + v2*(hCoeffs->rho63v2 + hCoeffs->rho63v3*v);
          break;
        case 2:
          deltalm = hCoeffs->delta62vh3 * vh3;
          rholm  = 1. + v2*(hCoeffs->rho62v2 + v*(hCoeffs->rho62v3
            + hCoeffs->rho62v4 * v));
          break;
        case 1:
          deltalm = hCoeffs->delta61vh3 * vh3;
          rholm  = 1. + v2*(hCoeffs->rho61v2 + hCoeffs->rho61v3*v);
          break;
        default:
          XLAL_ERROR( XLAL_EINVAL );
          break;
      }
      break;
    case 7:
      switch (m)
      {
        case 7:
          deltalm = hCoeffs->delta77vh3 * vh3;
          rholm   = 1. + v2*(hCoeffs->rho77v2 + hCoeffs->rho77v3 * v);
          break;
        case 6:
          deltalm = 0.0;
          rholm   = 1. + hCoeffs->rho76v2 * v2;
          break;
        case 5:
          deltalm = hCoeffs->delta75vh3 * vh3;
          rholm   = 1. + v2*(hCoeffs->rho75v2 + hCoeffs->rho75v3*v);
          break;
        case 4:
          deltalm = 0.0;
          rholm   = 1. + hCoeffs->rho74v2 * v2;
          break;
        case 3:
          deltalm = hCoeffs->delta73vh3 *vh3;
          rholm   = 1. + v2*(hCoeffs->rho73v2 + hCoeffs->rho73v3 * v);
          break;
        case 2:
          deltalm = 0.0;
          rholm   = 1. + hCoeffs->rho72v2 * v2;
          break;
        case 1:
          deltalm = hCoeffs->delta71vh3 * vh3;
          rholm   = 1. + v2*(hCoeffs->rho71v2 +hCoeffs->rho71v3 * v);
          break;
        default:
          XLAL_ERROR( XLAL_EINVAL );
          break;
      }
      break;
    case 8:
      switch (m)
      {
        case 8:
          deltalm = 0.0;
          rholm   = 1. + hCoeffs->rho88v2 * v2;
          break;
        case 7:
          deltalm = 0.0;
          rholm   = 1. + hCoeffs->rho87v2 * v2;
          break;
        case 6:
          deltalm = 0.0;
          rholm   = 1. + hCoeffs->rho86v2 * v2;
          break;
        case 5:
          deltalm = 0.0;
          rholm   = 1. + hCoeffs->rho85v2 * v2;
          break;
        case 4:
          deltalm = 0.0;
          rholm  = 1. + hCoeffs->rho84v2 * v2;
          break;
        case 3:
          deltalm = 0.0;
          rholm  = 1. + hCoeffs->rho83v2 * v2;
          break;
        case 2:
          deltalm = 0.0;
          rholm  = 1. + hCoeffs->rho82v2 * v2;
          break;
        case 1:
          deltalm = 0.0;
          rholm  = 1. + hCoeffs->rho81v2 * v2;
          break;
        default:
          XLAL_ERROR( XLAL_EINVAL );
          break;
      }
      break;
    default:
      XLAL_ERROR( XLAL_EINVAL );
      break;
  }

  /* Raise rholm to the lth power */
  rholmPwrl = 1.0;
  i = l;
  while ( i-- )
  {
    rholmPwrl *= rholm;
  }

  *hlm = Tlm * cexp(I * deltalm) * Slm * rholmPwrl;
  *hlm *= hNewton;

  return XLAL_SUCCESS;
} 

UNUSED static int GetGammaFuncPolyFitCoeffs(INT4 l, REAL8 *rc, REAL8 *ic, INT4 *n)
{
  switch (l)
  {
    case 2:
      *n=14;
      rc[0]=1.99999999999999955591e+00;  ic[0]=3.22608382017814060065e-17;
      rc[1]=2.02673319889760151182e-15;  ic[1]=-1.84556867019704529120e+00;
      rc[2]=-1.24646499595122994819e+00;  ic[2]=2.38602658260296928909e-15;
      rc[3]=-2.28093933212646085710e-14;  ic[3]=5.74994168928716131717e-01;
      rc[4]=2.30074940748511563848e-01;  ic[4]=-1.73122756583947918203e-13;
      rc[5]=-3.29706756505951927363e-13;  ic[5]=-7.37150467865764214004e-02;
      rc[6]=-2.20411092587177968871e-02;  ic[6]=2.57904962704984623933e-12;
      rc[7]=4.32434587204944421890e-12;  ic[7]=5.44875565760795897707e-03;
      rc[8]=1.35521988214179120967e-03;  ic[8]=-1.73999040425954351098e-11;
      rc[9]=-1.78850460906003159254e-11;  ic[9]=-2.64793232535916558436e-04;
      rc[10]=-6.11984972775516891273e-05;  ic[10]=5.82371495862477667561e-11;
      rc[11]=3.18539790579877081589e-11;  ic[11]=8.52522046525751256310e-06;
      rc[12]=2.39464425437062448761e-06;  ic[12]=-9.34804896897351069155e-11;
      rc[13]=-2.09042472516339025793e-11;  ic[13]=-1.14787420736747069446e-07;
      rc[14]=-9.90579705968170287577e-08;  ic[14]=5.73200976747826012988e-11;
      break;
    case 3:
      *n=17;
      rc[0]=5.99999999999999911182e+00;  ic[0]=4.99011064204401041956e-16;
      rc[1]=2.04973195880253699388e-14;  ic[1]=-7.53670601059080969009e+00;
      rc[2]=-5.58496365805052796816e+00;  ic[2]=1.22502192765571836485e-14;
      rc[3]=-7.42700086806577689228e-13;  ic[3]=2.97144750271350321924e+00;
      rc[4]=1.26521899117147063052e+00;  ic[4]=-1.55083512390841304998e-13;
      rc[5]=8.17204799869636230066e-12;  ic[5]=-4.51220080607706464093e-01;
      rc[6]=-1.39838374609512011704e-01;  ic[6]=6.50103640247139042316e-14;
      rc[7]=-4.14643020260693596399e-11;  ic[7]=3.83873716297114669915e-02;
      rc[8]=9.51441614463429015391e-03;  ic[8]=2.37695598241956503058e-12;
      rc[9]=1.12888026491461963262e-10;  ic[9]=-2.14957794673803622917e-03;
      rc[10]=-4.48393514099263866199e-04;  ic[10]=-7.42378423204585146718e-12;
      rc[11]=-1.75485886545179745541e-10;  ic[11]=8.67199225616744690448e-05;
      rc[12]=1.57220716535193093125e-05;  ic[12]=9.36684473171511441478e-12;
      rc[13]=1.56080462458914544178e-10;  ic[13]=-2.67028180292294481241e-06;
      rc[14]=-4.28885901884305392672e-07;  ic[14]=-5.46527490734194060957e-12;
      rc[15]=-7.39408789008707292106e-11;  ic[15]=6.52035336768981008158e-08;
      rc[16]=8.69039862896329577259e-09;  ic[16]=1.22242301716048307193e-12;
      rc[17]=1.44724959490529485391e-11;  ic[17]=-1.22581445662376487153e-09;
      break;
    case 4:
      *n=20;
      rc[0]=2.39999999999999893419e+01;  ic[0]=-1.79375063153381370776e-15;
      rc[1]=-1.10634121955375132991e-14;  ic[1]=-3.61468240423631002045e+01;
      rc[2]=-2.98765606427938372747e+01;  ic[2]=-1.02956644246267533580e-14;
      rc[3]=4.16495096194576585590e-13;  ic[3]=1.74707536688975011430e+01;
      rc[4]=8.03232346743089387076e+00;  ic[4]=-6.50524652014410466111e-13;
      rc[5]=-3.21364284008900387692e-12;  ic[5]=-3.07009931351610676487e+00;
      rc[6]=-1.01057357939570446881e+00;  ic[6]=1.07986553076934566939e-11;
      rc[7]=1.15512689194906326964e-11;  ic[7]=2.93387860735316030603e-01;
      rc[8]=7.64450379841131455461e-02;  ic[8]=-5.37265077091955424350e-11;
      rc[9]=-2.30248588550429423912e-11;  ic[9]=-1.81127271509467918653e-02;
      rc[10]=-3.94315684339589887092e-03;  ic[10]=1.28096431040525740113e-10;
      rc[11]=2.71804456474515377956e-11;  ic[11]=7.95272739854185808234e-04;
      rc[12]=1.49615791810396651625e-04;  ic[12]=-1.70295302298553896128e-10;
      rc[13]=-1.93839455267960327252e-11;  ic[13]=-2.64039266842940732901e-05;
      rc[14]=-4.39269713472772780060e-06;  ic[14]=1.33288162708960345280e-10;
      rc[15]=8.17893440051253035883e-12;  ic[15]=6.91140877456926402928e-07;
      rc[16]=1.03355362735877082546e-07;  ic[16]=-6.10931289175075353347e-11;
      rc[17]=-1.87572482321099490209e-12;  ic[17]=-1.45266181029833499685e-08;
      rc[18]=-1.96960645654392952387e-09;  ic[18]=1.51781964311192669405e-11;
      rc[19]=1.79893252013204647053e-13;  ic[19]=2.17581554143030092482e-10;
      rc[20]=2.63824506041024536198e-11;  ic[20]=-1.57870373502115787380e-12;
      break;
    case 5:
      *n=23;
      rc[0]=1.20000000000000042633e+02;  ic[0]=1.13235902068101571398e-14;
      rc[1]=-2.68580711873896078995e-13;  ic[1]=-2.04734120211816104984e+02;
      rc[2]=-1.85529627256336453911e+02;  ic[2]=-5.53675374019923124530e-13;
      rc[3]=7.24667905456982677761e-13;  ic[3]=1.17230328987316767098e+02;
      rc[4]=5.76323710061169549590e+01;  ic[4]=7.75349120885641420585e-12;
      rc[5]=1.34575358450624204433e-11;  ic[5]=-2.33828200354248423309e+01;
      rc[6]=-8.12296721090095630302e+00;  ic[6]=-4.23779389886580993942e-11;
      rc[7]=-8.24939441139052192424e-11;  ic[7]=2.47751288500960553662e+00;
      rc[8]=6.75613051984486956414e-01;  ic[8]=1.15396407193526087731e-10;
      rc[9]=2.12134515756304026917e-10;  ic[9]=-1.67008678498365176202e-01;
      rc[10]=-3.78285139396031602765e-02;  ic[10]=-1.78445067276060120749e-10;
      rc[11]=-3.08991596658912069653e-10;  ic[11]=7.91952744043218637149e-03;
      rc[12]=1.54335482581086395545e-03;  ic[12]=1.68968667767644648533e-10;
      rc[13]=2.78872035500946583284e-10;  ic[13]=-2.81641648383578878812e-04;
      rc[14]=-4.83698589492209104828e-05;  ic[14]=-1.01650540405357294398e-10;
      rc[15]=-1.61079253002597143154e-10;  ic[15]=7.85194741833700138257e-06;
      rc[16]=1.20914466986599434474e-06;  ic[16]=3.90377908649444825910e-11;
      rc[17]=5.95525990196317601460e-11;  ic[17]=-1.77236964938025159766e-07;
      rc[18]=-2.47545400967414738043e-08;  ic[18]=-9.27450680491957227862e-12;
      rc[19]=-1.36205633283653990394e-11;  ic[19]=3.30979269000267447385e-09;
      rc[20]=4.15354598784708039012e-10;  ic[20]=1.24222723694405735239e-12;
      rc[21]=1.75424613031773040268e-12;  ic[21]=-5.04810718827920047821e-11;
      rc[22]=-4.86699686780117796637e-12;  ic[22]=-7.17283789467363096761e-14;
      rc[23]=-9.72556415738544699394e-14;  ic[23]=5.10300941835457752682e-13;
      break;
    case 6:
      *n=25;
      rc[0]=7.20000000000000909495e+02;  ic[0]=5.28001673162517695448e-14;
      rc[1]=-3.13117587423797120456e-12;  ic[1]=-1.34840472127089833521e+03;
      rc[2]=-1.31791188374984949405e+03;  ic[2]=2.82191999372021256390e-12;
      rc[3]=2.56133477013332854903e-11;  ic[3]=8.88911601180224579366e+02;
      rc[4]=4.63024555024105040957e+02;  ic[4]=-2.14364464387596312510e-11;
      rc[5]=-3.88188227108362569922e-11;  ic[5]=-1.97929291218466545388e+02;
      rc[6]=-7.21206233010232864444e+01;  ic[6]=7.15892162472940014404e-11;
      rc[7]=-8.86158692265621330951e-11;  ic[7]=2.29880445203284793365e+01;
      rc[8]=6.53119119716107388030e+00;  ic[8]=-1.28255396081459799855e-10;
      rc[9]=3.57618590696209153211e-10;  ic[9]=-1.67766512216693586268e+00;
      rc[10]=-3.93979762513431308601e-01;  ic[10]=1.36576966122158473835e-10;
      rc[11]=-5.14983712264189376159e-10;  ic[11]=8.53456782471355007713e-02;
      rc[12]=1.71796570418033117678e-02;  ic[12]=-9.26124889845462365197e-11;
      rc[13]=4.20439675082866951318e-10;  ic[13]=-3.23320499834506454306e-03;
      rc[14]=-5.71861522475227042339e-04;  ic[14]=4.16097573782205617803e-11;
      rc[15]=-2.18136316826453884443e-10;  ic[15]=9.54819930997131852222e-05;
      rc[16]=1.51072928253370741044e-05;  ic[16]=-1.25617330598483915473e-11;
      rc[17]=7.47920567474547143920e-11;  ic[17]=-2.27283697745025969398e-06;
      rc[18]=-3.25954486265552767647e-07;  ic[18]=2.52045083638393292766e-12;
      rc[19]=-1.69506708544041468650e-11;  ic[19]=4.47040599146159426583e-08;
      rc[20]=5.84687893901174975494e-09;  ic[20]=-3.22020693288727630391e-13;
      rc[21]=2.44734623404294798720e-12;  ic[21]=-7.35958068995855980999e-10;
      rc[22]=-8.54942249491666587752e-11;  ic[22]=2.36676663190571902709e-14;
      rc[23]=-2.04094199432149993866e-13;  ic[23]=9.80344437954053310009e-12;
      rc[24]=8.27182905990709651727e-13;  ic[24]=-7.59773537376196240740e-16;
      rc[25]=7.48395603797179222147e-15;  ic[25]=-8.26328772225895778281e-14;
      break;
    case 7:
      *n=27;
      rc[0]=5.03999999999999363354e+03;  ic[0]=8.08505199756529414456e-13;
      rc[1]=-3.77555270085723605256e-11;  ic[1]=-1.01588330488962874369e+04;
      rc[2]=-1.05737879075195214682e+04;  ic[2]=-1.55380067535075754677e-11;
      rc[3]=3.26736549770790617111e-10;  ic[3]=7.54029309201148589636e+03;
      rc[4]=4.13008348634754020168e+03;  ic[4]=1.50904740928656868130e-10;
      rc[5]=-5.36622897935008170351e-10;  ic[5]=-1.84852959355544953723e+03;
      rc[6]=-7.02773654323887853934e+02;  ic[6]=-5.68940944683668321681e-10;
      rc[7]=-1.09961821899887395950e-09;  ic[7]=2.33036934951907966251e+02;
      rc[8]=6.87063829025273236084e+01;  ic[8]=9.05230524171791561356e-10;
      rc[9]=4.67915925420714685660e-09;  ic[9]=-1.82748470684823551835e+01;
      rc[10]=-4.43552346796885466063e+00;  ic[10]=-6.60334515983326137740e-10;
      rc[11]=-6.77453537365099480457e-09;  ic[11]=9.91399528116006512057e-01;
      rc[12]=2.05603288032340264513e-01;  ic[12]=1.53957368532390822723e-10;
      rc[13]=5.49039550567194623009e-09;  ic[13]=-3.98121049895816259134e-02;
      rc[14]=-7.23624325793481484176e-03;  ic[14]=9.38764418795747466811e-11;
      rc[15]=-2.81593725713054053619e-09;  ic[15]=1.24024195067216208566e-03;
      rc[16]=2.01236550240883546332e-04;  ic[16]=-8.69159357515885519897e-11;
      rc[17]=9.59415854516773518540e-10;  ic[17]=-3.10194237365626725232e-05;
      rc[18]=-4.55558697726970093924e-06;  ic[18]=3.19243207735721524406e-11;
      rc[19]=-2.20418589833414019046e-10;  ic[19]=6.39440106420022013871e-07;
      rc[20]=8.58468282014304558677e-08;  ic[20]=-6.63724082793533377908e-12;
      rc[21]=3.37828974379720595193e-11;  ic[21]=-1.10917231464300088345e-08;
      rc[22]=-1.36177764977618633107e-09;  ic[22]=8.13772830508736910580e-13;
      rc[23]=-3.31203507016232796813e-12;  ic[23]=1.64170668100052469187e-10;
      rc[24]=1.76010123251890346719e-11;  ic[24]=-5.49976193603900745597e-14;
      rc[25]=1.87897280002782972092e-13;  ic[25]=-2.03254376059447423920e-12;
      rc[26]=-1.47177541379656430421e-13;  ic[26]=1.58418430009589240547e-15;
      rc[27]=-4.69227664960393937304e-15;  ic[27]=1.69985088194030314866e-14;
      break;
    case 8:
      *n=27;
      rc[0]=4.03199999999986757757e+04;  ic[0]=4.88731363414744927849e-11;
      rc[1]=6.46076945267778727814e-11;  ic[1]=-8.63106643911703431513e+04;
      rc[2]=-9.47491363089909573318e+04;  ic[2]=-1.52970789362209320124e-09;
      rc[3]=-1.58122098723718233334e-09;  ic[3]=7.08961326436075760284e+04;
      rc[4]=4.05809609823038117611e+04;  ic[4]=1.04359137648829443483e-08;
      rc[5]=9.70002255358402196823e-09;  ic[5]=-1.89183202347473161353e+04;
      rc[6]=-7.47071882663170345040e+03;  ic[6]=-2.95757470428532679090e-08;
      rc[7]=-2.67823888838303836603e-08;  ic[7]=2.56706913378432864192e+03;
      rc[8]=7.82687995724326697200e+02;  ic[8]=4.45972347404129073509e-08;
      rc[9]=3.88830717768661186016e-08;  ic[9]=-2.14905159184369580316e+02;
      rc[10]=-5.37590324475326539755e+01;  ic[10]=-4.04193198643397612717e-08;
      rc[11]=-3.33454142334203613940e-08;  ic[11]=1.23667194279741572416e+01;
      rc[12]=2.63622435590602899325e+00;  ic[12]=2.35851399919971407868e-08;
      rc[13]=1.82217700450291106865e-08;  ic[13]=-5.24099960594630598365e-01;
      rc[14]=-9.77014273288818180241e-02;  ic[14]=-9.21283212051612456881e-09;
      rc[15]=-6.63862089441470429748e-09;  ic[15]=1.71581083838916494122e-02;
      rc[16]=2.84995156319164987327e-03;  ic[16]=2.45499385559987295848e-09;
      rc[17]=1.64975537354939059854e-09;  ic[17]=-4.49371654184225186902e-04;
      rc[18]=-6.74266253890703787938e-05;  ic[18]=-4.47079994109099182613e-10;
      rc[19]=-2.80764805981749258606e-10;  ic[19]=9.66709619156163409062e-06;
      rc[20]=1.32086745118445763675e-06;  ic[20]=5.46892983949244612727e-11;
      rc[21]=3.21964401878084589915e-11;  ic[21]=-1.74043804406358150694e-07;
      rc[22]=-2.14691855604613240520e-08;  ic[22]=-4.29375868850307455410e-12;
      rc[23]=-2.37801071127838256268e-12;  ic[23]=2.62864757976509641197e-09;
      rc[24]=2.72964040537827676704e-10;  ic[24]=1.95363606422740140723e-13;
      rc[25]=1.02144420892311827844e-13;  ic[25]=-3.15063928310940049923e-11;
      rc[26]=-2.07747165362108670628e-12;  ic[26]=-3.91443363089512690265e-15;
      rc[27]=-1.93856193825590045796e-15;  ic[27]=2.30153327964431122982e-13;
      break;
    default:
      XLAL_ERROR( XLAL_EINVAL );
      break;
  }
  return 0;
}

UNUSED static int GammaFunctionPolyFit(INT4 l, REAL8 hathatk, REAL8 *lnr, REAL8 *arg)
{
  REAL8 x=-2.0*hathatk,rc[30],ic[30],xn=1.,rsum=0.,isum=0.;
  COMPLEX16 sum=0.;
  INT4 i=0,n=-1;

  /* Set the polynomial coefficients */
  GetGammaFuncPolyFitCoeffs(l,rc,ic,&n);

  /* Perform interpolation of the Gamma(z) */
  for ( i=0; i<=n; i++)
  {
    rsum += rc[i]*xn;
    isum += ic[i]*xn*I;
    xn *= x;
  }

  *lnr = log(cabs(sum));
  *arg = carg(sum);

  return 0;
}

#endif /*_LALSIMIMRFACTORIZEDWAVEFORM_C*/
