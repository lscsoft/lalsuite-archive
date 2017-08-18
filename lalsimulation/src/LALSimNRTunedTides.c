/*
* Copyright (C) 2017 Tim Dietrich, Sebastiano Bernuzzi, Nathan Johnson-McDaniel, Shasvath J Kapadia, Francesco Pannarale and Sebastian Khan
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


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/Date.h>
#include <lal/Units.h>
#include <lal/LALSimIMR.h>

#include "LALSimNRTunedTides.h"

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

/**
 * function to swap masses and lambda to enforece m1 >= m2
 */
static int EnforcePrimaryMassIsm1(REAL8 *m1, REAL8 *m2, REAL8 *lambda1, REAL8 *lambda2){
    REAL8 lambda1_tmp, lambda2_tmp, m1_tmp, m2_tmp;
    if (*m1>*m2) {
       lambda1_tmp = *lambda1;
       lambda2_tmp = *lambda2;
       m1_tmp   = *m1;
       m2_tmp   = *m2;
   } else { /* swap spins and masses */
       lambda1_tmp = *lambda2;
       lambda2_tmp = *lambda1;
       m1_tmp   = *m2;
       m2_tmp   = *m1;
    }
    *m1 = m1_tmp;
    *m2 = m2_tmp;
    *lambda1 = lambda1_tmp;
    *lambda2 = lambda2_tmp;

    if (*m1 < *m2)
        XLAL_ERROR(XLAL_EDOM, "XLAL_ERROR in EnforcePrimaryMassIsm1. When trying\
 to enfore that m1 should be the larger mass.\
 After trying to enforce this m1 = %f and m2 = %f\n", *m1, *m2);

    return XLAL_SUCCESS;
}


/**
 * convienient function to compute tidal coupling constant. Eq. 2 in arXiv:1706.02969
 * given masses and lambda numbers
 */
double XLALSimNRTunedTidesComputeKappa2T(
    REAL8 m1_SI, /**< Mass of companion 1 (kg) */
    REAL8 m2_SI, /**< Mass of companion 2 (kg) */
    REAL8 lambda1, /**< (tidal deformability of mass 1) / m1^5 (dimensionless) */
    REAL8 lambda2 /**< (tidal deformability of mass 2) / m2^5 (dimensionless) */
)

{
    int errcode = EnforcePrimaryMassIsm1(&m1_SI, &m2_SI, &lambda1, &lambda2);
    XLAL_CHECK(XLAL_SUCCESS == errcode, errcode, "EnforcePrimaryMassIsm1 failed");

    const REAL8 m1 = m1_SI / LAL_MSUN_SI;
    const REAL8 m2 = m2_SI / LAL_MSUN_SI;
    const REAL8 mtot = m1 + m2;

    /* Xa and Xb are the masses normalised for a total mass = 1 */
    /* not the masses appear symmetrically so we don't need to switch them. */
    const REAL8 Xa = m1 / mtot;
    const REAL8 Xb = m2 / mtot;

    /**< tidal coupling constant. Eq. 2 in arXiv:1706.02969 */
    const REAL8 term1 = (Xa / Xb) * pow(Xa, 5.0) * lambda1;
    const REAL8 term2 = (Xb / Xa) * pow(Xb, 5.0) * lambda2;
    const REAL8 kappa2T = 3.0 * ( term1 + term2 );
    return kappa2T;
}

/**
 * compute the merger frequency of a BNS system.
 * Equation 6.15a from TD thesis
 * https://www.db-thueringen.de/servlets/MCRFileNodeServlet/dbt_derivate_00035321/Dietrich_Thesis_PDFA.pdf
 * Note these should agree with Equation 2 and Table 2 from arXiv:1504.01764
 * but there is a typo in the power of the exponents
 */
double XLALSimNRTunedTidesMergerFrequency(
    const REAL8 mtot_MSUN, /**< total mass of system (solar masses) */
    const REAL8 kappa2T /**< tidal coupling constant. Eq. 2 in arXiv:1706.02969 */
)
{

    //NOTE: check that m1 >= m2
    const REAL8 Q_0 = 0.3596;
    const REAL8 n_1 = 2.4384e-2;
    const REAL8 n_2 = -1.7167e-5;
    const REAL8 d_1 = 6.8865e-2;

    const REAL8 num = 1.0 + (n_1*kappa2T) + (n_2 * kappa2T*kappa2T);
    const REAL8 den = 1.0 + (d_1 * kappa2T);

    /* dimensionless angular frequency of merger */
    const REAL8 Momega_merger = Q_0 * (num / den);

    /* convert from angular frequency to frequency (divide by 2*pi)
     * and then convert from
     * dimensionless frequency to Hz (divide by mtot * LAL_MTSUN_SI)
     */
    const REAL8 fHz_merger = Momega_merger / (LAL_TWOPI) / (mtot_MSUN * LAL_MTSUN_SI);

    return fHz_merger;
}

/**
 * Internal function only
 * Function to call the frequency domain tidal correction
 * Equation (7) in arXiv:1706.02969
 */
double XLALSimNRTunedTidesFDTidalPhase(
    const REAL8 PN_x, /**< PN frequency parameter: PN_x = orb_freq^(2./3.) */
    const REAL8 PN_x_2, /**< PN frequency parameter: PN_x**2 */
    const REAL8 PN_x_3over2, /**< PN frequency parameter: PN_x**(3./2.) */
    const REAL8 PN_x_5over2, /**< PN frequency parameter: PN_x**(5./2.) */
    const REAL8 Xa, /**< Mass of companion 1 divided by total mass */
    const REAL8 Xb, /**< Mass of companion 2 divided by total mass */
    const REAL8 kappa2T /**< tidal coupling constant. Eq. 2 in arXiv:1706.02969 */
    )
{

    /* model parameters */
    const REAL8 c_Newt = 2.4375; /* 39.0 / 16.0 */

    const REAL8 n_1 = -17.428;
    const REAL8 n_3over2 = 31.867;
    const REAL8 n_2 = -26.414;
    const REAL8 n_5over2 = 62.362;

    const REAL8 d_1 = n_1 - 2.496; /* 3115.0/1248.0 */
    const REAL8 d_3over2 = 36.089;

    REAL8 tidal_phase = - kappa2T * c_Newt / (Xa * Xb) * PN_x_5over2;

    REAL8 num = 1.0 + (n_1 * PN_x) + (n_3over2 * PN_x_3over2) + (n_2 * PN_x_2) + (n_5over2 * PN_x_5over2) ;
    REAL8 den = 1.0 + (d_1 * PN_x) + (d_3over2 * PN_x_3over2) ;

    REAL8 ratio = num / den;

    tidal_phase *= ratio;

    return tidal_phase;
}

/**
 * Function to call the frequency domain tidal correction
 * Equation (7) in arXiv:1706.02969
 * over an array of input frequencies
 * Note internally we use m1>=m2 - this is enforced in the code.
 * So any can be supplied
 */
int XLALSimNRTunedTidesFDTidalPhaseFrequencySeries(
    const REAL8Sequence *phi_tidal, /**< [out] tidal phase frequency series */
    const REAL8Sequence *fHz, /**< list of input Gravitational wave Frequency in Hz to evaluate */
    REAL8 m1_SI, /**< Mass of companion 1 (kg) */
    REAL8 m2_SI, /**< Mass of companion 2 (kg) */
    REAL8 lambda1, /**< (tidal deformability of mass 1) / m1^5 (dimensionless) */
    REAL8 lambda2 /**< (tidal deformability of mass 2) / m2^5 (dimensionless) */
    )
{
    /* NOTE: internally m1 >= m2
     * This is enforced in the code below and we swap the lambda's
     * accordingly.
     */
    int errcode = EnforcePrimaryMassIsm1(&m1_SI, &m2_SI, &lambda1, &lambda2);
    XLAL_CHECK(XLAL_SUCCESS == errcode, errcode, "EnforcePrimaryMassIsm1 failed");

    const REAL8 m1 = m1_SI / LAL_MSUN_SI;
    const REAL8 m2 = m2_SI / LAL_MSUN_SI;
    const REAL8 mtot = m1 + m2;

    /* Xa and Xb are the masses normalised for a total mass = 1 */
    const REAL8 Xa = m1 / mtot;
    const REAL8 Xb = m2 / mtot;

    /**< tidal coupling constant.*/
    const REAL8 kappa2T = XLALSimNRTunedTidesComputeKappa2T(m1_SI, m2_SI, lambda1, lambda2);

    /* initialise */
    REAL8 PN_x = 0.0;
    REAL8 PN_x_2 = 0.0;
    REAL8 PN_x_3over2 = 0.0;
    REAL8 PN_x_5over2 = 0.0;

    for(UINT4 i = 0; i < (*fHz).length; i++){
        /* NRTunedTidesFDTidalPhase is Eq 7 in arXiv:1706.02969
         * and is a function of x = orb_freq^(2./3.)
         */
        PN_x = pow((*fHz).data[i] / 2.0, 2.0/3.0);
        PN_x_2 = PN_x * PN_x;
        PN_x_3over2 = pow(PN_x, 3.0/2.0);
        PN_x_5over2 = pow(PN_x, 5.0/2.0);

        (*phi_tidal).data[i] = XLALSimNRTunedTidesFDTidalPhase(PN_x, PN_x_2, PN_x_3over2, PN_x_5over2, Xa, Xb, kappa2T);
    }

    return XLAL_SUCCESS;
}
