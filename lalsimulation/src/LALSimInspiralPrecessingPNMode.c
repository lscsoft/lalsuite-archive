/*
 * Copyright (C) 2013 Evan Ochsner
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

#include <lal/LALSimInspiralPrecessingPNMode.h>
#include <lal/LALConstants.h>
#include "check_series_macros.h"

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

/**
 * FIXME - document this function!
 */
int XLALSimInspiralRadiationFrameToJFrame(
        REAL8TimeSeries **S1x_J,       /**< [returned] J-frame time series */
        REAL8TimeSeries **S1y_J,       /**< [returned] J-frame time series */
        REAL8TimeSeries **S1z_J,       /**< [returned] J-frame time series */
        REAL8TimeSeries **S2x_J,       /**< [returned] J-frame time series */
        REAL8TimeSeries **S2y_J,       /**< [returned] J-frame time series */
        REAL8TimeSeries **S2z_J,       /**< [returned] J-frame time series */
        REAL8TimeSeries **LNhatx_J,    /**< [returned] J-frame time series */
        REAL8TimeSeries **LNhaty_J,    /**< [returned] J-frame time series */
        REAL8TimeSeries **LNhatz_J,    /**< [returned] J-frame time series */
        const REAL8TimeSeries *V,      /**< post-Newtonian parameter */
        const REAL8TimeSeries *S1x,    /**< Spin1 vector x component */
        const REAL8TimeSeries *S1y,    /**< Spin1 vector y component */
        const REAL8TimeSeries *S1z,    /**< Spin1 vector z component */
        const REAL8TimeSeries *S2x,    /**< Spin2 vector x component */
        const REAL8TimeSeries *S2y,    /**< Spin2 vector y component */
        const REAL8TimeSeries *S2z,    /**< Spin2 vector z component */
        const REAL8TimeSeries *LNhatx, /**< unit orbital ang. mom. x comp. */
        const REAL8TimeSeries *LNhaty, /**< unit orbital ang. mom. y comp. */
        const REAL8TimeSeries *LNhatz, /**< unit orbital ang. mom. z comp. */
        const REAL8 m1,                /**< mass of companion 1 (kg) */
        const REAL8 m2,                /**< mass of companion 2 (kg) */
        const REAL8 fref               /**< frame is set by J at this GW freq.*/
        )
{
    LAL_CHECK_VALID_SERIES(V, XLAL_FAILURE);
    LAL_CHECK_VALID_SERIES(S1x, XLAL_FAILURE);
    LAL_CHECK_VALID_SERIES(S1y, XLAL_FAILURE);
    LAL_CHECK_VALID_SERIES(S1z, XLAL_FAILURE);
    LAL_CHECK_VALID_SERIES(S2x, XLAL_FAILURE);
    LAL_CHECK_VALID_SERIES(S2y, XLAL_FAILURE);
    LAL_CHECK_VALID_SERIES(S2z, XLAL_FAILURE);
    LAL_CHECK_VALID_SERIES(LNhatx, XLAL_FAILURE);
    LAL_CHECK_VALID_SERIES(LNhaty, XLAL_FAILURE);
    LAL_CHECK_VALID_SERIES(LNhatz, XLAL_FAILURE);
    LAL_CHECK_CONSISTENT_TIME_SERIES(V, S1x, XLAL_FAILURE);
    LAL_CHECK_CONSISTENT_TIME_SERIES(V, S1y, XLAL_FAILURE);
    LAL_CHECK_CONSISTENT_TIME_SERIES(V, S1z, XLAL_FAILURE);
    LAL_CHECK_CONSISTENT_TIME_SERIES(V, S2x, XLAL_FAILURE);
    LAL_CHECK_CONSISTENT_TIME_SERIES(V, S2y, XLAL_FAILURE);
    LAL_CHECK_CONSISTENT_TIME_SERIES(V, S2z, XLAL_FAILURE);
    LAL_CHECK_CONSISTENT_TIME_SERIES(V, LNhatx, XLAL_FAILURE);
    LAL_CHECK_CONSISTENT_TIME_SERIES(V, LNhaty, XLAL_FAILURE);
    LAL_CHECK_CONSISTENT_TIME_SERIES(V, LNhatz, XLAL_FAILURE);

    REAL8 LNx, LNy, LNz, LNmag, s1x, s1y, s1z, s2x, s2y, s2z, Jx, Jy, Jz;
    REAL8 E1x, E1y, E1z, E2x, E2y, E2z, norm, JdotN;
    REAL8 eta = m1 * m2 / (m1+m2) / (m1+m2);
    REAL8 m1_sec = LAL_G_SI * m1 / LAL_C_SI / LAL_C_SI / LAL_C_SI;
    REAL8 m2_sec = LAL_G_SI * m2 / LAL_C_SI / LAL_C_SI / LAL_C_SI;
    REAL8 M_sec = m1_sec + m2_sec;
    REAL8 M_sec_sq = M_sec * M_sec;
    REAL8 Vref = pow( LAL_PI * M_sec * fref, 1./3.);
    unsigned int j, refIdx = 0;

    *S1x_J = XLALCreateREAL8TimeSeries( S1x->name, &S1x->epoch, S1x->f0,
            S1x->deltaT, &S1x->sampleUnits, S1x->data->length);
    *S1y_J = XLALCreateREAL8TimeSeries( S1y->name, &S1y->epoch, S1y->f0,
            S1y->deltaT, &S1y->sampleUnits, S1y->data->length);
    *S1z_J = XLALCreateREAL8TimeSeries( S1z->name, &S1z->epoch, S1z->f0,
            S1z->deltaT, &S1z->sampleUnits, S1z->data->length);
    *S2x_J = XLALCreateREAL8TimeSeries( S2x->name, &S2x->epoch, S2x->f0,
            S2x->deltaT, &S2x->sampleUnits, S2x->data->length);
    *S2y_J = XLALCreateREAL8TimeSeries( S2y->name, &S2y->epoch, S2y->f0,
            S2y->deltaT, &S2y->sampleUnits, S2y->data->length);
    *S2z_J = XLALCreateREAL8TimeSeries( S2z->name, &S2z->epoch, S2z->f0,
            S2z->deltaT, &S2z->sampleUnits, S2z->data->length);
    *LNhatx_J = XLALCreateREAL8TimeSeries( LNhatx->name, &LNhatx->epoch,
            LNhatx->f0, LNhatx->deltaT, &LNhatx->sampleUnits,
            LNhatx->data->length);
    *LNhaty_J = XLALCreateREAL8TimeSeries( LNhaty->name, &LNhaty->epoch,
            LNhaty->f0, LNhaty->deltaT, &LNhaty->sampleUnits,
            LNhaty->data->length);
    *LNhatz_J = XLALCreateREAL8TimeSeries( LNhatz->name, &LNhatz->epoch,
            LNhatz->f0, LNhatz->deltaT, &LNhatz->sampleUnits,
            LNhatz->data->length);
    if( !S1x_J || !S1y_J || !S1z_J || !S2x_J || !S2y_J || !S2z_J
            || !LNhatx_J || !LNhaty_J || !LNhatz_J )
        XLAL_ERROR(XLAL_ENOMEM);

    // Find index at which the reference freq. is reached
    if( fref == 0 ) // fref==0 is a special case, set frame from J(fmin)
        Vref = V->data->data[0];
    while( V->data->data[refIdx] < Vref )
        refIdx++;

    // Find J = L_N + S_1 + S_2 at refIdx sample, normalize - this sets E3-axis
    LNx = LNhatx->data->data[refIdx];
    LNy = LNhaty->data->data[refIdx];
    LNz = LNhatz->data->data[refIdx];
    LNmag = M_sec_sq * eta / V->data->data[refIdx];
     /* EvolveOrbit functions return spin components with mass prefactors.
      * Specifically, S1z = (m1_sec/M_sec)^2 * \chi_1 * \hat{S}_1
      * and similarly for other components.
      * We want m1_sec^2 * \chi_1 * \hat{S}_1, so multiply by M_sec^2 */
    s1x = S1x->data->data[refIdx] * M_sec_sq;
    s1y = S1y->data->data[refIdx] * M_sec_sq;
    s1z = S1z->data->data[refIdx] * M_sec_sq;
    s2x = S2x->data->data[refIdx] * M_sec_sq;
    s2y = S2y->data->data[refIdx] * M_sec_sq;
    s2z = S2z->data->data[refIdx] * M_sec_sq;
    Jx = LNmag*LNx + s1x + s2x;
    Jy = LNmag*LNy + s1y + s2y;
    Jz = LNmag*LNz + s1z + s2z;
    norm = sqrt(Jx*Jx + Jy*Jy + Jz*Jz);
    Jx /= norm;
    Jy /= norm;
    Jz /= norm;
    // Now find E1 and E2 basis vectors
    // FIXME: What to do if J || N ?
    JdotN = Jz; // b/c N = {0, 0, 1}
    // E1 \propto N - (J \cdot N) J
    E1x = - JdotN * Jx;
    E1y = - JdotN * Jy;
    E1z = 1. - JdotN * Jx;
    norm = sqrt(E1x*E1x + E1y*E1y + E1z*E1z);
    E1x /= norm;
    E1y /= norm;
    E1z /= norm;
    // E2 = J x E1
    E2x = Jy * E1z - Jz * E1y;
    E2y = Jz * E1x - Jx * E1z;
    E2z = Jx * E1y - Jy * E1x;
    norm = sqrt(E2x*E2x + E2y*E2y + E2z*E2z);
    E2x /= norm;
    E2y /= norm;
    E2z /= norm;

    /* E1, E2, J are the {x,y,z} basis vectors of the J-frame.  Compute the
     * spin and L_N components in this frame and copy them to output series.
     * Explicitly, S1x_J = S1 \cdot E1, S1y_J = S1 \cdot E2, S1z_J = S1 \cdot J
     * and similarly for the other vectors in J-frame */
    for(j = 0; j < V->data->length; j++)
    {
        (*S1x_J)->data->data[j] = S1x->data->data[j]*E1x
            + S1y->data->data[j]*E1y + S1z->data->data[j]*E1z;
        (*S1y_J)->data->data[j] = S1x->data->data[j]*E2x
            + S1y->data->data[j]*E2y + S1z->data->data[j]*E2z;
        (*S1z_J)->data->data[j] = S1x->data->data[j]*Jx
            + S1y->data->data[j]*Jy + S1z->data->data[j]*Jz;
        (*S2x_J)->data->data[j] = S2x->data->data[j]*E1x
            + S2y->data->data[j]*E1y + S2z->data->data[j]*E1z;
        (*S2y_J)->data->data[j] = S2x->data->data[j]*E2x
            + S2y->data->data[j]*E2y + S2z->data->data[j]*E2z;
        (*S2z_J)->data->data[j] = S2x->data->data[j]*Jx
            + S2y->data->data[j]*Jy + S2z->data->data[j]*Jz;
        (*LNhatx_J)->data->data[j] = LNhatx->data->data[j]*E1x
            + LNhaty->data->data[j]*E1y + LNhatz->data->data[j]*E1z;
        (*LNhaty_J)->data->data[j] = LNhatx->data->data[j]*E2x
            + LNhaty->data->data[j]*E2y + LNhatz->data->data[j]*E2z;
        (*LNhatz_J)->data->data[j] = LNhatx->data->data[j]*Jx
            + LNhaty->data->data[j]*Jy + LNhatz->data->data[j]*Jz;
    }

    return XLAL_SUCCESS;
}

/**
 * Computes h(2,2) mode of spherical harmonic decomposition of
 * precessing post-Newtonian inspiral waveforms.
 *
 * This expression is derived to 1.5PN order in arXiv:0810.5336,
 * with the expressions provided in the Mathematica file ABFO09.m.
 * The expressions are further simplified in a NB found here:
 * FIXME - add URL of wiki page with NB deriving simpler form!
 */
COMPLEX16TimeSeries *XLALSimInspiralPrecessingPNMode22(
        const REAL8TimeSeries *V,      /**< post-Newtonian parameter */
        const REAL8TimeSeries *Phi,    /**< orbital phase */
        const REAL8TimeSeries *S1x,    /**< Spin1 x component in J frame */
        const REAL8TimeSeries *S1y,    /**< Spin1 y component in J frame */
        const REAL8TimeSeries *S1z,    /**< Spin1 z component in J frame */
        const REAL8TimeSeries *S2x,    /**< Spin2 x component in J frame */
        const REAL8TimeSeries *S2y,    /**< Spin2 y component in J frame */
        const REAL8TimeSeries *S2z,    /**< Spin2 z component in J frame */
        const REAL8TimeSeries *LNhatx, /**< unit orbital ang. mom. x comp. in J frame */
        const REAL8TimeSeries *LNhaty, /**< unit orbital ang. mom. y comp. in J frame */
        const REAL8TimeSeries *LNhatz, /**< unit orbital ang. mom. z comp. in J frame */
        const REAL8 m1,                /**< mass of companion 1 (kg) */
        const REAL8 m2,                /**< mass of companion 2 (kg) */
        const REAL8 r,                 /**< distance of source (m) */
        const int O                    /**< twice post-Newtonian order */
        )
{
    // NB: other series will be checked in XLALSimInspiralRadiationFrameToJFrame
    LAL_CHECK_VALID_SERIES(Phi, NULL);
    LAL_CHECK_CONSISTENT_TIME_SERIES(V, Phi, NULL);

    COMPLEX16TimeSeries *hlm = XLALCreateCOMPLEX16TimeSeries("h_22 mode",
            &V->epoch, 0., V->deltaT, &lalStrainUnit, V->data->length);
    if( !hlm ) XLAL_ERROR_NULL(XLAL_ENOMEM);

    size_t j;
    REAL8 M = m1+m2;
    REAL8 m1_sec = LAL_G_SI * m1 / LAL_C_SI / LAL_C_SI / LAL_C_SI;
    REAL8 m2_sec = LAL_G_SI * m2 / LAL_C_SI / LAL_C_SI / LAL_C_SI;
    REAL8 M_sec = m1_sec + m2_sec;
    REAL8 eta = m1*m2/M/M;
    REAL8 dm = (m1-m2)/M;
    REAL8 m1oM = m1/M;
    REAL8 m2oM = m2/M;
    REAL8 lalpi = LAL_PI;
    REAL8 prefac = -8.0*sqrt(LAL_PI/5.0)*M_sec*eta/r;
    COMPLEX16 term0, term1, term2, term3;
    COMPLEX16 expmalpha, expm2alpha, expm3alpha, expphi, expmphi;
    COMPLEX16 exp2phi, expm2phi;
    COMPLEX16 S1perp, S2perp, conjS1perp, conjS2perp;
    REAL8 v, v2, iota, alpha, s1x, s1y, s1z, s2x, s2y, s2z;
    REAL8 cosiota, cosiota2, cos2iotao2, cos4iotao2, cos5iotao2;
    REAL8 siniota, sin2iota, sin3iota, siniotao2, sin2iotao2, sin4iotao2;

    // Convert vectors from radiation frame to J frame
    /*ret = XLALSimInspiralRadiationFrameToJFrame(&S1x_J, &S1y_J, &S1z_J,
            &S2x_J, &S2y_J, &S2z_J, &LNhatx_J, &LNhaty_J, &LNhatz_J, V,
            S1x, S1y, S1z, S2x, S2y, S2z, LNhatx, LNhaty, LNhatz, m1, m2, fref);
    if(ret != XLAL_SUCCESS) XLAL_ERROR_NULL(XLAL_EFUNC);*/

    for(j=0; j < V->data->length; j++) {
        // Compute the dynamic quantities in the terms below
        v = V->data->data[j];
        v2 = v*v;
        expphi = cpolar(1., Phi->data->data[j]);
        expmphi = cpolar(1., - Phi->data->data[j]);
        exp2phi = cpolar(1., 2.*Phi->data->data[j]);
        expm2phi = cpolar(1., - 2.*Phi->data->data[j]);
        /* EvolveOrbit functions return spin components with mass prefactors.
         * Specifically, S1z = (m1/M)^2 * \chi_1 * \hat{S}_1
         * and similarly for other components.
         * Thus, we divide out the mass prefactors so s1z = \chi_1 * \hat{S}_1z
         */
        s1x = S1x->data->data[j] / m1oM / m1oM;
        s1y = S1y->data->data[j] / m1oM / m1oM;
        s1z = S1z->data->data[j] / m1oM / m1oM;
        s2x = S2x->data->data[j] / m2oM / m2oM;
        s2y = S2y->data->data[j] / m2oM / m2oM;
        s2z = S2z->data->data[j] / m2oM / m2oM;
        S1perp = crect(s1x, s1y);
        conjS1perp = crect(s1x, -s1y);
        S2perp = crect(s2x, s2y);
        conjS2perp = crect(s2x, -s2y);

        iota = acos(LNhatz->data->data[j]);
        alpha = atan2(LNhaty->data->data[j], LNhatx->data->data[j]);
        expmalpha = cpolar(1., - alpha);
        expm2alpha = cpolar(1., - 2.*alpha);
        expm3alpha = cpolar(1., - 3.*alpha);
        cosiota = cos(iota);
        cosiota2 = cos(2.*iota);
        cos2iotao2 = cos(iota/2.)*cos(iota/2.);
        cos4iotao2 = cos2iotao2*cos2iotao2;
        cos5iotao2 = cos4iotao2*cos(iota/2.);
        siniota = sin(iota);
        sin2iota = siniota*siniota;
        sin3iota = sin2iota*siniota;
        siniotao2 = sin(iota/2.);
        sin2iotao2 = siniotao2*siniotao2;
        sin4iotao2 = sin2iotao2*sin2iotao2;

        // zero the terms
        term0 = term1 = term2 = term3 = 0.;

        // Compute the terms up to the requested PN order
        switch(O) {
            default: // unsupported PN order
                XLALPrintError("XLAL Error - %s: PN order %d%s not supported\n", __func__, O/2, O%2?".5":"" );
                XLAL_ERROR_NULL(XLAL_EINVAL);
            case -1: // Use highest available PN order
            case 3:
                term3 = (expmalpha*((-3*conjS1perp*(3 + cosiota2) - 3*conjS2perp
                        *(3 + cosiota2))*eta + 12*conjS1perp*(3 + cosiota2)*m1oM
                        + 12*conjS2perp*(3 + cosiota2)*m2oM
                        + cos2iotao2*expm2phi*((-2*conjS1perp*(-9 + cosiota)
                        - 2*conjS2perp*(-9 + cosiota))*eta + conjS1perp
                        *(8 - 40*cosiota)*m1oM + conjS2perp
                        *(8 - 40*cosiota)*m2oM) + exp2phi*((2*conjS1perp
                        *(9 + cosiota) + 2*conjS2perp*(9 + cosiota))
                        *eta + conjS1perp*(8 + 40*cosiota)*m1oM + conjS2perp
                        *(8 + 40*cosiota)*m2oM)*sin2iotao2)*siniota)/48.
                        + (expm2alpha*(cos4iotao2*expm2phi*(168*lalpi - 56
                        *(-3 + 5*cosiota)*m1oM*s1z - 56*(-3 + 5*cosiota)
                        *m2oM*s2z - 14*(-5 + cosiota)*eta*(s1z + s2z))
                        + (-84*cosiota*m1oM*s1z - 84*cosiota*m2oM*s2z
                        + 21*cosiota*eta*(s1z + s2z))*sin2iota + exp2phi
                        *(168*lalpi - 56*(3 + 5*cosiota)*m1oM*s1z
                        - 56*(3 + 5*cosiota)*m2oM*s2z - 14*(5 + cosiota)*eta
                        *(s1z + s2z))*sin4iotao2 + dm*(cos2iotao2
                        *(-17 + 20*eta)*expmphi + (-17 + 20*eta)*expphi
                        *sin2iotao2)*siniota))/84. + (expm3alpha
                        *(3*((eta - 4*m1oM)*S1perp + (eta - 4*m2oM)*S2perp)
                        *sin3iota - 2*exp2phi*((eta + 20*m1oM)*S1perp
                        + (eta + 20*m2oM)*S2perp)*sin4iotao2*siniota
                        - 4*cos5iotao2*expm2phi*((eta + 20*m1oM)*S1perp
                        + (eta + 20*m2oM)*S2perp)*siniotao2))/24.;
            case 2:
                term2 = (expmalpha*(-(conjS1perp*m1oM) + conjS2perp*m2oM)
                        *(cos2iotao2*expmphi + expphi*sin2iotao2))/2.
                        + ((-107 + 55*eta)*expm2alpha*(cos4iotao2*expm2phi
                        + exp2phi*sin4iotao2))/42.;
            case 1:
                term1 = (dm*expm2alpha*(cos2iotao2*expmphi + expphi*sin2iotao2)
                        *siniota)/3.;
            case 0:
                term0 = expm2alpha*(cos4iotao2*expm2phi + exp2phi*sin4iotao2);
                break;
        }

        v = V->data->data[j];
        v2 = v*v;
        hlm->data->data[j] = prefac * v2
                * (term0 + v*(term1 + v*(term2 + v*(term3))));
    }

    return hlm;
}

/**
 * Computes h(2,-2) mode of spherical harmonic decomposition of
 * precessing post-Newtonian inspiral waveforms.
 *
 * This expression is derived to 1.5PN order in arXiv:0810.5336,
 * with the expressions provided in the Mathematica file ABFO09.m.
 * The expressions are further simplified in a NB found here:
 * FIXME - add URL of wiki page with NB deriving simpler form!
 *
 * We get negative modes via the symmetry (see text below Eq. 4.17 of 0810.5336)
 * h_{l-m} = (-1)^l h_{lm}^*(Phi + pi)
 */
COMPLEX16TimeSeries *XLALSimInspiralPrecessingPNMode2m2(
        const REAL8TimeSeries *V,      /**< post-Newtonian parameter */
        const REAL8TimeSeries *Phi,    /**< orbital phase */
        const REAL8TimeSeries *S1x,    /**< Spin1 x component in J frame */
        const REAL8TimeSeries *S1y,    /**< Spin1 y component in J frame */
        const REAL8TimeSeries *S1z,    /**< Spin1 z component in J frame */
        const REAL8TimeSeries *S2x,    /**< Spin2 x component in J frame */
        const REAL8TimeSeries *S2y,    /**< Spin2 y component in J frame */
        const REAL8TimeSeries *S2z,    /**< Spin2 z component in J frame */
        const REAL8TimeSeries *LNhatx, /**< unit orbital ang. mom. x comp. in J frame */
        const REAL8TimeSeries *LNhaty, /**< unit orbital ang. mom. y comp. in J frame */
        const REAL8TimeSeries *LNhatz, /**< unit orbital ang. mom. z comp. in J frame */
        const REAL8 m1,                /**< mass of companion 1 (kg) */
        const REAL8 m2,                /**< mass of companion 2 (kg) */
        const REAL8 r,                 /**< distance of source (m) */
        const int O                    /**< twice post-Newtonian order */
        )
{
    unsigned int i;
    // Copy Phi TimeSeries and apply a shift of Pi
    REAL8TimeSeries *PhiNew = XLALCutREAL8TimeSeries(Phi, 0, Phi->data->length);
    for(i=0; i < Phi->data->length; i++)
        PhiNew->data->data[i] += LAL_PI;
    // Generate m > 0 mode
    COMPLEX16TimeSeries *hlm = XLALSimInspiralPrecessingPNMode22(V, PhiNew,
            S1x, S1y, S1z, S2x, S2y, S2z, LNhatx, LNhaty, LNhatz, m1, m2, r, O);
    if( !hlm ) XLAL_ERROR_NULL(XLAL_EFUNC);
    // Apply conjugation and (-1)^l
    for(i=0; i < hlm->data->length; i++)
        hlm->data->data[i] = conj(hlm->data->data[i]);

    return hlm;
}

/**
 * Computes h(2,1) mode of spherical harmonic decomposition of
 * precessing post-Newtonian inspiral waveforms.
 *
 * This expression is derived to 1.5PN order in arXiv:0810.5336,
 * with the expressions provided in the Mathematica file ABFO09.m.
 * The expressions are further simplified in a NB found here:
 * FIXME - add URL of wiki page with NB deriving simpler form!
 */
COMPLEX16TimeSeries *XLALSimInspiralPrecessingPNMode21(
        const REAL8TimeSeries *V,      /**< post-Newtonian parameter */
        const REAL8TimeSeries *Phi,    /**< orbital phase */
        const REAL8TimeSeries *S1x,    /**< Spin1 x component in J frame */
        const REAL8TimeSeries *S1y,    /**< Spin1 y component in J frame */
        const REAL8TimeSeries *S1z,    /**< Spin1 z component in J frame */
        const REAL8TimeSeries *S2x,    /**< Spin2 x component in J frame */
        const REAL8TimeSeries *S2y,    /**< Spin2 y component in J frame */
        const REAL8TimeSeries *S2z,    /**< Spin2 z component in J frame */
        const REAL8TimeSeries *LNhatx, /**< unit orbital ang. mom. x comp. in J frame */
        const REAL8TimeSeries *LNhaty, /**< unit orbital ang. mom. y comp. in J frame */
        const REAL8TimeSeries *LNhatz, /**< unit orbital ang. mom. z comp. in J frame */
        const REAL8 m1,                /**< mass of companion 1 (kg) */
        const REAL8 m2,                /**< mass of companion 2 (kg) */
        const REAL8 r,                 /**< distance of source (m) */
        const int O                    /**< twice post-Newtonian order */
        )
{
    // NB: other series will be checked in XLALSimInspiralRadiationFrameToJFrame
    LAL_CHECK_VALID_SERIES(Phi, NULL);
    LAL_CHECK_CONSISTENT_TIME_SERIES(V, Phi, NULL);

    COMPLEX16TimeSeries *hlm = XLALCreateCOMPLEX16TimeSeries("h_21 mode",
            &V->epoch, 0., V->deltaT, &lalStrainUnit, V->data->length);
    if( !hlm ) XLAL_ERROR_NULL(XLAL_ENOMEM);

    size_t j;
    REAL8 M = m1+m2;
    REAL8 m1_sec = LAL_G_SI * m1 / LAL_C_SI / LAL_C_SI / LAL_C_SI;
    REAL8 m2_sec = LAL_G_SI * m2 / LAL_C_SI / LAL_C_SI / LAL_C_SI;
    REAL8 M_sec = m1_sec + m2_sec;
    REAL8 eta = m1*m2/M/M;
    REAL8 dm = (m1-m2)/M;
    REAL8 m1oM = m1/M;
    REAL8 m2oM = m2/M;
    REAL8 lalpi = LAL_PI;
    REAL8 prefac = -8.0*sqrt(LAL_PI/5.0)*M_sec*eta/r;
    COMPLEX16 term0, term1, term2, term3;
    COMPLEX16 expmalpha, expm2alpha, expphi, expmphi;
    COMPLEX16 exp2phi, expm2phi;
    COMPLEX16 S1perp, S2perp, conjS1perp, conjS2perp;
    REAL8 v, v2, iota, alpha, s1x, s1y, s1z, s2x, s2y, s2z;
    REAL8 cosiota, cosiota2, cosiota3, cos2iotao2;
    REAL8 siniota, siniota2, siniota3, sin2iota, siniotao2, sin2iotao2;

    for(j=0; j < V->data->length; j++) {
        // Compute the dynamic quantities in the terms below
        v = V->data->data[j];
        v2 = v*v;
        expphi = cpolar(1., Phi->data->data[j]);
        expmphi = cpolar(1., - Phi->data->data[j]);
        exp2phi = cpolar(1., 2.*Phi->data->data[j]);
        expm2phi = cpolar(1., - 2.*Phi->data->data[j]);
        /* EvolveOrbit functions return spin components with mass prefactors.
         * Specifically, S1z = (m1/M)^2 * \chi_1 * \hat{S}_1
         * and similarly for other components.
         * Thus, we divide out the mass prefactors so s1z = \chi_1 * \hat{S}_1z
         */
        s1x = S1x->data->data[j] / m1oM / m1oM;
        s1y = S1y->data->data[j] / m1oM / m1oM;
        s1z = S1z->data->data[j] / m1oM / m1oM;
        s2x = S2x->data->data[j] / m2oM / m2oM;
        s2y = S2y->data->data[j] / m2oM / m2oM;
        s2z = S2z->data->data[j] / m2oM / m2oM;
        S1perp = crect(s1x, s1y);
        conjS1perp = crect(s1x, -s1y);
        S2perp = crect(s2x, s2y);
        conjS2perp = crect(s2x, -s2y);

        iota = acos(LNhatz->data->data[j]);
        alpha = atan2(LNhaty->data->data[j], LNhatx->data->data[j]);
        expmalpha = cpolar(1., - alpha);
        expm2alpha = cpolar(1., - 2.*alpha);
        cosiota = cos(iota);
        cosiota2 = cos(2.*iota);
        cosiota3 = cos(3.*iota);
        cos2iotao2 = cos(iota/2.)*cos(iota/2.);
        siniota = sin(iota);
        siniota2 = sin(2.*iota);
        siniota3 = sin(3.*iota);
        sin2iota = siniota*siniota;
        siniotao2 = sin(iota/2.);
        sin2iotao2 = siniotao2*siniotao2;

        // zero the terms
        term0 = term1 = term2 = term3 = 0.;

        // Compute the terms up to the requested PN order
        switch(O) {
            default: // unsupported PN order
                XLALPrintError("XLAL Error - %s: PN order %d%s not supported\n", __func__, O/2, O%2?".5":"" );
                XLAL_ERROR_NULL(XLAL_EINVAL);
            case -1: // Use highest available PN order
            case 3:
                term3 = (cos2iotao2*expm2phi*((conjS1perp + conjS2perp)
                        *(-14 + 15*cosiota - cosiota2)*eta + 4*conjS1perp
                        *(-4 + 9*cosiota - 5*cosiota2)*m1oM + 4*conjS2perp
                        *(-4 + 9*cosiota - 5*cosiota2)*m2oM))/8.
                        + (9*cosiota*(conjS1perp*(-eta + 4*m1oM)
                        + conjS2perp*(-eta + 4*m2oM)))/16.
                        + (3*cosiota3*(conjS1perp*(-eta + 4*m1oM)
                        + conjS2perp*(-eta + 4*m2oM)))/16.
                        + (exp2phi*((conjS1perp + conjS2perp)
                        *(14 + 15*cosiota + cosiota2)*eta
                        + 4*conjS1perp*(4 + 9*cosiota + 5*cosiota2)*m1oM
                        + 4*conjS2perp*(4 + 9*cosiota + 5*cosiota2)*m2oM)
                        *sin2iotao2)/8. + (expm2alpha*(cos2iotao2*expm2phi
                        *(-8*cos2iotao2*(-7 + 10*cosiota)*m1oM*S1perp
                        - 8*cos2iotao2*(-7 + 10*cosiota)*m2oM*S2perp
                        - (4 + 5*cosiota + cosiota2)*eta*(S1perp + S2perp))
                        + (-24*cosiota*m1oM*S1perp - 24*cosiota*m2oM*S2perp
                        + 6*cosiota*eta*(S1perp + S2perp))*sin2iota
                        + exp2phi*sin2iotao2*(-((-4 + 5*cosiota - cosiota2)
                        *eta*(S1perp + S2perp))
                        - 8*(7 + 10*cosiota)*m1oM*S1perp*sin2iotao2
                        - 8*(7 + 10*cosiota)*m2oM*S2perp*sin2iotao2)))/8.
                        + (expmalpha*(dm*(cos2iotao2*(34 - 40*eta + cosiota
                        *(-68 + 80*eta))*expmphi + (-34 + 40*eta + cosiota
                        *(-68 + 80*eta))*expphi*sin2iotao2)
                        + (84*m1oM*s1z + 84*m2oM*s2z - 21*eta*(s1z + s2z)
                        + cos2iotao2*expm2phi*(-336*lalpi - 168*m1oM*s1z
                        - 168*m2oM*s2z - 70*eta*(s1z + s2z)) + exp2phi
                        *(336*lalpi - 168*m1oM*s1z - 168*m2oM*s2z - 70*eta
                        *(s1z + s2z))*sin2iotao2)*siniota
                        + (cos2iotao2*expm2phi*(280*m1oM*s1z + 280*m2oM*s2z
                        + 14*eta*(s1z + s2z)) + exp2phi*(-280*m1oM*s1z
                        - 280*m2oM*s2z - 14*eta*(s1z + s2z))*sin2iotao2)
                        *siniota2 + (-84*m1oM*s1z - 84*m2oM*s2z
                        + 21*eta*(s1z + s2z))*siniota3))/56.;
            case 2:
                term2 = (3*(expmphi - expphi)
                        *(conjS1perp*m1oM - conjS2perp*m2oM)*siniota)/4.
                        + expmalpha*((3*cos2iotao2*expmphi*(-(m1oM*s1z)
                        + m2oM*s2z))/2. + (3*expphi*(-(m1oM*s1z) + m2oM*s2z)
                        *sin2iotao2)/2. + (7.642857142857143 - (55*eta)/14.)
                        *(cos2iotao2*expm2phi - exp2phi*sin2iotao2)*siniota);
            case 1:
                term1 = dm*expmalpha*(cos2iotao2*(-1 + 2*cosiota)*expmphi
                        + (1 + 2*cosiota)*expphi*sin2iotao2);
            case 0:
                term0 = 3*expmalpha*(-(cos2iotao2*expm2phi)
                        + exp2phi*sin2iotao2)*siniota;
                break;
        }

        v = V->data->data[j];
        v2 = v*v;
        hlm->data->data[j] = prefac * v2
                * (term0 + v*(term1 + v*(term2 + v*(term3))));
    }

    return hlm;
}

/**
 * Computes h(2,-1) mode of spherical harmonic decomposition of
 * precessing post-Newtonian inspiral waveforms.
 *
 * This expression is derived to 1.5PN order in arXiv:0810.5336,
 * with the expressions provided in the Mathematica file ABFO09.m.
 * The expressions are further simplified in a NB found here:
 * FIXME - add URL of wiki page with NB deriving simpler form!
 *
 * We get negative modes via the symmetry (see text below Eq. 4.17 of 0810.5336)
 * h_{l-m} = (-1)^l h_{lm}^*(Phi + pi)
 */
COMPLEX16TimeSeries *XLALSimInspiralPrecessingPNMode2m1(
        const REAL8TimeSeries *V,      /**< post-Newtonian parameter */
        const REAL8TimeSeries *Phi,    /**< orbital phase */
        const REAL8TimeSeries *S1x,    /**< Spin1 x component in J frame */
        const REAL8TimeSeries *S1y,    /**< Spin1 y component in J frame */
        const REAL8TimeSeries *S1z,    /**< Spin1 z component in J frame */
        const REAL8TimeSeries *S2x,    /**< Spin2 x component in J frame */
        const REAL8TimeSeries *S2y,    /**< Spin2 y component in J frame */
        const REAL8TimeSeries *S2z,    /**< Spin2 z component in J frame */
        const REAL8TimeSeries *LNhatx, /**< unit orbital ang. mom. x comp. in J frame */
        const REAL8TimeSeries *LNhaty, /**< unit orbital ang. mom. y comp. in J frame */
        const REAL8TimeSeries *LNhatz, /**< unit orbital ang. mom. z comp. in J frame */
        const REAL8 m1,                /**< mass of companion 1 (kg) */
        const REAL8 m2,                /**< mass of companion 2 (kg) */
        const REAL8 r,                 /**< distance of source (m) */
        const int O                    /**< twice post-Newtonian order */
        )
{
    unsigned int i;
    // Copy Phi TimeSeries and apply a shift of Pi
    REAL8TimeSeries *PhiNew = XLALCutREAL8TimeSeries(Phi, 0, Phi->data->length);
    for(i=0; i < Phi->data->length; i++)
        PhiNew->data->data[i] += LAL_PI;
    // Generate m > 0 mode
    COMPLEX16TimeSeries *hlm = XLALSimInspiralPrecessingPNMode21(V, PhiNew,
            S1x, S1y, S1z, S2x, S2y, S2z, LNhatx, LNhaty, LNhatz, m1, m2, r, O);
    if( !hlm ) XLAL_ERROR_NULL(XLAL_EFUNC);
    // Apply conjugation and (-1)^l
    for(i=0; i < hlm->data->length; i++)
        hlm->data->data[i] = conj(hlm->data->data[i]);

    return hlm;
}

/**
 * Computes h(2,0) mode of spherical harmonic decomposition of
 * precessing post-Newtonian inspiral waveforms.
 *
 * This expression is derived to 1.5PN order in arXiv:0810.5336,
 * with the expressions provided in the Mathematica file ABFO09.m.
 * The expressions are further simplified in a NB found here:
 * FIXME - add URL of wiki page with NB deriving simpler form!
 */
COMPLEX16TimeSeries *XLALSimInspiralPrecessingPNMode20(
        const REAL8TimeSeries *V,      /**< post-Newtonian parameter */
        const REAL8TimeSeries *Phi,    /**< orbital phase */
        const REAL8TimeSeries *S1x,    /**< Spin1 x component in J frame */
        const REAL8TimeSeries *S1y,    /**< Spin1 y component in J frame */
        const REAL8TimeSeries *S1z,    /**< Spin1 z component in J frame */
        const REAL8TimeSeries *S2x,    /**< Spin2 x component in J frame */
        const REAL8TimeSeries *S2y,    /**< Spin2 y component in J frame */
        const REAL8TimeSeries *S2z,    /**< Spin2 z component in J frame */
        const REAL8TimeSeries *LNhatx, /**< unit orbital ang. mom. x comp. in J frame */
        const REAL8TimeSeries *LNhaty, /**< unit orbital ang. mom. y comp. in J frame */
        const REAL8TimeSeries *LNhatz, /**< unit orbital ang. mom. z comp. in J frame */
        const REAL8 m1,                /**< mass of companion 1 (kg) */
        const REAL8 m2,                /**< mass of companion 2 (kg) */
        const REAL8 r,                 /**< distance of source (m) */
        const int O                    /**< twice post-Newtonian order */
        )
{
    // NB: other series will be checked in XLALSimInspiralRadiationFrameToJFrame
    LAL_CHECK_VALID_SERIES(Phi, NULL);
    LAL_CHECK_CONSISTENT_TIME_SERIES(V, Phi, NULL);

    COMPLEX16TimeSeries *hlm = XLALCreateCOMPLEX16TimeSeries("h_20 mode",
            &V->epoch, 0., V->deltaT, &lalStrainUnit, V->data->length);
    if( !hlm ) XLAL_ERROR_NULL(XLAL_ENOMEM);

    size_t j;
    REAL8 M = m1+m2;
    REAL8 m1_sec = LAL_G_SI * m1 / LAL_C_SI / LAL_C_SI / LAL_C_SI;
    REAL8 m2_sec = LAL_G_SI * m2 / LAL_C_SI / LAL_C_SI / LAL_C_SI;
    REAL8 M_sec = m1_sec + m2_sec;
    REAL8 eta = m1*m2/M/M;
    REAL8 dm = (m1-m2)/M;
    REAL8 m1oM = m1/M;
    REAL8 m2oM = m2/M;
    REAL8 lalpi = LAL_PI;
    REAL8 prefac = -8.0*sqrt(LAL_PI/5.0)*M_sec*eta/r;
    COMPLEX16 term0, term1, term2, term3;
    COMPLEX16 expalpha, expmalpha, expphi, expmphi;
    COMPLEX16 exp2phi, expm2phi;
    COMPLEX16 S1perp, S2perp, conjS1perp, conjS2perp;
    REAL8 v, v2, iota, alpha, s1x, s1y, s1z, s2x, s2y, s2z;
    REAL8 cosiota, cosiota2, cos2iota;
    REAL8 siniota, siniota2, sin2iota;

    for(j=0; j < V->data->length; j++) {
        // Compute the dynamic quantities in the terms below
        v = V->data->data[j];
        v2 = v*v;
        expphi = cpolar(1., Phi->data->data[j]);
        expmphi = cpolar(1., - Phi->data->data[j]);
        exp2phi = cpolar(1., 2.*Phi->data->data[j]);
        expm2phi = cpolar(1., - 2.*Phi->data->data[j]);
        /* EvolveOrbit functions return spin components with mass prefactors.
         * Specifically, S1z = (m1/M)^2 * \chi_1 * \hat{S}_1
         * and similarly for other components.
         * Thus, we divide out the mass prefactors so s1z = \chi_1 * \hat{S}_1z
         */
        s1x = S1x->data->data[j] / m1oM / m1oM;
        s1y = S1y->data->data[j] / m1oM / m1oM;
        s1z = S1z->data->data[j] / m1oM / m1oM;
        s2x = S2x->data->data[j] / m2oM / m2oM;
        s2y = S2y->data->data[j] / m2oM / m2oM;
        s2z = S2z->data->data[j] / m2oM / m2oM;
        S1perp = crect(s1x, s1y);
        conjS1perp = crect(s1x, -s1y);
        S2perp = crect(s2x, s2y);
        conjS2perp = crect(s2x, -s2y);

        iota = acos(LNhatz->data->data[j]);
        alpha = atan2(LNhaty->data->data[j], LNhatx->data->data[j]);
        expalpha = cpolar(1., alpha);
        expmalpha = cpolar(1., - alpha);
        cosiota = cos(iota);
        cosiota2 = cos(2.*iota);
        cos2iota = cosiota*cosiota;
        siniota = sin(iota);
        siniota2 = sin(2.*iota);
        sin2iota = siniota*siniota;

        // zero the terms
        term0 = term1 = term2 = term3 = 0.;

        // Compute the terms up to the requested PN order
        switch(O) {
            default: // unsupported PN order
                XLALPrintError("XLAL Error - %s: PN order %d%s not supported\n", __func__, O/2, O%2?".5":"" );
                XLAL_ERROR_NULL(XLAL_EINVAL);
            case -1: // Use highest available PN order
            case 3:
                term3 = (42*(exp2phi + expm2phi)*lalpi*sin2iota)/5.
                        + cosiota*dm*(1.7 - 2*eta)*(expmphi - expphi)*siniota
                        + (7*expalpha*((conjS1perp - conjS2perp)*dm
                        *(-24*cos2iota + (2 + 12*cosiota + 10*cosiota2)*exp2phi
                        + (2 - 12*cosiota + 10*cosiota2)*expm2phi)*siniota
                        + (conjS1perp + conjS2perp)*(cos2iota*(-24 + 12*eta)
                        + (2 + 9*eta + cosiota2*(10 + eta) + cosiota
                        *(12 + 10*eta))*exp2phi + (2 + cosiota*(-12 - 10*eta)
                        + 9*eta + cosiota2*(10 + eta))*expm2phi)*siniota))/40.
                        + (7*expmalpha*(dm*(-24*cos2iota
                        + (2 - 12*cosiota + 10*cosiota2)*exp2phi
                        + (2 + 12*cosiota + 10*cosiota2)*expm2phi)
                        *(S1perp - S2perp)*siniota + (cos2iota*(-24 + 12*eta)
                        + (2 + cosiota*(-12 - 10*eta) + 9*eta
                        + cosiota2*(10 + eta))*exp2phi + (2 + 9*eta + cosiota2
                        *(10 + eta) + cosiota*(12 + 10*eta))*expm2phi)
                        *(S1perp + S2perp)*siniota))/40.
                        + (7*dm*(6 - 5*exp2phi - 5*expm2phi)*(s1z - s2z)
                        *siniota*siniota2)/10. + (7*(6 - 3*eta
                        + (-5 - eta/2.)*exp2phi + (-5 - eta/2.)*expm2phi)
                        *(s1z + s2z)*siniota*siniota2)/10.;
            case 2:
                term2 = (7*expalpha*((conjS1perp - conjS2perp)*((-1 + cosiota)
                        *expmphi + (-1 - cosiota)*expphi)
                        + (conjS1perp + conjS2perp)*dm*((-1 + cosiota)*expmphi
                        + (-1 - cosiota)*expphi)))/20.
                        + (7*expmalpha*(((1 + cosiota)*expmphi
                        + (1 - cosiota)*expphi)*(S1perp - S2perp)
                        + dm*((1 + cosiota)*expmphi + (1 - cosiota)*expphi)
                        *(S1perp + S2perp)))/20. + ((-107 + 55*eta)
                        *(exp2phi + expm2phi)*sin2iota)/10.
                        + (7*(expmphi - expphi)*(s1z - s2z)*siniota)/5.
                        + (7*dm*(expmphi - expphi)*(s1z + s2z)*siniota)/5.;
            case 1:
                term1 = (14*cosiota*dm*(-expmphi + expphi)*siniota)/5.;
            case 0:
                term0 = (21*(exp2phi + expm2phi)*sin2iota)/5.;
                break;
        }

        v = V->data->data[j];
        v2 = v*v;
        hlm->data->data[j] = prefac * v2
                * (term0 + v*(term1 + v*(term2 + v*(term3))));
    }

    return hlm;
}

/**
 * Computes h(3,3) mode of spherical harmonic decomposition of
 * precessing post-Newtonian inspiral waveforms.
 *
 * This expression is derived to 1.5PN order in arXiv:0810.5336,
 * with the expressions provided in the Mathematica file ABFO09.m.
 * The expressions are further simplified in a NB found here:
 * FIXME - add URL of wiki page with NB deriving simpler form!
 */
COMPLEX16TimeSeries *XLALSimInspiralPrecessingPNMode33(
        const REAL8TimeSeries *V,      /**< post-Newtonian parameter */
        const REAL8TimeSeries *Phi,    /**< orbital phase */
        const REAL8TimeSeries *S1x,    /**< Spin1 x component in J frame */
        const REAL8TimeSeries *S1y,    /**< Spin1 y component in J frame */
        const UNUSED REAL8TimeSeries *S1z, /**< Spin1 z component in J frame */
        const REAL8TimeSeries *S2x,    /**< Spin2 x component in J frame */
        const REAL8TimeSeries *S2y,    /**< Spin2 y component in J frame */
        const UNUSED REAL8TimeSeries *S2z, /**< Spin2 z component in J frame */
        const REAL8TimeSeries *LNhatx, /**< unit orbital ang. mom. x comp. in J frame */
        const REAL8TimeSeries *LNhaty, /**< unit orbital ang. mom. y comp. in J frame */
        const REAL8TimeSeries *LNhatz, /**< unit orbital ang. mom. z comp. in J frame */
        const REAL8 m1,                /**< mass of companion 1 (kg) */
        const REAL8 m2,                /**< mass of companion 2 (kg) */
        const REAL8 r,                 /**< distance of source (m) */
        const int O                    /**< twice post-Newtonian order */
        )
{
    // NB: other series will be checked in XLALSimInspiralRadiationFrameToJFrame
    LAL_CHECK_VALID_SERIES(Phi, NULL);
    LAL_CHECK_CONSISTENT_TIME_SERIES(V, Phi, NULL);

    COMPLEX16TimeSeries *hlm = XLALCreateCOMPLEX16TimeSeries("h_33 mode",
            &V->epoch, 0., V->deltaT, &lalStrainUnit, V->data->length);
    if( !hlm ) XLAL_ERROR_NULL(XLAL_ENOMEM);

    size_t j;
    REAL8 M = m1+m2;
    REAL8 m1_sec = LAL_G_SI * m1 / LAL_C_SI / LAL_C_SI / LAL_C_SI;
    REAL8 m2_sec = LAL_G_SI * m2 / LAL_C_SI / LAL_C_SI / LAL_C_SI;
    REAL8 M_sec = m1_sec + m2_sec;
    REAL8 eta = m1*m2/M/M;
    REAL8 dm = (m1-m2)/M;
    REAL8 m1oM = m1/M;
    REAL8 m2oM = m2/M;
    REAL8 prefac = -8.0*sqrt(LAL_PI/5.0)*M_sec*eta/r;
    COMPLEX16 term0, term1, term2, term3;
    COMPLEX16 expm2alpha, expm3alpha, expphi, expmphi;
    COMPLEX16 exp2phi, expm2phi, exp3phi, expm3phi;
    COMPLEX16 conjS1perp, conjS2perp;
    REAL8 v, v2, iota, alpha, s1x, s1y, s2x, s2y;
    REAL8 cos2iotao2, cos4iotao2, cos6iotao2;
    REAL8 siniota, sin2iota, siniotao2, sin2iotao2, sin4iotao2, sin6iotao2;

    for(j=0; j < V->data->length; j++) {
        // Compute the dynamic quantities in the terms below
        v = V->data->data[j];
        v2 = v*v;
        expphi = cpolar(1., Phi->data->data[j]);
        expmphi = cpolar(1., - Phi->data->data[j]);
        exp2phi = cpolar(1., 2.*Phi->data->data[j]);
        expm2phi = cpolar(1., - 2.*Phi->data->data[j]);
        exp3phi = cpolar(1., 3.*Phi->data->data[j]);
        expm3phi = cpolar(1., - 3.*Phi->data->data[j]);
        /* EvolveOrbit functions return spin components with mass prefactors.
         * Specifically, S1z = (m1/M)^2 * \chi_1 * \hat{S}_1
         * and similarly for other components.
         * Thus, we divide out the mass prefactors so s1z = \chi_1 * \hat{S}_1z
         */
        s1x = S1x->data->data[j] / m1oM / m1oM;
        s1y = S1y->data->data[j] / m1oM / m1oM;
        s2x = S2x->data->data[j] / m2oM / m2oM;
        s2y = S2y->data->data[j] / m2oM / m2oM;
        conjS1perp = crect(s1x, -s1y);
        conjS2perp = crect(s2x, -s2y);

        iota = acos(LNhatz->data->data[j]);
        alpha = atan2(LNhaty->data->data[j], LNhatx->data->data[j]);
        expm2alpha = cpolar(1., - 2.*alpha);
        expm3alpha = cpolar(1., - 3.*alpha);
        cos2iotao2 = cos(iota/2.)*cos(iota/2.);
        cos4iotao2 = cos2iotao2*cos2iotao2;
        cos6iotao2 = cos4iotao2*cos2iotao2;
        siniota = sin(iota);
        sin2iota = siniota*siniota;
        siniotao2 = sin(iota/2.);
        sin2iotao2 = siniotao2*siniotao2;
        sin4iotao2 = sin2iotao2*sin2iotao2;
        sin6iotao2 = sin4iotao2*sin2iotao2;

        // zero the terms
        term0 = term1 = term2 = term3 = 0.;

        // Compute the terms up to the requested PN order
        switch(O) {
            default: // unsupported PN order
                XLALPrintError("XLAL Error - %s: PN order %d%s not supported\n", __func__, O/2, O%2?".5":"" );
                XLAL_ERROR_NULL(XLAL_EINVAL);
            case -1: // Use highest available PN order
            case 3:
                term3 = (8*(conjS1perp + conjS2perp)*eta*expm2alpha*(cos4iotao2
                        *expm2phi - exp2phi*sin4iotao2))/9. 
                        + dm*expm3alpha*(cos6iotao2*(-4 + 2*eta)*expm3phi
                        + (0.07407407407407407 + eta/54.)*sin2iota
                        *(-(cos2iotao2*expmphi) + expphi*sin2iotao2)
                        + (4 - 2*eta)*exp3phi*sin6iotao2);
            case 2:
                term2 = (4*(1 - 3*eta)*expm3alpha*(cos4iotao2*expm2phi
                        - exp2phi*sin4iotao2)*siniota)/9.;
            case 1:
                term1 = dm*expm3alpha*(cos6iotao2*expm3phi + (sin2iota
                        *(cos2iotao2*expmphi - expphi*sin2iotao2))/36.
                        - exp3phi*sin6iotao2);
            case 0:
                term0 = 0.;
                break;
        }

        v = V->data->data[j];
        v2 = v*v;
        hlm->data->data[j] = prefac * v2
                * (term0 + v*(term1 + v*(term2 + v*(term3))));
    }

    return hlm;
}

/**
 * Computes h(3,-3) mode of spherical harmonic decomposition of
 * precessing post-Newtonian inspiral waveforms.
 *
 * This expression is derived to 1.5PN order in arXiv:0810.5336,
 * with the expressions provided in the Mathematica file ABFO09.m.
 * The expressions are further simplified in a NB found here:
 * FIXME - add URL of wiki page with NB deriving simpler form!
 *
 * We get negative modes via the symmetry (see text below Eq. 4.17 of 0810.5336)
 * h_{l-m} = (-1)^l h_{lm}^*(Phi + pi)
 */
COMPLEX16TimeSeries *XLALSimInspiralPrecessingPNMode3m3(
        const REAL8TimeSeries *V,      /**< post-Newtonian parameter */
        const REAL8TimeSeries *Phi,    /**< orbital phase */
        const REAL8TimeSeries *S1x,    /**< Spin1 x component in J frame */
        const REAL8TimeSeries *S1y,    /**< Spin1 y component in J frame */
        const REAL8TimeSeries *S1z,    /**< Spin1 z component in J frame */
        const REAL8TimeSeries *S2x,    /**< Spin2 x component in J frame */
        const REAL8TimeSeries *S2y,    /**< Spin2 y component in J frame */
        const REAL8TimeSeries *S2z,    /**< Spin2 z component in J frame */
        const REAL8TimeSeries *LNhatx, /**< unit orbital ang. mom. x comp. in J frame */
        const REAL8TimeSeries *LNhaty, /**< unit orbital ang. mom. y comp. in J frame */
        const REAL8TimeSeries *LNhatz, /**< unit orbital ang. mom. z comp. in J frame */
        const REAL8 m1,                /**< mass of companion 1 (kg) */
        const REAL8 m2,                /**< mass of companion 2 (kg) */
        const REAL8 r,                 /**< distance of source (m) */
        const int O                    /**< twice post-Newtonian order */
        )
{
    unsigned int i;
    // Copy Phi TimeSeries and apply a shift of Pi
    REAL8TimeSeries *PhiNew = XLALCutREAL8TimeSeries(Phi, 0, Phi->data->length);
    for(i=0; i < Phi->data->length; i++)
        PhiNew->data->data[i] += LAL_PI;
    // Generate m > 0 mode
    COMPLEX16TimeSeries *hlm = XLALSimInspiralPrecessingPNMode33(V, PhiNew,
            S1x, S1y, S1z, S2x, S2y, S2z, LNhatx, LNhaty, LNhatz, m1, m2, r, O);
    if( !hlm ) XLAL_ERROR_NULL(XLAL_EFUNC);
    // Apply conjugation and (-1)^l
    for(i=0; i < hlm->data->length; i++)
        hlm->data->data[i] = - conj(hlm->data->data[i]);

    return hlm;
}

/**
 * Computes h(3,2) mode of spherical harmonic decomposition of
 * precessing post-Newtonian inspiral waveforms.
 *
 * This expression is derived to 1.5PN order in arXiv:0810.5336,
 * with the expressions provided in the Mathematica file ABFO09.m.
 * The expressions are further simplified in a NB found here:
 * FIXME - add URL of wiki page with NB deriving simpler form!
 */
COMPLEX16TimeSeries *XLALSimInspiralPrecessingPNMode32(
        const REAL8TimeSeries *V,      /**< post-Newtonian parameter */
        const REAL8TimeSeries *Phi,    /**< orbital phase */
        const REAL8TimeSeries *S1x,    /**< Spin1 x component in J frame */
        const REAL8TimeSeries *S1y,    /**< Spin1 y component in J frame */
        const REAL8TimeSeries *S1z,    /**< Spin1 z component in J frame */
        const REAL8TimeSeries *S2x,    /**< Spin2 x component in J frame */
        const REAL8TimeSeries *S2y,    /**< Spin2 y component in J frame */
        const REAL8TimeSeries *S2z,    /**< Spin2 z component in J frame */
        const REAL8TimeSeries *LNhatx, /**< unit orbital ang. mom. x comp. in J frame */
        const REAL8TimeSeries *LNhaty, /**< unit orbital ang. mom. y comp. in J frame */
        const REAL8TimeSeries *LNhatz, /**< unit orbital ang. mom. z comp. in J frame */
        const REAL8 m1,                /**< mass of companion 1 (kg) */
        const REAL8 m2,                /**< mass of companion 2 (kg) */
        const REAL8 r,                 /**< distance of source (m) */
        const int O                    /**< twice post-Newtonian order */
        )
{
    // NB: other series will be checked in XLALSimInspiralRadiationFrameToJFrame
    LAL_CHECK_VALID_SERIES(Phi, NULL);
    LAL_CHECK_CONSISTENT_TIME_SERIES(V, Phi, NULL);

    COMPLEX16TimeSeries *hlm = XLALCreateCOMPLEX16TimeSeries("h_32 mode",
            &V->epoch, 0., V->deltaT, &lalStrainUnit, V->data->length);
    if( !hlm ) XLAL_ERROR_NULL(XLAL_ENOMEM);

    size_t j;
    REAL8 M = m1+m2;
    REAL8 m1_sec = LAL_G_SI * m1 / LAL_C_SI / LAL_C_SI / LAL_C_SI;
    REAL8 m2_sec = LAL_G_SI * m2 / LAL_C_SI / LAL_C_SI / LAL_C_SI;
    REAL8 M_sec = m1_sec + m2_sec;
    REAL8 eta = m1*m2/M/M;
    REAL8 dm = (m1-m2)/M;
    REAL8 m1oM = m1/M;
    REAL8 m2oM = m2/M;
    REAL8 prefac = -8.0*sqrt(LAL_PI/5.0)*M_sec*eta/r;
    COMPLEX16 term0, term1, term2, term3;
    COMPLEX16 expmalpha, expm2alpha, expphi, expmphi;
    COMPLEX16 exp2phi, expm2phi, exp3phi, expm3phi;
    COMPLEX16 conjS1perp, conjS2perp;
    REAL8 v, v2, iota, alpha, s1x, s1y, s1z, s2x, s2y, s2z;
    REAL8 cosiota, cos2iotao2, cos4iotao2;
    REAL8 siniota, siniotao2, sin2iotao2, sin4iotao2;

    for(j=0; j < V->data->length; j++) {
        // Compute the dynamic quantities in the terms below
        v = V->data->data[j];
        v2 = v*v;
        expphi = cpolar(1., Phi->data->data[j]);
        expmphi = cpolar(1., - Phi->data->data[j]);
        exp2phi = cpolar(1., 2.*Phi->data->data[j]);
        expm2phi = cpolar(1., - 2.*Phi->data->data[j]);
        exp3phi = cpolar(1., 3.*Phi->data->data[j]);
        expm3phi = cpolar(1., - 3.*Phi->data->data[j]);
        /* EvolveOrbit functions return spin components with mass prefactors.
         * Specifically, S1z = (m1/M)^2 * \chi_1 * \hat{S}_1
         * and similarly for other components.
         * Thus, we divide out the mass prefactors so s1z = \chi_1 * \hat{S}_1z
         */
        s1x = S1x->data->data[j] / m1oM / m1oM;
        s1y = S1y->data->data[j] / m1oM / m1oM;
        s1z = S1z->data->data[j] / m1oM / m1oM;
        s2x = S2x->data->data[j] / m2oM / m2oM;
        s2y = S2y->data->data[j] / m2oM / m2oM;
        s2z = S2z->data->data[j] / m2oM / m2oM;
        conjS1perp = crect(s1x, -s1y);
        conjS2perp = crect(s2x, -s2y);

        iota = acos(LNhatz->data->data[j]);
        alpha = atan2(LNhaty->data->data[j], LNhatx->data->data[j]);
        expmalpha = cpolar(1., - alpha);
        expm2alpha = cpolar(1., - 2.*alpha);
        cosiota = cos(iota);
        cos2iotao2 = cos(iota/2.)*cos(iota/2.);
        cos4iotao2 = cos2iotao2*cos2iotao2;
        siniota = sin(iota);
        siniotao2 = sin(iota/2.);
        sin2iotao2 = siniotao2*siniotao2;
        sin4iotao2 = sin2iotao2*sin2iotao2;

        // zero the terms
        term0 = term1 = term2 = term3 = 0.;

        // Compute the terms up to the requested PN order
        switch(O) {
            default: // unsupported PN order
                XLALPrintError("XLAL Error - %s: PN order %d%s not supported\n", __func__, O/2, O%2?".5":"" );
                XLAL_ERROR_NULL(XLAL_EINVAL);
            case -1: // Use highest available PN order
            case 3:
                term3 = 2*(conjS1perp + conjS2perp)*eta*expmalpha
                        *(cos2iotao2*expm2phi + exp2phi*sin2iotao2)*siniota
                        + expm2alpha*(-2*cos4iotao2*eta*expm2phi*(s1z + s2z)
                        + 2*eta*exp2phi*(s1z + s2z)*sin4iotao2 + dm
                        *(cos4iotao2*(-13.5 + (27*eta)/4.)*expm3phi
                        + cos2iotao2*((-1 + 3*cosiota)/6. + ((-1 + 3*cosiota)
                        *eta)/24.)*expmphi + ((-1 - 3*cosiota)/6.
                        + ((-1 - 3*cosiota)*eta)/24.)*expphi*sin2iotao2
                        + (-13.5 + (27*eta)/4.)*exp3phi*sin4iotao2)*siniota);
            case 2:
                term2 = (1 - 3*eta)*expm2alpha*(cos4iotao2*(2 - 3*cosiota)
                        *expm2phi + (2 + 3*cosiota)*exp2phi*sin4iotao2);
            case 1:
                term1 = (dm*expm2alpha*(54*cos4iotao2*expm3phi + cos2iotao2
                        *(1 - 3*cosiota)*expmphi + (1 + 3*cosiota)*expphi
                        *sin2iotao2 + 54*exp3phi*sin4iotao2)*siniota)/16.;
            case 0:
                term0 = 0.;
                break;
        }

        v = V->data->data[j];
        v2 = v*v;
        hlm->data->data[j] = prefac * v2
                * (term0 + v*(term1 + v*(term2 + v*(term3))));
    }

    return hlm;
}

/**
 * Computes h(3,-2) mode of spherical harmonic decomposition of
 * precessing post-Newtonian inspiral waveforms.
 *
 * This expression is derived to 1.5PN order in arXiv:0810.5336,
 * with the expressions provided in the Mathematica file ABFO09.m.
 * The expressions are further simplified in a NB found here:
 * FIXME - add URL of wiki page with NB deriving simpler form!
 *
 * We get negative modes via the symmetry (see text below Eq. 4.17 of 0810.5336)
 * h_{l-m} = (-1)^l h_{lm}^*(Phi + pi)
 */
COMPLEX16TimeSeries *XLALSimInspiralPrecessingPNMode3m2(
        const REAL8TimeSeries *V,      /**< post-Newtonian parameter */
        const REAL8TimeSeries *Phi,    /**< orbital phase */
        const REAL8TimeSeries *S1x,    /**< Spin1 x component in J frame */
        const REAL8TimeSeries *S1y,    /**< Spin1 y component in J frame */
        const REAL8TimeSeries *S1z,    /**< Spin1 z component in J frame */
        const REAL8TimeSeries *S2x,    /**< Spin2 x component in J frame */
        const REAL8TimeSeries *S2y,    /**< Spin2 y component in J frame */
        const REAL8TimeSeries *S2z,    /**< Spin2 z component in J frame */
        const REAL8TimeSeries *LNhatx, /**< unit orbital ang. mom. x comp. in J frame */
        const REAL8TimeSeries *LNhaty, /**< unit orbital ang. mom. y comp. in J frame */
        const REAL8TimeSeries *LNhatz, /**< unit orbital ang. mom. z comp. in J frame */
        const REAL8 m1,                /**< mass of companion 1 (kg) */
        const REAL8 m2,                /**< mass of companion 2 (kg) */
        const REAL8 r,                 /**< distance of source (m) */
        const int O                    /**< twice post-Newtonian order */
        )
{
    unsigned int i;
    // Copy Phi TimeSeries and apply a shift of Pi
    REAL8TimeSeries *PhiNew = XLALCutREAL8TimeSeries(Phi, 0, Phi->data->length);
    for(i=0; i < Phi->data->length; i++)
        PhiNew->data->data[i] += LAL_PI;
    // Generate m > 0 mode
    COMPLEX16TimeSeries *hlm = XLALSimInspiralPrecessingPNMode32(V, PhiNew,
            S1x, S1y, S1z, S2x, S2y, S2z, LNhatx, LNhaty, LNhatz, m1, m2, r, O);
    if( !hlm ) XLAL_ERROR_NULL(XLAL_EFUNC);
    // Apply conjugation and (-1)^l
    for(i=0; i < hlm->data->length; i++)
        hlm->data->data[i] = - conj(hlm->data->data[i]);

    return hlm;
}

/**
 * Computes h(3,1) mode of spherical harmonic decomposition of
 * precessing post-Newtonian inspiral waveforms.
 *
 * This expression is derived to 1.5PN order in arXiv:0810.5336,
 * with the expressions provided in the Mathematica file ABFO09.m.
 * The expressions are further simplified in a NB found here:
 * FIXME - add URL of wiki page with NB deriving simpler form!
 */
COMPLEX16TimeSeries *XLALSimInspiralPrecessingPNMode31(
        const REAL8TimeSeries *V,      /**< post-Newtonian parameter */
        const REAL8TimeSeries *Phi,    /**< orbital phase */
        const REAL8TimeSeries *S1x,    /**< Spin1 x component in J frame */
        const REAL8TimeSeries *S1y,    /**< Spin1 y component in J frame */
        const REAL8TimeSeries *S1z,    /**< Spin1 z component in J frame */
        const REAL8TimeSeries *S2x,    /**< Spin2 x component in J frame */
        const REAL8TimeSeries *S2y,    /**< Spin2 y component in J frame */
        const REAL8TimeSeries *S2z,    /**< Spin2 z component in J frame */
        const REAL8TimeSeries *LNhatx, /**< unit orbital ang. mom. x comp. in J frame */
        const REAL8TimeSeries *LNhaty, /**< unit orbital ang. mom. y comp. in J frame */
        const REAL8TimeSeries *LNhatz, /**< unit orbital ang. mom. z comp. in J frame */
        const REAL8 m1,                /**< mass of companion 1 (kg) */
        const REAL8 m2,                /**< mass of companion 2 (kg) */
        const REAL8 r,                 /**< distance of source (m) */
        const int O                    /**< twice post-Newtonian order */
        )
{
    // NB: other series will be checked in XLALSimInspiralRadiationFrameToJFrame
    LAL_CHECK_VALID_SERIES(Phi, NULL);
    LAL_CHECK_CONSISTENT_TIME_SERIES(V, Phi, NULL);

    COMPLEX16TimeSeries *hlm = XLALCreateCOMPLEX16TimeSeries("h_31 mode",
            &V->epoch, 0., V->deltaT, &lalStrainUnit, V->data->length);
    if( !hlm ) XLAL_ERROR_NULL(XLAL_ENOMEM);

    size_t j;
    REAL8 M = m1+m2;
    REAL8 m1_sec = LAL_G_SI * m1 / LAL_C_SI / LAL_C_SI / LAL_C_SI;
    REAL8 m2_sec = LAL_G_SI * m2 / LAL_C_SI / LAL_C_SI / LAL_C_SI;
    REAL8 M_sec = m1_sec + m2_sec;
    REAL8 eta = m1*m2/M/M;
    REAL8 dm = (m1-m2)/M;
    REAL8 m1oM = m1/M;
    REAL8 m2oM = m2/M;
    REAL8 prefac = -8.0*sqrt(LAL_PI/5.0)*M_sec*eta/r;
    COMPLEX16 term0, term1, term2, term3;
    COMPLEX16 expmalpha, expm2alpha, expphi, expmphi;
    COMPLEX16 exp2phi, expm2phi, exp3phi, expm3phi;
    COMPLEX16 S1perp, S2perp, conjS1perp, conjS2perp;
    REAL8 v, v2, iota, alpha, s1x, s1y, s1z, s2x, s2y, s2z;
    REAL8 cosiota, cosiota2, cos2iotao2;
    REAL8 siniota, sin2iota, siniotao2, sin2iotao2;

    for(j=0; j < V->data->length; j++) {
        // Compute the dynamic quantities in the terms below
        v = V->data->data[j];
        v2 = v*v;
        expphi = cpolar(1., Phi->data->data[j]);
        expmphi = cpolar(1., - Phi->data->data[j]);
        exp2phi = cpolar(1., 2.*Phi->data->data[j]);
        expm2phi = cpolar(1., - 2.*Phi->data->data[j]);
        exp3phi = cpolar(1., 3.*Phi->data->data[j]);
        expm3phi = cpolar(1., - 3.*Phi->data->data[j]);
        /* EvolveOrbit functions return spin components with mass prefactors.
         * Specifically, S1z = (m1/M)^2 * \chi_1 * \hat{S}_1
         * and similarly for other components.
         * Thus, we divide out the mass prefactors so s1z = \chi_1 * \hat{S}_1z
         */
        s1x = S1x->data->data[j] / m1oM / m1oM;
        s1y = S1y->data->data[j] / m1oM / m1oM;
        s1z = S1z->data->data[j] / m1oM / m1oM;
        s2x = S2x->data->data[j] / m2oM / m2oM;
        s2y = S2y->data->data[j] / m2oM / m2oM;
        s2z = S2z->data->data[j] / m2oM / m2oM;
        S1perp = crect(s1x, s1y);
        conjS1perp = crect(s1x, -s1y);
        S2perp = crect(s2x, s2y);
        conjS2perp = crect(s2x, -s2y);

        iota = acos(LNhatz->data->data[j]);
        alpha = atan2(LNhaty->data->data[j], LNhatx->data->data[j]);
        expmalpha = cpolar(1., - alpha);
        expm2alpha = cpolar(1., - 2.*alpha);
        cosiota = cos(iota);
        cosiota2 = cos(2.*iota);
        cos2iotao2 = cos(iota/2.)*cos(iota/2.);
        siniota = sin(iota);
        sin2iota = siniota*siniota;
        siniotao2 = sin(iota/2.);
        sin2iotao2 = siniotao2*siniotao2;

        // zero the terms
        term0 = term1 = term2 = term3 = 0.;

        // Compute the terms up to the requested PN order
        switch(O) {
            default: // unsupported PN order
                XLALPrintError("XLAL Error - %s: PN order %d%s not supported\n", __func__, O/2, O%2?".5":"" );
                XLAL_ERROR_NULL(XLAL_EINVAL);
            case -1: // Use highest available PN order
            case 3:
                term3 = 24*(conjS1perp + conjS2perp)*cos2iotao2*(-1 + cosiota)
                        *eta*exp2phi + 24*(conjS1perp + conjS2perp)
                        *(1 + cosiota)*eta*expm2phi*sin2iotao2 + expm2alpha
                        *(4*cos2iotao2*(1 + cosiota)*eta*expm2phi
                        *(-S1perp - S2perp) + 4*(-1 + cosiota)*eta*exp2phi
                        *(-S1perp - S2perp)*sin2iotao2) + expmalpha
                        *(dm*(-(cos2iotao2*(13 - 20*cosiota + 15*cosiota2)
                        *(4 + eta)*expmphi)/12. + ((13 + 20*cosiota
                        + 15*cosiota2)*(4 + eta)*expphi*sin2iotao2)/12.
                        + 135*(1 - eta/2.)*sin2iota*(-(cos2iotao2*expm3phi)
                        + exp3phi*sin2iotao2)) + 32*eta*(s1z + s2z)
                        *(-(cos2iotao2*expm2phi) - exp2phi*sin2iotao2)*siniota);
            case 2:
                term2 = 10*(1 - 3*eta)*expmalpha*(cos2iotao2*(1 - 3*cosiota)
                        *expm2phi + (-1 -3*cosiota)*exp2phi*sin2iotao2)*siniota;
            case 1:
                term1 = (dm*expmalpha*(cos2iotao2*(13 - 20*cosiota
                        + 15*cosiota2)*expmphi + (-13 - 20*cosiota - 15
                        *cosiota2)*expphi*sin2iotao2 + sin2iota*(270
                        *cos2iotao2*expm3phi - 270*exp3phi*sin2iotao2)))/8.;
            case 0:
                term0 = 0.;
                break;
        }

        v = V->data->data[j];
        v2 = v*v;
        hlm->data->data[j] = prefac * v2
                * (term0 + v*(term1 + v*(term2 + v*(term3))));
    }

    return hlm;
}

/**
 * Computes h(3,-1) mode of spherical harmonic decomposition of
 * precessing post-Newtonian inspiral waveforms.
 *
 * This expression is derived to 1.5PN order in arXiv:0810.5336,
 * with the expressions provided in the Mathematica file ABFO09.m.
 * The expressions are further simplified in a NB found here:
 * FIXME - add URL of wiki page with NB deriving simpler form!
 *
 * We get negative modes via the symmetry (see text below Eq. 4.17 of 0810.5336)
 * h_{l-m} = (-1)^l h_{lm}^*(Phi + pi)
 */
COMPLEX16TimeSeries *XLALSimInspiralPrecessingPNMode3m1(
        const REAL8TimeSeries *V,      /**< post-Newtonian parameter */
        const REAL8TimeSeries *Phi,    /**< orbital phase */
        const REAL8TimeSeries *S1x,    /**< Spin1 x component in J frame */
        const REAL8TimeSeries *S1y,    /**< Spin1 y component in J frame */
        const REAL8TimeSeries *S1z,    /**< Spin1 z component in J frame */
        const REAL8TimeSeries *S2x,    /**< Spin2 x component in J frame */
        const REAL8TimeSeries *S2y,    /**< Spin2 y component in J frame */
        const REAL8TimeSeries *S2z,    /**< Spin2 z component in J frame */
        const REAL8TimeSeries *LNhatx, /**< unit orbital ang. mom. x comp. in J frame */
        const REAL8TimeSeries *LNhaty, /**< unit orbital ang. mom. y comp. in J frame */
        const REAL8TimeSeries *LNhatz, /**< unit orbital ang. mom. z comp. in J frame */
        const REAL8 m1,                /**< mass of companion 1 (kg) */
        const REAL8 m2,                /**< mass of companion 2 (kg) */
        const REAL8 r,                 /**< distance of source (m) */
        const int O                    /**< twice post-Newtonian order */
        )
{
    unsigned int i;
    // Copy Phi TimeSeries and apply a shift of Pi
    REAL8TimeSeries *PhiNew = XLALCutREAL8TimeSeries(Phi, 0, Phi->data->length);
    for(i=0; i < Phi->data->length; i++)
        PhiNew->data->data[i] += LAL_PI;
    // Generate m > 0 mode
    COMPLEX16TimeSeries *hlm = XLALSimInspiralPrecessingPNMode31(V, PhiNew,
            S1x, S1y, S1z, S2x, S2y, S2z, LNhatx, LNhaty, LNhatz, m1, m2, r, O);
    if( !hlm ) XLAL_ERROR_NULL(XLAL_EFUNC);
    // Apply conjugation and (-1)^l
    for(i=0; i < hlm->data->length; i++)
        hlm->data->data[i] = - conj(hlm->data->data[i]);

    return hlm;
}

/**
 * Computes h(3,0) mode of spherical harmonic decomposition of
 * precessing post-Newtonian inspiral waveforms.
 *
 * This expression is derived to 1.5PN order in arXiv:0810.5336,
 * with the expressions provided in the Mathematica file ABFO09.m.
 * The expressions are further simplified in a NB found here:
 * FIXME - add URL of wiki page with NB deriving simpler form!
 */
COMPLEX16TimeSeries *XLALSimInspiralPrecessingPNMode30(
        const REAL8TimeSeries *V,      /**< post-Newtonian parameter */
        const REAL8TimeSeries *Phi,    /**< orbital phase */
        const REAL8TimeSeries *S1x,    /**< Spin1 x component in J frame */
        const REAL8TimeSeries *S1y,    /**< Spin1 y component in J frame */
        const REAL8TimeSeries *S1z,    /**< Spin1 z component in J frame */
        const REAL8TimeSeries *S2x,    /**< Spin2 x component in J frame */
        const REAL8TimeSeries *S2y,    /**< Spin2 y component in J frame */
        const REAL8TimeSeries *S2z,    /**< Spin2 z component in J frame */
        const REAL8TimeSeries *LNhatx, /**< unit orbital ang. mom. x comp. in J frame */
        const REAL8TimeSeries *LNhaty, /**< unit orbital ang. mom. y comp. in J frame */
        const REAL8TimeSeries *LNhatz, /**< unit orbital ang. mom. z comp. in J frame */
        const REAL8 m1,                /**< mass of companion 1 (kg) */
        const REAL8 m2,                /**< mass of companion 2 (kg) */
        const REAL8 r,                 /**< distance of source (m) */
        const int O                    /**< twice post-Newtonian order */
        )
{
    // NB: other series will be checked in XLALSimInspiralRadiationFrameToJFrame
    LAL_CHECK_VALID_SERIES(Phi, NULL);
    LAL_CHECK_CONSISTENT_TIME_SERIES(V, Phi, NULL);

    COMPLEX16TimeSeries *hlm = XLALCreateCOMPLEX16TimeSeries("h_30 mode",
            &V->epoch, 0., V->deltaT, &lalStrainUnit, V->data->length);
    if( !hlm ) XLAL_ERROR_NULL(XLAL_ENOMEM);

    size_t j;
    REAL8 M = m1+m2;
    REAL8 m1_sec = LAL_G_SI * m1 / LAL_C_SI / LAL_C_SI / LAL_C_SI;
    REAL8 m2_sec = LAL_G_SI * m2 / LAL_C_SI / LAL_C_SI / LAL_C_SI;
    REAL8 M_sec = m1_sec + m2_sec;
    REAL8 eta = m1*m2/M/M;
    REAL8 dm = (m1-m2)/M;
    REAL8 m1oM = m1/M;
    REAL8 m2oM = m2/M;
    REAL8 prefac = -8.0*sqrt(LAL_PI/5.0)*M_sec*eta/r;
    COMPLEX16 term0, term1, term2, term3;
    COMPLEX16 expalpha, expmalpha, expphi, expmphi;
    COMPLEX16 exp2phi, expm2phi, exp3phi, expm3phi;
    COMPLEX16 S1perp, S2perp, conjS1perp, conjS2perp;
    REAL8 v, v2, iota, alpha, s1x, s1y, s1z, s2x, s2y, s2z;
    REAL8 cosiota, cos2iota, siniota, sin2iota, sin3iota;

    for(j=0; j < V->data->length; j++) {
        // Compute the dynamic quantities in the terms below
        v = V->data->data[j];
        v2 = v*v;
        expphi = cpolar(1., Phi->data->data[j]);
        expmphi = cpolar(1., - Phi->data->data[j]);
        exp2phi = cpolar(1., 2.*Phi->data->data[j]);
        expm2phi = cpolar(1., - 2.*Phi->data->data[j]);
        exp3phi = cpolar(1., 3.*Phi->data->data[j]);
        expm3phi = cpolar(1., - 3.*Phi->data->data[j]);
        /* EvolveOrbit functions return spin components with mass prefactors.
         * Specifically, S1z = (m1/M)^2 * \chi_1 * \hat{S}_1
         * and similarly for other components.
         * Thus, we divide out the mass prefactors so s1z = \chi_1 * \hat{S}_1z
         */
        s1x = S1x->data->data[j] / m1oM / m1oM;
        s1y = S1y->data->data[j] / m1oM / m1oM;
        s1z = S1z->data->data[j] / m1oM / m1oM;
        s2x = S2x->data->data[j] / m2oM / m2oM;
        s2y = S2y->data->data[j] / m2oM / m2oM;
        s2z = S2z->data->data[j] / m2oM / m2oM;
        S1perp = crect(s1x, s1y);
        conjS1perp = crect(s1x, -s1y);
        S2perp = crect(s2x, s2y);
        conjS2perp = crect(s2x, -s2y);

        iota = acos(LNhatz->data->data[j]);
        alpha = atan2(LNhaty->data->data[j], LNhatx->data->data[j]);
        expalpha = cpolar(1., alpha);
        expmalpha = cpolar(1., - alpha);
        cosiota = cos(iota);
        cos2iota = cosiota*cosiota;
        siniota = sin(iota);
        sin2iota = siniota*siniota;
        sin3iota = sin2iota*siniota;

        // zero the terms
        term0 = term1 = term2 = term3 = 0.;

        // Compute the terms up to the requested PN order
        switch(O) {
            default: // unsupported PN order
                XLALPrintError("XLAL Error - %s: PN order %d%s not supported\n", __func__, O/2, O%2?".5":"" );
                XLAL_ERROR_NULL(XLAL_EINVAL);
            case -1: // Use highest available PN order
            case 3:
                term3 = (5*eta*(exp2phi - expm2phi)*(s1z + s2z)*sin2iota)/4.
                        + (5*(conjS1perp + conjS2perp)*eta*expalpha
                        *((1 + cosiota)*exp2phi + (1 - cosiota)*expm2phi)
                        *siniota)/12. + (5*eta*((-1 + cosiota)*exp2phi
                        + (-1 - cosiota)*expm2phi)*expmalpha*(S1perp + S2perp)
                        *siniota)/12. + (5*dm*((-270*exp3phi - 270*expm3phi
                        + 4*expmphi + 4*expphi + eta*(135*exp3phi
                        + 135*expm3phi + expmphi + expphi))*sin3iota
                        + cos2iota*(-16*expmphi + eta*(-4*expmphi - 4*expphi)
                        - 16*expphi)*siniota))/288.;
            case 2:
                term2 = ((25*cosiota*exp2phi)/24. - (25*cosiota*expm2phi)/24.)
                        *sin2iota + eta*((-25*cosiota*exp2phi)/8.
                        + (25*cosiota*expm2phi)/8.)*sin2iota;
            case 1:
                term1 = dm*(((75*exp3phi)/64. + (75*expm3phi)/64. - (5*expmphi)
                        /192. - (5*expphi)/192.)*sin3iota
                        + ((5*cos2iota*expmphi)/48. + (5*cos2iota*expphi)/48.)
                        *siniota);
            case 0:
                term0 = 0.;
                break;
        }

        v = V->data->data[j];
        v2 = v*v;
        hlm->data->data[j] = prefac * v2
                * (term0 + v*(term1 + v*(term2 + v*(term3))));
    }

    return hlm;
}

/**
 * Computes h(4,4) mode of spherical harmonic decomposition of
 * precessing post-Newtonian inspiral waveforms.
 *
 * This expression is derived to 1.5PN order in arXiv:0810.5336,
 * with the expressions provided in the Mathematica file ABFO09.m.
 * The expressions are further simplified in a NB found here:
 * FIXME - add URL of wiki page with NB deriving simpler form!
 */
COMPLEX16TimeSeries *XLALSimInspiralPrecessingPNMode44(
        const REAL8TimeSeries *V,      /**< post-Newtonian parameter */
        const REAL8TimeSeries *Phi,    /**< orbital phase */
        const UNUSED REAL8TimeSeries *S1x, /**< Spin1 x component in J frame */
        const UNUSED REAL8TimeSeries *S1y, /**< Spin1 y component in J frame */
        const UNUSED REAL8TimeSeries *S1z, /**< Spin1 z component in J frame */
        const UNUSED REAL8TimeSeries *S2x, /**< Spin2 x component in J frame */
        const UNUSED REAL8TimeSeries *S2y, /**< Spin2 y component in J frame */
        const UNUSED REAL8TimeSeries *S2z, /**< Spin2 z component in J frame */
        const REAL8TimeSeries *LNhatx, /**< unit orbital ang. mom. x comp. in J frame */
        const REAL8TimeSeries *LNhaty, /**< unit orbital ang. mom. y comp. in J frame */
        const REAL8TimeSeries *LNhatz, /**< unit orbital ang. mom. z comp. in J frame */
        const REAL8 m1,                /**< mass of companion 1 (kg) */
        const REAL8 m2,                /**< mass of companion 2 (kg) */
        const REAL8 r,                 /**< distance of source (m) */
        const int O                    /**< twice post-Newtonian order */
        )
{
    // NB: other series will be checked in XLALSimInspiralRadiationFrameToJFrame
    LAL_CHECK_VALID_SERIES(Phi, NULL);
    LAL_CHECK_CONSISTENT_TIME_SERIES(V, Phi, NULL);

    COMPLEX16TimeSeries *hlm = XLALCreateCOMPLEX16TimeSeries("h_44 mode",
            &V->epoch, 0., V->deltaT, &lalStrainUnit, V->data->length);
    if( !hlm ) XLAL_ERROR_NULL(XLAL_ENOMEM);

    size_t j;
    REAL8 M = m1+m2;
    REAL8 m1_sec = LAL_G_SI * m1 / LAL_C_SI / LAL_C_SI / LAL_C_SI;
    REAL8 m2_sec = LAL_G_SI * m2 / LAL_C_SI / LAL_C_SI / LAL_C_SI;
    REAL8 M_sec = m1_sec + m2_sec;
    REAL8 eta = m1*m2/M/M;
    REAL8 dm = (m1-m2)/M;
    REAL8 prefac = -8.0*sqrt(LAL_PI/5.0)*M_sec*eta/r;
    COMPLEX16 term0, term1, term2, term3;
    COMPLEX16 expm4alpha, expphi, expmphi;
    COMPLEX16 exp2phi, expm2phi, exp3phi, expm3phi, exp4phi, expm4phi;
    REAL8 v, v2, iota, alpha;
    REAL8 cos2iotao2, cos4iotao2, cos6iotao2, cos8iotao2;
    REAL8 siniota, sin2iota, sin3iota, siniotao2, sin2iotao2;
    REAL8 sin4iotao2, sin6iotao2, sin8iotao2;

    for(j=0; j < V->data->length; j++) {
        // Compute the dynamic quantities in the terms below
        v = V->data->data[j];
        v2 = v*v;
        expphi = cpolar(1., Phi->data->data[j]);
        expmphi = cpolar(1., - Phi->data->data[j]);
        exp2phi = cpolar(1., 2.*Phi->data->data[j]);
        expm2phi = cpolar(1., - 2.*Phi->data->data[j]);
        exp3phi = cpolar(1., 3.*Phi->data->data[j]);
        expm3phi = cpolar(1., - 3.*Phi->data->data[j]);
        exp4phi = cpolar(1., 4.*Phi->data->data[j]);
        expm4phi = cpolar(1., - 4.*Phi->data->data[j]);

        iota = acos(LNhatz->data->data[j]);
        alpha = atan2(LNhaty->data->data[j], LNhatx->data->data[j]);
        expm4alpha = cpolar(1., - 4.*alpha);
        cos2iotao2 = cos(iota/2.)*cos(iota/2.);
        cos4iotao2 = cos2iotao2*cos2iotao2;
        cos6iotao2 = cos4iotao2*cos2iotao2;
        cos8iotao2 = cos6iotao2*cos2iotao2;
        siniota = sin(iota);
        sin2iota = siniota*siniota;
        sin3iota = sin2iota*siniota;
        siniotao2 = sin(iota/2.);
        sin2iotao2 = siniotao2*siniotao2;
        sin4iotao2 = sin2iotao2*sin2iotao2;
        sin6iotao2 = sin4iotao2*sin2iotao2;
        sin8iotao2 = sin6iotao2*sin2iotao2;

        // zero the terms
        term0 = term1 = term2 = term3 = 0.;

        // Compute the terms up to the requested PN order
        switch(O) {
            default: // unsupported PN order
                XLALPrintError("XLAL Error - %s: PN order %d%s not supported\n", __func__, O/2, O%2?".5":"" );
                XLAL_ERROR_NULL(XLAL_EINVAL);
            case -1: // Use highest available PN order
            case 3:
                term3 = dm*(-1 + 2*eta)*expm4alpha*((3*(cos2iotao2*expmphi
                        + expphi*sin2iotao2)*sin3iota)/640.
                        + (81*(cos6iotao2*expm3phi + exp3phi*sin6iotao2)
                        *siniota)/160.);
            case 2:
                term2 = (-1 + 3*eta)*expm4alpha*(cos8iotao2*expm4phi
                        + (sin2iota*(cos4iotao2*expm2phi + exp2phi*sin4iotao2))
                        /16. + exp4phi*sin8iotao2);
            case 1:
                term1 = 0.;
            case 0:
                term0 = 0.;
                break;
        }

        v = V->data->data[j];
        v2 = v*v;
        hlm->data->data[j] = prefac * v2
                * (term0 + v*(term1 + v*(term2 + v*(term3))));
    }

    return hlm;
}

/**
 * Computes h(4,-4) mode of spherical harmonic decomposition of
 * precessing post-Newtonian inspiral waveforms.
 *
 * This expression is derived to 1.5PN order in arXiv:0810.5336,
 * with the expressions provided in the Mathematica file ABFO09.m.
 * The expressions are further simplified in a NB found here:
 * FIXME - add URL of wiki page with NB deriving simpler form!
 *
 * We get negative modes via the symmetry (see text below Eq. 4.17 of 0810.5336)
 * h_{l-m} = (-1)^l h_{lm}^*(Phi + pi)
 */
COMPLEX16TimeSeries *XLALSimInspiralPrecessingPNMode4m4(
        const REAL8TimeSeries *V,      /**< post-Newtonian parameter */
        const REAL8TimeSeries *Phi,    /**< orbital phase */
        const REAL8TimeSeries *S1x,    /**< Spin1 x component in J frame */
        const REAL8TimeSeries *S1y,    /**< Spin1 y component in J frame */
        const REAL8TimeSeries *S1z,    /**< Spin1 z component in J frame */
        const REAL8TimeSeries *S2x,    /**< Spin2 x component in J frame */
        const REAL8TimeSeries *S2y,    /**< Spin2 y component in J frame */
        const REAL8TimeSeries *S2z,    /**< Spin2 z component in J frame */
        const REAL8TimeSeries *LNhatx, /**< unit orbital ang. mom. x comp. in J frame */
        const REAL8TimeSeries *LNhaty, /**< unit orbital ang. mom. y comp. in J frame */
        const REAL8TimeSeries *LNhatz, /**< unit orbital ang. mom. z comp. in J frame */
        const REAL8 m1,                /**< mass of companion 1 (kg) */
        const REAL8 m2,                /**< mass of companion 2 (kg) */
        const REAL8 r,                 /**< distance of source (m) */
        const int O                    /**< twice post-Newtonian order */
        )
{
    unsigned int i;
    // Copy Phi TimeSeries and apply a shift of Pi
    REAL8TimeSeries *PhiNew = XLALCutREAL8TimeSeries(Phi, 0, Phi->data->length);
    for(i=0; i < Phi->data->length; i++)
        PhiNew->data->data[i] += LAL_PI;
    // Generate m > 0 mode
    COMPLEX16TimeSeries *hlm = XLALSimInspiralPrecessingPNMode44(V, PhiNew,
            S1x, S1y, S1z, S2x, S2y, S2z, LNhatx, LNhaty, LNhatz, m1, m2, r, O);
    if( !hlm ) XLAL_ERROR_NULL(XLAL_EFUNC);
    // Apply conjugation and (-1)^l
    for(i=0; i < hlm->data->length; i++)
        hlm->data->data[i] = conj(hlm->data->data[i]);

    return hlm;
}

/**
 * Computes h(4,3) mode of spherical harmonic decomposition of
 * precessing post-Newtonian inspiral waveforms.
 *
 * This expression is derived to 1.5PN order in arXiv:0810.5336,
 * with the expressions provided in the Mathematica file ABFO09.m.
 * The expressions are further simplified in a NB found here:
 * FIXME - add URL of wiki page with NB deriving simpler form!
 */
COMPLEX16TimeSeries *XLALSimInspiralPrecessingPNMode43(
        const REAL8TimeSeries *V,      /**< post-Newtonian parameter */
        const REAL8TimeSeries *Phi,    /**< orbital phase */
        const UNUSED REAL8TimeSeries *S1x, /**< Spin1 x component in J frame */
        const UNUSED REAL8TimeSeries *S1y, /**< Spin1 y component in J frame */
        const UNUSED REAL8TimeSeries *S1z, /**< Spin1 z component in J frame */
        const UNUSED REAL8TimeSeries *S2x, /**< Spin2 x component in J frame */
        const UNUSED REAL8TimeSeries *S2y, /**< Spin2 y component in J frame */
        const UNUSED REAL8TimeSeries *S2z, /**< Spin2 z component in J frame */
        const REAL8TimeSeries *LNhatx, /**< unit orbital ang. mom. x comp. in J frame */
        const REAL8TimeSeries *LNhaty, /**< unit orbital ang. mom. y comp. in J frame */
        const REAL8TimeSeries *LNhatz, /**< unit orbital ang. mom. z comp. in J frame */
        const REAL8 m1,                /**< mass of companion 1 (kg) */
        const REAL8 m2,                /**< mass of companion 2 (kg) */
        const REAL8 r,                 /**< distance of source (m) */
        const int O                    /**< twice post-Newtonian order */
        )
{
    // NB: other series will be checked in XLALSimInspiralRadiationFrameToJFrame
    LAL_CHECK_VALID_SERIES(Phi, NULL);
    LAL_CHECK_CONSISTENT_TIME_SERIES(V, Phi, NULL);

    COMPLEX16TimeSeries *hlm = XLALCreateCOMPLEX16TimeSeries("h_43 mode",
            &V->epoch, 0., V->deltaT, &lalStrainUnit, V->data->length);
    if( !hlm ) XLAL_ERROR_NULL(XLAL_ENOMEM);

    size_t j;
    REAL8 M = m1+m2;
    REAL8 m1_sec = LAL_G_SI * m1 / LAL_C_SI / LAL_C_SI / LAL_C_SI;
    REAL8 m2_sec = LAL_G_SI * m2 / LAL_C_SI / LAL_C_SI / LAL_C_SI;
    REAL8 M_sec = m1_sec + m2_sec;
    REAL8 eta = m1*m2/M/M;
    REAL8 dm = (m1-m2)/M;
    REAL8 prefac = -8.0*sqrt(LAL_PI/5.0)*M_sec*eta/r;
    COMPLEX16 term0, term1, term2, term3;
    COMPLEX16 expm3alpha, expphi, expmphi;
    COMPLEX16 exp2phi, expm2phi, exp3phi, expm3phi, exp4phi, expm4phi;
    REAL8 v, v2, iota, alpha;
    REAL8 cosiota, cos2iotao2, cos4iotao2, cos6iotao2;
    REAL8 siniota, sin2iota, siniotao2, sin2iotao2, sin4iotao2, sin6iotao2;

    for(j=0; j < V->data->length; j++) {
        // Compute the dynamic quantities in the terms below
        v = V->data->data[j];
        v2 = v*v;
        expphi = cpolar(1., Phi->data->data[j]);
        expmphi = cpolar(1., - Phi->data->data[j]);
        exp2phi = cpolar(1., 2.*Phi->data->data[j]);
        expm2phi = cpolar(1., - 2.*Phi->data->data[j]);
        exp3phi = cpolar(1., 3.*Phi->data->data[j]);
        expm3phi = cpolar(1., - 3.*Phi->data->data[j]);
        exp4phi = cpolar(1., 4.*Phi->data->data[j]);
        expm4phi = cpolar(1., - 4.*Phi->data->data[j]);

        iota = acos(LNhatz->data->data[j]);
        alpha = atan2(LNhaty->data->data[j], LNhatx->data->data[j]);
        expm3alpha = cpolar(1., - 3.*alpha);
        cosiota = cos(iota);
        cos2iotao2 = cos(iota/2.)*cos(iota/2.);
        cos4iotao2 = cos2iotao2*cos2iotao2;
        cos6iotao2 = cos4iotao2*cos2iotao2;
        siniota = sin(iota);
        sin2iota = siniota*siniota;
        siniotao2 = sin(iota/2.);
        sin2iotao2 = siniotao2*siniotao2;
        sin4iotao2 = sin2iotao2*sin2iotao2;
        sin6iotao2 = sin4iotao2*sin2iotao2;

        // zero the terms
        term0 = term1 = term2 = term3 = 0.;

        // Compute the terms up to the requested PN order
        switch(O) {
            default: // unsupported PN order
                XLALPrintError("XLAL Error - %s: PN order %d%s not supported\n", __func__, O/2, O%2?".5":"" );
                XLAL_ERROR_NULL(XLAL_EINVAL);
            case -1: // Use highest available PN order
            case 3:
                term3 = dm*(-1 + 2*eta)*expm3alpha*(cos6iotao2*(3 - 4*cosiota)
                        *expm3phi + sin2iota*(cos2iotao2
                        *(0.009259259259259259 - cosiota/27.)*expmphi
                        + (-0.009259259259259259 - cosiota/27.)*expphi
                        *sin2iotao2) + (-3 - 4*cosiota)*exp3phi*sin6iotao2);
            case 2:
                term2 = (20*(-1 + 3*eta)*expm3alpha*(cos4iotao2*(1 - 2*cosiota)
                        *expm2phi + 16*cos6iotao2*expm4phi
                        + (-1 - 2*cosiota)*exp2phi*sin4iotao2
                        - 16*exp4phi*sin6iotao2)*siniota)/81.;
            case 1:
                term1 = 0.;
            case 0:
                term0 = 0.;
                break;
        }

        v = V->data->data[j];
        v2 = v*v;
        hlm->data->data[j] = prefac * v2
                * (term0 + v*(term1 + v*(term2 + v*(term3))));
    }

    return hlm;
}

/**
 * Computes h(4,-3) mode of spherical harmonic decomposition of
 * precessing post-Newtonian inspiral waveforms.
 *
 * This expression is derived to 1.5PN order in arXiv:0810.5336,
 * with the expressions provided in the Mathematica file ABFO09.m.
 * The expressions are further simplified in a NB found here:
 * FIXME - add URL of wiki page with NB deriving simpler form!
 *
 * We get negative modes via the symmetry (see text below Eq. 4.17 of 0810.5336)
 * h_{l-m} = (-1)^l h_{lm}^*(Phi + pi)
 */
COMPLEX16TimeSeries *XLALSimInspiralPrecessingPNMode4m3(
        const REAL8TimeSeries *V,      /**< post-Newtonian parameter */
        const REAL8TimeSeries *Phi,    /**< orbital phase */
        const REAL8TimeSeries *S1x,    /**< Spin1 x component in J frame */
        const REAL8TimeSeries *S1y,    /**< Spin1 y component in J frame */
        const REAL8TimeSeries *S1z,    /**< Spin1 z component in J frame */
        const REAL8TimeSeries *S2x,    /**< Spin2 x component in J frame */
        const REAL8TimeSeries *S2y,    /**< Spin2 y component in J frame */
        const REAL8TimeSeries *S2z,    /**< Spin2 z component in J frame */
        const REAL8TimeSeries *LNhatx, /**< unit orbital ang. mom. x comp. in J frame */
        const REAL8TimeSeries *LNhaty, /**< unit orbital ang. mom. y comp. in J frame */
        const REAL8TimeSeries *LNhatz, /**< unit orbital ang. mom. z comp. in J frame */
        const REAL8 m1,                /**< mass of companion 1 (kg) */
        const REAL8 m2,                /**< mass of companion 2 (kg) */
        const REAL8 r,                 /**< distance of source (m) */
        const int O                    /**< twice post-Newtonian order */
        )
{
    unsigned int i;
    // Copy Phi TimeSeries and apply a shift of Pi
    REAL8TimeSeries *PhiNew = XLALCutREAL8TimeSeries(Phi, 0, Phi->data->length);
    for(i=0; i < Phi->data->length; i++)
        PhiNew->data->data[i] += LAL_PI;
    // Generate m > 0 mode
    COMPLEX16TimeSeries *hlm = XLALSimInspiralPrecessingPNMode43(V, PhiNew,
            S1x, S1y, S1z, S2x, S2y, S2z, LNhatx, LNhaty, LNhatz, m1, m2, r, O);
    if( !hlm ) XLAL_ERROR_NULL(XLAL_EFUNC);
    // Apply conjugation and (-1)^l
    for(i=0; i < hlm->data->length; i++)
        hlm->data->data[i] = conj(hlm->data->data[i]);

    return hlm;
}

/**
 * Computes h(4,2) mode of spherical harmonic decomposition of
 * precessing post-Newtonian inspiral waveforms.
 *
 * This expression is derived to 1.5PN order in arXiv:0810.5336,
 * with the expressions provided in the Mathematica file ABFO09.m.
 * The expressions are further simplified in a NB found here:
 * FIXME - add URL of wiki page with NB deriving simpler form!
 */
COMPLEX16TimeSeries *XLALSimInspiralPrecessingPNMode42(
        const REAL8TimeSeries *V,      /**< post-Newtonian parameter */
        const REAL8TimeSeries *Phi,    /**< orbital phase */
        const UNUSED REAL8TimeSeries *S1x, /**< Spin1 x component in J frame */
        const UNUSED REAL8TimeSeries *S1y, /**< Spin1 y component in J frame */
        const UNUSED REAL8TimeSeries *S1z, /**< Spin1 z component in J frame */
        const UNUSED REAL8TimeSeries *S2x, /**< Spin2 x component in J frame */
        const UNUSED REAL8TimeSeries *S2y, /**< Spin2 y component in J frame */
        const UNUSED REAL8TimeSeries *S2z, /**< Spin2 z component in J frame */
        const REAL8TimeSeries *LNhatx, /**< unit orbital ang. mom. x comp. in J frame */
        const REAL8TimeSeries *LNhaty, /**< unit orbital ang. mom. y comp. in J frame */
        const REAL8TimeSeries *LNhatz, /**< unit orbital ang. mom. z comp. in J frame */
        const REAL8 m1,                /**< mass of companion 1 (kg) */
        const REAL8 m2,                /**< mass of companion 2 (kg) */
        const REAL8 r,                 /**< distance of source (m) */
        const int O                    /**< twice post-Newtonian order */
        )
{
    // NB: other series will be checked in XLALSimInspiralRadiationFrameToJFrame
    LAL_CHECK_VALID_SERIES(Phi, NULL);
    LAL_CHECK_CONSISTENT_TIME_SERIES(V, Phi, NULL);

    COMPLEX16TimeSeries *hlm = XLALCreateCOMPLEX16TimeSeries("h_42 mode",
            &V->epoch, 0., V->deltaT, &lalStrainUnit, V->data->length);
    if( !hlm ) XLAL_ERROR_NULL(XLAL_ENOMEM);

    size_t j;
    REAL8 M = m1+m2;
    REAL8 m1_sec = LAL_G_SI * m1 / LAL_C_SI / LAL_C_SI / LAL_C_SI;
    REAL8 m2_sec = LAL_G_SI * m2 / LAL_C_SI / LAL_C_SI / LAL_C_SI;
    REAL8 M_sec = m1_sec + m2_sec;
    REAL8 eta = m1*m2/M/M;
    REAL8 dm = (m1-m2)/M;
    REAL8 prefac = -8.0*sqrt(LAL_PI/5.0)*M_sec*eta/r;
    COMPLEX16 term0, term1, term2, term3;
    COMPLEX16 expm2alpha, expphi, expmphi;
    COMPLEX16 exp2phi, expm2phi, exp3phi, expm3phi, exp4phi, expm4phi;
    REAL8 v, v2, iota, alpha;
    REAL8 cosiota, cosiota2, cos2iotao2, cos4iotao2;
    REAL8 siniota, sin2iota, siniotao2, sin2iotao2, sin4iotao2;

    for(j=0; j < V->data->length; j++) {
        // Compute the dynamic quantities in the terms below
        v = V->data->data[j];
        v2 = v*v;
        expphi = cpolar(1., Phi->data->data[j]);
        expmphi = cpolar(1., - Phi->data->data[j]);
        exp2phi = cpolar(1., 2.*Phi->data->data[j]);
        expm2phi = cpolar(1., - 2.*Phi->data->data[j]);
        exp3phi = cpolar(1., 3.*Phi->data->data[j]);
        expm3phi = cpolar(1., - 3.*Phi->data->data[j]);
        exp4phi = cpolar(1., 4.*Phi->data->data[j]);
        expm4phi = cpolar(1., - 4.*Phi->data->data[j]);

        iota = acos(LNhatz->data->data[j]);
        alpha = atan2(LNhaty->data->data[j], LNhatx->data->data[j]);
        expm2alpha = cpolar(1., - 2.*alpha);
        cosiota = cos(iota);
        cosiota2 = cos(2.*iota);
        cos2iotao2 = cos(iota/2.)*cos(iota/2.);
        cos4iotao2 = cos2iotao2*cos2iotao2;
        siniota = sin(iota);
        sin2iota = siniota*siniota;
        siniotao2 = sin(iota/2.);
        sin2iotao2 = siniotao2*siniotao2;
        sin4iotao2 = sin2iotao2*sin2iotao2;

        // zero the terms
        term0 = term1 = term2 = term3 = 0.;

        // Compute the terms up to the requested PN order
        switch(O) {
            default: // unsupported PN order
                XLALPrintError("XLAL Error - %s: PN order %d%s not supported\n", __func__, O/2, O%2?".5":"" );
                XLAL_ERROR_NULL(XLAL_EINVAL);
            case -1: // Use highest available PN order
            case 3:
                term3 = (dm*(1 - 2*eta)*expm2alpha*(1134*cos4iotao2
                        *(-1 + 2*cosiota)*expm3phi - 3*cos2iotao2
                        *(6 - 7*cosiota + 7*cosiota2)*expmphi
                        - 3*(6 + 7*cosiota + 7*cosiota2)*expphi*sin2iotao2
                        - 1134*(1 + 2*cosiota)*exp3phi*sin4iotao2)*siniota)/80.;
            case 2:
                term2 = (-1 + 3*eta)*expm2alpha*(cos4iotao2*(4.5 - 7*cosiota
                        + (7*cosiota2)/2.)*expm2phi
                        + (4.5 + 7*cosiota + (7*cosiota2)/2.)*exp2phi
                        *sin4iotao2 + sin2iota*(28*cos4iotao2*expm4phi
                        + 28*exp4phi*sin4iotao2));
            case 1:
                term1 = 0.;
            case 0:
                term0 = 0.;
                break;
        }

        v = V->data->data[j];
        v2 = v*v;
        hlm->data->data[j] = prefac * v2
                * (term0 + v*(term1 + v*(term2 + v*(term3))));
    }

    return hlm;
}

/**
 * Computes h(4,-2) mode of spherical harmonic decomposition of
 * precessing post-Newtonian inspiral waveforms.
 *
 * This expression is derived to 1.5PN order in arXiv:0810.5336,
 * with the expressions provided in the Mathematica file ABFO09.m.
 * The expressions are further simplified in a NB found here:
 * FIXME - add URL of wiki page with NB deriving simpler form!
 *
 * We get negative modes via the symmetry (see text below Eq. 4.17 of 0810.5336)
 * h_{l-m} = (-1)^l h_{lm}^*(Phi + pi)
 */
COMPLEX16TimeSeries *XLALSimInspiralPrecessingPNMode4m2(
        const REAL8TimeSeries *V,      /**< post-Newtonian parameter */
        const REAL8TimeSeries *Phi,    /**< orbital phase */
        const REAL8TimeSeries *S1x,    /**< Spin1 x component in J frame */
        const REAL8TimeSeries *S1y,    /**< Spin1 y component in J frame */
        const REAL8TimeSeries *S1z,    /**< Spin1 z component in J frame */
        const REAL8TimeSeries *S2x,    /**< Spin2 x component in J frame */
        const REAL8TimeSeries *S2y,    /**< Spin2 y component in J frame */
        const REAL8TimeSeries *S2z,    /**< Spin2 z component in J frame */
        const REAL8TimeSeries *LNhatx, /**< unit orbital ang. mom. x comp. in J frame */
        const REAL8TimeSeries *LNhaty, /**< unit orbital ang. mom. y comp. in J frame */
        const REAL8TimeSeries *LNhatz, /**< unit orbital ang. mom. z comp. in J frame */
        const REAL8 m1,                /**< mass of companion 1 (kg) */
        const REAL8 m2,                /**< mass of companion 2 (kg) */
        const REAL8 r,                 /**< distance of source (m) */
        const int O                    /**< twice post-Newtonian order */
        )
{
    unsigned int i;
    // Copy Phi TimeSeries and apply a shift of Pi
    REAL8TimeSeries *PhiNew = XLALCutREAL8TimeSeries(Phi, 0, Phi->data->length);
    for(i=0; i < Phi->data->length; i++)
        PhiNew->data->data[i] += LAL_PI;
    // Generate m > 0 mode
    COMPLEX16TimeSeries *hlm = XLALSimInspiralPrecessingPNMode42(V, PhiNew,
            S1x, S1y, S1z, S2x, S2y, S2z, LNhatx, LNhaty, LNhatz, m1, m2, r, O);
    if( !hlm ) XLAL_ERROR_NULL(XLAL_EFUNC);
    // Apply conjugation and (-1)^l
    for(i=0; i < hlm->data->length; i++)
        hlm->data->data[i] = conj(hlm->data->data[i]);

    return hlm;
}

/**
 * Computes h(4,1) mode of spherical harmonic decomposition of
 * precessing post-Newtonian inspiral waveforms.
 *
 * This expression is derived to 1.5PN order in arXiv:0810.5336,
 * with the expressions provided in the Mathematica file ABFO09.m.
 * The expressions are further simplified in a NB found here:
 * FIXME - add URL of wiki page with NB deriving simpler form!
 */
COMPLEX16TimeSeries *XLALSimInspiralPrecessingPNMode41(
        const REAL8TimeSeries *V,      /**< post-Newtonian parameter */
        const REAL8TimeSeries *Phi,    /**< orbital phase */
        const UNUSED REAL8TimeSeries *S1x, /**< Spin1 x component in J frame */
        const UNUSED REAL8TimeSeries *S1y, /**< Spin1 y component in J frame */
        const UNUSED REAL8TimeSeries *S1z, /**< Spin1 z component in J frame */
        const UNUSED REAL8TimeSeries *S2x, /**< Spin2 x component in J frame */
        const UNUSED REAL8TimeSeries *S2y, /**< Spin2 y component in J frame */
        const UNUSED REAL8TimeSeries *S2z, /**< Spin2 z component in J frame */
        const REAL8TimeSeries *LNhatx, /**< unit orbital ang. mom. x comp. in J frame */
        const REAL8TimeSeries *LNhaty, /**< unit orbital ang. mom. y comp. in J frame */
        const REAL8TimeSeries *LNhatz, /**< unit orbital ang. mom. z comp. in J frame */
        const REAL8 m1,                /**< mass of companion 1 (kg) */
        const REAL8 m2,                /**< mass of companion 2 (kg) */
        const REAL8 r,                 /**< distance of source (m) */
        const int O                    /**< twice post-Newtonian order */
        )
{
    // NB: other series will be checked in XLALSimInspiralRadiationFrameToJFrame
    LAL_CHECK_VALID_SERIES(Phi, NULL);
    LAL_CHECK_CONSISTENT_TIME_SERIES(V, Phi, NULL);

    COMPLEX16TimeSeries *hlm = XLALCreateCOMPLEX16TimeSeries("h_41 mode",
            &V->epoch, 0., V->deltaT, &lalStrainUnit, V->data->length);
    if( !hlm ) XLAL_ERROR_NULL(XLAL_ENOMEM);

    size_t j;
    REAL8 M = m1+m2;
    REAL8 m1_sec = LAL_G_SI * m1 / LAL_C_SI / LAL_C_SI / LAL_C_SI;
    REAL8 m2_sec = LAL_G_SI * m2 / LAL_C_SI / LAL_C_SI / LAL_C_SI;
    REAL8 M_sec = m1_sec + m2_sec;
    REAL8 eta = m1*m2/M/M;
    REAL8 dm = (m1-m2)/M;
    REAL8 prefac = -8.0*sqrt(LAL_PI/5.0)*M_sec*eta/r;
    COMPLEX16 term0, term1, term2, term3;
    COMPLEX16 expmalpha, expphi, expmphi;
    COMPLEX16 exp2phi, expm2phi, exp3phi, expm3phi, exp4phi, expm4phi;
    REAL8 v, v2, iota, alpha;
    REAL8 cosiota, cosiota2, cosiota3, cos2iotao2;
    REAL8 siniota, sin2iota, sin3iota, siniotao2, sin2iotao2;

    for(j=0; j < V->data->length; j++) {
        // Compute the dynamic quantities in the terms below
        v = V->data->data[j];
        v2 = v*v;
        expphi = cpolar(1., Phi->data->data[j]);
        expmphi = cpolar(1., - Phi->data->data[j]);
        exp2phi = cpolar(1., 2.*Phi->data->data[j]);
        expm2phi = cpolar(1., - 2.*Phi->data->data[j]);
        exp3phi = cpolar(1., 3.*Phi->data->data[j]);
        expm3phi = cpolar(1., - 3.*Phi->data->data[j]);
        exp4phi = cpolar(1., 4.*Phi->data->data[j]);
        expm4phi = cpolar(1., - 4.*Phi->data->data[j]);

        iota = acos(LNhatz->data->data[j]);
        alpha = atan2(LNhaty->data->data[j], LNhatx->data->data[j]);
        expmalpha = cpolar(1., - alpha);
        cosiota = cos(iota);
        cosiota2 = cos(2.*iota);
        cosiota3 = cos(3.*iota);
        cos2iotao2 = cos(iota/2.)*cos(iota/2.);
        siniota = sin(iota);
        sin2iota = siniota*siniota;
        sin3iota = sin2iota*siniota;
        siniotao2 = sin(iota/2.);
        sin2iotao2 = siniotao2*siniotao2;

        // zero the terms
        term0 = term1 = term2 = term3 = 0.;

        // Compute the terms up to the requested PN order
        switch(O) {
            default: // unsupported PN order
                XLALPrintError("XLAL Error - %s: PN order %d%s not supported\n", __func__, O/2, O%2?".5":"" );
                XLAL_ERROR_NULL(XLAL_EINVAL);
            case -1: // Use highest available PN order
            case 3:
                term3 = (dm*(1 - 2*eta)*expmalpha*(cos2iotao2*(-15 + 30*cosiota
                        - 21*cosiota2 + 14*cosiota3)*expmphi + (15 + 30*cosiota
                        + 21*cosiota2 + 14*cosiota3)*expphi*sin2iotao2
                        + sin2iota*(378*cos2iotao2*(-1 + 4*cosiota)*expm3phi
                        + 378*(1 + 4*cosiota)*exp3phi*sin2iotao2)))/8.;
            case 2:
                term2 = ((1 - 3*eta)*expmalpha*((-560*cos2iotao2*expm4phi
                        + 560*exp4phi*sin2iotao2)*sin3iota + (-10*cos2iotao2
                        *(6 - 7*cosiota + 7*cosiota2)*expm2phi
                        + 10*(6 + 7*cosiota + 7*cosiota2)*exp2phi*sin2iotao2)
                        *siniota))/3.;
            case 1:
                term1 = 0.;
            case 0:
                term0 = 0.;
                break;
        }

        v = V->data->data[j];
        v2 = v*v;
        hlm->data->data[j] = prefac * v2
                * (term0 + v*(term1 + v*(term2 + v*(term3))));
    }

    return hlm;
}

/**
 * Computes h(4,-1) mode of spherical harmonic decomposition of
 * precessing post-Newtonian inspiral waveforms.
 *
 * This expression is derived to 1.5PN order in arXiv:0810.5336,
 * with the expressions provided in the Mathematica file ABFO09.m.
 * The expressions are further simplified in a NB found here:
 * FIXME - add URL of wiki page with NB deriving simpler form!
 *
 * We get negative modes via the symmetry (see text below Eq. 4.17 of 0810.5336)
 * h_{l-m} = (-1)^l h_{lm}^*(Phi + pi)
 */
COMPLEX16TimeSeries *XLALSimInspiralPrecessingPNMode4m1(
        const REAL8TimeSeries *V,      /**< post-Newtonian parameter */
        const REAL8TimeSeries *Phi,    /**< orbital phase */
        const REAL8TimeSeries *S1x,    /**< Spin1 x component in J frame */
        const REAL8TimeSeries *S1y,    /**< Spin1 y component in J frame */
        const REAL8TimeSeries *S1z,    /**< Spin1 z component in J frame */
        const REAL8TimeSeries *S2x,    /**< Spin2 x component in J frame */
        const REAL8TimeSeries *S2y,    /**< Spin2 y component in J frame */
        const REAL8TimeSeries *S2z,    /**< Spin2 z component in J frame */
        const REAL8TimeSeries *LNhatx, /**< unit orbital ang. mom. x comp. in J frame */
        const REAL8TimeSeries *LNhaty, /**< unit orbital ang. mom. y comp. in J frame */
        const REAL8TimeSeries *LNhatz, /**< unit orbital ang. mom. z comp. in J frame */
        const REAL8 m1,                /**< mass of companion 1 (kg) */
        const REAL8 m2,                /**< mass of companion 2 (kg) */
        const REAL8 r,                 /**< distance of source (m) */
        const int O                    /**< twice post-Newtonian order */
        )
{
    unsigned int i;
    // Copy Phi TimeSeries and apply a shift of Pi
    REAL8TimeSeries *PhiNew = XLALCutREAL8TimeSeries(Phi, 0, Phi->data->length);
    for(i=0; i < Phi->data->length; i++)
        PhiNew->data->data[i] += LAL_PI;
    // Generate m > 0 mode
    COMPLEX16TimeSeries *hlm = XLALSimInspiralPrecessingPNMode41(V, PhiNew,
            S1x, S1y, S1z, S2x, S2y, S2z, LNhatx, LNhaty, LNhatz, m1, m2, r, O);
    if( !hlm ) XLAL_ERROR_NULL(XLAL_EFUNC);
    // Apply conjugation and (-1)^l
    for(i=0; i < hlm->data->length; i++)
        hlm->data->data[i] = conj(hlm->data->data[i]);

    return hlm;
}

/**
 * Computes h(4,0) mode of spherical harmonic decomposition of
 * precessing post-Newtonian inspiral waveforms.
 *
 * This expression is derived to 1.5PN order in arXiv:0810.5336,
 * with the expressions provided in the Mathematica file ABFO09.m.
 * The expressions are further simplified in a NB found here:
 * FIXME - add URL of wiki page with NB deriving simpler form!
 */
COMPLEX16TimeSeries *XLALSimInspiralPrecessingPNMode40(
        const REAL8TimeSeries *V,      /**< post-Newtonian parameter */
        const REAL8TimeSeries *Phi,    /**< orbital phase */
        const UNUSED REAL8TimeSeries *S1x, /**< Spin1 x component in J frame */
        const UNUSED REAL8TimeSeries *S1y, /**< Spin1 y component in J frame */
        const UNUSED REAL8TimeSeries *S1z, /**< Spin1 z component in J frame */
        const UNUSED REAL8TimeSeries *S2x, /**< Spin2 x component in J frame */
        const UNUSED REAL8TimeSeries *S2y, /**< Spin2 y component in J frame */
        const UNUSED REAL8TimeSeries *S2z, /**< Spin2 z component in J frame */
        const UNUSED REAL8TimeSeries *LNhatx, /**< unit orbital ang. mom. x comp. in J frame */
        const UNUSED REAL8TimeSeries *LNhaty, /**< unit orbital ang. mom. y comp. in J frame */
        const REAL8TimeSeries *LNhatz, /**< unit orbital ang. mom. z comp. in J frame */
        const REAL8 m1,                /**< mass of companion 1 (kg) */
        const REAL8 m2,                /**< mass of companion 2 (kg) */
        const REAL8 r,                 /**< distance of source (m) */
        const int O                    /**< twice post-Newtonian order */
        )
{
    // NB: other series will be checked in XLALSimInspiralRadiationFrameToJFrame
    LAL_CHECK_VALID_SERIES(Phi, NULL);
    LAL_CHECK_CONSISTENT_TIME_SERIES(V, Phi, NULL);

    COMPLEX16TimeSeries *hlm = XLALCreateCOMPLEX16TimeSeries("h_40 mode",
            &V->epoch, 0., V->deltaT, &lalStrainUnit, V->data->length);
    if( !hlm ) XLAL_ERROR_NULL(XLAL_ENOMEM);

    size_t j;
    REAL8 M = m1+m2;
    REAL8 m1_sec = LAL_G_SI * m1 / LAL_C_SI / LAL_C_SI / LAL_C_SI;
    REAL8 m2_sec = LAL_G_SI * m2 / LAL_C_SI / LAL_C_SI / LAL_C_SI;
    REAL8 M_sec = m1_sec + m2_sec;
    REAL8 eta = m1*m2/M/M;
    REAL8 dm = (m1-m2)/M;
    REAL8 prefac = -8.0*sqrt(LAL_PI/5.0)*M_sec*eta/r;
    COMPLEX16 term0, term1, term2, term3;
    COMPLEX16 expphi, expmphi;
    COMPLEX16 exp2phi, expm2phi, exp3phi, expm3phi, exp4phi, expm4phi;
    REAL8 v, v2, iota;
    REAL8 cosiota, cosiota2, siniota, sin2iota, sin3iota, sin4iota;

    for(j=0; j < V->data->length; j++) {
        // Compute the dynamic quantities in the terms below
        v = V->data->data[j];
        v2 = v*v;
        expphi = cpolar(1., Phi->data->data[j]);
        expmphi = cpolar(1., - Phi->data->data[j]);
        exp2phi = cpolar(1., 2.*Phi->data->data[j]);
        expm2phi = cpolar(1., - 2.*Phi->data->data[j]);
        exp3phi = cpolar(1., 3.*Phi->data->data[j]);
        expm3phi = cpolar(1., - 3.*Phi->data->data[j]);
        exp4phi = cpolar(1., 4.*Phi->data->data[j]);
        expm4phi = cpolar(1., - 4.*Phi->data->data[j]);

        iota = acos(LNhatz->data->data[j]);
        cosiota = cos(iota);
        cosiota2 = cos(2.*iota);
        siniota = sin(iota);
        sin2iota = siniota*siniota;
        sin3iota = sin2iota*siniota;
        sin4iota = sin2iota*sin2iota;

        // zero the terms
        term0 = term1 = term2 = term3 = 0.;

        // Compute the terms up to the requested PN order
        switch(O) {
            default: // unsupported PN order
                XLALPrintError("XLAL Error - %s: PN order %d%s not supported\n", __func__, O/2, O%2?".5":"" );
                XLAL_ERROR_NULL(XLAL_EINVAL);
            case -1: // Use highest available PN order
            case 3:
                term3 = (3*cosiota*dm*(1 - 2*eta)*(-378*(exp3phi - expm3phi)
                        *sin3iota + (expmphi - expphi - 7*cosiota2
                        *(-expmphi + expphi))*siniota))/4.;
            case 2:
                term2 = (1 - 3*eta)*((-25 - 35*cosiota2)*(exp2phi + expm2phi)
                        *sin2iota - 280*(exp4phi + expm4phi)*sin4iota);
            case 1:
                term1 = 0.;
            case 0:
                term0 = 0.;
                break;
        }

        v = V->data->data[j];
        v2 = v*v;
        hlm->data->data[j] = prefac * v2
                * (term0 + v*(term1 + v*(term2 + v*(term3))));
    }

    return hlm;
}
