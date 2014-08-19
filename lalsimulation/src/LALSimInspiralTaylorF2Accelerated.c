/*
 *  Copyright (C) 2014 Andrew Lundgren
 *  Based from code found in:
 *    - LALSimInspiralTaylorF2.c
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

#include <stdlib.h>
#include <math.h>
#include <lal/Date.h>
#include <lal/FrequencySeries.h>
#include <lal/LALConstants.h>
#include <lal/LALDatatypes.h>
#include <lal/LALSimInspiral.h>
#include <lal/Units.h>
#include <lal/XLALError.h>
#include "LALSimInspiralPNCoefficients.c"

#define MKL_Complex16 COMPLEX16
#include <mkl_vml_functions.h>


/**
 * NOTE: This is a work in progress. Right now it's just a proof of concept
 * that vectorized math functions can speed up the TaylorF2 waveform
 * Results are not guaranteed to be accurate or even correct.
 * There's also a lot of improvements that could be made in the use of
 * vectorized functions.
 */
int XLALSimInspiralTaylorF2Accelerated(
        COMPLEX16FrequencySeries **htilde_out, /**< FD waveform */
        const REAL8 phi_ref,                   /**< reference orbital phase (rad) */
        const REAL8 deltaF,                    /**< frequency resolution */
        const REAL8 m1_SI,                     /**< mass of companion 1 (kg) */
        const REAL8 m2_SI,                     /**< mass of companion 2 (kg) */
        const REAL8 S1z,                       /**<  z component of the spin of companion 1 */
        const REAL8 S2z,                       /**<  z component of the spin of companion 2  */
        const REAL8 fStart,                    /**< start GW frequency (Hz) */
        const REAL8 fEnd,                      /**< highest GW frequency (Hz) of waveform generation - if 0, end at Schwarzschild ISCO */
        const REAL8 f_ref,                     /**< Reference GW frequency (Hz) - if 0 reference point is coalescence */
        const REAL8 r,                         /**< distance of source (m) */
        const REAL8 quadparam1,                /**< quadrupole deformation parameter of body 1 (dimensionless, 1 for BH) */
        const REAL8 quadparam2,                /**< quadrupole deformation parameter of body 2 (dimensionless, 1 for BH) */
        const REAL8 lambda1,                   /**< (tidal deformation of body 1)/(mass of body 1)^5 */
        const REAL8 lambda2,                   /**< (tidal deformation of body 2)/(mass of body 2)^5 */
        const LALSimInspiralSpinOrder spinO,  /**< twice PN order of spin effects */
        const LALSimInspiralTidalOrder tideO,  /**< flag to control tidal effects */
        const INT4 phaseO,                     /**< twice PN phase order */
        const INT4 amplitudeO                  /**< twice PN amplitude order */
        )
{
    /* external: SI; internal: solar masses */
    const REAL8 m1 = m1_SI / LAL_MSUN_SI;
    const REAL8 m2 = m2_SI / LAL_MSUN_SI;
    const REAL8 m = m1 + m2;
    const REAL8 m_sec = m * LAL_MTSUN_SI;  /* total mass in seconds */
    const REAL8 eta = m1 * m2 / (m * m);
    const REAL8 piM = LAL_PI * m_sec;
    const REAL8 vISCO = 1. / sqrt(6.);
    const REAL8 fISCO = vISCO * vISCO * vISCO / piM;
    const REAL8 m1OverM = m1 / m;
    const REAL8 m2OverM = m2 / m;
    REAL8 shft, amp0, f_max;
    size_t i, n, iStart, wfLen;
    COMPLEX16 *data = NULL;
    LIGOTimeGPS tC = {0, 0};

    /* phasing coefficients */
    PNPhasingSeries pfa;
    XLALSimInspiralPNPhasing_F2(&pfa, m1, m2, S1z, S2z, S1z*S1z, S2z*S2z, S1z*S2z, quadparam1, quadparam2, spinO);

    REAL8 pfaN = 0.;
    REAL8 pfa2 = 0.; REAL8 pfa3 = 0.; REAL8 pfa4 = 0.;
    REAL8 pfa5 = 0.; REAL8 pfl5 = 0.;
    REAL8 pfa6 = 0.; REAL8 pfl6 = 0.;
    REAL8 pfa7 = 0.;

    switch (phaseO)
    {
        case -1:
        case 7:
            pfa7 = pfa.v[7];
        case 6:
            pfa6 = pfa.v[6];
            pfl6 = pfa.vlogv[6];
        case 5:
            pfa5 = pfa.v[5];
            pfl5 = pfa.vlogv[5];
        case 4:
            pfa4 = pfa.v[4];
        case 3:
            pfa3 = pfa.v[3];
        case 2:
            pfa2 = pfa.v[2];
	case 1:
	    XLALPrintWarning( "There is no 0.5PN phase coefficient, returning Newtonian-order phase.\n" );
        case 0:
            pfaN = pfa.v[0];
            break;
        default:
            XLAL_ERROR(XLAL_ETYPE, "Invalid phase PN order %s", phaseO);
    }

    if (!((0 == amplitudeO) || (-1 == amplitudeO)))
    {
        XLAL_ERROR(XLAL_ETYPE, "Invalid amplitude PN order %s", amplitudeO);
    }
    
    /* Generate tidal terms separately.
     * Enums specifying tidal order are in LALSimInspiralWaveformFlags.h
     */
    REAL8 pft10 = 0.;
    REAL8 pft12 = 0.;
    switch( tideO )
    {
	    case LAL_SIM_INSPIRAL_TIDAL_ORDER_ALL:
        case LAL_SIM_INSPIRAL_TIDAL_ORDER_6PN:
            pft12 = pfaN * ( XLALSimInspiralTaylorF2Phasing_12PNTidalCoeff(m1OverM, lambda1)
                            + XLALSimInspiralTaylorF2Phasing_12PNTidalCoeff(m2OverM, lambda2) );
        case LAL_SIM_INSPIRAL_TIDAL_ORDER_5PN:
            pft10 = pfaN * ( XLALSimInspiralTaylorF2Phasing_10PNTidalCoeff(m1OverM, lambda1)
                            + XLALSimInspiralTaylorF2Phasing_10PNTidalCoeff(m2OverM, lambda2) );
        case LAL_SIM_INSPIRAL_TIDAL_ORDER_0PN:
            break;
        default:
            XLAL_ERROR(XLAL_EINVAL, "Invalid tidal PN order %s", tideO);
    }

    /* The flux and energy coefficients below are used to compute SPA amplitude corrections */
    const REAL8 FTaN = XLALSimInspiralPNFlux_0PNCoeff(eta);
    const REAL8 dETaN = 2. * XLALSimInspiralPNEnergy_0PNCoeff(eta);

    COMPLEX16FrequencySeries *htilde;

    /* Perform some initial checks */
    if (!htilde_out) XLAL_ERROR(XLAL_EFAULT);
    if (*htilde_out) XLAL_ERROR(XLAL_EFAULT);
    if (m1_SI <= 0) XLAL_ERROR(XLAL_EDOM);
    if (m2_SI <= 0) XLAL_ERROR(XLAL_EDOM);
    if (fStart <= 0) XLAL_ERROR(XLAL_EDOM);
    if (f_ref < 0) XLAL_ERROR(XLAL_EDOM);
    if (r <= 0) XLAL_ERROR(XLAL_EDOM);

    /* allocate htilde */
    if ( fEnd == 0. ) // End at ISCO
        f_max = fISCO;
    else // End at user-specified freq.
        f_max = fEnd;
    n = (size_t) (f_max / deltaF + 1);
    XLALGPSAdd(&tC, -1 / deltaF);  /* coalesce at t=0 */
    htilde = XLALCreateCOMPLEX16FrequencySeries("htilde: FD waveform", &tC, 0.0, deltaF, &lalStrainUnit, n);
    if (!htilde) XLAL_ERROR(XLAL_EFUNC);
    memset(htilde->data->data, 0, n * sizeof(COMPLEX16));
    XLALUnitDivide(&htilde->sampleUnits, &htilde->sampleUnits, &lalSecondUnit);

    /* extrinsic parameters */
    amp0 = -4. * m1 * m2 / r * LAL_MRSUN_SI * LAL_MTSUN_SI * sqrt(LAL_PI/12.L);
    shft = LAL_TWOPI * (tC.gpsSeconds + 1e-9 * tC.gpsNanoSeconds);

    /* Fill with non-zero vals from fStart to f_max */
    iStart = (size_t) ceil(fStart / deltaF);
    data = htilde->data->data;

    /* Compute the SPA phase at the reference point
     * N.B. f_ref == 0 means we define the reference time/phase at "coalescence"
     * when the frequency approaches infinity. In that case,
     * the integrals Eq. 3.15 of arXiv:0907.0700 vanish when evaluated at
     * f_ref == infinity. If f_ref is finite, we must compute the SPA phase
     * evaluated at f_ref, store it as ref_phasing and subtract it off.
     */
    REAL8 ref_phasing = 0.;
    if( f_ref != 0. ) {
        const REAL8 vref = cbrt(piM*f_ref);
        const REAL8 logvref = log(vref);
        const REAL8 v2ref = vref * vref;
        const REAL8 v3ref = vref * v2ref;
        const REAL8 v4ref = vref * v3ref;
        const REAL8 v5ref = vref * v4ref;
        const REAL8 v6ref = vref * v5ref;
        const REAL8 v7ref = vref * v6ref;
        const REAL8 v8ref = vref * v7ref;
        const REAL8 v10ref = v2ref * v8ref;
        const REAL8 v12ref = v2ref * v10ref;
        ref_phasing += pfa7 * v7ref;
        ref_phasing += (pfa6 + pfl6 * logvref) * v6ref;
        ref_phasing += (pfa5 + pfl5 * logvref) * v5ref;
        ref_phasing += pfa4 * v4ref;
        ref_phasing += pfa3 * v3ref;
        ref_phasing += pfa2 * v2ref;
        ref_phasing += pfaN;

        /* Tidal terms in reference phasing */
        ref_phasing += pft12 * v12ref;
        ref_phasing += pft10 * v10ref;

        ref_phasing /= v5ref;
    } /* End of if(f_ref != 0) block */

    REAL8 *workspace;
    wfLen = n; /* With smart indexing, we can reduce this */
    workspace = (REAL8 *) malloc(4*wfLen*sizeof(REAL8));
    REAL8 *work_v3 = workspace;
    REAL8 *work_v = workspace + wfLen;
    REAL8 *work_logv = workspace + 2*wfLen;
    REAL8 *work_phasing = workspace + 3*wfLen;
    
    const REAL8 vStep = piM*deltaF;
    for (i = iStart; i < n; i++)
    {
        work_v3[i] = vStep*i;
    }

    vdCbrt(wfLen, work_v3, work_v);
    vdLn(wfLen, work_v, work_logv);

    for (i = iStart; i < n; i++)
    {
        const REAL8 f = i * deltaF;
        const REAL8 v = work_v[i];
        const REAL8 v2 = v * v;
        const REAL8 v3 = work_v3[i];
        const REAL8 v4 = v * v3;
        const REAL8 v5 = v * v4;
        const REAL8 v6 = v * v5;
        const REAL8 v7 = v * v6;
        const REAL8 v8 = v * v7;
        const REAL8 v10 = v2 * v8;
        const REAL8 v12 = v2 * v10;
        const REAL8 logv = work_logv[i];
        REAL8 phasing = 0.;
        REAL8 dEnergy = 0.;
        REAL8 flux = 0.;
        REAL8 amp;

        phasing += pfa7 * v7;
        phasing += (pfa6 + pfl6 * logv) * v6;
        phasing += (pfa5 + pfl5 * logv) * v5;
        phasing += pfa4 * v4;
        phasing += pfa3 * v3;
        phasing += pfa2 * v2;
        phasing += pfaN;

        /* Tidal terms in phasing */
        phasing += pft12 * v12;
        phasing += pft10 * v10;

        phasing /= v5;
        flux = FTaN * v10;
        dEnergy = dETaN * v;
        // Note the factor of 2 b/c phi_ref is orbital phase
        phasing += shft * f - 2.*phi_ref - ref_phasing - LAL_PI_4;
        work_phasing[i] = phasing;
        amp = amp0 * sqrt(-dEnergy/flux) * v;
        data[i] = amp;
        /* data[i] = amp*(cos(phasing) - 1.0j*sin(phasing)); */
    }

    /*vzCIS(wfLen, work_phasing, data);*/

    /* Note: The cos and sin could probably be done faster using vzCIS 
     * which does vzCIS z => exp(i*z), then doing a vetor mul (vzMul?) by amp
     * But I don't want to dig into the MKL complex structure right now
     */
    
    *htilde_out = htilde;
    free((void *) workspace);
    return XLAL_SUCCESS;
}
