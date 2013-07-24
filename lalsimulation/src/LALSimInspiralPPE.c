/*
 *  Copyright (C) 2007 Jolien Creighton, B.S. Sathyaprakash, Thomas Cokelaer
 *  Copyright (C) 2012 Leo Singer
 *  Copyright (C) 2012 Walter Del Pozzo
 *  Assembled from code found in:
 *    - LALInspiralStationaryPhaseApproximation2.c
 *    - LALInspiralChooseModel.c
 *    - LALInspiralSetup.c
 *    - LALSimInspiralTaylorF2ReducedSpin.c
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
#include <lal/LALSimInspiralTestGRParams.h>

/**
 * Find the least nonnegative integer power of 2 that is
 * greater than or equal to n.  Inspired by similar routine
 * in gstlal.
 */
static size_t CeilPow2(double n) {
    double signif;
    int exponent;
    signif = frexp(n, &exponent);
    if (signif < 0)
        return 1;
    if (signif == 0.5)
        exponent -= 1;
    return ((size_t) 1) << exponent;
}

/**
 * Computes the stationary phase approximation to the Fourier transform of
 * a chirp waveform with phase given by Eq.\eqref{eq_InspiralFourierPhase_f2}
 * and amplitude given by expanding \f$1/\sqrt{\dot{F}}\f$. If the PN order is
 * set to -1, then the highest implemented order is used.
 * The phasing is then modified for the testing parameters if present. See Li et al, 2011
 * \author W. Del Pozzo
 */
int XLALSimInspiralPPE(
        COMPLEX16FrequencySeries **htilde_out, /**< FD waveform */
        const REAL8 phic,                /**< coalescence GW phase (rad) */
        const REAL8 deltaF,              /**< frequency resolution */
        const REAL8 m1_SI,               /**< mass of companion 1 (kg) */
        const REAL8 m2_SI,               /**< mass of companion 2 (kg) */
        const REAL8 fStart,              /**< start GW frequency (Hz) */
        const REAL8 r,                   /**< distance of source (m) */
        const INT4 phaseO,               /**< twice PN phase order */
        const INT4 amplitudeO,            /**< twice PN amplitude order */
        const LALSimInspiralTestGRParam *extraParams /**< structure of testing parameters */
        )
{
    const int beta = 0., sigma = 0.;

    /* external: SI; internal: solar masses */
    const REAL8 m1 = m1_SI / LAL_MSUN_SI;
    const REAL8 m2 = m2_SI / LAL_MSUN_SI;
    const REAL8 m = m1 + m2;
    const REAL8 m_sec = m * LAL_MTSUN_SI;  /* total mass in seconds */
    const REAL8 eta = m1 * m2 / (m * m);
    const REAL8 piM = LAL_PI * m_sec;
    const REAL8 vISCO = 1. / sqrt(6.);
    const REAL8 fISCO = vISCO * vISCO * vISCO / piM;
    const REAL8 v0 = cbrt(piM * fStart);
    REAL8 aPPE=0.0,alphaPPE=0.0,bPPE=0.0,betaPPE=0.0;
    REAL8 shft, amp0, f_max;
    size_t i, n, iStart, iISCO;
    COMPLEX16 *data = NULL;
    LIGOTimeGPS tC = {0, 0};
    
    /* phasing coefficients */
    REAL8 pfaN = XLALSimInspiralTaylorF2_NewtCoeff(eta);
    REAL8 pfa0 = XLALSimInspiralTaylorF2_0PNCoeff();
    REAL8 pfa1 = XLALSimInspiralTaylorF2_05PNCoeff();
    REAL8 pfa2 = XLALSimInspiralTaylorF2_1PNCoeff(eta);
    REAL8 pfa3 = XLALSimInspiralTaylorF2_15PNCoeff(beta);
    REAL8 pfa4 = XLALSimInspiralTaylorF2_2PNCoeff(eta,sigma);
    REAL8 pfa5 = XLALSimInspiralTaylorF2_25PNCoeff(eta);
    REAL8 pfl5 = XLALSimInspiralTaylorF2_25PNLogCoeff(eta);
    REAL8 pfa6 = XLALSimInspiralTaylorF2_3PNCoeff(eta);
    REAL8 pfl6 = XLALSimInspiralTaylorF2_3PNLogCoeff();
    REAL8 pfa7 = XLALSimInspiralTaylorF2_35PNCoeff(eta);
    /* modify for the GR testing coefficients */

    if (extraParams!=NULL) 
    {
        if (XLALSimInspiralTestGRParamExists(extraParams,"aPPE")) aPPE = XLALSimInspiralGetTestGRParam(extraParams,"aPPE");
        if (XLALSimInspiralTestGRParamExists(extraParams,"alphaPPE")) aPPE = XLALSimInspiralGetTestGRParam(extraParams,"alphaPPE");
        if (XLALSimInspiralTestGRParamExists(extraParams,"bPPE")) aPPE = XLALSimInspiralGetTestGRParam(extraParams,"bPPE");
        if (XLALSimInspiralTestGRParamExists(extraParams,"betaPPE")) aPPE = XLALSimInspiralGetTestGRParam(extraParams,"betaPPE");
    }
    
    /* flux coefficients */
    const REAL8 FTaN = XLALSimInspiralPNFlux_0PNCoeff(eta);
    const REAL8 FTa2 = XLALSimInspiralPNFlux_2PNCoeff(eta);
    const REAL8 FTa3 = XLALSimInspiralPNFlux_3PNCoeff(eta);
    const REAL8 FTa4 = XLALSimInspiralPNFlux_4PNCoeff(eta);
    const REAL8 FTa5 = XLALSimInspiralPNFlux_5PNCoeff(eta);
    const REAL8 FTl6 = XLALSimInspiralPNFlux_6PNLogCoeff(eta);
    const REAL8 FTa6 = XLALSimInspiralPNFlux_6PNCoeff(eta);
    const REAL8 FTa7 = XLALSimInspiralPNFlux_7PNCoeff(eta);

    /* energy coefficients */
    const REAL8 dETaN = 2. * XLALSimInspiralPNEnergy_0PNCoeff(eta);
    const REAL8 dETa1 = 2. * XLALSimInspiralPNEnergy_2PNCoeff(eta);
    const REAL8 dETa2 = 3. * XLALSimInspiralPNEnergy_4PNCoeff(eta);
    const REAL8 dETa3 = 4. * XLALSimInspiralPNEnergy_6PNCoeff(eta);

    COMPLEX16FrequencySeries *htilde;

    /* Perform some initial checks */
    if (!htilde_out) XLAL_ERROR(XLAL_EFAULT);
    if (*htilde_out) XLAL_ERROR(XLAL_EFAULT);
    if (phic < 0) XLAL_ERROR(XLAL_EDOM);
    if (m1_SI <= 0) XLAL_ERROR(XLAL_EDOM);
    if (m2_SI <= 0) XLAL_ERROR(XLAL_EDOM);
    if (fStart <= 0) XLAL_ERROR(XLAL_EDOM);
    if (r <= 0) XLAL_ERROR(XLAL_EDOM);

    /* allocate htilde */
    f_max = CeilPow2(fISCO);
    n = f_max / deltaF + 1;
    XLALGPSAdd(&tC, -1 / deltaF);  /* coalesce at t=0 */
    htilde = XLALCreateCOMPLEX16FrequencySeries("htilde: FD waveform", &tC, 0.0, deltaF, &lalStrainUnit, n);
    if (!htilde) XLAL_ERROR(XLAL_EFUNC);
    memset(htilde->data->data, 0, n * sizeof(COMPLEX16));
    XLALUnitDivide(&htilde->sampleUnits, &htilde->sampleUnits, &lalSecondUnit);

    /* extrinsic parameters */
    //phi0 = phic;
    amp0 = -4. * m1 * m2 / r * LAL_MRSUN_SI * LAL_MTSUN_SI * sqrt(LAL_PI/12.L);
    shft = -LAL_TWOPI * (tC.gpsSeconds + 1e-9 * tC.gpsNanoSeconds);

    iStart = (size_t) ceil(fStart / deltaF);
    iISCO = (size_t) (fISCO / deltaF);
    iISCO = (iISCO < n) ? iISCO : n;  /* overflow protection; should we warn? */
    data = htilde->data->data;
    for (i = iStart; i < iISCO; i++) {
        const REAL8 f = i * deltaF;
        const REAL8 v = cbrt(piM*f);
        const REAL8 v2 = v * v;
        const REAL8 v3 = v * v2;
        const REAL8 v4 = v * v3;
        const REAL8 v5 = v * v4;
        const REAL8 v6 = v * v5;
        const REAL8 v7 = v * v6;
        const REAL8 v8 = v * v7;
        const REAL8 v9 = v * v8;
        const REAL8 v10 = v * v9;
        REAL8 phasing = 0.;
        REAL8 dEnergy = 0.;
        REAL8 flux = 0.;
        REAL8 amp;

        switch (phaseO)
        {
            case -1:
            case 7:
                phasing += pfa7 * v7;
            case 6:
                phasing += (pfa6 + pfl6 * log(4.*v) ) * v6;
            case 5:
                phasing += (pfa5 + pfl5 * log(v/v0)) * v5;
            case 4:
                phasing += pfa4 * v4;
            case 3:
                phasing += pfa3 * v3;
            case 2:
                phasing += pfa2 * v2;
			case 1:
				phasing += pfa1 * v;
            case 0:
                phasing += pfa0;
                break;
            default:
                XLALDestroyCOMPLEX16FrequencySeries(htilde);
                XLAL_ERROR(XLAL_ETYPE);
        }
        switch (amplitudeO)
        {
            case -1:
            case 7:
                flux += FTa7 * v7;
            case 6:
                flux += (FTa6 + FTl6*log(16.*v2)) * v6;
                dEnergy += dETa3 * v6;
            case 5:
                flux += FTa5 * v5;
            case 4:
                flux += FTa4 * v4;
                dEnergy += dETa2 * v4;
            case 3:
                flux += FTa3 * v3;
            case 2:
                flux += FTa2 * v2;
                dEnergy += dETa1 * v2;
            case 0:
                flux += 1.;
                dEnergy += 1.;
                break;
            default:
                XLALDestroyCOMPLEX16FrequencySeries(htilde);
                XLAL_ERROR(XLAL_ETYPE);
        }
        phasing *= pfaN / v5;
        flux *= FTaN * v10;
        dEnergy *= dETaN * v;

        phasing += shft * f -2.0* phic +betaPPE*pow(v,bPPE)*pow(eta,3.0*bPPE/5.0);
        amp = amp0 * sqrt(-dEnergy/flux) * v * (1.0+alphaPPE*pow(v,aPPE));
        data[i] = amp * cos(phasing - LAL_PI_4) - amp * sin(phasing- LAL_PI_4) * 1.0j;
    }

    *htilde_out = htilde;
    return XLAL_SUCCESS;
}

