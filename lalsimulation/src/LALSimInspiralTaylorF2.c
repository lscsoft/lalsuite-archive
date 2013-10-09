/*
 *  Copyright (C) 2007 Jolien Creighton, B.S. Sathyaprakash, Thomas Cokelaer
 *  Copyright (C) 2012 Leo Singer, Evan Ochsner, Les Wade, Alex Nitz
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
#include <lal/LALSimInspiralEOS.h>
#include <lal/Units.h>
#include <lal/XLALError.h>
#include "LALSimInspiralPNCoefficients.c"

/**
 * Computes the stationary phase approximation to the Fourier transform of
 * a chirp waveform with phase given by \eqref{eq_InspiralFourierPhase_f2}
 * and amplitude given by expanding \f$1/\sqrt{\dot{F}}\f$. If the PN order is
 * set to -1, then the highest implemented order is used.
 *
 * See arXiv:0810.5336 and arXiv:astro-ph/0504538 for spin corrections
 * to the phasing.
 * See arXiv:1303.7412 for spin-orbit phasing corrections at 3 and 3.5PN order
 */
int XLALSimInspiralTaylorF2(
        COMPLEX16FrequencySeries **htilde_out, /**< FD waveform */
        const REAL8 phic,                      /**< orbital coalescence phase (rad) */
        const REAL8 deltaF,                    /**< frequency resolution */
        const REAL8 m1_SI,                     /**< mass of companion 1 (kg) */
        const REAL8 m2_SI,                     /**< mass of companion 2 (kg) */
        const REAL8 S1z,                       /**<  z component of the spin of companion 1 */
        const REAL8 S2z,                       /**<  z component of the spin of companion 2  */
        const REAL8 fStart,                    /**< start GW frequency (Hz) */
        const REAL8 fEnd,                      /**< highest GW frequency (Hz) of waveform generation - if 0, end at Schwarzschild ISCO */
        const REAL8 r,                         /**< distance of source (m) */
        const REAL8 lambda1,                   /**< (tidal deformation of body 1)/(mass of body 1)^5 */
        const REAL8 lambda2,                   /**< (tidal deformation of body 2)/(mass of body 2)^5 */
        const LALSimInspiralSpinOrder spinO,  /**< twice PN order of spin effects */
        const LALSimInspiralTidalOrder tideO,  /**< flag to control tidal effects */
        const INT4 phaseO,                     /**< twice PN phase order */
        const INT4 amplitudeO,                 /**< twice PN amplitude order */
        const LALEquationOfState eos
        )
{
    const REAL8 lambda = -1987./3080.;
    const REAL8 theta = -11831./9240.;

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
    const REAL8 chi1 = m1 / m;
    const REAL8 chi2 = m2 / m;
    REAL8 fCONT = fEnd; 
    REAL8 lam1 = lambda1;
    REAL8 lam2 = lambda2;
    REAL8 shft, amp0, f_max;
    size_t i, n, iStart;
    COMPLEX16 *data = NULL;
    LIGOTimeGPS tC = {0, 0};

    /* phasing coefficients */
    const REAL8 pfaN = 3.L/(128.L * eta);
    const REAL8 pfa2 = 5.L*(743.L/84.L + 11.L * eta)/9.L;
    const REAL8 pfa3 = -16.L*LAL_PI;
    const REAL8 pfa4 = 5.L*(3058.673L/7.056L + 5429.L/7.L * eta
                     + 617.L * eta*eta)/72.L;
    const REAL8 pfa5 = 5.L/9.L * (7729.L/84.L - 13.L * eta) * LAL_PI;
    const REAL8 pfl5 = 5.L/3.L * (7729.L/84.L - 13.L * eta) * LAL_PI;
    const REAL8 pfa6 = (11583.231236531L/4.694215680L
                     - 640.L/3.L * LAL_PI * LAL_PI - 6848.L/21.L*LAL_GAMMA)
                     + eta * (-15335.597827L/3.048192L
                     + 2255./12. * LAL_PI * LAL_PI
                     - 1760./3.*theta +12320./9.*lambda)
                     + eta*eta * 76055.L/1728.L - eta*eta*eta * 127825.L/1296.L;
    const REAL8 pfl6 = -6848.L/21.L;
    const REAL8 pfa7 = LAL_PI * 5.L/756.L * ( 15419335.L/336.L
                     + 75703.L/2.L * eta - 14809.L * eta*eta);

    /* Spin coefficients */
    REAL8 pn_beta = 0;
    REAL8 pn_sigma = 0;
    REAL8 pn_gamma = 0;
    REAL8 psiSO3 = 0., psiSO35 = 0.; // 3PN and 3.5PN spin-orbit phasing terms

    REAL8 d = (m1 - m2) / (m1 + m2);
    REAL8 xs = .5 * (S1z + S2z);
    REAL8 xa = .5 * (S1z - S2z);
//    fprintf(stderr, "TaylorF2 Waveform Generator; EOS: %d\n", eos);
    lam1 = XLALSimInspiralEOSLambda(eos, m1);
    lam2 = XLALSimInspiralEOSLambda(eos, m2);
    lam1 = lam1/(m1*LAL_MTSUN_SI*m1*LAL_MTSUN_SI*m1*LAL_MTSUN_SI*m1*LAL_MTSUN_SI*m1*LAL_MTSUN_SI);
    lam2 = lam2/(m2*LAL_MTSUN_SI*m2*LAL_MTSUN_SI*m2*LAL_MTSUN_SI*m2*LAL_MTSUN_SI*m2*LAL_MTSUN_SI);
    fCONT = XLALSimInspiralContactFrequency(m1, lam1, m2, lam2);
    
//	fprintf(stderr, "Lambdas: %e, %e\n", lam1, lam2);

    REAL8 qm_def1 = XLALSimInspiralEOSQfromLambda(lam1); /* The QM deformability parameters */
    REAL8 qm_def2 = XLALSimInspiralEOSQfromLambda(lam2); /* This is 1 for black holes and larger for neutron stars */

//    fprintf(stderr,"Quadparams: %e, %e\n", qm_def1, qm_def2);

    switch( spinO )
    {
        case LAL_SIM_INSPIRAL_SPIN_ORDER_ALL:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_35PN:
            psiSO35 = (chi1 * (-8980424995./762048. + 6586595.*eta/756.
                    - 305.*eta*eta/36.) + d * (170978035./48384.
                    - 2876425.*eta/672. - 4735.*eta*eta/144.) ) * chi1 * S1z
                    + (chi2 * (-8980424995./762048. + 6586595.*eta/756.
                    - 305.*eta*eta/36.) - d * (170978035./48384.
                    - 2876425.*eta/672. - 4735.*eta*eta/144.) ) * chi2 * S2z;
        case LAL_SIM_INSPIRAL_SPIN_ORDER_3PN:
            psiSO3 = LAL_PI * ( (260.*chi1 + 1490./3.) * chi1 * S1z
                    + (260.*chi2 + 1490./3.) * chi2 * S2z);
        case LAL_SIM_INSPIRAL_SPIN_ORDER_25PN:
            /* Compute 2.5PN SO correction */
            // See Eq. (6.25) in arXiv:0810.5336
            pn_gamma = (732985.L/2268.L - 24260.L/81.L * eta - 340.L/9.L * eta * eta ) * xs;
            pn_gamma += (732985.L/2268.L +140.L/9.0L * eta) * xa * d;
        case LAL_SIM_INSPIRAL_SPIN_ORDER_2PN:
            /* Compute 2.0PN SS, QM, and self-spin */
            // See Eq. (6.24) in arXiv:0810.5336
            // 9b,c,d in arXiv:astro-ph/0504538
            pn_sigma = eta * (721.L/48.L *S1z*S2z-247.L/48.L*S1z*S2z);
            pn_sigma += (720*qm_def1 - 1)/96.0 * (chi1*chi1*S1z*S1z);
            pn_sigma += (720*qm_def2 - 1)/96.0 * (chi2*chi2*S2z*S2z);
            pn_sigma -= (240*qm_def1 - 7)/96.0 * (chi1*chi1*S1z*S1z);
            pn_sigma -= (240*qm_def2 - 7)/96.0 * (chi2*chi2*S2z*S2z);
        case LAL_SIM_INSPIRAL_SPIN_ORDER_15PN:
            /* Compute 1.5PN SO correction */
            // Eq. (6.23) in arXiv:0810.5336
            pn_beta = (113.L/12.L- 19.L/3.L * eta) * xs + 113.L/12.L * d * xa;
        case LAL_SIM_INSPIRAL_SPIN_ORDER_1PN:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_05PN:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_0PN:
            break;
        default:
            XLALPrintError("XLAL Error - %s: Invalid spin PN order %s\n",
                    __func__, spinO );
            XLAL_ERROR(XLAL_EINVAL);
            break;
    }

    /* Tidal coefficients for phasing, fluz, and energy */
    REAL8 pft10 = 0.;
    REAL8 pft12 = 0.;
    REAL8 pft13 = 0.;
    REAL8 pft14 = 0.;
    REAL8 pft15 = 0.;
    switch( tideO )
    {
        case LAL_SIM_INSPIRAL_TIDAL_ORDER_ALL:
        case LAL_SIM_INSPIRAL_TIDAL_ORDER_75PN:
            pft15 = lam2* chi2*chi2*chi2*chi2 * 1.L/28.L*LAL_PI*(27719.L - 22127.L*chi2 + 7022.L*chi2*chi2 - 10232.L*chi2*chi2*chi2)
                + lam1* chi1*chi1*chi1*chi1 * 1.L/28.L*LAL_PI*(27719.L - 22127.L*chi1 + 7022.L*chi1*chi1 - 10232.L*chi1*chi1*chi1);
        case LAL_SIM_INSPIRAL_TIDAL_ORDER_7PN:
            pft14 = - lam2 * chi2*chi2*chi2*chi2 * 24.L*(39927845.L/508032.L                           - 480043345.L/9144576.L*chi2 + 9860575.L/127008.L*chi2*chi2 - 421821905.L/2286144.L*chi2*chi2*chi2 + 4359700.L/35721.L*chi2*chi2*chi2*chi2 - 10578445.L/285768.L*chi2*chi2*chi2*chi2*chi2)
                   - lam1 * chi1*chi1*chi1*chi1 * 24.L*(39927845.L/508032.L 
                - 480043345.L/9144576.L*chi1 + 9860575.L/127008.L*chi1*chi1 - 421821905.L/2286144.L*chi1*chi1*chi1 + 4359700.L/35721.L*chi1*chi1*chi1*chi1 - 10578445.L/285768.L*chi1*chi1*chi1*chi1*chi1);
        case LAL_SIM_INSPIRAL_TIDAL_ORDER_65PN:
            pft13 = - lam2 * chi2*chi2*chi2*chi2 * 24.L*(12.L - 11.L*chi2)*LAL_PI
                    - lam1 * chi1*chi1*chi1*chi1 * 24.L*(12.L - 11.L*chi1)*LAL_PI;
        case LAL_SIM_INSPIRAL_TIDAL_ORDER_6PN:
            pft12 = - 5.L * lam2 * chi2*chi2*chi2*chi2 * (3179.L - 919.L*chi2
                    - 2286.L*chi2*chi2 + 260.L*chi2*chi2*chi2)/28.L
                    - 5.L * lam1 * chi1*chi1*chi1*chi1 * (3179.L - 919.L*chi1
                    - 2286.L*chi1*chi1 + 260.L*chi1*chi1*chi1)/28.L;
        case LAL_SIM_INSPIRAL_TIDAL_ORDER_5PN:
            pft10 = - 24.L * lam2 * chi2*chi2*chi2*chi2 * (1.L + 11.L*chi1)
                    - 24.L * lam1 * chi1*chi1*chi1*chi1 * (1.L + 11.L*chi2);
        case LAL_SIM_INSPIRAL_TIDAL_ORDER_0PN:
            break;
        default:
            XLALPrintError("XLAL Error - %s: Invalid tidal PN order %s\n",
                    __func__, tideO );
            XLAL_ERROR(XLAL_EINVAL);
            break;
    }

//    fprintf(stderr, "6PN Tidal: %e\n", pft12);
//    fprintf(stderr, "5PN Tidal: %e\n", pft10);

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
    if (m1_SI <= 0) XLAL_ERROR(XLAL_EDOM);
    if (m2_SI <= 0) XLAL_ERROR(XLAL_EDOM);
    if (fStart <= 0) XLAL_ERROR(XLAL_EDOM);
    if (r <= 0) XLAL_ERROR(XLAL_EDOM);

    /* allocate htilde */
    if ( fEnd == 0. ) // Use min{fISCO,fCONT} as termination frequency
        f_max = fISCO < fCONT ? fISCO : fCONT;
    else // End at user-specified freq.
        f_max = fEnd;
    // printf("fISCO = %e\n", fISCO);
    // printf("fCONT = %e\n", fCONT);
    // printf("f_max = %e\n", f_max);    
    n = (size_t) (f_max / deltaF + 1);
    XLALGPSAdd(&tC, -1 / deltaF);  /* coalesce at t=0 */
    htilde = XLALCreateCOMPLEX16FrequencySeries("htilde: FD waveform", &tC, 0.0, deltaF, &lalStrainUnit, n);
    if (!htilde) XLAL_ERROR(XLAL_EFUNC);
    memset(htilde->data->data, 0, n * sizeof(COMPLEX16));
    XLALUnitDivide(&htilde->sampleUnits, &htilde->sampleUnits, &lalSecondUnit);

    /* extrinsic parameters */
    amp0 = -4. * m1 * m2 / r * LAL_MRSUN_SI * LAL_MTSUN_SI * sqrt(LAL_PI/12.L);
    shft = LAL_TWOPI * (tC.gpsSeconds + 1e-9 * tC.gpsNanoSeconds);

    const REAL8 log4=log(4.0);
    const REAL8 logv0=log(v0);
    
    /* Fill with non-zero vals from fStart to f_max */
    iStart = (size_t) ceil(fStart / deltaF);
    data = htilde->data->data;

    FILE *phaseTF2;

    phaseTF2 = fopen("phaseTF2.out","w");

    for (i = iStart; i < n; i++) {
        const REAL8 f = i * deltaF;
        const REAL8 v = cbrt(piM*f);
        const REAL8 logv = log(v);
        const REAL8 v2 = v * v;
        const REAL8 v3 = v * v2;
        const REAL8 v4 = v * v3;
        const REAL8 v5 = v * v4;
        const REAL8 v6 = v * v5;
        const REAL8 v7 = v * v6;
        const REAL8 v8 = v * v7;
        const REAL8 v9 = v * v8;
        const REAL8 v10 = v * v9;
        const REAL8 v12 = v2 * v10;
        const REAL8 v13 = v * v12;
        const REAL8 v14 = v * v13;
        const REAL8 v15 = v * v14;
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
                phasing += (pfa6 + pfl6 * (log4+logv)) * v6;
            case 5:
                phasing += (pfa5 + pfl5 * (logv-logv0)) * v5;
            case 4:
                phasing += pfa4 * v4;
            case 3:
                phasing += pfa3 * v3;
            case 2:
                phasing += pfa2 * v2;
            case 0:
                phasing += 1.;
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
                flux += (FTa6 + FTl6*2.0*(log4+logv)) * v6;
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

        switch( spinO )
        {
            case LAL_SIM_INSPIRAL_SPIN_ORDER_ALL:
            case LAL_SIM_INSPIRAL_SPIN_ORDER_35PN:
                phasing += psiSO35 * v7;
            case LAL_SIM_INSPIRAL_SPIN_ORDER_3PN:
                phasing += psiSO3 * v6;
            case LAL_SIM_INSPIRAL_SPIN_ORDER_25PN:
                phasing += -pn_gamma * (1 + 3*(logv-logv0)) * v5;
            case LAL_SIM_INSPIRAL_SPIN_ORDER_2PN:
                phasing += -10.L*pn_sigma * v4;
            case LAL_SIM_INSPIRAL_SPIN_ORDER_15PN:
                phasing += 4.L*pn_beta * v3;
            case LAL_SIM_INSPIRAL_SPIN_ORDER_1PN:
            case LAL_SIM_INSPIRAL_SPIN_ORDER_05PN:
            case LAL_SIM_INSPIRAL_SPIN_ORDER_0PN:
                break;
            default:
                XLALPrintError("XLAL Error - %s: Invalid spin PN order %s\n",
                        __func__, spinO );
                XLAL_ERROR(XLAL_EINVAL);
                break;
        }
//	fprintf(stderr,"Phase before tidal: %e\n", phasing);
        switch( tideO )
        {
            case LAL_SIM_INSPIRAL_TIDAL_ORDER_ALL:
            case LAL_SIM_INSPIRAL_TIDAL_ORDER_75PN:
                phasing += pft15 * v15;
            case LAL_SIM_INSPIRAL_TIDAL_ORDER_7PN:
                phasing += pft14 * v14;
            case LAL_SIM_INSPIRAL_TIDAL_ORDER_65PN:
                phasing += pft13 * v13;
            case LAL_SIM_INSPIRAL_TIDAL_ORDER_6PN:
                phasing += pft12 * v12;
//              fprintf(stderr,"Added v12, phasing = %e\n",phasing);
            case LAL_SIM_INSPIRAL_TIDAL_ORDER_5PN:
                phasing += pft10 * v10;
//              fprintf(stderr,"Added v10, phasing = %e\n",phasing);
            case LAL_SIM_INSPIRAL_TIDAL_ORDER_0PN:
                break;
            default:
                XLALPrintError("XLAL Error - %s: Invalid tidal PN order %s\n",
                        __func__, tideO );
                XLAL_ERROR(XLAL_EINVAL);
                break;
        }
//	fprintf(stderr,"Phasing after tidal: %e\n", phasing);
        phasing *= pfaN / v5;
        flux *= FTaN * v10;
        dEnergy *= dETaN * v;
        // Note the factor of 2 b/c phic is orbital phase
        phasing += shft * f - 2.*phic;
        amp = amp0 * sqrt(-dEnergy/flux) * v;
        data[i] = amp * cos(phasing - LAL_PI_4)
                - amp * sin(phasing - LAL_PI_4) * 1.0j;

        fprintf(phaseTF2, "%e %e\n", f, phasing - LAL_PI_4);

    }
       
    fclose(phaseTF2);

    *htilde_out = htilde;
    return XLAL_SUCCESS;
}
