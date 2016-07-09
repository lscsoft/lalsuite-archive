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
#include <lal/Sequence.h>
#include <lal/LALDatatypes.h>
#include <lal/LALSimInspiral.h>
#include <lal/Units.h>
#include <lal/XLALError.h>
#include "LALSimInspiralPNCoefficients.c"

#ifndef _OPENMP
#define omp ignore
#endif

/**
 * @addtogroup LALSimInspiralTaylorXX_c
 * @{
 *
 * @review TaylorF2 routines reviewed by Frank Ohme, Andrew Lundgren, Alex Nitz,
 * Alex Nielsen, Salvatore Vitale, Jocelyn Read, Sebastian Khan.
 * The review concluded with git hash 6106138b2140ffb11bc38fc914e0a1de7082dc4d (Nov 2014)
 *
 * @name Routines for TaylorF2 Waveforms
 * @sa
 * Section IIIF of Alessandra Buonanno, Bala R Iyer, Evan
 * Ochsner, Yi Pan, and B S Sathyaprakash, "Comparison of post-Newtonian
 * templates for compact binary inspiral signals in gravitational-wave
 * detectors", Phys. Rev. D 80, 084043 (2009), arXiv:0907.0700v1
 *
 * @{
 */

/** \brief Returns structure containing TaylorF2 phasing coefficients for given
 *  physical parameters.
 */
int XLALSimInspiralTaylorF2AlignedPhasing(
        PNPhasingSeries **pn,   /**< phasing coefficients (output) */
        const REAL8 m1,         /**< mass of body 1 */
        const REAL8 m2,		/**< mass of body 2 */
        const REAL8 chi1,	/**< aligned spin parameter of body 1 */
        const REAL8 chi2,	/**< aligned spin parameter of body 2 */
        const REAL8 qm_def1,	/**< quadrupole-monopole parameter of body 1 (set 1 for BH) */
        const REAL8 qm_def2,	/**< quadrupole-monopole parameter of body 2 (set 1 for BH) */
        const LALSimInspiralSpinOrder spinO,  /**< PN order for spin contributions */
        const LALSimInspiralTestGRParam *p /**< Linked list containing the extra testing GR parameters >**/
	)
{
    PNPhasingSeries *pfa;

    if (!pn) XLAL_ERROR(XLAL_EFAULT);
    if (*pn) XLAL_ERROR(XLAL_EFAULT);


    pfa = (PNPhasingSeries *) LALMalloc(sizeof(PNPhasingSeries));

    XLALSimInspiralPNPhasing_F2(pfa, m1, m2, chi1, chi2, chi1*chi1, chi2*chi2, chi1*chi2, qm_def1, qm_def2, spinO, p);

    *pn = pfa;

    return XLAL_SUCCESS;
}

int XLALSimInspiralTaylorF2Core(
        COMPLEX16FrequencySeries **htilde_out, /**< FD waveform */
	const REAL8Sequence *freqs,            /**< frequency points at which to evaluate the waveform (Hz) */
        const REAL8 phi_ref,                   /**< reference orbital phase (rad) */
        const REAL8 m1_SI,                     /**< mass of companion 1 (kg) */
        const REAL8 m2_SI,                     /**< mass of companion 2 (kg) */
        const REAL8 S1z,                       /**<  z component of the spin of companion 1 */
        const REAL8 S2z,                       /**<  z component of the spin of companion 2  */
        const REAL8 f_ref,                     /**< Reference GW frequency (Hz) - if 0 reference point is coalescence */
	const REAL8 shft,		       /**< time shift to be applied to frequency-domain phase (sec)*/
        const REAL8 r,                         /**< distance of source (m) */
        const REAL8 quadparam1,                /**< quadrupole deformation parameter of body 1 (dimensionless, 1 for BH) */
        const REAL8 quadparam2,                /**< quadrupole deformation parameter of body 2 (dimensionless, 1 for BH) */
        const REAL8 lambda1,                   /**< (tidal deformation of body 1)/(mass of body 1)^5 */
        const REAL8 lambda2,                   /**< (tidal deformation of body 2)/(mass of body 2)^5 */
        const LALSimInspiralSpinOrder spinO,  /**< twice PN order of spin effects */
        const LALSimInspiralTidalOrder tideO,  /**< flag to control tidal effects */
        const INT4 phaseO,                     /**< twice PN phase order */
        const INT4 amplitudeO,                  /**< twice PN amplitude order */
        const LALSimInspiralTestGRParam *p /**< Linked list containing the extra testing GR parameters >**/
        )
{

    if (!htilde_out) XLAL_ERROR(XLAL_EFAULT);
    if (!freqs) XLAL_ERROR(XLAL_EFAULT);
    /* external: SI; internal: solar masses */
    //static int calls_debug = 0;
    const REAL8 m1 = m1_SI / LAL_MSUN_SI;
    const REAL8 m2 = m2_SI / LAL_MSUN_SI;
    const REAL8 m = m1 + m2;
    const REAL8 m_sec = m * LAL_MTSUN_SI;  /* total mass in seconds */
    const REAL8 eta = m1 * m2 / (m * m);
    const REAL8 piM = LAL_PI * m_sec;
    const REAL8 m1OverM = m1 / m;
    const REAL8 m2OverM = m2 / m;
    REAL8 amp0;
    size_t i;
    COMPLEX16 *data = NULL;
    LIGOTimeGPS tC = {0, 0};
    INT4 iStart = 0;

    COMPLEX16FrequencySeries *htilde = NULL;

    if (*htilde_out) { //case when htilde_out has been allocated in XLALSimInspiralTaylorF2
	    htilde = *htilde_out;
	    iStart = htilde->data->length - freqs->length; //index shift to fill pre-allocated data
	    if(iStart < 0) XLAL_ERROR(XLAL_EFAULT);
    }
    else { //otherwise allocate memory here
	    htilde = XLALCreateCOMPLEX16FrequencySeries("htilde: FD waveform", &tC, freqs->data[0], 0., &lalStrainUnit, freqs->length);
	    if (!htilde) XLAL_ERROR(XLAL_EFUNC);
	    XLALUnitMultiply(&htilde->sampleUnits, &htilde->sampleUnits, &lalSecondUnit);
    }

    /* phasing coefficients */
    PNPhasingSeries pfa;
    XLALSimInspiralPNPhasing_F2(&pfa, m1, m2, S1z, S2z, S1z*S1z, S2z*S2z, S1z*S2z, quadparam1, quadparam2, spinO, p);

    REAL8 pfaN = 0.; REAL8 pfa1 = 0.;
    REAL8 pfa2 = 0.; REAL8 pfa3 = 0.; REAL8 pfa4 = 0.;
    REAL8 pfa5 = 0.; REAL8 pfl5 = 0.;
    REAL8 pfa6 = 0.; REAL8 pfl6 = 0.;
    REAL8 pfa7 = 0.;

    //fprintf(stdout, "=========== DEBUG by Jeongcho Kim & Chunglee KIm ===========\n");
    //fprintf(stdout, "==Function : LALSimInspiralTaylorF2.c:XLALSimInspiralTaylorF2()\n");
    //fprintf(stdout, "Sequence of calls = %d\n", ++calls_debug);
    //fprintf(stdout, "(S1z, S2z) = (%f, %f)\n", S1z, S2z);
    //fprintf(stdout, "=========================================================\n");
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
            pfa1 = pfa.v[1];
        case 0:
            pfaN = pfa.v[0];
            break;
        default:
            XLAL_ERROR(XLAL_ETYPE, "Invalid phase PN order %d", phaseO);
    }

    /* Validate expansion order arguments.
     * This must be done here instead of in the OpenMP parallel loop
     * because when OpenMP parallelization is turned on, early exits
     * from loops (via return or break statements) are not permitted.
     */

    /* Validate amplitude PN order. */
    switch (amplitudeO)
    {
        case -1:
        case 7:
        case 6:
        case 5:
        case 4:
        case 3:
        case 2:
        case 0:
            break;
        default:
            XLAL_ERROR(XLAL_ETYPE, "Invalid amplitude PN order %d", amplitudeO);
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
	    pft12 = pfaN * (lambda1*XLALSimInspiralTaylorF2Phasing_12PNTidalCoeff(m1OverM) + lambda2*XLALSimInspiralTaylorF2Phasing_12PNTidalCoeff(m2OverM) );
        case LAL_SIM_INSPIRAL_TIDAL_ORDER_5PN:
            pft10 = pfaN * ( lambda1*XLALSimInspiralTaylorF2Phasing_10PNTidalCoeff(m1OverM) + lambda2*XLALSimInspiralTaylorF2Phasing_10PNTidalCoeff(m2OverM) );
        case LAL_SIM_INSPIRAL_TIDAL_ORDER_0PN:
            break;
        default:
            XLAL_ERROR(XLAL_EINVAL, "Invalid tidal PN order %d", tideO);
    }

    /* The flux and energy coefficients below are used to compute SPA amplitude corrections */

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


    /* Perform some initial checks */
    if (m1_SI <= 0) XLAL_ERROR(XLAL_EDOM);
    if (m2_SI <= 0) XLAL_ERROR(XLAL_EDOM);
    if (f_ref < 0) XLAL_ERROR(XLAL_EDOM);
    if (r <= 0) XLAL_ERROR(XLAL_EDOM);

    /* extrinsic parameters */
    amp0 = -4. * m1 * m2 / r * LAL_MRSUN_SI * LAL_MTSUN_SI * sqrt(LAL_PI/12.L);

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
        const REAL8 v9ref = vref * v8ref;
        const REAL8 v10ref = vref * v9ref;
        const REAL8 v12ref = v2ref * v10ref;
        ref_phasing += pfa7 * v7ref;
        ref_phasing += (pfa6 + pfl6 * logvref) * v6ref;
        ref_phasing += (pfa5 + pfl5 * logvref) * v5ref;
        ref_phasing += pfa4 * v4ref;
        ref_phasing += pfa3 * v3ref;
        ref_phasing += pfa2 * v2ref;
        ref_phasing += pfa1 * vref;
        ref_phasing += pfaN;

        /* Tidal terms in reference phasing */
        ref_phasing += pft12 * v12ref;
        ref_phasing += pft10 * v10ref;

        ref_phasing /= v5ref;
    } /* End of if(f_ref != 0) block */

    #pragma omp parallel for
    for (i = 0; i < freqs->length; i++) {
        const REAL8 f = freqs->data[i];
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
        phasing += pfa1 * v;
        phasing += pfaN;

        /* Tidal terms in phasing */
        phasing += pft12 * v12;
        phasing += pft10 * v10;

    /* WARNING! Amplitude orders beyond 0 have NOT been reviewed!
     * Use at your own risk. The default is to turn them off.
     * These do not currently include spin corrections.
     * Note that these are not higher PN corrections to the amplitude.
     * They are the corrections to the leading-order amplitude arising
     * from the stationary phase approximation. See for instance
     * Eq 6.9 of arXiv:0810.5336
     */
	switch (amplitudeO)
        {
            case 7:
                flux += FTa7 * v7;
            case 6:
                flux += (FTa6 + FTl6*logv) * v6;
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
            case -1: /* Default to no SPA amplitude corrections */
            case 0:
                flux += 1.;
                dEnergy += 1.;
        }

        phasing /= v5;
        flux *= FTaN * v10;
        dEnergy *= dETaN * v;
        // Note the factor of 2 b/c phi_ref is orbital phase
        phasing += shft * f - 2.*phi_ref - ref_phasing;
        amp = amp0 * sqrt(-dEnergy/flux) * v;
        data[i+iStart] = amp * cos(phasing - LAL_PI_4)
                - amp * sin(phasing - LAL_PI_4) * 1.0j;
    }

    *htilde_out = htilde;
    return XLAL_SUCCESS;
}

/**
 * Computes the stationary phase approximation to the Fourier transform of
 * a chirp waveform. The amplitude is given by expanding \f$1/\sqrt{\dot{F}}\f$.
 * If the PN order is set to -1, then the highest implemented order is used.
 *
 * @note f_ref is the GW frequency at which phi_ref is defined. The most common
 * choice in the literature is to choose the reference point as "coalescence",
 * when the frequency becomes infinite. This is the behavior of the code when
 * f_ref==0. If f_ref > 0, phi_ref sets the orbital phase at that GW frequency.
 *
 * See arXiv:0810.5336 and arXiv:astro-ph/0504538 for spin corrections
 * to the phasing.
 * See arXiv:1303.7412 for spin-orbit phasing corrections at 3 and 3.5PN order
 *
 * The spin and tidal order enums are defined in LALSimInspiralWaveformFlags.h
 */
int XLALSimInspiralTaylorF2(
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
        const INT4 amplitudeO,                  /**< twice PN amplitude order */
        const LALSimInspiralTestGRParam *p /**< Linked list containing the extra testing GR parameters >**/
        )
{
    /* external: SI; internal: solar masses */
    const REAL8 m1 = m1_SI / LAL_MSUN_SI;
    const REAL8 m2 = m2_SI / LAL_MSUN_SI;
    const REAL8 m = m1 + m2;
    const REAL8 m_sec = m * LAL_MTSUN_SI;  /* total mass in seconds */
    // const REAL8 eta = m1 * m2 / (m * m);
    const REAL8 piM = LAL_PI * m_sec;
    const REAL8 vISCO = 1. / sqrt(6.);
    const REAL8 fISCO = vISCO * vISCO * vISCO / piM;
    //const REAL8 m1OverM = m1 / m;
    // const REAL8 m2OverM = m2 / m;
    REAL8 shft, f_max;
    size_t i, n;
    INT4 iStart;
    REAL8Sequence *freqs = NULL;
    LIGOTimeGPS tC = {0, 0};
    int ret;

    COMPLEX16FrequencySeries *htilde = NULL;

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
    if (f_max <= fStart) XLAL_ERROR(XLAL_EDOM);

    n = (size_t) (f_max / deltaF + 1);
    XLALGPSAdd(&tC, -1 / deltaF);  /* coalesce at t=0 */
    htilde = XLALCreateCOMPLEX16FrequencySeries("htilde: FD waveform", &tC, 0.0, deltaF, &lalStrainUnit, n);
    if (!htilde) XLAL_ERROR(XLAL_EFUNC);
    memset(htilde->data->data, 0, n * sizeof(COMPLEX16));
    XLALUnitMultiply(&htilde->sampleUnits, &htilde->sampleUnits, &lalSecondUnit);

    /* Fill with non-zero vals from fStart to f_max */
    iStart = (INT4) ceil(fStart / deltaF);

    /* Sequence of frequencies where waveform model is to be evaluated */
    freqs = XLALCreateREAL8Sequence(n - iStart);

    /* extrinsic parameters */
    shft = LAL_TWOPI * (tC.gpsSeconds + 1e-9 * tC.gpsNanoSeconds);

    #pragma omp parallel for
    for (i = iStart; i < n; i++) {
        freqs->data[i-iStart] = i * deltaF;
    }
    ret = XLALSimInspiralTaylorF2Core(&htilde, freqs, phi_ref, m1_SI, m2_SI,
                                      S1z, S2z, f_ref, shft, r, quadparam1, quadparam2,
                                      lambda1, lambda2, spinO, tideO, phaseO, amplitudeO, p);

    XLALDestroyREAL8Sequence(freqs);

    *htilde_out = htilde;

    return ret;
}

/*
 *  Copyright (C) 2014 Jeongcho Kim, Hyung Won Lee, Chunglee Kim
 *  All equations and Appendix D cited here are in Evan Ochsner's thesis.
 *  Adapting features in:
 *    - LALSimInspiral.c
 *    - LALSimInspiral.h
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

#define MAX_HARMONICS    7 // maximum hamoincs to integrate
#define MAX_PHASE_ORDER    7 // maximum PN order for phase
#define MAX_AMP_ORDER    5 // maximum PN order for amplitude

/**
 * Find the least nonnegative integer power of 2 that is
 * greater than or equal to n.  Inspired by similar routine in gstlal.
 */
/*static size_t CeilPow2(double n) {
    double signif;
    int exponent;
    signif = frexp(n, &exponent);
    if (signif < 0)
        return 1;
    if (signif == 0.5)
        exponent -= 1;
    return ((size_t) 1) << exponent;
}*/

/* Calculate the spin corrections for TaylorF2 
        reference -> <http://arxiv.org/pdf/0810.5336v3.pdf>
*/

typedef struct tag_sf2_spin_corr_amp {
    REAL8 beta;
    REAL8 sigma;
    REAL8 gamma;
}sf2_spin_corr_amp;
typedef struct tag_sf2_amp_corr_param {
    REAL8 costh;      // cos(theta), cosine of inclination of the source
    REAL8 sinth;      // sin(theta), sine of inclination of the source
    REAL8 f_cross;    // F_cross
    REAL8 f_plus;     // F_plus
    REAL8 m;          // m1+m1 total mass in solar mass
    REAL8 eta;        // eta = m1*m1/m^2
    REAL8 delta;      // delta = (m1-m1)/m
    REAL8 s1x;        // x component of spin1 
    REAL8 s1y;        // y component of spin1 
    REAL8 s1z;        // z component of spin1 
    REAL8 s2x;        // x component of spin2 
    REAL8 s2y;        // y component of spin2 
    REAL8 s2z;        // z component of spin2 
    REAL8 LNhatx;     // x component of LNhat 
    REAL8 LNhaty;     // y component of LNhat 
    REAL8 LNhatz;     // z component of LNhat 
    REAL8 Xs_x;
    REAL8 Xs_y;
    REAL8 Xs_z;
    REAL8 Xa_x;
    REAL8 Xa_y;
    REAL8 Xa_z;
    REAL8 Xs2;        // Xs*Xs (=chi_s dot chi_s)
    REAL8 Xa2;        // Xa*Xa (=chi_a dot chi_a)
    REAL8 XsXa;       // Xs dot Xa (=xchi_s dot chi_a)
    REAL8 XsLn;       // chi_s dot Ln_hat
    REAL8 XaLn;       // chi_a dot Ln_hat
    REAL8 costh2;
    REAL8 costh3;
    REAL8 costh4;
    REAL8 costh5;
    REAL8 costh6;
    REAL8 sinth2;
    REAL8 sinth3;
    REAL8 sinth4;
    REAL8 sinth5;
} sf2_amp_corr_param;

// add new parameters for spin of the campanion 2 to calculate 
// Xs(=chi_s) and Xa(=chi_a), also it need to calculate XsLn(=chi_s dot Ln_hat)
// XaLn(=chi_a dot Ln_hat)
sf2_spin_corr_amp sf2_spin_corrections_amp(
        const REAL8 m1,            /** mass of companion 1 (solar masses) */
        const REAL8 m2,            /** mass of companion 2 (solar masses) */
        const REAL8 S1x,           /** x component of the spin of companion 1 */
        const REAL8 S1y,           /** y component of the spin of companion 1 */
        const REAL8 S1z,           /** z component of the spin of companion 1 */
        const REAL8 S2x,           /** x component of the spin of companion 2 */
        const REAL8 S2y,           /** y component of the spin of companion 2 */
        const REAL8 S2z,           /** z component of the spin of companion 2 */
        const REAL8 lnhatx,        /** x component of the vector Ln hat */
        const REAL8 lnhaty,        /** y component of the vector Ln hat */
        const REAL8 lnhatz    ) ;  /** z component of the vector Ln hat */

sf2_spin_corr_amp sf2_spin_corrections_amp(
        const REAL8 m1,            /** mass of companion 1 (solar masses) */
        const REAL8 m2,            /** mass of companion 2 (solar masses) */ 
        const REAL8 S1x,           /** x component of the spin of companion 1 */
        const REAL8 S1y,           /** y component of the spin of companion 1 */
        const REAL8 S1z,           /** z component of the spin of companion 1 */
        const REAL8 S2x,           /** x component of the spin of companion 2 */
        const REAL8 S2y,           /** y component of the spin of companion 2 */ 
        const REAL8 S2z,           /** z component of the spin of companion 2 */
        const REAL8 lnhatx,        /** x component of the vector Ln hat */
        const REAL8 lnhaty,        /** y component of the vector Ln hat */
        const REAL8 lnhatz         /** z component of the vector Ln hat */
        )                          /** z component of the spin of companion 2  */
{
    sf2_spin_corr_amp spin_corrections;
    REAL8 M = m1 + m2;
    REAL8 eta = m1 * m2 / (M * M);
    REAL8 delta = (m1 - m2) / (m1 + m2);
    REAL8 Xs_x = 0.5*(S1x+S2x);
    REAL8 Xs_y = 0.5*(S1y+S2y);
    REAL8 Xs_z = 0.5*(S1z+S2z);
    REAL8 Xa_x = 0.5*(S1x-S2x);
    REAL8 Xa_y = 0.5*(S1y-S2y);
    REAL8 Xa_z = 0.5*(S1z-S2z);
    REAL8 Xs2 = Xs_x*Xs_x+Xs_y*Xs_y+Xs_z*Xs_z;  // Xs*Xs (=chi_s dot chi_s)
    REAL8 Xa2 = Xa_x*Xa_x+Xa_y*Xa_y+Xa_z*Xa_z;  // Xa*Xa (=chi_a dot chi_a)
    REAL8 XsXa = Xs_x*Xa_x + Xs_y*Xa_y + Xs_z*Xa_z;  // Xs dot Xa (=xchi_s dot chi_a)
    REAL8 XsLn = (Xs_x*lnhatx)+(Xs_y*lnhaty)+(Xs_z*lnhatz); // chi_s dot Ln_hat
    REAL8 XaLn = (Xa_x*lnhatx)+(Xa_y*lnhaty)+(Xa_z*lnhatz); // chi_a dot Ln_hat
    

    // correcte beta, gamma, sigma based on the Evans thesis
    REAL8 sf2_beta = (113.L/12.L- 19.L/3.L * eta) * (XsLn) + 113.L/12.L * delta * (XaLn);
    
    REAL8 sf2_sigma = eta * (721.L/48.L * ((XsLn)*(XsLn) - (XaLn)*(XaLn))-247.L/48.L*(Xs2 - Xa2));
    sf2_sigma += (1-2*eta)* (719/96.0 * ((XsLn)*(XsLn) + (XaLn)*(XaLn)) - 233.L/96.L * (Xs2 + Xa2));
    sf2_sigma += delta * (719/48.0 *(XsLn)*(XaLn) - 233.L/48.L * XsXa);
    
    // should be check the sign of gamma with respect to the Evans thesis Eq. (4.84)
    REAL8 sf2_gamma = (732985.L/2268.L - 24260.L/81.L * eta - 340.L/9.L * (eta*eta) ) * (XsLn);
    sf2_gamma += (732985.L/2268.L +140.L/9.0L * eta) * delta * XaLn;
    

    spin_corrections.beta = sf2_beta;
    spin_corrections.sigma = sf2_sigma;
    spin_corrections.gamma = -sf2_gamma;
    return spin_corrections;
}

/**
 calculate coefficients for each SPA phase PN order
 assume the maximum twice PN order to be 9, array size is 10
*/
int sf2_psi_SPA_coeffs_PN_order(
    REAL8 *PN_coeffs, /** coeffs for each PN order*/
    const REAL8 m1, /**< Mass of body 1, in Msol */
    const REAL8 m2, /**< Mass of body 2, in Msol */
    const REAL8 chi1L, /**< Component of dimensionless spin 1 along Lhat */
    const REAL8 chi2L, /**< Component of dimensionless spin 2 along Lhat */
    const REAL8 chi1sq, /**< Magnitude of dimensionless spin 1 */
    const REAL8 chi2sq, /**< Magnitude of dimensionless spin 2 */
    const REAL8 chi1dotchi2, /**< Dot product of dimensionles spin 1 and spin 2 */
    const REAL8 qm_def1, /**< Quadrupole deformation parameter of body 1 (dimensionless) */
    const REAL8 qm_def2, /**< Quadrupole deformation parameter of body 2 (dimensionless) */
    const REAL8 lambda1,                   /**< (tidal deformation of body 1)/(mass of body 1)^5 */
    const REAL8 lambda2,                   /**< (tidal deformation of body 2)/(mass of body 2)^5 */
    const LALSimInspiralTidalOrder tideO,  /**< flag to control tidal effects */
    const LALSimInspiralSpinOrder spinO /**< Enums specifying spin order are in LALSimInspiralWaveformFlags.h */
    );
/**
 calculate coefficients for each SPA phase PN order
 assume the maximum twice PN order to be PN_PHASING_SERIES_MAX_ORDER
 PN coefficients are different, since ammplitude correction requires 
 higher harmonics not only second harmonic.
*/
int sf2_psi_SPA_coeffs_PN_order(
    REAL8 *PN_coeffs, /** coeffs for each PN order*/
    const REAL8 m1, /**< Mass of body 1, in Msol */
    const REAL8 m2, /**< Mass of body 2, in Msol */
    const REAL8 chi1L, /**< Component of dimensionless spin 1 along Lhat */
    const REAL8 chi2L, /**< Component of dimensionless spin 2 along Lhat */
    const REAL8 chi1sq, /**< Magnitude of dimensionless spin 1 */
    const REAL8 chi2sq, /**< Magnitude of dimensionless spin 2 */
    const REAL8 chi1dotchi2, /**< Dot product of dimensionles spin 1 and spin 2 */
    const REAL8 qm_def1, /**< Quadrupole deformation parameter of body 1 (dimensionless) */
    const REAL8 qm_def2, /**< Quadrupole deformation parameter of body 2 (dimensionless) */
    const REAL8 lambda1,                   /**< (tidal deformation of body 1)/(mass of body 1)^5 */
    const REAL8 lambda2,                   /**< (tidal deformation of body 2)/(mass of body 2)^5 */
    const LALSimInspiralTidalOrder tideO,  /**< flag to control tidal effects */
    const LALSimInspiralSpinOrder spinO /**< Enums specifying spin order are in LALSimInspiralWaveformFlags.h */
    )
{
    PNPhasingSeries pfa;
    memset(&pfa, 0x00, sizeof(pfa));
    XLALSimInspiralPNPhasing_F2(
        &pfa, /**< \todo UNDOCUMENTED */
        m1, /**< Mass of body 1, in Msol */
        m2, /**< Mass of body 2, in Msol */
        chi1L, /**< Component of dimensionless spin 1 along Lhat */
        chi2L, /**< Component of dimensionless spin 2 along Lhat */
        chi1sq, /**< Magnitude of dimensionless spin 1 */
        chi2sq, /**< Magnitude of dimensionless spin 2 */
        chi1dotchi2, /**< Dot product of dimensionles spin 1 and spin 2 */
        qm_def1, /**< Quadrupole deformation parameter of body 1 (dimensionless) */
        qm_def2, /**< Quadrupole deformation parameter of body 2 (dimensionless) */
        spinO, /**< Enums specifying spin order are in LALSimInspiralWaveformFlags.h */
        NULL /**< nonGRParams pointer */
        );
    REAL8 eta = m1*m2/(m1+m2)/(m1+m2);
    REAL8 m1OverM = m1/(m1+m2);
    REAL8 m2OverM = m2/(m1+m2);
    // recover PN coeffs without prefactor for higher harminics 
    for(int i=0; i<=PN_PHASING_SERIES_MAX_ORDER; i++)
    {
      PN_coeffs[i] = 128.0*eta/3.0*pfa.v[i];
    }
    //PN_coeffs[0] = 1.0; /** 0th oprder newtonian */
    //PN_coeffs[1] = 0.0; /** 0.5 PN order v^1*/
    //PN_coeffs[2] = 3715.0/756.0 + 55.0*eta/9.0; /** 1.0 PN order v^2*/
    //PN_coeffs[3] = 4.0*spin_corrections->beta - 16.0*LAL_PI ; /** 1.5 PN order v^3*/
    //PN_coeffs[4] = 15293365.0/508032.0 + 27145.0*eta/504.0 + 3085.0*eta*eta/72.0 
    //                - 10.0*spin_corrections->sigma; /** 2.0 PN order v^4*/
    //PN_coeffs[5] = 38645.0*LAL_PI/756.0 - 65.0*LAL_PI*eta/9.0 - spin_corrections->gamma; /** 2.5 PN order v^5*/
    //PN_coeffs[6] = 11583231236531.0/4694215680.0 - 6848.0*euler_number/21.0 - 640.0*LAL_PI*LAL_PI/3.0
    //                + (2255.0*LAL_PI*LAL_PI/12.0 - 15737765635.0/3048192.0)*eta 
    //                + 76055.0*eta*eta/1728.0 - 127825.0*eta*eta*eta/1296.0; /** 3.0 PN order v^6 from Evans nb file*/
    //PN_coeffs[7] = LAL_PI*(77096675.0/254016.0 + 378515.0*eta/1512.0 
    //                - 74045.0*eta*eta/756.0); /** 3.5 PN order v^7 from Evans nb file*/
    //PN_coeffs[8] = 0.0; /** 4.0 PN order v^8*/
    //PN_coeffs[9] = 0.0; /** 4.5 PN order v^9*/
    /* Generate tidal terms separately.
     * Enums specifying tidal order are in LALSimInspiralWaveformFlags.h
     */
    REAL8 pft10 = 0.;
    REAL8 pft12 = 0.;
    switch( tideO )
    {
	case LAL_SIM_INSPIRAL_TIDAL_ORDER_ALL:
        case LAL_SIM_INSPIRAL_TIDAL_ORDER_6PN:
	    pft12 = (lambda1*XLALSimInspiralTaylorF2Phasing_12PNTidalCoeff(m1OverM) + lambda2*XLALSimInspiralTaylorF2Phasing_12PNTidalCoeff(m2OverM) );
        case LAL_SIM_INSPIRAL_TIDAL_ORDER_5PN:
            pft10 = ( lambda1*XLALSimInspiralTaylorF2Phasing_10PNTidalCoeff(m1OverM) + lambda2*XLALSimInspiralTaylorF2Phasing_10PNTidalCoeff(m2OverM) );
        case LAL_SIM_INSPIRAL_TIDAL_ORDER_0PN:
            break;
        default:
            XLAL_ERROR(XLAL_EINVAL, "Invalid tidal PN order %d", tideO);
    }
    PN_coeffs[12] += pft12;
    PN_coeffs[10] += pft10;
    return 0;
}

// prototype for new function for psi_SPA, SPA phase calculation in Eq. (4.81)
REAL8 sf2_psi_SPA(
    const REAL8 f,          /** frequency of GW  */
    const int k,            /** harmonics number  */
    const REAL8 shft,       /** shift factor 2*Pi*t_c */
    const REAL8 phic,       /** orbital coalescence phase (rad) */
    const LALSimInspiralTidalOrder tideO,  /**< flag to control tidal effects */
    const int phaseOrder,   /** twice PN phase order */
    const REAL8 *PN_coeffs, /** coeffs for each PN order*/
    const REAL8 m,          /** total mass */
    const REAL8 eta         /** symmetric mass ratio Eq (4.12) */
    );
// definition of psi SPA calculationi in Eq. (4.81)
REAL8 sf2_psi_SPA(
    const REAL8 f,          /** frequency of GW  */
    const int k,            /** harmonics number  */
    const REAL8 shft,       /** shift factor 2*Pi*t_c */
    const REAL8 phic,       /** orbital coalescence phase (rad) */
    const LALSimInspiralTidalOrder tideO,  /**< flag to control tidal effects */
    const int phaseOrder,   /** twice PN phase order */
    const REAL8 *PN_coeffs, /** coeffs for each PN order this should be calculated before calling this function*/
    const REAL8 m,          /** total mass */
    const REAL8 eta         /** symmetric mass ratio Eq (4.12) */
    )
{
    REAL8 psi;
    REAL8 f_k = f/k;
    const REAL8 m_sec = m * LAL_MTSUN_SI;  /* total mass in seconds */
    const REAL8 two_piM = LAL_TWOPI * m_sec;
    const REAL8 v = cbrt(two_piM*f_k);
    REAL8 vs[PN_PHASING_SERIES_MAX_ORDER+1]={0};
    REAL8 pre_factor;
    int i;
    // set values for power of v
    vs[0] = 1.0;
    vs[5] = v*v*v*v*v;
    //for(i=1; i<=phaseOrder; i++)  
    for(i=1; i<=PN_PHASING_SERIES_MAX_ORDER; i++)  
    {
        vs[i] = vs[i-1]*v;
    }
    pre_factor = 0.0117187500000/vs[5]/eta;    //3.0/256.0/vs[5]/eta;
    // calculate phase value upto the given order
    psi = 0;
    for(i=0; i<=phaseOrder; i++)
    {
        if(i==5)
        {
            psi += PN_coeffs[i]*(1.0+3.0*log(vs[1]))*vs[i];
        }
        else if(i==6)
        {
            psi += (PN_coeffs[i] - 6848.0/21.0*log(4.0*vs[1]))*vs[i];
        }
        else
        {
            psi += PN_coeffs[i]*vs[i];
        }
    }
    switch( tideO )
    {
	case LAL_SIM_INSPIRAL_TIDAL_ORDER_ALL:
        case LAL_SIM_INSPIRAL_TIDAL_ORDER_6PN:
            psi += PN_coeffs[12] * vs[12];
        case LAL_SIM_INSPIRAL_TIDAL_ORDER_5PN:
            psi += PN_coeffs[10] * vs[10];
        case LAL_SIM_INSPIRAL_TIDAL_ORDER_0PN:
            break;
        default:
            XLAL_ERROR(XLAL_EINVAL, "Invalid tidal PN order %d", tideO);
    }
    return shft*f_k - phic + pre_factor*psi;
}

// prototype for the function calculating amplitude corrections for
// given frequency, PN order, and harmonics.
// returns amplitude value of COMPLEX16 data
COMPLEX16 sf2_amp_SPA_plus(
    REAL8 f,       /** given ferquency */
    INT4 n,        /** twice order of PN for amplitude correction*/
    INT4 k,        /** harmonic number */
    REAL8 fStart,  /** initial frequency */
    REAL8 fISCO,   /** the final cut ferquency corresponding to f_LSO*/
    sf2_amp_corr_param *amp_corr_param /** various parameters*/
    );
// prototype for the function calculating amplitude corrections for
// given frequency, PN order, and harmonics
// returns amplitude value of COMPLEX16 data
COMPLEX16 sf2_amp_SPA_plus(
    REAL8 f,      /** given ferquency */
    INT4 n,       /** twice order of PN for amplitude correction*/
    INT4 k,       /** harmonic number */
    REAL8 fStart, /** initial frequency */
    REAL8 fISCO,  /** the final cut ferquency corresponding to f_LSO*/
    sf2_amp_corr_param *param /** various parameters*/
    )
{
    COMPLEX16 amp;
    REAL8 f_k = f/k;
    const REAL8 m_sec = param->m * LAL_MTSUN_SI;  /* total mass in seconds */
    const REAL8 two_piM = LAL_TWOPI * m_sec;
    const REAL8 v = cbrt(two_piM*f_k);
    REAL8 pre_factor;
    amp = 0.0 + 0.0j;
    if(f < fStart || f > k*fISCO) return amp; // return 0 if out of range
    //if(f < k*Fin || f > k*fISCO) return amp; // return 0 if out of range, Fin is the orbital frequency at start of observation Eq. (4.87)
    pre_factor = pow(v, n-3.5);
    switch(n)
    {
        case 0:  // Newtonian case
            if(k==2)
            {
                amp = -(1+param->costh2);
                amp = amp/sqrt(2.0);
            }
            else
            {
                amp = 0.0 + 0.0j;
            }
            break;
        case 1:  // 0.5 PN order
            if(k==1)
            {
                amp = -(5.0/8.0+param->costh2/8.0);
                amp = param->sinth*param->delta*amp;
                }
            else if(k==3)
            {
                amp = (1+param->costh2);
                amp = 9.0/8.0*param->sinth*param->delta*amp/sqrt(3.0);
            }
            else
            {
                amp = 0.0 + 0.0j;
            }
            break;
        case 2:  // 1.0 PN order
            if(k==1)
            {
                amp = 1;
                amp = param->sinth*(param->delta*param->XsLn + param->XaLn)*amp;
            }
            else if(k==2)
            {
                amp = (1385.0/672.0 - 109.0*param->eta/24.0 + 
                    (265.0/672.0+11.0*param->eta/24.0)*param->costh2 +
                    (-1.0/3.0 + param->eta)*param->costh4);
                amp = amp/sqrt(2.0);
            }
            else if(k==4)
            {
                amp = (1.0+param->costh2);
                amp = param->sinth2*amp*2/3.0*(3.0*param->eta-1.0);
            }
            else
            {
                amp = 0.0 + 0.0j;
            }
            break;
        case 3: // 1.5 PN order
            if(k==1)
            {
                amp = (-2119.0/5376.0 -263.0*param->eta/192.0 +
                    (937.0/5376.0 -3.0*param->eta/64.0)*param->costh2 +
                    (-1.0/192.0 + param->eta/96.0)*param->costh4);
                amp = param->sinth*param->delta*amp;
            }
            else if(k==2)
            {
                amp = (-27.0/8.0*(1.0+param->costh2)*param->delta*param->XaLn +
                    (-27.0/8.0*(1.0+param->costh2)+
                    0.5*(9.0-7.0*param->costh2)*param->eta)*param->XsLn);
                amp = amp/sqrt(2.0);
            }
            else if(k==3)
            {
                amp = (-6969.0/1792.0 + 81.0*param->eta/16.0 +
                    (-2811.0/1792.0+27.0*param->eta/64.0)*param->costh2 +
                    81.0/128.0 * (1.0-2.0*param->eta)*param->costh4);
                amp = param->sinth*param->delta*amp/sqrt(3.0);
            }
            else if(k==5)
            {
                amp = (1.0+param->costh2);
                amp = param->sinth3*param->delta*amp*625.0/384.0/sqrt(5.0)*(1.0-2.0*param->eta);
            }
            else
            {
                amp = 0.0 + 0.0j;
            }
            break;
        case 4: // 2.0 PN order
            if(k==1)
            {
                amp = (5.0*LAL_PI/8.0 + LAL_PI/8.0*param->costh2) 
                    + (((11.0/40.0+5.0*log(2.0)/4.0) + (7.0/40.0 + log(2.0)/4.0)*param->costh2))*1.0j;
                amp = param->sinth*param->delta*amp;
            }
            else if(k==2)
            {
                amp = (113419241.0/40642560.0 + 152987.0*param->eta/16128.0 - 
                    11099.0*param->eta*param->eta/1152.0 + (165194153.0/40642560.0 -
                    149.0*param->eta/1792.0 + 6709.0*param->eta*param->eta/1152.0)*param->costh2 +
                    (1693.0/2016.0 -5723.0*param->eta/2016.0 +13.0*param->eta*param->eta/12.0)*param->costh4 +
                    (-1.0/24.0 +5.0*param->eta/24.0 - 5.0*param->eta*param->eta/24.0)*param->costh6 + 
                    param->eta*(721.0/96.0*(param->XsLn*param->XsLn - param->XaLn*param->XaLn) - 
                    439.0/96.0*(param->Xs2 - param->Xa2))*(1.0 + param->costh2));
                amp = amp/sqrt(2.0);
            }
            else if(k==3)
            {
                amp = (9.0*LAL_PI/8.0*(1.0+param->costh*param->costh))
                    + ((-189.0/40.0 + 27.0*log(1.5)/4.0)*(1.0+param->costh*param->costh))*1.0j;
                amp = param->sinth*param->delta*amp/sqrt(3.0);
            }
            else if(k==4)
            {
                amp = (16109.0/2520.0 - 13367.0*param->eta/504.0 + 
                    39.0*param->eta*param->eta*0.5 +
                    (16.0/15.0 - 16.0*param->eta/3.0 + 16.0*param->eta*param->eta/3.0)*param->costh2*(param->costh4 - 3.0) + 
                    (-10733.0/2520.0 + 7991.0*param->eta/504.0 - 
                    53.0*param->eta*param->eta/6.0)*param->costh4);
                amp = 0.5*amp;
            }
            else if(k==6)
            {
                amp = (1.0 + param->costh2);
                amp = (-81.0/40.0 + 81.0*param->eta/8.0 - 81.0*param->eta*param->eta/8.0)*param->sinth4*amp/sqrt(6.0);
            }
            else
            {
                amp = 0.0 + 0.0j;
            }
            break;
        case 5: // 2.5 PN order
            if(k==1)
            {
                amp = (37533839.0/325140480.0 + 76171.0*param->eta/43008.0 - 
                    8407.0*param->eta*param->eta/4608.0 + 
                    (-29850823.0/325140480.0 + 56543.0*param->eta/129024.0 +
                    139.0*param->eta*param->eta/576.0)*param->costh2 + 
                    (255.0/14336.0 - 2659.0*param->eta/64512.0 +
                    127.0*param->eta*param->eta/9216.0)*param->costh4 + 
                    (-1.0/9216.0 + param->eta/2304.0 - param->eta*param->eta/3072.0)*param->costh6);
                amp = param->sinth*param->delta*amp;
            }
            else if(k==2)
            {
                amp = (85.0*LAL_PI/64.0*(1.0-4.0*param->eta)*(1.0+param->costh2))+
                    ((-9.0/5.0 + 32.0*param->eta + 
                    (14.0/5.0 *(1.0+4.0*param->eta))*param->costh2 + 
                    (7.0/5.0*(1.0-4.0*param->eta))*param->costh4))*1.0j; 
                amp = amp/sqrt(2.0);
            }
            else if(k==3)
            {
                amp = (-8781361.0/7225344.0 -366781.0*param->eta/17920.0 +
                    15193.0*param->eta*param->eta/1280.0 + 
                    (-238136057.0/36126720.0 + 37829.0*param->eta/71680.0 -
                    7073.0*param->eta*param->eta/1280.0)*param->costh2 + 
                    (-328347.0/143360.0 + 396009.0*param->eta/71680.0 -
                    10557.0*param->eta*param->eta/5120.0)*param->costh4 + 
                    (729.0/5120.0 - 729.0*param->eta/1280.0 +
                    2187.0*param->eta*param->eta/5120.0)*param->costh6);
                amp = param->sinth*param->delta*amp/sqrt(3.0);
            }
            else if(k==4)
            {
                amp = (8.0*LAL_PI/3.0*(3.0*param->eta - 1.0)*(1.0+param->costh2))+
                    ((56.0/5.0 - 1193.0*param->eta/30.0 +
                    32.0*log(2.0)/3.0*(3.0*param->eta - 1.0))*(1.0+param->costh2))*1.0j;
                amp = 0.5*param->sinth2*amp;
            }
            else if(k==5)
            {
                amp = (-854375.0/86016.0 + 3919375.0*param->eta/129024.0 -
                    160625.0*param->eta*param->eta/9216.0 + 
                    (40625.0/9216.0 - 40625.0*param->eta/2304.0 +
                    40625.0*param->eta*param->eta/3072.0)*param->costh2 + 
                    (1863125.0/258048.0 - 2519375.0*param->eta/129024.0 +
                    85625.0*param->eta*param->eta/9216.0)*param->costh4 + 
                    (-15625.0/9216.0 +15625.0*param->eta/2304.0 -
                    15625.0*param->eta*param->eta/3072.0)*param->costh6);
                amp = param->sinth*param->delta*amp/sqrt(5.0);
            }
            else if(k==7)
            {
                amp =(1.0+param->costh2);
                amp = param->sinth5*param->delta*amp/sqrt(7.0)*
                      (117649.0/46080.0 - 117649.0*param->eta/11520.0 + 117649.0*param->eta*param->eta/15360.0);
            }
            else
            {
                amp = 0.0 + 0.0j;
            }
            break;
        case 6:  // 3.0 PN order
            amp = 0.0 + 0.0j;
            break;
        case 7:  // 3.5 PN order
            amp = 0.0 + 0.0j;
            break;
        case 8:  // 4.0 PN order
            amp = 0.0 + 0.0j;
            break;
        case 9:  // 4.5 PN order
            amp = 0.0 + 0.0j;
            break;
        default:
            amp = 0.0 + 0.0j;
            break;
    }
    return pre_factor*amp;
}

COMPLEX16 sf2_amp_SPA_cross(
    REAL8 f,      /** given ferquency */
    INT4 n,       /** twice order of PN for amplitude correction*/
    INT4 k,       /** harmonic number */
    REAL8 fStart, /** initial frequency */
    REAL8 fISCO,  /** the final cut ferquency corresponding to f_LSO*/
    sf2_amp_corr_param *param /** various parameters*/
    );
COMPLEX16 sf2_amp_SPA_cross(
    REAL8 f,      /** given ferquency */
    INT4 n,       /** twice order of PN for amplitude correction*/
    INT4 k,       /** harmonic number */
    REAL8 fStart, /** initial frequency */
    REAL8 fISCO,  /** the final cut ferquency corresponding to f_LSO*/
    sf2_amp_corr_param *param /** various parameters*/
    )
{
    COMPLEX16 amp;
    REAL8 f_k = f/k;
    const REAL8 m_sec = param->m * LAL_MTSUN_SI;  /* total mass in seconds */
    const REAL8 two_piM = LAL_TWOPI * m_sec;
    const REAL8 v = cbrt(two_piM*f_k);
    REAL8 pre_factor;
    amp = 0.0 + 0.0j;
    if(f < fStart || f > k*fISCO) return amp; // return 0 if out of range
    //if(f < k*Fin || f > k*fISCO) return amp; // return 0 if out of range, Fin is the orbital frequency at the start of observation Eq. (4.87)
    pre_factor = pow(v, n-3.5);
    switch(n)
    {
        case 0:  // Newtonian case
            if(k==2)
            {
                amp = - 2.0*param->costh*1.0j;
                amp = amp/sqrt(2.0);
            }
            else
            {
                amp = 0.0 + 0.0j;
            }
            break;
        case 1:  // 0.5 PN order
            if(k==1)
            {
                amp = - 0.75*param->costh*1.0j;
                amp = param->sinth*param->delta*amp;
        }
            else if(k==3)
            {
                amp = 2.0*param->costh*1.0j;
                amp = 9.0/8.0*param->sinth*param->delta*amp/sqrt(3.0);
            }
            else
            {
                amp = 0.0 + 0.0j;
            }
            break;
        case 2:  // 1.0 PN order
            if(k==1)
            {
                amp = param->costh*1.0j;
                amp = param->sinth*(param->delta*param->XsLn + param->XaLn)*amp;
            }
            else if(k==2)
            {
                amp = ((387.0/112.0-85.0*param->eta/12.0)*param->costh +
                    (-4.0/3.0 + 4.0*param->eta)*param->costh3)*1.0j;
                amp = amp/sqrt(2.0);
            }
            else if(k==4)
            {
                amp = (2.0*param->costh)*1.0j;
                amp = param->sinth2*amp*2/3.0*(3.0*param->eta-1.0);
            }
            else
            {
                amp = 0.0 + 0.0j;
            }
            break;
        case 3: // 1.5 PN order
            if(k==1)
            {
                amp = (-(155.0/896.0 + 145.0*param->eta/96.0)*param->costh +
                    5.0/96.0*(2.0*param->eta-1.0)*param->costh3)*1.0j;
                amp = param->sinth*param->delta*amp;
            }
            else if(k==2)
            {
                amp = (-27.0/4.0*param->delta*param->XaLn +
                    (-27.0/4.0 +(5.0-4.0*param->costh2)*param->eta)*param->XsLn)*param->costh*1.0j;
                amp = amp/sqrt(2.0);
            }
            else if(k==3)
            {
                amp = ((-6213.0/896.0 + 135.0*param->eta/16.0)*param->costh + 
                    135.0/64.0*(1.0-2.0*param->eta)*param->costh3)*1.0j;
                amp = param->sinth*param->delta*amp/sqrt(3.0);
            }
            else if(k==5)
            {
                amp = (2.0*param->costh)*1.0j;
                amp = param->sinth3*param->delta*amp*625.0/384.0/sqrt(5.0)*(1.0-2.0*param->eta);
            }
            else
            {
                amp = 0.0 + 0.0j;
            }
            break;
        case 4: // 2.0 PN order
            if(k==1)
            {
                amp = -(9.0/20.0 + 3.0*log(2.0)/2.0)*param->costh
                    + ((0.75*LAL_PI)*param->costh)*1.0j;
                amp = param->sinth*param->delta*amp;
            }
            else if(k==2)
            {
                amp = ((114020009.0/20321280.0 + 133411*param->eta/8064.0 -
                    7499.0*param->eta*param->eta/576.0)*param->costh + 
                    (5777.0/2520.0 - 5555.0*param->eta/504.0 + 34.0*param->eta*param->eta/3.0)*param->costh3 + 
                    (-0.25 + 5.0*param->eta/4.0 - 5.0*param->eta*param->eta/4.0)*param->costh5 +
                    param->eta*(721.0/96.0*(param->XsLn*param->XsLn - param->XaLn*param->XaLn) - 
                    439.0/96.0*(param->Xs2 - param->Xa2))*2.0*param->costh)*1.0j;
                amp = amp/sqrt(2.0);
            }
            else if(k==3)
            {
                amp = (-2.0*(-189.0/40.0 + 27.0*log(1.5)/4.0)*param->costh)
                    + (9.0*LAL_PI/4.0*param->costh)*1.0j;
                amp = param->sinth*param->delta*amp/sqrt(3.0);
            }
            else if(k==4)
            {
                amp = ((2953.0/252.0 - 12023.0*param->eta/252.0 +
                    101.0*param->eta*param->eta/3.0)*param->costh + 
                    (-18797.0/1260.0 + 16055.0*param->eta/252.0 - 
                    149.0*param->eta*param->eta/3.0)*param->costh3 + 
                    (16.0/5.0 -16.0*param->eta + 16.0*param->eta*param->eta)*param->costh5)*1.0j;
                amp = 0.5*amp;
            }
            else if(k==6)
            {
                amp = (2.0*param->costh)*1.0j;
                amp = (-81.0/40.0 + 81.0*param->eta/8.0 - 81.0*param->eta*param->eta/8.0)*param->sinth4*amp/sqrt(6.0);
            }
            else
            {
                amp = 0.0 + 0.0j;
            }
            break;
        case 5: // 2.5 PN order
            if(k==1)
            {
                amp = ((-3453823.0/54190080.0 +163015.0*param->eta/64512.0 -
                    4237.0*param->eta*param->eta/2304.0)*param->costh + 
                    (34373.0/322560.0 - 11755.0*param->eta/32256.0 +
                    631.0*param->eta*param->eta/2304.0)*param->costh3 + 
                    (-7.0/4608.0 + 7.0*param->eta/1152.0 - 7.0*param->eta*param->eta/1536.0)*param->costh5)*1.0j;
                amp = param->sinth*param->delta*amp;
            }
            else if(k==2)
            {
                amp = ((2.0 - 282.0*param->eta/5.0)*param->costh + (-22.0/5.0 + 94.0*param->eta/5.0)*param->costh3)+
                    ((85.0*LAL_PI/32.0*(1.0 - 4.0*param->eta))*param->costh)*1.0j; 
                amp = amp/sqrt(2.0);
            }
            else if(k==3)
            {
                amp = ((-63633869.0/18063360.0 - 89609.0*param->eta/2560.0 +
                    697.0*param->eta*param->eta/40.0)*param->costh + 
                    (-508689.0/71680.0 + 812727.0*param->eta/35840.0 -
                    4707.0*param->eta*param->eta/320.0)*param->costh3 + 
                    (1701.0/2560.0 - 1701.0*param->eta/640.0 +
                    5103.0*param->eta*param->eta/2560.0)*param->costh5)*1.0j;
                amp = param->sinth*param->delta*amp/sqrt(3.0);
            }
            else if(k==4)
            {
                amp = -    (56.0/5.0 - 1193.0*param->eta/30.0 +
                    32.0*log(2.0)/3.0*(3.0*param->eta - 1.0))*2.0*param->costh +
                    ((8.0*LAL_PI/3.0*(3.0*param->eta - 1.0))*2.0*param->costh)*1.0j;
                amp = 0.5*param->sinth2*amp;
            }
            else if(k==5)
            {
                amp = ((-2388125.0/129024.0 + 3569375.0*param->eta/64512.0 -
                    141875.0*param->eta*param->eta/4608.0)*param->costh + 
                    (3000625.0/129024.0 - 1598125.0*param->eta/21504.0 +
                    51875.0*param->eta*param->eta/1152.0)*param->costh3 + 
                    (-21875.0/4608.0 + 21875.0*param->eta/1152.0 -
                    21875.0*param->eta*param->eta/1536.0)*param->costh5)*1.0j;
                amp = param->sinth*param->delta*amp/sqrt(5.0);
            }
            else if(k==7)
            {
                amp =(2.0*param->costh)*1.0j;
                amp = param->sinth5*param->delta*amp/sqrt(7.0)*
                        (117649.0/46080.0 - 117649.0*param->eta/11520.0 + 
                         117649.0*param->eta*param->eta/15360.0);
            }
            else
            {
                amp = 0.0 + 0.0j;
            }
            break;
        case 6:  // 3.0 PN order
            amp = 0.0 + 0.0j;
            break;
        case 7:  // 3.5 PN order
            amp = 0.0 + 0.0j;
            break;
        case 8:  // 4.0 PN order
            amp = 0.0 + 0.0j;
            break;
        case 9:  // 4.5 PN order
            amp = 0.0 + 0.0j;
            break;
        default:
            amp = 0.0 + 0.0j;
            break;
    }
    return pre_factor*amp;
}
/**
 * Computing the stationary phase approximation to the Fourier transform of
 * a chirp waveform with phase given by Eq.(4.82) 
 * and amplitude given by Eq. (4.72) with C_k^(n) given in Appendix D  
 * PN order is the highest order being used to calculate Eq. (4.72)
 */
int XLALSimInspiralTaylorF2AmpPlus(
    COMPLEX16FrequencySeries **htilde_out, /* frequency-domain waveform */
    REAL8 phic,                     /* orbital coalescence phase (rad) */
    REAL8 deltaF,                   /* sampling frequency (Hz) */
    REAL8 inclination,              /* inclination of source (rad), corresponds to theta of N_hat */
    REAL8 m1_SI,                    /* mass1 (kg) */
    REAL8 m2_SI,                    /* mass2 (kg) */
    REAL8 S1x,                      /* initial value of S1x */
    REAL8 S1y,                      /* initial value of S1y */
    REAL8 S1z,                      /* initial value of S1z */
    REAL8 S2x,                      /* initial value of S2x */
    REAL8 S2y,                      /* initial value of S2y */    
    REAL8 S2z,                      /* initial value of S2z */
    REAL8 LNhatx,                   /* initial value of LNhatx */
    REAL8 LNhaty,                   /* initial value of LNhaty */
    REAL8 LNhatz,                   /* initial value of LNhatz */
    const REAL8 f_ref,              /**< Reference GW frequency (Hz) - if 0 reference point is coalescence */
    REAL8 fStart,                   /* start GW frequency (Hz) */
    REAL8 f_max0,                   /* ending GW frequency (Hz) */
    REAL8 r,                        /* distance of source (m) */
    const REAL8 quadparam1,                /**< quadrupole deformation parameter of body 1 (dimensionless, 1 for BH) */
    const REAL8 quadparam2,                /**< quadrupole deformation parameter of body 2 (dimensionless, 1 for BH) */
    const REAL8 lambda1,                   /**< (tidal deformation of body 1)/(mass of body 1)^5 */
    const REAL8 lambda2,                   /**< (tidal deformation of body 2)/(mass of body 2)^5 */
    const LALSimInspiralSpinOrder spinO,  /**< twice PN order of spin effects */
    const LALSimInspiralTidalOrder tideO,  /**< flag to control tidal effects */
    int phaseO,                     /* twice PN phase order */
    int amplitudeO                  /* twice PN amplitude order */
   )
{
    static int calls_debug = 0;
    const REAL8 lambda = -1987./3080.;
    const REAL8 theta = -11831./9240.;

    /* external: SI; internal: solar masses */
    const REAL8 m1 = m1_SI / LAL_MSUN_SI;
    const REAL8 m2 = m2_SI / LAL_MSUN_SI;
    const REAL8 m = m1 + m2;
    const REAL8 m_sec = m * LAL_MTSUN_SI;  /* total mass in seconds */
    const REAL8 eta = m1 * m2 / (m * m);
    const REAL8 delta = (m1 - m2)/ m;
    const REAL8 piM = LAL_PI * m_sec;
    const REAL8 vISCO = 1. / sqrt(6.);
    const REAL8 fISCO = vISCO * vISCO * vISCO / piM;
    const REAL8 v0 = cbrt(piM * fStart);
    REAL8 shft, amp0, f_max, f;
    size_t i, n, k, iStart, iEnd;
    int mm;
    COMPLEX16 *data = NULL; /** actual storage for amplitude data pointer */
    COMPLEX16FrequencySeries *htilde; /** waveform data */
    LIGOTimeGPS tC = {0, 0}; 
    REAL8 overall_factor; 
    REAL8 PN_coeffs_SPA[PN_PHASING_SERIES_MAX_ORDER+1]; /** PN order phasing calculation coefficient series*/
    sf2_spin_corr_amp spin_corrections_SPA; /** spin correction coeffs beta, sigma, gamma Eqs. (4.82,83,84)*/
    REAL8 f_plus, f_cross, costh, sinth;
    COMPLEX16 amp;
    sf2_amp_corr_param amp_corr_param;
    const REAL8 m1OverM = m1 / m;
    const REAL8 m2OverM = m2 / m;

    REAL8 phasing;
    //COMPLEX16 prec_fac;

    //fprintf(stdout, "=========== DEBUG by Jeongcho Kim & Chunglee KIm ===========\n");
    //fprintf(stdout, "==Function : LALSimInspiralTaylorF2Amp.c:XLALSimInspiralTaylorF2AmpPlus()\n");
    //fprintf(stdout, "Sequence of calls = %d\n", ++calls_debug);
    //fprintf(stdout, "(S1x, S1y, S1z) = (%f, %f, %f)\n", S1x, S1y, S1z);
    //fprintf(stdout, "(S2x, S2y, S2z) = (%f, %f, %f)\n", S2x, S2y, S2z);
    //fprintf(stdout, "=========================================================\n");

    /** calculate frequency independent coeffcients */
    overall_factor = m*m/r*LAL_MTSUN_SI*LAL_MRSUN_SI*sqrt(5.0*LAL_PI*eta/48.0); /** overall factor of Eq. (4.72) */
    // calculate spin correction coeffs
    //spin_corrections_SPA = 
    //    sf2_spin_corrections_amp(m1, m2, S1x, S1y, S1z, S2x, S2y, S2z, LNhatx, LNhaty, LNhatz);
    // calculate PN_coeffs_SPA 
    REAL8 chi1L = S1x*LNhatz + S1y*LNhaty + S1z*LNhatz;
    REAL8 chi2L = S2x*LNhatz + S2y*LNhaty + S2z*LNhatz;
    REAL8 chi1sq = S1x*S1x + S1y*S1y + S1z*S1z;
    REAL8 chi2sq = S2x*S2x + S2y*S2y + S2z*S2z;
    REAL8 chi1dotchi2 = S1x*S2x + S1y*S2y + S1z*S2z;
    sf2_psi_SPA_coeffs_PN_order(PN_coeffs_SPA, 
        m1, /**< Mass of body 1, in Msol */
        m2, /**< Mass of body 2, in Msol */
        chi1L, /**< Component of dimensionless spin 1 along Lhat */
        chi2L, /**< Component of dimensionless spin 2 along Lhat */
        chi1sq, /**< Magnitude of dimensionless spin 1 */
        chi2sq, /**< Magnitude of dimensionless spin 2 */
        chi1dotchi2, /**< Dot product of dimensionles spin 1 and spin 2 */
        quadparam1, /**< Quadrupole deformation parameter of body 1 (dimensionless) */
        quadparam2, /**< Quadrupole deformation parameter of body 2 (dimensionless) */
        lambda1,                   /**< (tidal deformation of body 1)/(mass of body 1)^5 */
        lambda2,                   /**< (tidal deformation of body 2)/(mass of body 2)^5 */
        tideO,  /**< flag to control tidal effects */
        spinO /**< Enums specifying spin order are in LALSimInspiralWaveformFlags.h */
        );

    f_cross = 1.0; //0.5*(1+costh*costh)*cos(2*ra)*sin(2*psi) + costh*sin(2*ra)*cos(2*psi);
    f_plus  = 1.0; //0.5*(1+costh*costh)*cos(2*ra)*cos(2*psi) - costh*sin(2*ra)*sin(2*psi);
    amp_corr_param.costh = cos(inclination); // cosine of inclination of the source
    amp_corr_param.sinth = sin(inclination); // sine of inclination of the source
    amp_corr_param.costh2 = amp_corr_param.costh*amp_corr_param.costh;
    amp_corr_param.costh3 = amp_corr_param.costh2*amp_corr_param.costh;
    amp_corr_param.costh4 = amp_corr_param.costh3*amp_corr_param.costh;
    amp_corr_param.costh5 = amp_corr_param.costh4*amp_corr_param.costh;
    amp_corr_param.costh6 = amp_corr_param.costh5*amp_corr_param.costh;
    amp_corr_param.sinth2 = amp_corr_param.sinth*amp_corr_param.sinth;
    amp_corr_param.sinth3 = amp_corr_param.sinth2*amp_corr_param.sinth;
    amp_corr_param.sinth4 = amp_corr_param.sinth3*amp_corr_param.sinth;
    amp_corr_param.sinth5 = amp_corr_param.sinth4*amp_corr_param.sinth;
    amp_corr_param.f_cross = f_cross;
    amp_corr_param.f_plus = f_plus;
    amp_corr_param.m = m;
    amp_corr_param.eta = eta;
    amp_corr_param.delta = delta;
    amp_corr_param.s1x = S1x;
    amp_corr_param.s1y = S1y;
    amp_corr_param.s1z = S1z;
    amp_corr_param.s2x = S2x;
    amp_corr_param.s2y = S2y;
    amp_corr_param.s2z = S2z;
    amp_corr_param.LNhatx = LNhatx;
    amp_corr_param.LNhaty = LNhaty;
    amp_corr_param.LNhatz = LNhatz;
    amp_corr_param.Xs_x = 0.5*(amp_corr_param.s1x + amp_corr_param.s2x);
    amp_corr_param.Xs_y = 0.5*(amp_corr_param.s1y + amp_corr_param.s2y);
    amp_corr_param.Xs_z = 0.5*(amp_corr_param.s1z + amp_corr_param.s2z);
    amp_corr_param.Xa_x = 0.5*(amp_corr_param.s1x - amp_corr_param.s2x);
    amp_corr_param.Xa_y = 0.5*(amp_corr_param.s1y - amp_corr_param.s2y);
    amp_corr_param.Xa_z = 0.5*(amp_corr_param.s1z - amp_corr_param.s2z);
    amp_corr_param.Xs2 = amp_corr_param.Xs_x*amp_corr_param.Xs_x+amp_corr_param.Xs_y*amp_corr_param.Xs_y+amp_corr_param.Xs_z*amp_corr_param.Xs_z;  // Xs*Xs (=chi_s dot chi_s)
    amp_corr_param.Xa2 = amp_corr_param.Xa_x*amp_corr_param.Xa_x+amp_corr_param.Xa_y*amp_corr_param.Xa_y+amp_corr_param.Xa_z*amp_corr_param.Xa_z;  // Xa*Xa (=chi_a dot chi_a)
    amp_corr_param.XsXa = amp_corr_param.Xs_x*amp_corr_param.Xa_x + amp_corr_param.Xs_y*amp_corr_param.Xa_y + amp_corr_param.Xs_z*amp_corr_param.Xa_z;  // Xs dot Xa (=xchi_s dot chi_a)
    amp_corr_param.XsLn = (amp_corr_param.Xs_x*amp_corr_param.LNhatx)+(amp_corr_param.Xs_y*amp_corr_param.LNhaty)+(amp_corr_param.Xs_z*amp_corr_param.LNhatz); // chi_s dot Ln_hat
    amp_corr_param.XaLn = (amp_corr_param.Xa_x*amp_corr_param.LNhatx)+(amp_corr_param.Xa_y*amp_corr_param.LNhaty)+(amp_corr_param.Xa_z*amp_corr_param.LNhatz); // chi_a dot Ln_hat

    /* initial checks */
    if (!htilde_out) {
        printf("htilde_out is NULL.\n");
        XLAL_ERROR(XLAL_EFAULT);
    }
    if (*htilde_out) {
        printf("*htilde_out is NOT NULL.\n");
        XLAL_ERROR(XLAL_EFAULT);
    }
    if (m1_SI <= 0) XLAL_ERROR(XLAL_EDOM);
    if (m2_SI <= 0) XLAL_ERROR(XLAL_EDOM);
    if (fStart <= 0) XLAL_ERROR(XLAL_EDOM);
    if (r <= 0) XLAL_ERROR(XLAL_EDOM);

    /* allocate htilde */
    if ( f_max0 == 0. ) // End at ISCO
        //f_max = MAX_HARMONICS*0.5*fISCO; // experiment for upper limit
        f_max = fISCO;
    else // End at user-specified freq.
        f_max = f_max0;
    n = (size_t) (f_max / deltaF + 1);

    if (amplitudeO < 0)
        amplitudeO = MAX_AMP_ORDER;
    if (phaseO < 0)
        amplitudeO = MAX_PHASE_ORDER;


    XLALGPSAdd(&tC, -1 / deltaF);  /* coalesce at t=0 */
    htilde = XLALCreateCOMPLEX16FrequencySeries("hplustilde: FD waveform", &tC, 0.0, deltaF, &lalStrainUnit, n);
    if (!htilde) {
        XLAL_ERROR(XLAL_EFUNC);
    }
    else {
    }
    memset(htilde->data->data, 0, n * sizeof(COMPLEX16));
    XLALUnitDivide(&htilde->sampleUnits, &htilde->sampleUnits, &lalSecondUnit);

    /* extrinsic parameters */
    shft = LAL_TWOPI * (tC.gpsSeconds + 1e-9 * tC.gpsNanoSeconds); // 2*Pi*t_c

    iStart = (size_t) ceil(fStart / deltaF);
    iEnd = n;
    //iEnd = (size_t) (f_max / deltaF);
    //iEnd = (iEnd < n) ? iEnd : n;  /* overflow protection; should we warn? */

    // add f_ref effect
    REAL8 phasing_ref[MAX_HARMONICS+1]={0};
    if(f_ref > 0.0)
    {
      for(k=1; k<=MAX_HARMONICS; k++)
      {
        phasing_ref[k] = sf2_psi_SPA(f_ref, k, shft, phic, tideO, phaseO, PN_coeffs_SPA, m, eta)
                         - shft*f_ref/k + phic; // we need this since we want only Psi_SPA^(k)(f_ref/k)
      }
    }

    data = htilde->data->data;
    for (i = iStart; i < iEnd; i++) {
      f = i * deltaF;
      data[i] = 0.0 + 0.0j;
      for (k = 1; k <= MAX_HARMONICS; k++) // up to 7th harmonics
      {            
        phasing = sf2_psi_SPA(f, k, shft, phic, tideO, phaseO, PN_coeffs_SPA, m, eta) - phasing_ref[k];
        for (n = 0; n <= (size_t)amplitudeO; n++)
        {
          amp = sf2_amp_SPA_plus(f, n, k, fStart, fISCO, &amp_corr_param);
          data[i] += amp*(cos(k*phasing - LAL_PI_4) - sin(k*phasing - LAL_PI_4)*1.0j); 
        }
      }
      data[i] = overall_factor*data[i];
    }

    *htilde_out = htilde;
    return XLAL_SUCCESS;
}
/**
 * Computes the stationary phase approximation to the Fourier transform of
 * a chirp waveform with phase given by Eq.(4.82)
 * and amplitude given by Eq. (4.72) with C_k^(n) given in Appendix D. 
 * PN order is the highest order being used by Eq (4.72)
 */
int XLALSimInspiralTaylorF2AmpCross(
    COMPLEX16FrequencySeries **htilde_out, /** frequency-domain waveform */
    REAL8 phic,                            /** orbital coalescence phase (rad) */
    REAL8 deltaF,                          /** sampling frequency (Hz) */
    REAL8 inclination,                     /** inclination of source (rad), corresponds to theta of N_hat */
    REAL8 m1_SI,                           /** mass 1 (kg) */
    REAL8 m2_SI,                           /** mass 2 (kg) */
    REAL8 S1x,                             /** initial value of S1x */
    REAL8 S1y,                             /** initial value of S1y */
    REAL8 S1z,                             /** initial value of S1z */
    REAL8 S2x,                             /** initial value of S2x */
    REAL8 S2y,                             /** initial value of S2y */    
    REAL8 S2z,                             /** initial value of S2z */
    REAL8 LNhatx,                          /** initial value of LNhatx */
    REAL8 LNhaty,                          /** initial value of LNhaty */
    REAL8 LNhatz,                          /** initial value of LNhatz */
    const REAL8 f_ref,                     /**< Reference GW frequency (Hz) - if 0 reference point is coalescence */
    REAL8 fStart,                          /** start GW frequency (Hz) */
    REAL8 f_max0,                          /** ending GW frequency (Hz) */
    REAL8 r,                               /** distance of source (m) */
    const REAL8 quadparam1,                /**< quadrupole deformation parameter of body 1 (dimensionless, 1 for BH) */
    const REAL8 quadparam2,                /**< quadrupole deformation parameter of body 2 (dimensionless, 1 for BH) */
    const REAL8 lambda1,                   /**< (tidal deformation of body 1)/(mass of body 1)^5 */
    const REAL8 lambda2,                   /**< (tidal deformation of body 2)/(mass of body 2)^5 */
    const LALSimInspiralSpinOrder spinO,  /**< twice PN order of spin effects */
    const LALSimInspiralTidalOrder tideO,  /**< flag to control tidal effects */
    int phaseO,                            /** twice PN phase order */
    int amplitudeO                         /** twice PN amplitude order */
   )
{
    //static int calls_debug = 0;
    //const REAL8 lambda = -1987./3080.;
    //const REAL8 theta = -11831./9240.;

    /* external: SI; internal: solar masses */
    const REAL8 m1 = m1_SI / LAL_MSUN_SI;
    const REAL8 m2 = m2_SI / LAL_MSUN_SI;
    const REAL8 m = m1 + m2;
    const REAL8 m_sec = m * LAL_MTSUN_SI;  /* total mass in seconds */
    const REAL8 eta = m1 * m2 / (m * m);
    const REAL8 delta = (m1 - m2)/ m;
    const REAL8 piM = LAL_PI * m_sec;
    const REAL8 vISCO = 1. / sqrt(6.);
    const REAL8 fISCO = vISCO * vISCO * vISCO / piM;
    //const REAL8 v0 = cbrt(piM * fStart);
    REAL8 shft, f_max, f;
    size_t i, n, k, iStart, iEnd;
    //int mm;
    COMPLEX16 *data = NULL; /** actual storage for amplitude data pointer */
    COMPLEX16FrequencySeries *htilde; /** waveform data */
    LIGOTimeGPS tC = {0, 0}; 
    REAL8 overall_factor; 
    REAL8 PN_coeffs_SPA[PN_PHASING_SERIES_MAX_ORDER+1]; /** PN order phasing calculation coefficient series*/
    //sf2_spin_corr_amp spin_corrections_SPA; /** spin correction coeffs beta, sigma, gamma Eqs. (4.82,83,84)*/
    REAL8 f_plus, f_cross;
    COMPLEX16 amp;
    sf2_amp_corr_param amp_corr_param;
    //const REAL8 m1OverM = m1 / m;
    //const REAL8 m2OverM = m2 / m;

    REAL8 alpha, alpha_ref, zeta, zeta_ref, beta, phasing;
    COMPLEX16 prec_fac;

    //fprintf(stdout, "=========== DEBUG by Jeongcho Kim & Chunglee KIm ===========\n");
    //fprintf(stdout, "==Function : LALSimInspiralTaylorF2Amp.c:XLALSimInspiralTaylorF2AmpCross()\n");
    //fprintf(stdout, "Sequence of calls = %d\n", ++calls_debug);
    //fprintf(stdout, "(S1x, S1y, S1z) = (%f, %f, %f)\n", S1x, S1y, S1z);
    //fprintf(stdout, "(S2x, S2y, S2z) = (%f, %f, %f)\n", S2x, S2y, S2z);
    //fprintf(stdout, "=========================================================\n");

    /** calculate frequency independent coeffcients */
    overall_factor = m*m/r*LAL_MTSUN_SI*LAL_MRSUN_SI*sqrt(5.0*LAL_PI*eta/48.0); /** overall factor of Eq. (4.72) */
    // calculate spin correction coeffs
    //spin_corrections_SPA = 
    //    sf2_spin_corrections_amp(m1, m2, S1x, S1y, S1z, S2x, S2y, S2z, LNhatx, LNhaty, LNhatz);
    // calculate PN_coeffs_SPA 
    REAL8 chi1L = S1x*LNhatz + S1y*LNhaty + S1z*LNhatz;
    REAL8 chi2L = S2x*LNhatz + S2y*LNhaty + S2z*LNhatz;
    REAL8 chi1sq = S1x*S1x + S1y*S1y + S1z*S1z;
    REAL8 chi2sq = S2x*S2x + S2y*S2y + S2z*S2z;
    REAL8 chi1dotchi2 = S1x*S2x + S1y*S2y + S1z*S2z;
    sf2_psi_SPA_coeffs_PN_order(PN_coeffs_SPA, 
        m1, /**< Mass of body 1, in Msol */
        m2, /**< Mass of body 2, in Msol */
        chi1L, /**< Component of dimensionless spin 1 along Lhat */
        chi2L, /**< Component of dimensionless spin 2 along Lhat */
        chi1sq, /**< Magnitude of dimensionless spin 1 */
        chi2sq, /**< Magnitude of dimensionless spin 2 */
        chi1dotchi2, /**< Dot product of dimensionles spin 1 and spin 2 */
        quadparam1, /**< Quadrupole deformation parameter of body 1 (dimensionless) */
        quadparam2, /**< Quadrupole deformation parameter of body 2 (dimensionless) */
        lambda1,                   /**< (tidal deformation of body 1)/(mass of body 1)^5 */
        lambda2,                   /**< (tidal deformation of body 2)/(mass of body 2)^5 */
        tideO,  /**< flag to control tidal effects */
        spinO /**< Enums specifying spin order are in LALSimInspiralWaveformFlags.h */
        );

    f_cross = 1.0; //0.5*(1+costh*costh)*cos(2*ra)*sin(2*psi) + costh*sin(2*ra)*cos(2*psi);
    f_plus  = 1.0; //0.5*(1+costh*costh)*cos(2*ra)*cos(2*psi) - costh*sin(2*ra)*sin(2*psi);
    amp_corr_param.costh = cos(inclination); // cosine of inclination of the source
    amp_corr_param.sinth = sin(inclination); // sine of inclination of the source
    amp_corr_param.costh2 = amp_corr_param.costh*amp_corr_param.costh;
    amp_corr_param.costh3 = amp_corr_param.costh2*amp_corr_param.costh;
    amp_corr_param.costh4 = amp_corr_param.costh3*amp_corr_param.costh;
    amp_corr_param.costh5 = amp_corr_param.costh4*amp_corr_param.costh;
    amp_corr_param.costh6 = amp_corr_param.costh5*amp_corr_param.costh;
    amp_corr_param.sinth2 = amp_corr_param.sinth*amp_corr_param.sinth;
    amp_corr_param.sinth3 = amp_corr_param.sinth2*amp_corr_param.sinth;
    amp_corr_param.sinth4 = amp_corr_param.sinth3*amp_corr_param.sinth;
    amp_corr_param.sinth5 = amp_corr_param.sinth4*amp_corr_param.sinth;
    amp_corr_param.f_cross = f_cross;
    amp_corr_param.f_plus = f_plus;
    amp_corr_param.m = m;
    amp_corr_param.eta = eta;
    amp_corr_param.delta = delta;
    amp_corr_param.s1x = S1x;
    amp_corr_param.s1y = S1y;
    amp_corr_param.s1z = S1z;
    amp_corr_param.s2x = S2x;
    amp_corr_param.s2y = S2y;
    amp_corr_param.s2z = S2z;
    amp_corr_param.LNhatx = LNhatx;
    amp_corr_param.LNhaty = LNhaty;
    amp_corr_param.LNhatz = LNhatz;
    amp_corr_param.Xs_x = 0.5*(amp_corr_param.s1x + amp_corr_param.s2x);
    amp_corr_param.Xs_y = 0.5*(amp_corr_param.s1y + amp_corr_param.s2y);
    amp_corr_param.Xs_z = 0.5*(amp_corr_param.s1z + amp_corr_param.s2z);
    amp_corr_param.Xa_x = 0.5*(amp_corr_param.s1x - amp_corr_param.s2x);
    amp_corr_param.Xa_y = 0.5*(amp_corr_param.s1y - amp_corr_param.s2y);
    amp_corr_param.Xa_z = 0.5*(amp_corr_param.s1z - amp_corr_param.s2z);
    amp_corr_param.Xs2 = amp_corr_param.Xs_x*amp_corr_param.Xs_x+amp_corr_param.Xs_y*amp_corr_param.Xs_y+amp_corr_param.Xs_z*amp_corr_param.Xs_z;  // Xs*Xs (=chi_s dot chi_s)
    amp_corr_param.Xa2 = amp_corr_param.Xa_x*amp_corr_param.Xa_x+amp_corr_param.Xa_y*amp_corr_param.Xa_y+amp_corr_param.Xa_z*amp_corr_param.Xa_z;  // Xa*Xa (=chi_a dot chi_a)
    amp_corr_param.XsXa = amp_corr_param.Xs_x*amp_corr_param.Xa_x + amp_corr_param.Xs_y*amp_corr_param.Xa_y + amp_corr_param.Xs_z*amp_corr_param.Xa_z;  // Xs dot Xa (=xchi_s dot chi_a)
    amp_corr_param.XsLn = (amp_corr_param.Xs_x*amp_corr_param.LNhatx)+(amp_corr_param.Xs_y*amp_corr_param.LNhaty)+(amp_corr_param.Xs_z*amp_corr_param.LNhatz); // chi_s dot Ln_hat
    amp_corr_param.XaLn = (amp_corr_param.Xa_x*amp_corr_param.LNhatx)+(amp_corr_param.Xa_y*amp_corr_param.LNhaty)+(amp_corr_param.Xa_z*amp_corr_param.LNhatz); // chi_a dot Ln_hat

    /* Perform some initial checks */
    if (!htilde_out) {
        XLAL_ERROR(XLAL_EFAULT);
    }
    if (*htilde_out) {
        XLAL_ERROR(XLAL_EFAULT);
    }
    if (m1_SI <= 0) XLAL_ERROR(XLAL_EDOM);
    if (m2_SI <= 0) XLAL_ERROR(XLAL_EDOM);
    if (fStart <= 0) XLAL_ERROR(XLAL_EDOM);
    if (r <= 0) XLAL_ERROR(XLAL_EDOM);

    /* allocate htilde */
    if ( f_max0 == 0. ) // End at ISCO
        //f_max = MAX_HARMONICS*0.5*fISCO; // expperiment for upper limit for integration
        f_max = fISCO; 
    else // End at user-specified freq.
        f_max = f_max0;
    n = (size_t)(f_max / deltaF + 1);

    if (amplitudeO < 0)
        amplitudeO = MAX_AMP_ORDER;
    if (phaseO < 0)
        amplitudeO = MAX_PHASE_ORDER;

    XLALGPSAdd(&tC, -1 / deltaF);  /* coalesce at t=0 */
    htilde = XLALCreateCOMPLEX16FrequencySeries("hcrosstilde: FD waveform", &tC, 0.0, deltaF, &lalStrainUnit, n);
    if (!htilde) {
        XLAL_ERROR(XLAL_EFUNC);
    }
    else {
    }
    memset(htilde->data->data, 0, n * sizeof(COMPLEX16));
    XLALUnitDivide(&htilde->sampleUnits, &htilde->sampleUnits, &lalSecondUnit);

    /* extrinsic parameters */
    shft = LAL_TWOPI * (tC.gpsSeconds + 1e-9 * tC.gpsNanoSeconds); // 2*Pi*t_c

    iStart = (size_t) ceil(fStart / deltaF);
    iEnd = n;
    // add f_ref effect
    REAL8 phasing_ref[MAX_HARMONICS+1]={0};
    if(f_ref > 0.0)
    {
      for(k=1; k<=MAX_HARMONICS; k++)
      {
        phasing_ref[k] = sf2_psi_SPA(f_ref, k, shft, phic, tideO, phaseO, PN_coeffs_SPA, m, eta)
                         - shft*f_ref/k + phic; // we need this since we want only Psi_SPA^(k)(f_ref/k)
      }
    }

    data = htilde->data->data;
    for (i = iStart; i < iEnd; i++) {
      f = i * deltaF;
      data[i] = 0.0 + 0.0j;
      for (k = 1; k <= MAX_HARMONICS; k++) // up to 7th harmonics
      {
        phasing = sf2_psi_SPA(f, k, shft, phic, tideO, phaseO, PN_coeffs_SPA, m, eta) - phasing_ref[k];
        for (n = 0; n <= amplitudeO; n++)
        {
            amp = sf2_amp_SPA_cross(f, n, k, fStart, fISCO, &amp_corr_param);
            data[i] += amp*(cos(k*phasing - LAL_PI_4) - sin(k*phasing - LAL_PI_4)*1.0j); // changed sign for time reversal
        }
      }
      data[i] = overall_factor*data[i];
    }

    *htilde_out = htilde;
    return XLAL_SUCCESS;
}

/** @} */
/** @} */
