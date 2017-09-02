/*
 * Copyright (C) 2015 Sebastian Khan, Michael Puerrer
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

#include <math.h>
#include <lal/LALSimIMR.h>
#include <lal/FrequencySeries.h>
#include <lal/Sequence.h>
#include <lal/Units.h>
#include <lal/LALConstants.h>


#ifndef _OPENMP
#define omp ignore
#endif


// FIXME: limits below
/**
 * @addtogroup LALSimIMRTIDAL_c
 *
 * @{
 *
 * @name IMRPhenomD_NRTidal
 *
 * @author Sebastian Khan, Michael Puerrer
 *
 * @brief C code for IMRPhenomD Husa:2015iqa, Khan:2015jqa with added
 * tidal phase correction from arXiv:1706.02969.
 *
 * This is a frequency domain model that adds tidal modifications of
 * the phasing to the IMRPhenomD model.
 *
 * @note Parameter ranges:
 *   * ? <= eta <= 0.25
 *   * 0 <= Lambda_i <= ?
 *   * -1 <= chi_i <= 1
 *
 *  Aligned component spin on neutron stars.
 *  Symmetric mass-ratio eta = m1*m2/(m1+m2)^2.
 *  Total mass Mtot.
 *
 * @{
 */


/**
 * Driver routine to compute the spin-aligned, inspiral-merger-ringdown
 * phenomenological waveform IMRPhenomD in the frequency domain.
 *
 * Reference:
 * - Waveform: Eq. 35 and 36 in arXiv:1508.07253
 * - Coefficients: Eq. 31 and Table V in arXiv:1508.07253
 *
 *  All input parameters should be in SI units. Angles should be in radians.
 */
int XLALSimIMRPhenomD_NRTidal_GenerateFD(
    COMPLEX16FrequencySeries **htilde, /**< [out] FD waveform */
    const REAL8 phi0,                  /**< Orbital phase at fRef (rad) */
    const REAL8 fRef_in,               /**< reference frequency (Hz) */
    const REAL8 deltaF,                /**< Sampling frequency (Hz) */
    const REAL8 m1_SI_in,                 /**< Mass of companion 1 (kg) */
    const REAL8 m2_SI_in,                 /**< Mass of companion 2 (kg) */
    const REAL8 chi1_in,                  /**< Aligned-spin parameter of companion 1 */
    const REAL8 chi2_in,                  /**< Aligned-spin parameter of companion 2 */
    const REAL8 f_min,                 /**< Starting GW frequency (Hz) */
    const REAL8 f_max,                 /**< End frequency; 0 defaults to Mf = \ref f_CUT for IMRPhenomD */
    const REAL8 distance,               /**< Distance of source (m) */
    const REAL8 lambda1_in,               /**< (tidal deformability of mass 1) / m1^5 (dimensionless) */
    const REAL8 lambda2_in,               /**< (tidal deformability of mass 2) / m2^5 (dimensionless) */
    const LALSimInspiralTestGRParam *extraParams /**< linked list containing the extra testing GR parameters */
) {

  /* swap masses, spins and lambdas to enforce that m1>=m2 */
  REAL8 chi1, chi2, m1_SI, m2_SI, lambda1, lambda2;
  if (m1_SI_in>m2_SI_in) {
     lambda1 = lambda1_in;
     lambda2 = lambda2_in;
     chi1 = chi1_in;
     chi2 = chi2_in;
     m1_SI   = m1_SI_in;
     m2_SI   = m2_SI_in;
  } else { // swap spins, lambdas and masses
     lambda1 = lambda2_in;
     lambda2 = lambda1_in;
     chi1 = chi2_in;
     chi2 = chi1_in;
     m1_SI   = m2_SI_in;
     m2_SI   = m1_SI_in;
  }


  /* external: SI; internal: solar masses */
  const REAL8 m1 = m1_SI / LAL_MSUN_SI;
  const REAL8 m2 = m2_SI / LAL_MSUN_SI;

  /* check inputs for sanity */
  XLAL_CHECK(0 != htilde, XLAL_EFAULT, "htilde is null");
  if (*htilde) XLAL_ERROR(XLAL_EFAULT);
  if (fRef_in < 0) XLAL_ERROR(XLAL_EDOM, "fRef_in must be positive (or 0 for 'ignore')\n");
  if (deltaF <= 0) XLAL_ERROR(XLAL_EDOM, "deltaF must be positive\n");
  if (m1 <= 0) XLAL_ERROR(XLAL_EDOM, "m1 must be positive\n");
  if (m2 <= 0) XLAL_ERROR(XLAL_EDOM, "m2 must be positive\n");
  if (f_min <= 0) XLAL_ERROR(XLAL_EDOM, "f_min must be positive\n");
  if (f_max < 0) XLAL_ERROR(XLAL_EDOM, "f_max must be greater than 0\n");
  if (distance <= 0) XLAL_ERROR(XLAL_EDOM, "distance must be positive\n");

  if (chi1 > 1.0 || chi1 < -1.0 || chi2 > 1.0 || chi2 < -1.0)
    XLAL_ERROR(XLAL_EDOM, "Spins outside the range [-1,1] are not supported\n");

  // if no reference frequency given, set it to the starting GW frequency
  REAL8 fRef = (fRef_in == 0.0) ? f_min : fRef_in;

  int status = XLALSimIMRPhenomDGenerateFD(
      htilde,
      phi0,
      fRef,
      deltaF,
      m1_SI,
      m2_SI,
      chi1,
      chi2,
      f_min,
      f_max,
      distance,
      extraParams
  );
  XLAL_CHECK(XLAL_SUCCESS == status, status, "Failed to generate IMRPhenomD waveform.");

  size_t ind_min = (size_t) (f_min / deltaF);
  // size_t ind_max = (size_t) (f_max / deltaF);
  size_t ind_max = (size_t) ((*htilde)->data->length);

  // Get FD tidal phase correction from arXiv:1706.02969
  REAL8Sequence *freqs = XLALCreateREAL8Sequence(ind_max);
  // populate frequencies
  for (size_t i = ind_min; i < ind_max; i++)
  {
    REAL8 fHz = i * deltaF; // geometric frequency
    freqs->data[i] = fHz;
  }

  REAL8Sequence *phi_tidal = XLALCreateREAL8Sequence(freqs->length);
  REAL8Sequence *amp_tidal = XLALCreateREAL8Sequence(freqs->length);
  status = XLALSimNRTunedTidesFDTidalPhaseFrequencySeries(
    phi_tidal, amp_tidal, freqs,
    m1_SI, m2_SI, lambda1, lambda2
  );
  XLAL_CHECK(XLAL_SUCCESS == status, status, "XLALSimNRTunedTidesFDTidalPhaseFrequencySeries Failed.");

  /* loop over htilde and apply taper and phase shift */
  for (size_t i=ind_min; i<ind_max; i++) { // loop over frequency points in sequence
    // Apply tidal phase correction and amplitude taper
    COMPLEX16 Corr = amp_tidal->data[i] * cexp(-I*phi_tidal->data[i]);
    (*htilde)->data->data[i] *= Corr;
  }

  XLALDestroyREAL8Sequence(freqs);
  XLALDestroyREAL8Sequence(phi_tidal);
  XLALDestroyREAL8Sequence(amp_tidal);

  return XLAL_SUCCESS;
}

/** @} */

/** @} */
