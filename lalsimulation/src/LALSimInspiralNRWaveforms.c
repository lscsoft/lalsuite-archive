#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <time.h>
#include <stdbool.h>
#include <alloca.h>
#include <string.h>
#include <libgen.h>

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

#include <gsl/gsl_errno.h>
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_spline.h>
#include <lal/Units.h>
#include <lal/SeqFactories.h>
#include <lal/LALConstants.h>
#include <lal/XLALError.h>
#include <lal/FrequencySeries.h>
#include <lal/Date.h>
#include <lal/StringInput.h>
#include <lal/Sequence.h>
#include <lal/LALStdio.h>
#include <lal/FileIO.h>

#include <lal/LALSimInspiral.h>
#include <lal/LALSimIMR.h>
#include <lal/LALConfig.h>
#include <lal/SphericalHarmonics.h>

#include <hdf5.h>
#include <lal/H5FileIO.h>

#include "LALSimIMRSEOBNRROMUtilities.c"

static const INT4 ROMDataHDF5_VERSION_MAJOR = 1;
static const INT4 ROMDataHDF5_VERSION_MINOR = 0;
static const INT4 ROMDataHDF5_VERSION_MICRO = 0;

static UINT4 XLALSimInspiralNRWaveformGetDataFromHDF5File(
  REAL8Vector** output,                    /**< Returned vector uncompressed */
  LALH5File* pointer,                     /**< Pointer to HDF5 file */
  REAL8 totalMass,                        /**< Total mass of system for scaling */
  REAL8 startTime,                        /**< Start time of veturn vector */
  size_t length,                          /**< Length of returned vector */
  REAL8 deltaT,                           /**< Sample rate of returned vector */
  const char *keyName                     /**< Name of vector to uncompress */
  )
{
  UINT4 idx;
  size_t comp_data_length;
  REAL8 massTime;
  gsl_interp_accel *acc;
  gsl_spline *spline;
  gsl_vector *knotsVector, *dataVector;
  LALH5File *group = XLALH5GroupOpen(pointer, keyName);
  knotsVector=dataVector=NULL;

  ReadHDF5RealVectorDataset(group, "knots", &knotsVector);
  ReadHDF5RealVectorDataset(group, "data", &dataVector);

  *output = XLALCreateREAL8Vector(length);

  comp_data_length = dataVector->size;
  /* SPLINE STUFF */
  acc = gsl_interp_accel_alloc();
  spline = gsl_spline_alloc(gsl_interp_cspline, comp_data_length);
  gsl_spline_init(spline, knotsVector->data, dataVector->data,
                  comp_data_length);

  for (idx = 0; idx < length; idx++)
  {
    massTime = (startTime + idx*deltaT) / (totalMass * LAL_MTSUN_SI);
    (*output)->data[idx] = gsl_spline_eval(spline, massTime, acc);
  }

  gsl_vector_free(knotsVector);
  gsl_vector_free(dataVector);
  gsl_spline_free (spline);
  gsl_interp_accel_free (acc);

  return XLAL_SUCCESS;
}


int XLALSimInspiralNRWaveformGetHplusHcross(
        REAL8TimeSeries **hplus,               /**< Output h_+ vector */
        REAL8TimeSeries **hcross,              /**< Output h_x vector */
        REAL8 phiRef,                          /**< orbital phase at reference pt. */
        REAL8 inclination,                     /**< inclination angle */
        REAL8 deltaT,                          /**< sampling interval (s) */
        REAL8 m1,                              /**< mass of companion 1 (kg) */
        REAL8 m2,                              /**< mass of companion 2 (kg) */
        REAL8 r,                               /**< distance of source (m) */
        REAL8 fStart,                          /**< start GW frequency (Hz) */
        UNUSED REAL8 fRef,                     /**< reference GW frequency (Hz) */
        UNUSED REAL8 s1x,                      /**< initial value of S1x */
        UNUSED REAL8 s1y,                      /**< initial value of S1y */
        UNUSED REAL8 s1z,                      /**< initial value of S1z */
        UNUSED REAL8 s2x,                      /**< initial value of S2x */
        UNUSED REAL8 s2y,                      /**< initial value of S2y */
        UNUSED REAL8 s2z,                      /**< initial value of S2z */
        const char *NRDataFile                 /**< Location of NR HDF file */
        )
{
  /* Declarations */
  UINT4 curr_idx;
  INT4 model, modem;
  size_t array_length;
  REAL8 Mflower, time_start_M, time_start_s, time_end_M, time_end_s;
  REAL8 chi, est_start_time, curr_h_real, curr_h_imag;
  REAL8 distance_scale_fac;
  COMPLEX16 curr_ylm;
  /* These keys follow a strict formulation and cannot be longer than 11
   * characters */
  char *amp_key = XLALMalloc(20);
  char *phase_key = XLALMalloc(20);
  gsl_vector *tmpVector=NULL;
  LALH5File *file, *group;
  LIGOTimeGPS tmpEpoch;
  tmpEpoch.gpsSeconds = 0;
  tmpEpoch.gpsNanoSeconds = 0;
  REAL8Vector *curr_amp, *curr_phase;

  /* I can't understand mass in seconds, so use solar masses. NR files will use
   * solar masses as well, so easier for that conversion
   */
  m1 = m1 / LAL_MSUN_SI;
  m2 = m2 / LAL_MSUN_SI;

  file = XLALH5FileOpen(NRDataFile, "r");

  /* FIXME: Add here some sanity checks that input is sane.
   *        Examples, are m1, m2 s1, s2 consistent with file?
   *        Something like
   * XLALH5FileQueryScalarAttributeValue(&nrMass1, file, "mass1')
   * if fabs((mass1 - m1)/mass1) < 1E-4
   * {
   *    XLALERROR()
   * }
   */

  /* First figure out time series that is needed.
   * Demand that 22 mode that is present and use that to figure this out
   */
  XLALH5FileQueryScalarAttributeValue(&Mflower, file, "f_lower_at_1MSUN");
  /* Figure out start time of data */
  group = XLALH5GroupOpen(file, "amp_l2_m2");
  ReadHDF5RealVectorDataset(group, "knots", &tmpVector);
  time_start_M = (REAL8)(gsl_vector_get(tmpVector, 0));
  time_end_M = (REAL8)(gsl_vector_get(tmpVector, tmpVector->size - 1));
  gsl_vector_free(tmpVector);
  time_start_s = time_start_M * (m1 + m2) * LAL_MTSUN_SI;
  time_end_s = time_end_M * (m1 + m2) * LAL_MTSUN_SI;

  /* We don't want to return the *entire* waveform if it will be much longer
   * than the specified f_lower. Therefore guess waveform length using
   * the SEOBNR_ROM function and add 10% for safety.
   * FIXME: Is this relevant for precessing waveforms?
   */
  chi = XLALSimIMRPhenomBComputeChi(m1, m2, s1z, s2z);
  est_start_time = XLALSimIMRSEOBNRv2ChirpTimeSingleSpin(m1 * LAL_MSUN_SI,
                                                m2 * LAL_MSUN_SI, chi, fStart);
  est_start_time = (-est_start_time) * 1.1;
  if (est_start_time > time_start_s)
  {
    /* Restrict start time of waveform */
    time_start_s = est_start_time;
    time_start_M = time_start_s / ((m1 + m2) * LAL_MTSUN_SI);
  }
  else if (fStart < Mflower / (m1 + m2) )
  {
     XLAL_ERROR(XLAL_EDOM, "WAVEFORM IS NOT LONG ENOUGH TO REACH f_low. %e %e",
                fStart, Mflower, Mflower / (m1 + m2));
  }

  array_length = (UINT4)( (time_end_s - time_start_s) / deltaT);
  
  /* Create the return time series, use arbitrary epoch here. We set this
   * properly later. */
  *hplus  = XLALCreateREAL8TimeSeries("H_PLUS", &tmpEpoch, 0.0, deltaT,
                                      &lalStrainUnit, array_length );
  *hcross = XLALCreateREAL8TimeSeries("H_CROSS", &tmpEpoch, 0.0, deltaT,
                                      &lalStrainUnit, array_length );
  for (curr_idx = 0; curr_idx < array_length; curr_idx++)
  {
    (*hplus)->data->data[curr_idx] = 0.;
    (*hcross)->data->data[curr_idx] = 0.;
  }

  /* Create the distance scale factor */
  distance_scale_fac = (m1 + m2) * LAL_MRSUN_SI / r;

  /* Generate the waveform */
  for (model=2; model < 9 ; model++)
  {
    for (modem=-model; modem < (model+1); modem++)
    {
      sprintf(amp_key, "amp_l%d_m%d", model, modem);
      sprintf(phase_key, "phase_l%d_m%d", model, modem);

      /* Check that both groups exist */
      if (XLALH5CheckGroupExists(file, amp_key))
      {
        continue;
      }
      if (XLALH5CheckGroupExists(file, phase_key))
      {
        continue;
      }

      /* Get amplitude and phase from file */
      XLALSimInspiralNRWaveformGetDataFromHDF5File(&curr_amp, file, (m1 + m2),
                                  time_start_s, array_length, deltaT, amp_key);
      XLALSimInspiralNRWaveformGetDataFromHDF5File(&curr_phase, file, (m1 + m2),
                                time_start_s, array_length, deltaT, phase_key);

      /* FIXME: Why is this -2? It is also in the python code */
      curr_ylm = XLALSpinWeightedSphericalHarmonic(inclination, phiRef, -2,
                                                   model, modem);
      for (curr_idx = 0; curr_idx < array_length; curr_idx++)
      {
        curr_h_real = curr_amp->data[curr_idx]
                    * cos(curr_phase->data[curr_idx]) * distance_scale_fac;
        curr_h_imag = curr_amp->data[curr_idx] 
                    * sin(curr_phase->data[curr_idx]) * distance_scale_fac;
        (*hplus)->data->data[curr_idx] = (*hplus)->data->data[curr_idx]
               + curr_h_real * creal(curr_ylm) + curr_h_imag * cimag(curr_ylm);
        /* FIXME: No idea whether these should be minus or plus, guessing
         *        minus for now. Only affects some polarization phase.
         */
        (*hcross)->data->data[curr_idx] = (*hcross)->data->data[curr_idx]
               + curr_h_real * cimag(curr_ylm) - curr_h_imag * creal(curr_ylm);
      }
      XLALDestroyREAL8Vector(curr_amp);
      XLALDestroyREAL8Vector(curr_phase);
    }
  }

  XLALFree(phase_key);
  XLALFree(amp_key);
  XLALH5FileClose(file);

  return XLAL_SUCCESS;
}
