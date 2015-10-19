#ifdef LAL_HDF5_ENABLED
#include <lal/H5FileIO.h>
static const INT4 ROMDataHDF5_VERSION_MAJOR = 1;
static const INT4 ROMDataHDF5_VERSION_MINOR = 0;
static const INT4 ROMDataHDF5_VERSION_MICRO = 0;

int XLALSimInspiralNRWaveformGetDataFromHDF5File(
        REAL8TimeSeries **output,               /**< Returned vector uncompressed */
        LALH5File pointer,                      /**< Pointer to HDF5 file */
        LIGOTimeGPS epoch,                  /**< Start time of returned vector */
        size_t length,                          /**< Length of returned vector */
        REAL8 deltaT,                           /**< Sample rate of returned vector */
        const char *keyName                     /**< Name of vector to uncompress */
        )
{
  gsl_vector* degVector, knotsVector, dataVector;
  LALH5File *group = XLALH5GroupOpen(file, keyName);

  ReadHDF5RealVectorDataset(group, "deg", &degVector);
  ReadHDF5RealVectorDataset(group, "knots", &knotsVector);
  ReadHDF5RealVectorDataset(group, "data", &dataVector);

  *output = XLALCreateREAL8TimeSeries(keyName, &epoch, 0, deltaT,
                                      &lalStrainUnit, length)

  /* FIXME: The univariate spline bit goes in here 
   *        The 3 gsl_vectors are the input and output->data->data the Vector
   *        to write to. 
   */

  gsl_vector_free(degVector);
  gsl_vector_free(knotsVector);
  gsl_vector_free(dataVector);

  return output
}

int XLALSimInspiralNRWaveformGetHplusHcross(
        REAL8TimeSeries **hplus,                /**< OUTPUT h_+ vector */
        REAL8TimeSeries **hcross,               /**< OUTPUT h_x vector */
        REAL8 phiRef,                   /**< orbital phase at reference pt. */
        REAL8 inclination,              /**< inclination angle */
        REAL8 deltaT,                   /**< sampling interval (s) */
        REAL8 m1,                       /**< mass of companion 1 (kg) */
        REAL8 m2,                       /**< mass of companion 2 (kg) */
        REAL8 fStart,                   /**< start GW frequency (Hz) */
        REAL8 fRef,                     /**< reference GW frequency (Hz) */
        REAL8 s1x,                      /**< initial value of S1x */
        REAL8 s1y,                      /**< initial value of S1y */
        REAL8 s1z,                      /**< initial value of S1z */
        REAL8 s2x,                      /**< initial value of S2x */
        REAL8 s2y,                      /**< initial value of S2y */
        REAL8 s2z,                      /**< initial value of S2z */
        const char *NRDataFile,         /**< Location of NR HDF file */
        )
{
  /* Declarations */
  UINT4 model, modem, test_int, curr_idx;
  size_t array_length;
  REAL8 Mflower, time_start_M, time_start_s, time_end_m, time_end_s;
  REAL8 chi, est_start_time, curr_h_real, curr_h_imag, curr_ylm;
  /* These keys follow a strict formulation and cannot be longer than 11
   * characters */
  char *amp_key = XLALMalloc(20);
  char *phase_key = XLALMalloc(20);
  gsl_vector *tmpVector;
  LALH5File *file, *group;
  LIGOTimeGPS tmpEpoch;
  tmpEpoch.gpsSeconds = 0;
  tmpEpoch.gpsNanoSeconds = 0;
  REAL8TimeSeries curr_amp, curr_phase;

  *file = XLALH5FileOpen(NRDataFile, "r");

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
  *group = XLALH5GroupOpen(file, "amp_l1_m2");
  ReadHDF5RealVectorDataset(group, "knots", &tmpVector);
  time_start_M = (REAL8)(gsl_vector_get(tmpVector, 0));
  time_end_M = (REAL8)(gsl_vector_get(tmpVector, tmpVector->size - 1));
  gsl_vector_free(tmpVector);
  time_start_s = time_start_M * (m1 + m2);
  time_end_s = time_end_M * (m1 + m2);
  fprintf(stderr, "hybrid start time in M: %f", time_start_M);
  fprintf(stderr, "hybrdid start time in s: %f", time_start_s);

  /* We don't want to return the *entire* waveform if it will be much longer
   * than the specified f_lower. Therefore guess waveform length using
   * the SEOBNR_ROM function and add 10% for safety.
   * FIXME: Is this relevant for precessing waveforms?
   */
  chi = XLALSimIMRPhenomBComputeChi(m1, m2, s1z, s2z);
  est_start_time = XLALSimIMRSEOBNRv2ChirpTimeSingleSpin(m1, m2, chi, fStart);
  est_start_time = (-est_start_time) * 1.1;
  if (est_start_time > time_start_s)
  {
    /* Restrict start time of waveform */
    time_start_s = est_start_time;
    time_start_M = time_start_s / (m1 + m2);
  else if (fStart < Mflower / (m1 + m2))
  {
     XLAL_ERROR(XLAL_EDOM, "WAVEFORM IS NOT LONG ENOUGH TO REACH f_low. %e %e",
                fStart, Mflower / (m1 + m2));
  }

  array_length = UINT4( (time_end_s - time_start_s) / delta_t);
  
  /* Create the return time series, use arbitrary epoch here. We set this
   * properly later. */
  *hplus  = XLALCreateREAL8TimeSeries("H_PLUS", &tmpEpoch, 0.0, deltaT,
                                      &lalStrainUnit, array_length );
  *hcross = XLALCreateREAL8TimeSeries("H_CROSS", &tmpEpoch, 0.0, deltaT,
                                      &lalStrainUnit, array_length );

  /* Create the distance scale factor */
  distance_scale_fac = (m1 + m2) * MRSUN_SI / (MTSUN_SI * PC_SI * 1.0e6);

  /* Generate the waveform */
  for (model=2; model < 9 ; model++)
  {
    for (modem=-model; model < (model+1); modem++)
    {
      sprintf(amp_key, "amp_l%d_m%d", model, modem);
      sprintf(phase_key, "phase_l%d_m%d", model, modem);

      /* Check that both groups exist */
      if H5Gget_objinfo(file, amp_key, 0, NULL)
      {
        continue;
      }
      if H5Gget_objinfo(file, phase_key, 0, NULL)
      {
        continue;
      }
      fprintf(stderr, "Using %d, %d mode.\n", model, modem);

      /* Get amplitude and phase from file */
      XLALSimInspiralNRWaveformGetDataFromHDF5File(&curr_amp, file,
                                                array_length, deltaT, amp_key);
      XLALSimInspiralNRWaveformGetDataFromHDF5File(&curr_phase, file, 
                                              array_length, deltaT, phase_key);

      for (curr_idx = 0; curr_idx < array_length; curr_idx++)
      {
        curr_h_real = curr_amp.data.data[curr_idx]
                    * cos(curr_phase.data.data[curr_idx]) * distance_scale_fac;
        curr_h_imag = curr_amp.data.data[curr_idx] 
                    * sin(curr_phase.data.data[curr_idx]) * distance_scale_fac;
        hplus.data.data[curr_idx] = hplus.data.data[curr_idx]
                  + curr_h_real * curr_ylm.real - curr_h_imag * curr_ylm.imag;
        /* FIXME: No idea whether these should be minus or plus, guessing
         *        minus for now. Only affects some polarization phase.
         */
        hcross.data.data[curr_idx] = hcross.data.data[curr_idx]
                  - curr_h_real * curr_ylm.imag - curr_h_imag * curr_ylm.real;
        XLALDestroyREAL8TimeSeries(curr_amp);
        XLALDestroyREAL8TimeSeries(curr_phase);
      }
    }
  }

  XLALFree(phase_key);
  XLALFree(amp_key);
  XLALH5FileClose(file);
}


#endif

