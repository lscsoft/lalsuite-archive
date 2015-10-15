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
  gsl_vector* degVector, knotsVector, dataVector
  LALH5File *group = XLALH5GroupOpen(file, keyName)

  ReadHDF5RealVectorDataset(group, "deg", &degVector)
  ReadHDF5RealVectorDataset(group, "knots", &knotsVector)
  ReadHDF5RealVectorDataset(group, "data", &dataVector)

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
  REAL8 Mflower, time_start_M, time_start_s, time_end_m, time_end_s;
  REAL8 chi, est_start_time;
  size_t tmp_length;
  gsl_vector *tmpVector;
  LALH5File *file, *group;

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
  XLALH5FileQueryScalarAttributeValue(&Mflower, file, "f_lower_at_1MSUN")
  /* Figure out start time of data */
  *group = XLALH5GroupOpen(file, "amp_l1_m2")
  ReadHDF5RealVectorDataset(group, "knots", &tmpVector)
  time_start_M = (REAL8)(gsl_vector_get(tmpVector, 0));
  time_end_M = (REAL8)(gsl_vector_get(tmpVector, tmpVector->size - 1));
  gsl_vector_free(tmpVector);
  time_start_s = time_start_M * (m1 + m2)
  time_end_s = time_end_M * (m1 + m2)
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
  
    
  /* Create the return time series */
  *hplus  = XLALCreateREAL8TimeSeries("H_PLUS", &tc, 0.0, deltaT,
                                      &lalStrainUnit, sigReVec->length );
  *hcross = XLALCreateREAL8TimeSeries("H_CROSS", &tc, 0.0, deltaT,
                                      &lalStrainUnit, sigImVec->length );


  XLALH5FileClose(file);
}


#endif

