#include <pycbc.h>

/* more includes... */

/* basic function for returning data starting with a frame cache */
real_vector_t * get_data_from_cache(
    const char        *cache_name,
    const char        *channel_name,
    unsigned long int t_start,
    double            duration
    )
{
  real_vector_t   *output;
  FrCache         *cache;
  FrStream        *stream;
  REAL4TimeSeries *lal_timeseries;

  if (cache_name == NULL)
    return NULL;
  if (channel_name == NULL)
    return NULL;
  if (duration < 0) 
    return NULL;
  
  /* open the data cache and use it to get a frame stream */
  cache  = XLALFrImportCache( cache_name );
  stream = XLALFrCacheOpen( cache );
  XLALFrDestroyCache( cache );

  /* set the mode of the frame stream */
  XLALFrSetMode( stream, LAL_FR_VERBOSE_MODE );

  /* read the series */
  lal_timeseries =  get_data_from_frame_stream( channel_name, t_start, duration, stream );

  /* convert to pycbc vector */
  output = convert_timeseries_lal2pycbc ( lal_timeseries);

  /* close the frame stream */
  XLALFrClose( stream );

  return output;
}


/* low-level routine to read single-precision frame data 
   mostly taken from ring in lalapps
   -- should not used at the pycb top level */
REAL4TimeSeries * get_data_from_frame_stream( const char *channel_name,
					    unsigned_long_int t_start_nsec, 
					    double duration, 
					    FrStream *stream )
{
  double srate;
  int npoints;
  ldiv_t t_start_sec;
  REAL4TimeSeries *output_timeseries;
  int vector_size;
  LIGOTimeGPS t_start_gps;

  /** some error checking */
  if ( channel_name == NULL)
    return NULL;
  if (duration < 0)
    return NULL;
  if (stream == NULL)
    return NULL;

  output_timeseries = calloc( 1, sizeof( *output_timeseries ) );

  /* use glibc for long integer division */
  /* make this into a function or use lal?*/
  t_start_sec = ldiv( t_start_nsec, PYLAL_BILLION_INT);
  t_start_gps.gpsSeconds = t_start_sec.quot;
  t_start_gps.gpsNanoSeconds = t_start_sec.rem;
  output_timeseries->epoch = t_start_gps;  

  XLALFrSeek( stream, &t_start_gps );
 
  strncpy( output_timeseries->name, channel_name, sizeof( output_timeseries->name ) - 1 );

  /* call first time to get sample rate */
  XLALFrGetREAL4TimeSeriesMetadata( output_timeseries, stream );

  /* compute sample rateand number of points to request */
  srate   = 1.0/output_timeseries->deltaT;
  npoints = floor( duration*srate + 0.5 ); /* round to nearest integer */

  /* create the output data */
  output_timeseries->data = XLALCreateREAL4Vector( npoints );

  /* now get the data */
  XLALFrGetREAL4TimeSeries( &series, stream );
  
  return output_timeseries;
}

/* return value is number of segments if it worked */
static int get_data_segments(
    real_vector_t           **output_segments,
    real_vector_t           *input_timeseries,
    real_vector_t           *invspec,
    REAL4FFTPlan            *fwdplan,
    int                     num_overlap_segments,
    double                  segment_duration,
    double                  stride_duration
    )
{
  real_vector_t *segments=NULL;
  UINT4  index;
  REAL4TimeSeries *input_timeseries_lal;
  REAL4FrequencySeries *invspec_lal;

  /* some error checking */
  if (input_data == NULL) 
    return -1;
  if (fwdplan == NULL)
    return -1;
  /* ...more */

  /* memory for segments */
  segments = calloc( num_overlap_segments, sizeof(*segments) );

  /* convert to lal */
  input_timeseries_lal =  convert_timeseries_pycbc2lal (input_timeseries);
  invspec_lal = convert_freqseries_pycbc2lal (invspec);

  for ( index = 0; index < num_overlap_segments; ++index )
    compute_data_segment( segments[index], input_timeseries_lal, invspec_lal,
			  segment_duration, stride_duration, fwdplan );

  /* free the converted data */
  /* note-- don't free the data contained in these structures because these are 
     not owned by these lal structures */
  free(input_timeseries_lal);
  free(invspec_lal);

  output_segments = &segments;
  return num_overlap_segments;
}



/* routine to compute a single overwhitened data segment */
int compute_data_segment(
    COMPLEX8FrequencySeries  *segment,
    UINT4                     segmentNumber,
    REAL4TimeSeries          *series,
    REAL4FrequencySeries     *invspec,
    REAL8                     segmentDuration,
    REAL8                     strideDuration,
    REAL4FFTPlan             *fwdPlan
    )
{
  REAL4TimeSeries seg;
  REAL4Vector     vec;
  INT8            ns;
  UINT4           segmentLength;
  UINT4           segmentStride;

  segmentLength  = floor( segmentDuration/series->deltaT + 0.5 );
  segmentStride  = floor( strideDuration/series->deltaT + 0.5 );

  /* allocate memory for the data */
  segment->data = XLALCreateCOMPLEX8Vector( segmentLength/2 + 1 );

  /* create a time series that contains only the relevant data */
  seg        = *series;
  vec.length = segmentLength;
  vec.data   = seg.data->data + segmentNumber * segmentStride;
  seg.data   = &vec;
  ns  = epoch_to_ns( &seg.epoch );
  ns += sec_to_ns( segmentNumber * strideDuration );
  ns_to_epoch( &seg.epoch, ns );

  /* fft the data */
  XLALREAL4TimeFreqFFT( segment, &seg, fwdPlan );

  /* multiply by the inverse spectrum */
  XLALSCVectorMultiply( segment->data, invspec->data, segment->data );
  XLALUnitMultiply( &segment->sampleUnits, &segment->sampleUnits,
      &invspec->sampleUnits );

  return 0;
}

/* computes the inverse power spectrum */
static real_vector_t *get_invspec(
    real_vector_t         *channel,
    REAL4FFTPlan          *fwdplan,
    REAL4FFTPlan          *revplan,
    double                low_cutoff_frequency,
    double                segment_duration,
    double                stride_duration
    double                truncate_duration
    )
{
  real_vector_t *invspec = NULL;

  invspec = compute_average_spectrum( channel, segment_duration,
				      stride_duration, fwdplan);   
  
  /* invert spectrum */
  invert_spectrum( invspec, sample_rate, stride_duration,
		   truncate_duration, low_cutoff_frequency, fwdplan,
		   revplan );

  return invspec;
}

