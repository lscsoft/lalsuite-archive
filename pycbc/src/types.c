#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <cufft.h>
#include <lal/LALDatatypes.h>
#include "pycbc.h"

/* real vector manipulation */
real_vector_t* new_real_vector_t( 
    int length, enum cbc_memory_meta_types_t memory_location )
{
  real_vector_t* c = (real_vector_t*) malloc( sizeof(real_vector_t) );

  c->t_start = 0;
  c->dx = 1;
  c->vector_length = length;
  c->owner = 1;
  c->memory_type = memory_location;

  if ( memory_location == gpu_cuda_global_memory )
  {
    c->element_size_bytes = sizeof(cufftReal);
    cudaMalloc( (void**) &(c->data), length * c->element_size_bytes );
  }
  else if ( memory_location == cpu_generic_memory )
  {
    c->element_size_bytes = sizeof(REAL4);
    c->data = malloc( length * c->element_size_bytes );
  }

#ifdef PYCBC_MEM_DEBUG
  fprintf( stderr, "created real_vector_t at %p\n", c );
#endif
  return c;
}
  
void delete_real_vector_t( real_vector_t* p )
{
#ifdef PYCBC_MEM_DEBUG
  fprintf( stderr, "deleting real_vector_t at %p\n", p );
#endif
  if ( p->owner )
  {
    if ( p->memory_type == gpu_cuda_global_memory )
    {
      cudaFree( p->data );
    }
    else if ( p->memory_type == cpu_generic_memory )
    {
      free( p->data );
    }
  }
  free( p );
}

/* complex vector manipulation */
complex_vector_t* new_complex_vector_t( 
    int length, enum cbc_memory_meta_types_t memory_location )
{
  complex_vector_t* c = (complex_vector_t*) malloc( sizeof(complex_vector_t) );

  c->t_start = 0;
  c->dx = 1;
  c->vector_length = length;
  c->owner = 1;
  c->memory_type = memory_location;

  if ( memory_location == gpu_cuda_global_memory )
  {
    c->element_size_bytes = sizeof(cufftComplex);
    cudaMalloc( (void**) &(c->data), length * c->element_size_bytes );
  }
  else if ( memory_location == cpu_generic_memory )
  {
    c->element_size_bytes = sizeof(COMPLEX8);
    c->data = malloc( length * c->element_size_bytes );
  }

#ifdef PYCBC_MEM_DEBUG
  fprintf( stderr, "created complex_vector_t at %p\n", c );
#endif
  return c;
}
  
void delete_complex_vector_t( complex_vector_t* p )
{
#ifdef PYCBC_MEM_DEBUG
  fprintf( stderr, "deleting complex_vector_t at %p\n", p );
#endif
  if ( p->owner )
  {
    if ( p->memory_type == gpu_cuda_global_memory )
    {
      cudaFree( p->data );
    }
    else if ( p->memory_type == cpu_generic_memory )
    {
      free( p->data );
    }
  }
  free( p );
}


/* get a pycbc time series from a lal real4 time series */
real_vector_t *convert_timeseries_lal2pycbc ( REAL4TimeSeries *lal_timeseries)
{
  real_vector_t *pycbc_timeseries = NULL;

  if (lal_timeseries == NULL)
    return pycbc_timeseries;

  pycbc_timeseries = (real_vector_t *) malloc( sizeof(*pycbc_timeseries) );

  pycbc_timeseries->t_start = lal_timeseries->epoch.gpsSeconds * XLAL_BILLION_INT8 + lal_timeseries->epoch.gpsNanoSeconds;
  pycbc_timeseries->dx = lal_timeseries->deltaT;
  pycbc_timeseries->vector_length = lal_timeseries->data->length;

  pycbc_timeseries->element_size_bytes = sizeof(REAL4);
  pycbc_timeseries->memory_type = cpu_generic_memory;
  pycbc_timeseries->owner = 0; /* data owned by lal type */

  pycbc_timeseries->data = lal_timeseries->data->data;
  
  return pycbc_timeseries;
}

/* get a pycbc freq series from a lal real4 freq series */
real_vector_t *convert_freqseries_lal2pycbc ( REAL4FrequencySeries *lal_freqseries)
{
  real_vector_t *pycbc_freqseries = NULL;

  if (lal_freqseries == NULL)
    return pycbc_freqseries;

  pycbc_freqseries = (real_vector_t *) malloc( sizeof(*pycbc_freqseries) );

  pycbc_freqseries->t_start = lal_freqseries->epoch.gpsSeconds * XLAL_BILLION_INT8 + lal_freqseries->epoch.gpsNanoSeconds;
  pycbc_freqseries->dx = lal_freqseries->deltaF;
  pycbc_freqseries->vector_length = lal_freqseries->data->length;

  pycbc_freqseries->element_size_bytes = sizeof(REAL4);
  pycbc_freqseries->memory_type = cpu_generic_memory;
  pycbc_freqseries->owner = 0; /* data owned by lal type */

  pycbc_freqseries->data = lal_freqseries->data->data;
  
  return pycbc_freqseries;
}


/* get a lal time series from a pycbc vector */
REAL4TimeSeries *convert_timeseries_pycbc2lal ( real_vector_t *pycbc_timeseries)
{
  REAL4TimeSeries *lal_timeseries = NULL;
  ldiv_t t_start_sec;

  if (pycbc_timeseries == NULL) 
    return lal_timeseries;

  lal_timeseries = (REAL4TimeSeries *)malloc( sizeof(*timeseries) );

  /* use glibc for long integer division */
  /* make this into a function or use lal?*/
  t_start_sec = ldiv( pycbc_timeseries->t_start, XLAL_BILLION_INT8);
  lal_timeseries->epoch.gpsSeconds = t_start_sec.quot;
  lal_timeseries->epoch.gpsNanoSeconds = t_start_sec.rem;

  lal_timeseries->deltaT = pycbc_timeseries->dx;
  lal_timeseries->f0 = 0; /* not heterodyned */
  lal_timeseries->data->length = pycbc_timeseries->vector_length;
  lal_timeseries->data->data = pycbc_timeseries->data;
  
  return lal_timeseries;
}

/* get a lal freq series from a pycbc vector */
REAL4FrequencySeries *convert_freqseries_pycbc2lal ( real_vector_t *pycbc_freqseries)
{
  REAL4FrequencySeries *lal_freqseries = NULL;
  ldiv_t t_start_sec;

  if (pycbc_freqseries == NULL) 
    return lal_freqseries;

  lal_freqseries = (REAL4FrequencySeries *)malloc( sizeof(*freqseries) );

  /* use glibc for long integer division */
  /* make this into a function or use lal?*/
  t_start_sec = ldiv( pycbc_freqseries->t_start, XLAL_BILLION_INT8);
  lal_freqseries->epoch.gpsSeconds = t_start_sec.quot;
  lal_freqseries->epoch.gpsNanoSeconds = t_start_sec.rem;

  lal_freqseries->deltaT = pycbc_freqseries->dx;
  lal_freqseries->data->length = pycbc_freqseries->vector_length;
  lal_freqseries->data->length = pycbc_freqseries->vector_length;
  lal_freqseries->data->data = pycbc_freqseries->data;
  
  return lal_freqseries;
}
