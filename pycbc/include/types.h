#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pycbc.h>

#ifndef PYCBC_TYPES_H
#define PYCBC_TYPES_H

real_vector_t* new_real_vector_t( int length, 
    enum cbc_memory_meta_types_t memory_location );
void delete_real_vector_t( real_vector_t *p );

complex_vector_t* new_complex_vector_t( int length, 
    enum cbc_memory_meta_types_t memory_location );
void delete_complex_vector_t( complex_vector_t *p );

/* lal <-> pycbc conversion functions */
real_vector_t *convert_timeseries_lal2pycbc ( REAL4TimeSeries *lal_timeseries);

real_vector_t *convert_freqseries_lal2pycbc ( REAL4FrequencySeries *lal_freqseries);

REAL4TimeSeries *convert_timeseries_pycbc2lal ( real_vector_t *pycbc_timeseries);

REAL4FrequencySeries *convert_freqseries_pycbc2lal ( real_vector_t *pycbc_freqseries);

#endif /* PYCBC_TYPES_H */
