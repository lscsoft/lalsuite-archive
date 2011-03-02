#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <cuda.h>
#include <cufft.h>
#include <lal/LALDatatypes.h>

#ifndef PYCBC_H
#define PYCBC_H

enum cbc_memory_meta_types_t {
  gpu_cuda_global_memory,
  gpu_cuda_constant_memory,
  gpu_cuda_texture_memory,
  gpu_opencl_global_memory,
  gpu_opencl_constant_memory,
  gpu_opencl_local_memory,
  gpu_opencl_private_memory,
  cpugpu_cuda_zero_latency_memory,
  cpu_generic_memory,
  cpu_pinned_memory,
  cpu_fftw_aligned_memory,
  non_specified_memory
};

typedef struct
{
  unsigned long int t_start;
  double dx;
  unsigned int vector_length;
  size_t element_size_bytes;
  enum cbc_memory_meta_types_t memory_type;
  int owner;
  void *data;
}
real_vector_t;

typedef real_vector_t complex_vector_t;

typedef struct
{
  int sign;
  int size;
  enum cbc_memory_meta_types_t memory_type;
  void *plan;
}
real_fft_plan_t;

typedef real_fft_plan_t complex_fft_plan_t;

/* Note about python memory usage:                              */
/* If a function returns allocated memory, it must have a       */
/* corresponding %newtype swig declaration, so that the memory  */
/* can be garbage collected by the python interpreter. Certain  */
/* functions may or may not allocate memory depending on the    */
/* command line arguments. In this case there should be an      */
/* output typemap that provides the correct garbage collection  */
/* information for swig. The default typemaps are listed below  */
/* for the most common operations. These can be applied to a    */
/* particular function with the %apply directive                */

/* fft functions */
int forward_real_fft( complex_vector_t* output, real_vector_t* input, 
    real_fft_plan_t* plan );
int reverse_real_fft( real_vector_t* output, complex_vector_t* input, 
    real_fft_plan_t* plan );
int complex_vec_fft( complex_vector_t* output, complex_vector_t* input, 
    complex_fft_plan_t* plan );

/* demonstration of a swig-wrapped psd function */
#ifdef SWIG
%newobject pwelch;
#endif
real_vector_t* welch_power_spectrum( 
    real_vector_t* self, int length, int overlap );

#endif /* PYCBC_H */
