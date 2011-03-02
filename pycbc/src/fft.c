#include <pycbc.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <cufft.h>
#include <fftw3.h>

complex_fft_plan_t* new_complex_fft_plan_t( int size, int fwdflg, int level, 
    enum cbc_memory_meta_types_t memory_location )
{
  complex_fft_plan_t* plan = 
    (complex_fft_plan_t*) malloc( sizeof( complex_fft_plan_t ) );

  /* create the plan */
  if ( memory_location == gpu_cuda_global_memory )
  {
    cufftHandle* cuplan_ptr = (cufftHandle*) malloc( sizeof(cufftHandle) );
    cufftPlan1d( cuplan_ptr, size, CUFFT_C2C, 1 );
    plan->plan = cuplan_ptr;
    plan->memory_type = gpu_cuda_global_memory;
  }
  else if ( memory_location == cpu_generic_memory )
  {
    float* tmp1;
    float* tmp2;
    int flags = FFTW_UNALIGNED;

    switch ( level )
    {
      case 0: /* estimate */
        flags |= FFTW_ESTIMATE;
        break;
      default: /* exhaustive measurement */
        flags |= FFTW_EXHAUSTIVE;
        /* fall-through */
      case 2: /* lengthy measurement */
        flags |= FFTW_PATIENT;
        /* fall-through */
      case 1: /* measure the best plan */
        flags |= FFTW_MEASURE;
        break;
    }

    /* alloc plan an tmp arrays */
    tmp1 = malloc( size * sizeof( *tmp1 ) );
    tmp2 = malloc( size * sizeof( *tmp2 ) );

    if ( fwdflg ) /* forward */
      plan->plan = fftwf_plan_r2r_1d( size, tmp1, tmp2, FFTW_R2HC, flags );
    else /* reverse */
      plan->plan = fftwf_plan_r2r_1d( size, tmp1, tmp2, FFTW_HC2R, flags );

    plan->memory_type = cpu_generic_memory;

    /* free tmp arrays */
    free( tmp2 );
    free( tmp1 );
  }

  /* set remaining plan fields */
  plan->size = size;
  plan->sign = ( fwdflg ? -1 : 1 );

  return plan;
}

void delete_complex_fft_plan_t( complex_fft_plan_t* plan )
{
  if ( plan->memory_type == gpu_cuda_global_memory )
    cufftDestroy( plan->plan );
  else if ( plan->memory_type == gpu_cuda_global_memory )
    fftw_destroy_plan( plan->plan );

  free( plan );
}

#if 0
int complex_vector_fft( complex_vector_t* output, complex_vector_t* input, 
    complex_fft_plan_t* plan )
{
  cufftHandle cuda_plan = *( (cufftHandle*) plan->plan );

  if ( plan->sign = -1 )
  {
    cufftExecC2C( cuda_plan, input, output, CUFFT_FORWARD );
  }
  else
  {
    cufftExecC2C( cuda_plan, input, output, CUFFT_INVERSE );
  }

  /* vector step */
  output->dx = 1.0 / (input->length * input->dx);

  return 0;
}
#endif
