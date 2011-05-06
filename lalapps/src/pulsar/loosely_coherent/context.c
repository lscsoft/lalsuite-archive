#include <stdio.h>
#include <math.h>
#include "global.h"
#include "context.h"
#include "util.h"
#include "cmdline.h"
#include "dataset.h"

#define restrict

#include <lal/LALDatatypes.h>
#include <lal/AVFactories.h>
#include <lal/ComplexFFT.h>

#define XALLOC(var, call)	while((var=call)==NULL) { \
	fprintf(stderr, "*** Could not allocate " #var "=" #call "\n"); \
	condor_safe_sleep(10); \
	}

struct gengetopt_args_info args_info;

LOOSE_CONTEXT * create_context(void)
{
LOOSE_CONTEXT *ctx;
int wing_step;
int day_samples=round(2.0*SIDEREAL_DAY/args_info.coherence_length_arg);

ctx=do_alloc(1, sizeof(*ctx));
	
ctx->nsamples=1+ceil(2.0*(max_gps()-min_gps())/args_info.coherence_length_arg);
wing_step=round(ctx->nsamples*args_info.coherence_length_arg/SIDEREAL_DAY);
ctx->nsamples=day_samples*floor(ctx->nsamples/day_samples);
ctx->nsamples=round235up_int(ctx->nsamples);

TODO("increase plan optimization level to 1")
fprintf(stderr, "Creating FFT plan of length %d\n", ctx->nsamples);
XALLOC(ctx->fft_plan, XLALCreateForwardCOMPLEX16FFTPlan(ctx->nsamples, 0));

XALLOC(ctx->plus_samples, XLALCreateCOMPLEX16Vector(ctx->nsamples));
XALLOC(ctx->cross_samples, XLALCreateCOMPLEX16Vector(ctx->nsamples));
XALLOC(ctx->plus_fft, XLALCreateCOMPLEX16Vector(ctx->nsamples));
XALLOC(ctx->cross_fft, XLALCreateCOMPLEX16Vector(ctx->nsamples));
XALLOC(ctx->plus_te_fft, XLALCreateCOMPLEX16Vector(ctx->nsamples));
XALLOC(ctx->cross_te_fft, XLALCreateCOMPLEX16Vector(ctx->nsamples));

return(ctx);
}