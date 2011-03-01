# include <math.h>
# include <stdio.h>
# include <stdlib.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>
#include <lal/LALGSL.h>

# include <lal/LALStdlib.h>
# include <lal/LALConstants.h>
# include <lal/SimulateCoherentGW.h>
# include <lal/GeneratePPNInspiral.h>
# include <lal/LIGOMetadataTables.h>
# include <lal/LALDatatypes.h>
# include <lal/LALComplex.h>


#ifdef  __cplusplus
extern "C" {
#endif


void LALInspiralTaylorF2Consistency(
   LALStatus        *status,
   REAL4Vector      *signalvec,
   InspiralTemplate *params
   );