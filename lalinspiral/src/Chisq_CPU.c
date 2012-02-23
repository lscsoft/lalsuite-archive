/*
*  Copyright (C) 2010 Karsten Wiesner, Stas Babak, Duncan Brown, Eirini Messaritaki, Jolien Creighton, Reinhard Prix, Craig Robinson
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with with program; see the file COPYING. If not, write to the
*  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
*  MA  02111-1307  USA
*/

/*-----------------------------------------------------------------------
 *
 * File Name: Chisq_CPU.c
 *
 * Author: Wiesner, K., Anderson, W. G., and Brown, D. A.
 *
 * Revision: $Id$
 *
 *-----------------------------------------------------------------------
 */
#include <stdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/AVFactories.h>
#include <lal/ComplexFFT.h>
#include <lal/FindChirp.h>
#include <lal/FindChirpChisq.h>
#include "Chisq_CPU.h"

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

void Chisq_CPU (REAL4* chisq, COMPLEX8* q, COMPLEX8* qtilde, FindChirpChisqParams* params, 
    UINT4 numPoints, UINT4 numChisqBins, REAL4 chisqNorm, LALStatus UNUSED *status)
{
  
  UINT4* chisqBin     = params->chisqBinVec->data;
  COMPLEX8* qtildeBin = params->qtildeBinVec->data;
  int rc;

  for ( UINT4 l = 0; l < numChisqBins; ++l )
    {
      memset( qtildeBin, 0, numPoints * sizeof(COMPLEX8) );
      
      memcpy( qtildeBin + chisqBin[l], qtilde + chisqBin[l],
	      (chisqBin[l+1] - chisqBin[l]) * sizeof(COMPLEX8) );
      
      rc = XLALCOMPLEX8VectorFFT(params->qBinVecPtr[l], \
          params->qtildeBinVec, params->plan);
      if (rc != 0)
        XLAL_ERROR_VOID(rc);
    }
  
  memset( chisq, 0, numPoints * sizeof(REAL4) );
  
  for ( UINT4 j = 0; j < numPoints; ++j )
    {
      for ( UINT4 l = 0; l < numChisqBins; ++l )
	{
	  REAL4 Xl = params->qBinVecPtr[l]->data[j].re;
	  REAL4 Yl = params->qBinVecPtr[l]->data[j].im;
	  REAL4 deltaXl = chisqNorm * Xl -
	    (chisqNorm * q[j].re / (REAL4) (numChisqBins));
	  REAL4 deltaYl = chisqNorm * Yl -
	    (chisqNorm * q[j].im / (REAL4) (numChisqBins));
	  
	  chisq[j] += deltaXl * deltaXl + deltaYl * deltaYl;
	}
    }
}
  
void Chisq_CPU_PTF (REAL4* chisq, FindChirpChisqInput* input, FindChirpChisqParams* params, 
		LALStatus UNUSED *status)
{
	UINT4                 i, j, l, m;
	UINT4                 numPoints, len;

	UINT4                 numChisqPts;
	UINT4                 numChisqBins;
	UINT4                *chisqBin;
	REAL4                 chisqNorm;

	COMPLEX8             *qtildeBin;

	/* PTF workspace */
	UINT4                 bin_lo, bin_hi;
	REAL8                 norm_fac, bin_norm, hc_l, hs_l, delta_hc, delta_hs;
	REAL8                 v1[5], v2[5];
	REAL4                *hc;
	REAL4                *hs;
	REAL4                *PTFB;
	REAL4                *PTFsegNorm;
	REAL4                *PTFP;
	COMPLEX8             *PTFq;
        COMPLEX8             *PTFqtilde;

	/*
	 *
	 * point local pointers to structure pointers
	 *
	 */

	numPoints    = input->qVec->length;
	PTFq         = input->PTFqVec->data;
	PTFqtilde    = input->PTFqtildeVec->data;
	PTFP         = input->PTFPVec->data;
	PTFsegNorm   = input->PTFsegNormVec->data;

	PTFB         = params->PTFB->data;
	hc           = params->PTFhc->data;
	hs           = params->PTFhs->data;
	qtildeBin    = params->qtildeBinVec->data;
	numChisqPts  = params->chisqBinVec->length;
	numChisqBins = numChisqPts - 1;
	chisqBin     = params->chisqBinVec->data;
	chisqNorm    = sqrt( params->norm );
	len          = numPoints / 2 + 1;

	/* 
	 *
	 * fill the numBins time series vectors for the chi-squared statistic
	 *
	 */


	/* Construct B_ij for every bin */

	memset( params->PTFB->data, 0, numChisqBins * 25 * sizeof(REAL4) );
	for ( l = 0; l < numChisqBins; ++l ) 
	{
		bin_lo = chisqBin[l];
		bin_hi = chisqBin[l+1] > (UINT4) (input->kmax -1) ? (UINT4) (input->kmax -1) : 
			(UINT4) chisqBin[l+1];
		if ( bin_hi < bin_lo )
		{
			fprintf(stderr, "Bin boundaries error for PTF\n");
			ABORT( status, FINDCHIRPH_ECRUP, FINDCHIRPH_MSGECRUP );
		}
		for (i = 0; i < 5; i++)
		{ 
			for ( m = 0; m < 5; m++ )
			{  
				PTFB[l * 25 + i * 5 + m] = PTFsegNorm[(i + 5 * m) * len + bin_hi] - 
					PTFsegNorm[(i + 5 * m) * len + bin_lo]; 
				PTFB[l * 25 + i * 5 + m] *= 4.0 * input->deltaF;
			}
		}
	} 

	/* 
	 *
	 * calculate the chi-squared value at each time
	 *
	 */

	memset( params->PTFhc->data, 0, numPoints * sizeof(REAL4) );
	memset( params->PTFhs->data, 0, numPoints * sizeof(REAL4) );

	for ( j = 0; j < numPoints; j++ )
	{ 
		for ( i =0; i < 5; i++ )
		{
			hc[j] += PTFP[i * numPoints + j] * PTFq[i * numPoints + j].re;
			hs[j] += PTFP[i * numPoints + j] * PTFq[i * numPoints + j].im;
		}
	}

	memset( chisq, 0, numPoints * sizeof(REAL4) );
	norm_fac = 4.0 / ((REAL4) (numPoints) * (REAL4) (numPoints));

	for ( l = 0; l < numChisqBins; ++l )
	{  
		for (i = 0; i < 5; i++)
		{ 
			memset( params->qtildeBinVec->data, 0, numPoints * sizeof(COMPLEX8) );
			memcpy( qtildeBin + chisqBin[l], PTFqtilde + i * len + 
					chisqBin[l], (chisqBin[l+1] - chisqBin[l]) * sizeof(COMPLEX8) );

			LALCOMPLEX8VectorFFT( status->statusPtr, params->PTFqBinVecPtr[i],
					params->qtildeBinVec, params->plan );
			CHECKSTATUSPTR( status );
		}
		for ( j = 0; j < numPoints; j++ ) /* beginning of main loop over time */
		{ 

			bin_norm = 0.0;
			/* compute normalization factors c_l*/
			for (i = 0; i < 5; i++)
			{ 
				for ( m = 0; m < 5; m++ )
				{  
					bin_norm += PTFP[i * numPoints + j] * PTFP[m * numPoints + j] * 
						PTFB[ l * 25 + i * 5 + m];
				}
			}

			hc_l = hs_l = 0.0;
			/* Compute Chisq */ 
			for (i = 0; i < 5; i++)
			{ 
				v1[i] = params->PTFqBinVecPtr[i]->data[j].re;
				v2[i] = params->PTFqBinVecPtr[i]->data[j].im;
				hc_l += PTFP[i * numPoints + j] * v1[i];
				hs_l += PTFP[i * numPoints + j] * v2[i];
			}
			delta_hc = hc_l - bin_norm * (REAL8) hc[j];
			delta_hs = hs_l - bin_norm * (REAL8) hs[j];
			chisq[j] += (REAL4) (( delta_hc * delta_hc + delta_hs * delta_hs ) * norm_fac / bin_norm);
			if( isnan( chisq[j] ) ) fprintf(stderr,"Chisq value is nan for bin %d at point=%d\n",l,j);
		}
		/* fprintf(stderr,"%e\n",chisq[j]); */
	} /* End of main loop over time */
}
