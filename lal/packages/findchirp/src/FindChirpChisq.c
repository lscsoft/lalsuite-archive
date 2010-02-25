/*
*  Copyright (C) 2007 Stas Babak, Duncan Brown, Eirini Messaritaki, Jolien Creighton, Reinhard Prix, Craig Robinson
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
 * File Name: FindChirpChisq.c
 *
 * Author: Anderson, W. G., and Brown, D. A.
 *
 * Revision: $Id$
 *
 *-----------------------------------------------------------------------
 */

#if 0
<lalVerbatim file="FindChirpChisqCV">
Author: Anderson, W. G., and Brown D. A.
$Id$
</lalVerbatim>

<lalLaTeX>
\subsection{Module \texttt{FindChirpChisq.c}}
\label{ss:FindChirpChisq.c}

Module to implement the $\chi^2$ veto for the stationary phase chirp.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{FindChirpChisqCP}
\idx{LALFindChirpChisqVeto()}

\subsubsection*{Description}

The function \texttt{LALFindChirpChisqVeto()} perfoms a $\chi^2$ veto on
an entire data segment using the algorithm described below. On exit the
vector \texttt{chisqVec} contains the value $\chi^2(t_j)$ for the data
segment.

\subsubsection*{Algorithm}

chisq algorithm here

\subsubsection*{Uses}
\begin{verbatim}
LALCreateReverseComplexFFTPlan()
LALDestroyComplexFFTPlan()
LALCCreateVector()
LALCDestroyVector()
LALCOMPLEX8VectorFFT()
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{FindChirpChisqCV}}
</lalLaTeX>
#endif

#include <stdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/AVFactories.h>
#include <lal/ComplexFFT.h>
#include <lal/FindChirp.h>
#include <lal/FindChirpChisq.h>
#include <lal/FindChirpPTF.h>
#include <lal/MatrixUtils.h>

#ifndef isnan
int isnan(double);
#define isnan(x) ((isnan)((double)(x)))
#endif

NRCSID (FINDCHIRPCHISQC, "$Id$");


/* <lalVerbatim file="FindChirpChisqCP"> */
void
LALFindChirpComputeChisqBins(
    LALStatus                  *status,
    UINT4Vector                *chisqBinVec,
    FindChirpSegment           *fcSeg,
    UINT4                       kmax
    )
/* </lalVerbatim> */
{
  UINT4         k, incIdx;
  REAL4        *tmpltPower;
  UINT4        *chisqBin = NULL;
  UINT4         numChisqBins;
  UINT4         chisqPt;
  REAL4         increment;
  REAL4         nextBin;
  REAL4         partSum;

  INITSTATUS( status, "LALFindChirpComputeChisqBins", FINDCHIRPCHISQC );
  ATTATCHSTATUSPTR( status );


  ASSERT( chisqBinVec, status,
      FINDCHIRPCHISQH_ENULL, FINDCHIRPCHISQH_MSGENULL );
  ASSERT( chisqBinVec->length > 1, status,
      FINDCHIRPCHISQH_ECHIZ, FINDCHIRPCHISQH_MSGECHIZ );
  ASSERT( ! chisqBinVec->data, status,
      FINDCHIRPCHISQH_ENNUL, FINDCHIRPCHISQH_MSGENNUL );

  ASSERT( fcSeg, status,
      FINDCHIRPCHISQH_ENULL, FINDCHIRPCHISQH_MSGENULL );
  ASSERT( fcSeg->segNorm, status,
      FINDCHIRPCHISQH_ENULL, FINDCHIRPCHISQH_MSGENULL );
  ASSERT( fcSeg->segNorm->data, status,
      FINDCHIRPCHISQH_ENULL, FINDCHIRPCHISQH_MSGENULL );
  ASSERT( fcSeg->segNorm->length, status,
      FINDCHIRPCHISQH_ENUMZ, FINDCHIRPCHISQH_MSGENUMZ );
  ASSERT( fcSeg->tmpltPowerVec, status,
      FINDCHIRPCHISQH_ENULL, FINDCHIRPCHISQH_MSGENULL );
  ASSERT( fcSeg->tmpltPowerVec->length, status,
      FINDCHIRPCHISQH_ENUMZ, FINDCHIRPCHISQH_MSGENUMZ );
  ASSERT( fcSeg->tmpltPowerVec->data, status,
      FINDCHIRPCHISQH_ENULL, FINDCHIRPCHISQH_MSGENULL );


  /*
   *
   * calculate the chisq bins for the segment and template
   *
   */


  /* number of chisq bins is one less than the number of bin boundaries */
  numChisqBins = chisqBinVec->length - 1;

  /* vector containing k^(-7/3) / S_h(f) */
  tmpltPower = fcSeg->tmpltPowerVec->data;

  /* compute the amount of power in a chisq bin */
  incIdx = kmax > fcSeg->segNorm->length-1 ? fcSeg->segNorm->length-1 : kmax;
  increment = fcSeg->segNorm->data[incIdx] / (REAL4) numChisqBins;

  /* allocate memory for the bin boundaries */
  chisqBin = chisqBinVec->data = (UINT4 *)
    LALCalloc( chisqBinVec->length, sizeof(UINT4) );
  if ( ! chisqBinVec->data )
  {
    ABORT( status, FINDCHIRPCHISQH_EALOC, FINDCHIRPCHISQH_MSGEALOC );
  }

  /* initalize the first bin boundary */
  nextBin   = increment;
  chisqPt   = 0;
  partSum   = 0.0;
  chisqBin[chisqPt++] = 0;

  /* calculate the frequencies of the chi-squared bin boundaries */
  for ( k = 1; k < incIdx; ++k )
  {
    partSum += tmpltPower[k];
    if ( partSum >= nextBin )
    {
      chisqBin[chisqPt++] = k;
      nextBin += increment;
      if ( chisqPt == numChisqBins ) break;
    }
  }

  /* check that we have sucessfully allocated all the bins */
  if ( k == fcSeg->tmpltPowerVec->length && chisqPt != numChisqBins )
  {
    /* if we have reaced the end of the template power vec and not */
    /* allocated all the bin boundaries then there is a problem    */
    ABORT( status, FINDCHIRPCHISQH_EBINS, FINDCHIRPCHISQH_MSGEBINS );
  }

  /* the last bin boundary is at can be at Nyquist since   */
  /* qtilde is zero above the ISCO of the current template */
  chisqBin[numChisqBins] = fcSeg->data->data->length;


  /* normal exit */
  DETATCHSTATUSPTR( status );
  RETURN( status );
}


/* <lalVerbatim file="FindChirpChisqCP"> */
void
LALFindChirpChisqVeto (
    LALStatus                  *status,
    REAL4Vector                *chisqVec,
    FindChirpChisqInput        *input,
    FindChirpChisqParams       *params
    )
/* </lalVerbatim> */
{
  UINT4                 i, j, l, m;
  UINT4                 numPoints, len;

  REAL4                *chisq;

  COMPLEX8             *q;
  COMPLEX8             *qtilde;

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

  INITSTATUS( status, "LALFindChirpChisqVeto", FINDCHIRPCHISQC );
  ATTATCHSTATUSPTR( status );


  /*
   *
   * check that the arguments are reasonable
   *
   */


  /* check that the output pointer is non-null and has room to store data */
  ASSERT( chisqVec, status,
      FINDCHIRPCHISQH_ENULL, FINDCHIRPCHISQH_MSGENULL );
  ASSERT( chisqVec->data, status,
      FINDCHIRPCHISQH_ENULL, FINDCHIRPCHISQH_MSGENULL );

  /* check that the parameter structure exists */
  ASSERT( params, status, FINDCHIRPCHISQH_ENULL, FINDCHIRPCHISQH_MSGENULL );

  /* check that the chisq bin vector is reasonable */
  ASSERT( params->chisqBinVec, status,
      FINDCHIRPCHISQH_ENULL, FINDCHIRPCHISQH_MSGENULL );
  ASSERT( params->chisqBinVec->data, status,
      FINDCHIRPCHISQH_ENULL, FINDCHIRPCHISQH_MSGENULL );
  ASSERT( params->chisqBinVec->length > 0, status,
      FINDCHIRPCHISQH_ECHIZ, FINDCHIRPCHISQH_MSGECHIZ );

  /* check that the fft plan exists */
  ASSERT( params->plan, status,
      FINDCHIRPCHISQH_ENULL, FINDCHIRPCHISQH_MSGENULL );

  /* check that the input exists */
  ASSERT( input, status, FINDCHIRPCHISQH_ENULL, FINDCHIRPCHISQH_MSGENULL );



  /* check that we are using the correct approximant */
  switch ( params->approximant )
  {
    case TaylorT1:
    case TaylorT2:
    case TaylorT3:
    case TaylorF2:
    case GeneratePPN:
    case PadeT1:
    case EOB:
    case EOBNR:
    case FindChirpSP:
      /* check that the input contains some data */
      ASSERT( input->qVec, status,
          FINDCHIRPCHISQH_ENULL, FINDCHIRPCHISQH_MSGENULL );
      ASSERT( input->qVec->data, status,
          FINDCHIRPCHISQH_ENULL, FINDCHIRPCHISQH_MSGENULL );
      ASSERT( input->qtildeVec, status,
          FINDCHIRPCHISQH_ENULL, FINDCHIRPCHISQH_MSGENULL );
      ASSERT( input->qtildeVec->data, status,
          FINDCHIRPCHISQH_ENULL, FINDCHIRPCHISQH_MSGENULL );

      /* check that the workspace vectors exist */
      ASSERT( params->qtildeBinVec, status,
          FINDCHIRPCHISQH_ENULL, FINDCHIRPCHISQH_MSGENULL );
      ASSERT( params->qtildeBinVec->data, status,
          FINDCHIRPCHISQH_ENULL, FINDCHIRPCHISQH_MSGENULL );
      ASSERT( params->qtildeBinVec->length > 0, status,
          FINDCHIRPCHISQH_ECHIZ, FINDCHIRPCHISQH_MSGECHIZ );
      /*
       *
       * point local pointers to structure pointers
       *
       */


      chisq     = chisqVec->data;

      numPoints = input->qVec->length;
      q         = input->qVec->data;
      qtilde    = input->qtildeVec->data;

      numChisqPts  = params->chisqBinVec->length;
      numChisqBins = numChisqPts - 1;
      chisqBin     = params->chisqBinVec->data;
      chisqNorm    = sqrt( params->norm );

      qtildeBin = params->qtildeBinVec->data;


      /*
       *
       * fill the numBins time series vectors for the chi-squared statistic
       *
       */


      for ( l = 0; l < numChisqBins; ++l )
      {
        memset( qtildeBin, 0, numPoints * sizeof(COMPLEX8) );

        memcpy( qtildeBin + chisqBin[l], qtilde + chisqBin[l],
            (chisqBin[l+1] - chisqBin[l]) * sizeof(COMPLEX8) );

        LALCOMPLEX8VectorFFT( status->statusPtr, params->qBinVecPtr[l],
            params->qtildeBinVec, params->plan );
        CHECKSTATUSPTR( status );
      }


      /*
       *
       * calculate the chi-squared value at each time
       *
       */


      memset( chisq, 0, numPoints * sizeof(REAL4) );

      for ( j = 0; j < numPoints; ++j )
      {
        for ( l = 0; l < numChisqBins; ++l )
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
      break;

    case FindChirpPTF:

      /* check that the input contains some data */
      ASSERT( input->PTFqVec, status, 
          FINDCHIRPCHISQH_ENULL, FINDCHIRPCHISQH_MSGENULL );
      ASSERT( input->PTFqVec->data, status, 
          FINDCHIRPCHISQH_ENULL, FINDCHIRPCHISQH_MSGENULL );
      ASSERT( input->PTFqtildeVec, status, 
          FINDCHIRPCHISQH_ENULL, FINDCHIRPCHISQH_MSGENULL );
      ASSERT( input->PTFqtildeVec->data, status, 
          FINDCHIRPCHISQH_ENULL, FINDCHIRPCHISQH_MSGENULL );
      ASSERT( input->PTFPVec, status, 
          FINDCHIRPCHISQH_ENULL, FINDCHIRPCHISQH_MSGENULL );
      ASSERT( input->PTFPVec->data, status, 
          FINDCHIRPCHISQH_ENULL, FINDCHIRPCHISQH_MSGENULL );
      ASSERT( input->PTFsegNormVec, status, 
          FINDCHIRPCHISQH_ENULL, FINDCHIRPCHISQH_MSGENULL );
      ASSERT( input->PTFsegNormVec->data, status, 
          FINDCHIRPCHISQH_ENULL, FINDCHIRPCHISQH_MSGENULL );

      /* check that the workspace vectors exist */

      ASSERT( params->PTFB, status, 
          FINDCHIRPCHISQH_ENULL, FINDCHIRPCHISQH_MSGENULL );
      ASSERT( params->PTFB->data, status, 
          FINDCHIRPCHISQH_ENULL, FINDCHIRPCHISQH_MSGENULL );

      /*
       *
       * point local pointers to structure pointers
       *
       */


      chisq        = chisqVec->data;
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

      memset( chisqVec->data, 0, numPoints * sizeof(REAL4) );
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
      break;

    default:
      ABORT( status, FINDCHIRPCHISQH_EIAPX, FINDCHIRPCHISQH_MSGEIAPX );
      break;
  }




  /* normal exit */
  DETATCHSTATUSPTR( status );
  RETURN( status );
}
