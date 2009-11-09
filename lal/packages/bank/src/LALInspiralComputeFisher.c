/*
*  Copyright (C) 2009 Tjonnie Li, Chris Van Den Broeck
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

#if 0
<lalVerbatim file="LALInspiralComputeFisherCV">
Author: Tjonnie Li, Chris Van Den Broeck
$Id$
</lalVerbatim>

<lalLaTeX>
\subsection{Module \texttt{LALInspiralComputeFisher.c}}

Module to compute the Fisher matrix for a variety of time domain waveforms.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{XLALInspiralComputeFisherMatrixCP}
\idx{XLALInspiralComputeFisherMatrix()}
\begin{itemize}
   \item \texttt{fisher,} Output, the Fisher matrix at the point defined by \texttt{params}
   \item \texttt{psd,} Input, the power spectral density of the data
   \item \texttt{params,} Input, the parameters where Fisher must be computed.
\end{itemize}

\input{XLALInspiralComputeTDWDerivCP}
\idx{XLALInspiralComputeTDWDeriv()}
\begin{itemize}
   \item \texttt{Wderiv,} Output, the time derivative of waveform at the point
   point defined by \texttt{params}
   \item \texttt{psd,}  Input, the power spectral density of the data
   \item \texttt{params,} Input, the parameters where Fisher must be computed
   \item \texttt{paramid,} Input, id of the parameter to take derivative on
   \item \texttt{initdelta,} Input, initial difference in parameters
   \item \texttt{tolerance,} Input, stop iteration when difference between two
   bisections is smaller than tolerance.
\end{itemize}


\subsubsection*{Description}
We calculate the usual Fisher matrix.

\subsubsection*{Algorithm}

\subsubsection*{Uses}

\begin{verbatim}
LALMalloc
LALFree
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{LALInspiralComputeFisherCV}}

</lalLaTeX>
#endif

#include <stdlib.h>
#include <lal/LALStdlib.h>
#include <stdio.h>
#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>
#include <lal/LALInspiralBank.h>
#include <lal/LALNoiseModels.h>
#include <lal/MatrixUtils.h>
#include <lal/FindChirpPTF.h>


/* <lalVerbatim file="XLALInspiralComputeFisherMatrixCP">  */
INT4 XLALInspiralComputeFisherMatrix (
    REAL8Vector				   *fisher,
    REAL8FrequencySeries                   *psd,
    InspiralTemplate                       *params
    )
/* </lalVerbatim> */
{
  /* XLAL error handling */
  INT4 errcode = XLAL_SUCCESS;
  static const char* func = "XLALInspiralComputeFisherMatrix";

  /* number of points in a time-domain segment */
  UINT4 N = 2 * (psd->data->length - 1);
  UINT4 i, j, k;

  /* get a local copy of the extrinstic parameters */
  REAL8 phi0 = params->startPhase;

  /* bounds on the power spectrum integration for the moments */
  REAL8 deltaT = params->tSampling;

  /* local pointers to the Q and Qtilde data */

  REAL8VectorSequence			*PTFQ;
  COMPLEX16VectorSequence		*PTFQtilde;
  REAL4Vector				*PTFphi;
  REAL4Vector				*PTFomega_2_3;
  REAL4VectorSequence			*PTFe1;
  REAL4VectorSequence			*PTFe2;
  REAL8FFTPlan				*fwdPlan;
  REAL8FFTPlan				*revPlan;
  COMPLEX16Vector			*intrinsicderiv;
  COMPLEX16VectorSequence	        *derivs;
  REAL8		                         initdelta = 0.001;
  REAL8		                         tolerance = 0.001;
  REAL8Vector	                 	*invpsd;
  REAL8Vector	                	*tempnorm;
  COMPLEX16Vector               	*tempnormtilde;
  REAL8                                 *det;
  LALStatus                              status;

  PTFQ           = XLALCreateREAL8VectorSequence( 5, N );
  PTFQtilde      = XLALCreateCOMPLEX16VectorSequence( 5, N / 2 + 1 );
  PTFphi         = XLALCreateVector( N );
  PTFomega_2_3   = XLALCreateVector( N );
  PTFe1          = XLALCreateVectorSequence( 3, N );
  PTFe2          = XLALCreateVectorSequence( 3, N );
  fwdPlan        = XLALCreateForwardREAL8FFTPlan( N, 0 );
  revPlan        = XLALCreateReverseREAL8FFTPlan( N, 0 );
  det            = NULL;

  memset( &status, 0, sizeof(LALStatus) );

  /* call the PTF waveform code */
  errcode = XLALFindChirpPTFWaveform( PTFphi, PTFomega_2_3, PTFe1, PTFe2,
      params, deltaT);

  if ( errcode != XLAL_SUCCESS ) XLAL_ERROR( func, errcode );



  /* Waveform derivatives */
  intrinsicderiv = XLALCreateCOMPLEX16Vector(N / 2 + 1);
  derivs = XLALCreateCOMPLEX16VectorSequence(4, N / 2 + 1);


  /* Compute intrinsic derivatives */
  for (i = 0; i < 4; ++i)
  {
    errcode = XLALInspiralComputePTFWDeriv(intrinsicderiv, psd, params, i + 1, initdelta, tolerance);
    if ( errcode != XLAL_SUCCESS )
    {
      fprintf( stderr, "XLALInspiralComputePTFWDeriv failed\n" );
      exit( 1 );
    }
    for (k = 0; k < N / 2 + 1; ++k)
    {
      derivs->data[i * (N / 2 + 1) + k].re = cos(phi0) * intrinsicderiv->data[k].re
        + sin(phi0) * intrinsicderiv->data[k].im;
      derivs->data[i * (N / 2 + 1) + k].im = cos(phi0) * intrinsicderiv->data[k].im
        - sin(phi0) * intrinsicderiv->data[k].re;
    }
    derivs->data[i * (N / 2 + 1)].re = (cos(phi0) + sin(phi0)) * intrinsicderiv->data[0].re;
    derivs->data[i * (N / 2 + 1)].im = (cos(phi0) + sin(phi0)) * intrinsicderiv->data[0].im;
    derivs->data[i * (N / 2 + 1) + N / 2].re = (cos(phi0) + sin(phi0)) * intrinsicderiv->data[N / 2].re;
    derivs->data[i * (N / 2 + 1) + N / 2].im = (cos(phi0) + sin(phi0)) * intrinsicderiv->data[N / 2].im;
  }

  /* Compute full metric */
  invpsd		= XLALCreateREAL8Vector(N / 2 + 1);
  tempnorm		= XLALCreateREAL8Vector(N);
  tempnormtilde = XLALCreateCOMPLEX16Vector(N / 2 + 1);

  for (k = 0; k < N / 2 + 1; ++k)
  {
    if (psd->data->data[k] == 0.)
      invpsd->data[k] = 0.;
    else
      invpsd->data[k] = 1.0 / psd->data->data[k];
  }

  for (i = 0; i < 4; ++i)
  {
    for (j = 0; j < i + 1; ++j)
    {
      for (k = 0; k < N / 2 + 1; ++k)
      {
        tempnormtilde->data[k].re =
          ( derivs->data[i * (N/2+1) + k].re * derivs->data[j * (N/2+1) + k].re
            + derivs->data[i * (N/2+1) + k].im * derivs->data[j * (N/2+1) + k].im)
          * invpsd->data[k];
        tempnormtilde->data[k].im =
          ( derivs->data[i * (N/2+1) + k].im * derivs->data[j * (N/2+1) + k].re
            - derivs->data[i * (N/2+1) + k].re * derivs->data[j * (N/2+1) + k].im)
          * invpsd->data[k];
      }
 
      /* Inverse Fourier of tempnorm */
      XLALREAL8ReverseFFT(tempnorm, tempnormtilde, revPlan);

      fisher->data[i * (i + 1) / 2 + j] = 4.0 * tempnorm->data[0];

    }
  }

 


  /* this memory deallocation code should be moved to a separate function */
  XLALDestroyCOMPLEX16Vector(intrinsicderiv);
  XLALDestroyCOMPLEX16VectorSequence(derivs);
  XLALDestroyCOMPLEX16Vector(tempnormtilde);
  XLALDestroyREAL8Vector(tempnorm);
  XLALDestroyREAL8VectorSequence( PTFQ );
  XLALDestroyCOMPLEX16VectorSequence( PTFQtilde );
  XLALDestroyVector( PTFphi );
  XLALDestroyVector( PTFomega_2_3 );
  XLALDestroyVectorSequence( PTFe1 );
  XLALDestroyVectorSequence( PTFe2 );
  XLALDestroyREAL8FFTPlan( fwdPlan );
  XLALDestroyREAL8FFTPlan( revPlan );

  /* normal exit */
  return errcode;
}





