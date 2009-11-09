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
<lalVerbatim file="FisherTestCV">
Author: Tjonnie Li, Chris Van Den Broeck
$Id$
</lalVerbatim>

<lalLaTeX>

\subsection{Program \texttt{FisherTest.c}}
\label{ss:FisherTest.c}

Calculates the Fisher matrix at a single point in parameter space.

\subsubsection*{Usage}

\subsubsection*{Description}

\subsubsection*{Exit codes}

\subsubsection*{Algorithm}

\subsubsection*{Uses}
\begin{verbatim}
lalDebugLevel
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{FisherTestCV}}

</lalLaTeX>
#endif

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <lal/LALInspiralBank.h>
#include <lal/LALStdio.h>
#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>

NRCSID(METRICTESTC,"$Id$");

int lalDebugLevel = 1;

int main( void )
{
  UINT4 i;
  LALStatus status;
  INT4 errcode;
  /* create memory for the Fisher matrix */
  REAL8Vector *fisher;
  InspiralTemplate tmplt;
  REAL8FrequencySeries psd;
  void (*noisemodel)(LALStatus*,REAL8*,REAL8) = LALLIGOIPsd;
  REAL8 deltaT = 1.0/16384.0;
  UINT4 N = 16384;
  REAL8 deltaF = 1.0/((REAL8)N * deltaT);

  fisher = XLALCreateREAL8Vector( 16 );

  memset( &status, 0, sizeof(LALStatus) );
  memset( &tmplt, 0, sizeof(InspiralTemplate) );

  tmplt.massChoice = m1Andm2;
  tmplt.mass1 = 5.0;
  tmplt.mass2 = 1.4;
  tmplt.chi = 0.8;
  tmplt.kappa = 0.5;

  tmplt.fLower=100.0;
  tmplt.fCutoff = 1.0 / (2.0 * deltaT);
  tmplt.tSampling = deltaT;

  tmplt.sourceTheta = LAL_PI / 3.;
  tmplt.sourcePhi = LAL_PI / 6.;
  tmplt.polarisationAngle = LAL_PI / 4.;
  tmplt.startPhase = 0.;
  tmplt.startTime = 0.;
  tmplt.signalAmplitude = 1.0;

  LALInspiralParameterCalc(&status, &tmplt);
  if ( status.statusCode )
  {
    REPORTSTATUS( &status );
    exit( 1 );
  }

  /* create memory for the PSD */
  memset( &psd, 0, sizeof(REAL8FrequencySeries) );
  LALDCreateVector( &status, &(psd.data), N / 2 + 1 );
  if ( status.statusCode )
  {
    REPORTSTATUS( &status );
    exit( 1 );
  }

  psd.deltaF = deltaF;

  /* create a LIGO PSD */
  LALNoiseSpectralDensity( &status, psd.data, noisemodel, psd.deltaF );
  if ( status.statusCode )
  {
    REPORTSTATUS( &status );
    exit( 1 );
  }

  errcode = XLALInspiralComputeFisherMatrix( fisher, &psd, &tmplt );
  if ( errcode != XLAL_SUCCESS )
  {
    fprintf( stderr, "XLALInspiralComputeFisherMatrix failed\n" );
    exit( 1 );
  }

/*
  fprintf( stderr, "Printing out components of the projected metric in the intrinsic parameter space\n");
  for ( i = 0; i < 10; ++i )
  {
    fprintf( stderr, "Gamma[%d] = %e\n", i, fisher.Gamma[i] );
  }
*/

  /* destory memory for fisher */
  XLALDestroyREAL8Vector( fisher );

  return 0;
}
