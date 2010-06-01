/*
*  Copyright (C) 2007 Diego Fazi, Duncan Brown
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
 * File Name: FindChirpPTFFilter.c
 *
 * Author: Brown, D. A. and Fazi, D.
 *
 * Revision: $Id$
 *
 *-----------------------------------------------------------------------
 */

#if 0
<lalVerbatim file="FindChirpPTFFilterCV">
Author: Brown, D. A. and Fazi, D.
$Id$
</lalVerbatim>

<lalLaTeX>
\input{FindChirpPTFFilterCDoc}

\vfill{\footnotesize\input{FindChirpPTFFilterCV}}
#endif

#include <config.h>
#include <stdlib.h>
#include <math.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALError.h>
#include <lal/LALConstants.h>
#include <lal/Date.h>
#include <lal/AVFactories.h>
#include <lal/FindChirp.h>
#include <lal/LALInspiral.h>
#include <lal/FindChirpPTF.h>
#include <lal/MatrixUtils.h>

NRCSID (FINDCHIRPPTFFILTERC, "$Id$");

/* <lalVerbatim file="FindChirpPTFFilterCP"> */

double rint(double x);

void
LALFindChirpPTFFilterSegment (
    LALStatus                  *status,
    SnglInspiralTable         **eventList,
    FindChirpFilterInput       *input,
    FindChirpFilterParams      *params
    )
/* </lalVerbatim> */
{
  UINT4                 i, j, k, l, kmax, kmin;
  UINT4                 numPoints;
  UINT4                 deltaEventIndex;
  UINT4                 ignoreIndex;
  REAL4                 deltaT, deltaF, fFinal, f_min, length;
  REAL4                 N, r, s, x, y;
  REAL8                 u1[5], u2[5], v1[5], v2[5];
  REAL8                 max_eigen, c, alpha, beta;
  REAL8                 v1_dot_u1, v1_dot_u2, v2_dot_u1, v2_dot_u2, 
                        v1u1_minus_v2u2;
  REAL4                *Binv        = NULL;
  REAL4                *PTFP        = NULL;
  REAL8                *rho         = NULL;
  COMPLEX8             *q           = NULL;
  COMPLEX8             *PTFQtilde   = NULL;
  COMPLEX8             *qtilde      = NULL;
  COMPLEX8             *PTFqtilde   = NULL;
  COMPLEX8             *PTFq        = NULL;
  COMPLEX8             *inputData   = NULL;
  COMPLEX8             *tmpltSignal = NULL;
  COMPLEX8Vector        qVec;
  /* FindChirpBankVetoData clusterInput; */
  
  INITSTATUS( status, "LALFindChirpPTFFilter", FINDCHIRPPTFFILTERC );
  ATTATCHSTATUSPTR( status );


  /*
   *
   * check that the arguments are reasonable
   *
   */


  /* make sure the output handle exists, but points to a null pointer */
  ASSERT( eventList, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT( !*eventList, status, FINDCHIRPH_ENNUL, FINDCHIRPH_MSGENNUL );

  /* make sure that the parameter structure exists */
  ASSERT( params, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );

  /* check that the filter parameters are reasonable */
  ASSERT( params->deltaT > 0, status,
      FINDCHIRPH_EDTZO, FINDCHIRPH_MSGEDTZO );
  ASSERT( params->rhosqThresh >= 0, status,
      FINDCHIRPH_ERHOT, FINDCHIRPH_MSGERHOT );
  ASSERT( params->chisqThresh >= 0, status, 
      FINDCHIRPH_ECHIT, FINDCHIRPH_MSGECHIT );
  ASSERT( params->chisqDelta >= 0, status,
      FINDCHIRPH_ECHIT, FINDCHIRPH_MSGECHIT );

  /* check that the fft plan exists */
  ASSERT( params->invPlan, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );

  /* check that the workspace vectors exist */
  ASSERT( params->qtildeVec, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT( params->qtildeVec->data, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );

  /* check that the chisq parameter and input structures exist */
  ASSERT( params->chisqParams, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT( params->chisqInput, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  
  /* if a rhosqVec vector has been created, check we can store data in it */
  if ( params->rhosqVec ) 
  {
    ASSERT( params->rhosqVec->data->data, status, 
        FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
    ASSERT( params->rhosqVec->data, status, 
        FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  }
  if ( params->cVec )
  { 
    ASSERT( params->cVec->data->data, status,  
        FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL  );
    ASSERT( params->cVec->data, status,
        FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  }

  /* if we are doing a chisq, check we can store the data */
  if ( input->segment->chisqBinVec->length )
  { 
    ASSERT( params->chisqVec, status,
        FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
    ASSERT( params->chisqVec->data, status,
        FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  }

  /* make sure that the input structure exists */
  ASSERT( input, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );

  /* make sure that the input structure contains some input */
  ASSERT( input->fcTmplt, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT( input->segment, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );

  /* make sure that the filter has been initialized for the correct */
  /* approximant                                                    */
  if ( params->approximant != FindChirpPTF )
  {
    ABORT( status, FINDCHIRPH_EAPRX, FINDCHIRPH_MSGEAPRX );
  }

  /* make sure the approximant in the tmplt and segment agree */
  if ( params->approximant != input->fcTmplt->tmplt.approximant ||
      params->approximant != input->segment->approximant )
  {
    ABORT( status, FINDCHIRPH_EAPRX, FINDCHIRPH_MSGEAPRX );
  }
  
#if 0  
       ASSERT( input->fcTmplt->tmplt.approximant == FindChirpPTF, status,
       FINDCHIRPH_EAPRX, FINDCHIRPH_MSGEAPRX );
       ASSERT( input->segment->approximant == FindChirpPTF, status,
       FINDCHIRPH_EAPRX, FINDCHIRPH_MSGEAPRX );
#endif

  /*
   *
   * point local pointers to input and output pointers
   *
   */

  /* number of points in a segment */
  numPoints    = params->PTFqVec->vectorLength;
  N            = (REAL4) numPoints;
  /* workspace vectors */
  rho          = params->PTFsnrVec->data;
  q            = params->qVec->data;
  qtilde       = params->qtildeVec->data;
  PTFq         = params->PTFqVec->data;
  PTFqtilde    = params->PTFqtildeVec->data;
  PTFP         = params->PTFPVec->data;
  qVec.length  = numPoints;

  /* template and data */
  tmpltSignal  = input->fcTmplt->data->data;
  inputData    = input->segment->data->data->data;
  length       = input->segment->data->data->length;
  PTFQtilde    = input->fcTmplt->PTFQtilde->data;
  Binv         = input->fcTmplt->PTFBinverse->data;

  /* number of points and frequency cutoffs */
  deltaT       = (REAL4) params->deltaT;
  deltaF       = 1.0 / ( deltaT * N );
  fFinal       = (REAL4) input->fcTmplt->tmplt.fFinal;
  f_min        = (REAL4) input->fcTmplt->tmplt.fLower;
  kmax         = fFinal / deltaF < numPoints/2 ? fFinal / deltaF : numPoints/2;
  kmin         = f_min / deltaF > 1.0 ? f_min/ deltaF : 1;

  /*
   *
   * compute viable search regions in the snrsq vector
   *
   */

  if ( input->fcTmplt->tmplt.tC <= 0 )
  {
    ABORT( status, FINDCHIRPH_ECHTZ, FINDCHIRPH_MSGECHTZ );
  }

  deltaEventIndex = (UINT4) rint( (input->fcTmplt->tmplt.tC / deltaT) + 1.0 );

  /* ignore corrupted data at start and end */
  params->ignoreIndex = ignoreIndex = 
    ( input->segment->invSpecTrunc / 2 ) + deltaEventIndex;

  if ( lalDebugLevel & LALINFO )
  {
    CHAR infomsg[256];

    LALSnprintf( infomsg, sizeof(infomsg) / sizeof(*infomsg),
        "m1 = %e, m2 = %e, chi = %e, kappa = %e "
        "=> %e seconds => %d points\n"
        "invSpecTrunc = %d => ignoreIndex = %d\n", 
        input->fcTmplt->tmplt.mass1, input->fcTmplt->tmplt.mass2, 
        input->fcTmplt->tmplt.chi, input->fcTmplt->tmplt.kappa, 
        input->fcTmplt->tmplt.tC , deltaEventIndex, 
        input->segment->invSpecTrunc, ignoreIndex );
    LALInfo( status, infomsg );
  }

  /* XXX check that we are not filtering corrupted data XXX */
  /* XXX this is hardwired to 1/4 segment length        XXX */
  if ( ignoreIndex > numPoints / 4 )
  {
    ABORT( status, FINDCHIRPH_ECRUP, FINDCHIRPH_MSGECRUP );
  }
  /* XXX reset ignoreIndex to one quarter of a segment XXX */
  params->ignoreIndex = ignoreIndex = numPoints / 4;

  if ( lalDebugLevel & LALINFO )
  {
    CHAR infomsg[256];

    LALSnprintf( infomsg, sizeof(infomsg) / sizeof(*infomsg), 
        "filtering from %d to %d\n",
        ignoreIndex, numPoints - ignoreIndex );
    LALInfo( status, infomsg );
  }

  /*
   *
   * compute the PTF filter statistic
   *
   */

  /* clear the snr output vector and workspace*/
  memset( params->PTFsnrVec->data, 0, numPoints * sizeof(REAL8) );
  memset( params->qVec->data, 0, numPoints * sizeof(COMPLEX8) );
  memset( params->PTFPVec->data, 0, 5 * numPoints * sizeof(REAL4) );
  memset( params->PTFqVec->data, 0, 5 * numPoints * sizeof(COMPLEX8) );
  memset( params->PTFqtildeVec->data, 0, 5 * ( numPoints / 2 + 1 ) 
      * sizeof(COMPLEX8) );
  
  for ( i = 0; i < 5; ++i )
  {

    /* compute qtilde using data and Qtilde */

    memset( params->qtildeVec->data, 0, numPoints * sizeof(COMPLEX8) );

    /* qtilde positive frequency, not DC or nyquist */
    for ( k = kmin; k < kmax ; ++k )
    {
      r = inputData[k].re; /* Re[s_tilde] */
      s = inputData[k].im; /* Im[s_tilde] */
      x = PTFQtilde[i * (numPoints / 2 + 1) + k].re; /* Re[Q*_0_tilde] */
      y = 0 - PTFQtilde[i * (numPoints / 2 + 1) + k].im; /* Im[Q*_0_tilde] */
      qtilde[k].re = PTFqtilde[ i * (numPoints / 2 + 1) + k].re =  
                     2 * (r*x - s*y);
      qtilde[k].im = PTFqtilde[ i * (numPoints / 2 + 1) + k].im =  
                     2 * (r*y + s*x);
    }
    qVec.data = params->PTFqVec->data + (i * numPoints);

    /* inverse fft to get <s,Q_0>+i<s,Q_pi> */
    LALCOMPLEX8VectorFFT( status->statusPtr, &qVec, params->qtildeVec, 
        params->invPlan );
    CHECKSTATUSPTR( status );
  }

  /* now we have PTFqVec->data which contains <s|Q^I_0> + i <s|Q^I_\pi/2> */
   
  for ( j = 0; j < numPoints; ++j ) /* beginning of main loop over time */
  {  
    for (i = 0; i < 5; i++)
    { 
      v1[i] = PTFq[i * numPoints + j].re;
      v2[i] = PTFq[i * numPoints + j].im;
    }
    /* construct the vectors u[i] = B^(-1) v[i] */
    for (i = 0; i < 5; i++)
    { 
      u1[i] = 0.0;
      u2[i] = 0.0;
      for ( l = 0; l < 5; l++ )
      {  
        u1[i] = u1[i] + (REAL8) Binv[i * 5 + l] * v1[l];
        u2[i] = u2[i] + (REAL8) Binv[i * 5 + l] * v2[l];
      }
    }

    /* Compute SNR */ 
    v1_dot_u1 = v1_dot_u2 = v2_dot_u1 = v2_dot_u2 = max_eigen = 0.0;
    for (i = 0; i < 5; i++)
    {
      v1_dot_u1 = v1_dot_u1 + v1[i] * u1[i];
      v1_dot_u2 = v1_dot_u2 + v1[i] * u2[i];
      v2_dot_u1 = v2_dot_u1 + v2[i] * u1[i];
      v2_dot_u2 = v2_dot_u2 + v2[i] * u2[i];
    }
    v1u1_minus_v2u2 = v1_dot_u1 - v2_dot_u2;
    max_eigen = 0.5 * ( v1_dot_u1 + v2_dot_u2 + sqrt( v1u1_minus_v2u2 * 
          v1u1_minus_v2u2 + 4 * v1_dot_u2 * v2_dot_u1 ));
     q[j].re = (REAL4) (2.0 * sqrt(max_eigen) / N); /* snr */
     rho[j] = sqrt ( v1[0] * v1[0] + v2[0] * v2[0] ); 
    
    /* evaluate extrinsic parameters for every coalescence time */
    c = ( max_eigen - v1_dot_u1 ) / v1_dot_u2 ;
    alpha = v1_dot_u1 + c * (v2_dot_u1 + v1_dot_u2) + c * c * v2_dot_u2;
    alpha = sqrt (1.0 / alpha);
    beta = c * alpha;
    for (i=0; i<5; i++)
    {
      PTFP[i * numPoints + j] = (REAL4) (alpha * u1[i] + beta * u2[i]);
    }
   
  } /* End of main loop over time */
  
  

  /* 
   *
   * calculate signal to noise squared 
   *
   */


  /* if full snrsq vector is required, set it to zero */
  if ( params->rhosqVec )
    memset( params->rhosqVec->data->data, 0, numPoints * sizeof( REAL4 ) );
  
  if (params->cVec )
    memset( params->cVec->data->data, 0, numPoints * sizeof( COMPLEX8 ) );
  
  /* if full snrsq vector is required, store the snrsq */
  if ( params->rhosqVec ) 
  {
    memcpy( params->rhosqVec->name, input->segment->data->name,
        LALNameLength * sizeof(CHAR) );
    memcpy( &(params->rhosqVec->epoch), &(input->segment->data->epoch), 
        sizeof(LIGOTimeGPS) );
    params->rhosqVec->deltaT = input->segment->deltaT;

    for ( j = 0; j < numPoints; ++j )
    {
      params->rhosqVec->data->data[j] =  q[j].re * q[j].re;
      
    }
  }


  /*
   *
   * look for and cluster events in the snr vector
   *
   */

#if 0  
  clusterInput.length       = 1; /* do not do the bank veto */
  clusterInput.qVecArray    = NULL;
  clusterInput.fcInputArray = NULL;
  clusterInput.ccMat        = NULL;
  clusterInput.normMat      = NULL;
  clusterInput.spec         = NULL;
  clusterInput.resp         = NULL;
#endif
  

  /* cluster events searches the qVec for events */
  /* params->qVec = params->PTFsnrVec; */
  input->fcTmplt->norm = 1.0;
 
  if ( params->cVec )
  {
    memcpy( params->cVec->name, input->segment->data->name,
        LALNameLength * sizeof(CHAR) );
    memcpy( &(params->cVec->epoch), &(input->segment->data->epoch),
        sizeof(LIGOTimeGPS) );
    params->cVec->deltaT = input->segment->deltaT;

    for ( j = 0; j < numPoints; ++j )
    {
      params->cVec->data->data[j].re = q[j].re;
      params->cVec->data->data[j].im = 0.0 ;
    }
  }
  
#if 0
  LALFindChirpClusterEvents( status->statusPtr, eventList, 
     input, params, &clusterInput, 0 );
  CHECKSTATUSPTR( status );

  params->qVec = NULL;
#endif  
  
  /* normal exit */
  DETATCHSTATUSPTR( status );
  RETURN( status );
}

