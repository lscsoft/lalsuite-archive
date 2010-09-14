#if 0
<lalVerbatim file="FindChirpACTDFilterCV">
Author: Brown D. A., McKechan D. J. A., Robinson C.
$Id$
</lalVerbatim>

<lalLaTeX>
\subsection{Module \texttt{FindChirpACTDFilter.c}}
\label{ss:FindChirpACTDFilter.c}

This module provides the core of the matched filter for amplitude corrected
binary inspiral chirps.


\subsubsection*{Prototypes}
\vspace{0.1in}
\input{FindChirpACTDFilterCP}
\idx{LALFindChirpACTDFilterSegment()}

\subsubsection*{Description}

\subsubsection*{Algorithm}

\subsubsection*{Notes}

\vfill{\footnotesize\input{FindChirpACTDFilterCV}}
</lalLaTeX>
#endif

#include <math.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALError.h>
#include <lal/LALConstants.h>
#include <lal/Date.h>
#include <lal/AVFactories.h>
#include <lal/FindChirp.h>
#include <lal/FindChirpChisq.h>
#include <lal/FindChirpACTD.h>
#include <lal/NonCentralChisq.h>

double rint(double x);
/* debugging */
extern int vrbflg;                      /* verbocity of lal function    */

static REAL4
XLALFindHarmonicRhoSq( REAL8 expRhoSq, REAL8 cdfLevel );

NRCSID (FINDCHIRPACTDFILTERC, "$Id$");


/* <lalVerbatim file="FindChirpACTDFilterCP"> */
void
LALFindChirpACTDFilterSegment (
    LALStatus                  *status,
    SnglInspiralTable         **eventList,
    FindChirpFilterInput       *input,
    FindChirpFilterParams      *params
    )
/* </lalVerbatim> */
{
  UINT4                 i, j, k, kmax;
  UINT4                 numPoints;
  UINT4                 tmpltLength;
  UINT4                 deltaEventIndex;
  UINT4                 ignoreIndex;
  REAL8                 deltaT;
  REAL8                 deltaF;
  REAL4                 normFac;
  REAL4                 normFacSq;
  COMPLEX8Vector      **qtilde;
  REAL4Vector         **q;
  COMPLEX8             *inputData     = NULL;
  COMPLEX8Vector        tmpltSignal[NACTDTILDEVECS];

  /*
  SnglInspiralTable    *thisEvent     = NULL;
  REAL4                 modqsqThresh;
  BOOLEAN               haveEvent     = 0;
  */

  INITSTATUS( status, "LALFindChirpACTDFilter", FINDCHIRPACTDFILTERC );
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
  ASSERT( params->invPlanACTD, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );

  /* check that the workspace vectors exist */
  ASSERT( params->qVecACTD, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT( params->qtildeVecACTD, status, FINDCHIRPH_ENULL,
          FINDCHIRPH_MSGENULL );

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

  /* check the approximant */
  if ( params->approximant != AmpCorPPN )
  {
    ABORT( status, FINDCHIRPH_EUAPX, FINDCHIRPH_MSGEUAPX );
  }

  /* make sure the approximant in the tmplt and segment agree */
  if ( params->approximant != input->fcTmplt->tmplt.approximant ||
      params->approximant != input->segment->approximant )
  {
    ABORT( status, FINDCHIRPH_EAPRX, FINDCHIRPH_MSGEAPRX );
  }



  /*
   *
   * point local pointers to input and output pointers
   *
   */

  /* number of points in a segment */
  numPoints = params->qVecACTD[0]->length;
  tmpltLength = input->fcTmplt->ACTDtilde->vectorLength;


  /* workspace vectors */
  qtilde = params->qtildeVecACTD;
  q      = params->qVecACTD;

  for( i = 0; i < NACTDTILDEVECS; ++i )
  {
    /* filter */
	  memset( qtilde[i]->data, 0, qtilde[i]->length * sizeof( COMPLEX8 ) );
    memset( q[i]->data, 0, numPoints * sizeof( REAL4 ) );

		/* template */
    tmpltSignal[i].length = tmpltLength;
    tmpltSignal[i].data   = input->fcTmplt->ACTDtilde->data + i * tmpltLength;
  }


  /* data */
  inputData = input->segment->data->data->data;


  deltaT = params->deltaT;

  deltaF = 1.0 / ( deltaT * (REAL8)(numPoints)  );
  /*kmax = input->fcTmplt->tmplt.fFinal / deltaF < numPoints/2 ?
  input->fcTmplt->tmplt.fFinal / deltaF : numPoints/2;*/
  /* Craig: I suspect the limits of integration are incorrect. I will
     artificially set it to numPoints/2. */
  kmax = numPoints / 2;


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

    snprintf( infomsg, sizeof(infomsg) / sizeof(*infomsg),
        "m1 = %e m2 = %e => %e seconds => %d points\n"
        "invSpecTrunc = %d => ignoreIndex = %d\n",
        input->fcTmplt->tmplt.mass1, input->fcTmplt->tmplt.mass2,
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

    snprintf( infomsg, sizeof(infomsg) / sizeof(*infomsg),
        "filtering from %d to %d\n",
        ignoreIndex, numPoints - ignoreIndex );
    LALInfo( status, infomsg );
  }


  /*
   *
   * compute qtildes and qs
   *
   */


  /* qtilde positive frequency, not DC or nyquist */
  for( i = 0; i < NACTDTILDEVECS; ++i )
  {
    if ( (i % NACTDVECS >= input->fcTmplt->startVecACTD )
      && (i % NACTDVECS < input->fcTmplt->stopVecACTD ) )
    {
      for ( k = 1; k < kmax; ++k )
      {
        REAL4 r = inputData[k].re;
        REAL4 s = inputData[k].im;
        REAL4 x = tmpltSignal[i].data[k].re;
        REAL4 y = 0 - tmpltSignal[i].data[k].im; /* NB: Complex conj. */

        qtilde[i]->data[k].re = r*x - s*y;
        qtilde[i]->data[k].im = r*y + s*x;
      }

      /* inverse fft to get q */
      if ( XLALREAL4ReverseFFT( q[i], qtilde[i], params->invPlanACTD ) == 
                                                                XLAL_FAILURE )
      {
        ABORTXLAL( status );
      }
    }
  }




  /*
   *
   * calculate signal to noise squared
   *
   */

  if (params->cVec )
    memset( params->cVec->data->data, 0, numPoints * sizeof( COMPLEX8 ) );

  /* if full snrsq vector is required, store the snrsq */
  if ( params->rhosqVec )
  {
    memset( params->rhosqVec->data->data, 0, numPoints * sizeof( REAL4 ) );
    /* normalisation */

    normFac = 2.0 * deltaT / (REAL4)(numPoints);
	  normFacSq = normFac * normFac;

    memcpy( params->rhosqVec->name, input->segment->data->name,
        LALNameLength * sizeof(CHAR) );
    memcpy( &(params->rhosqVec->epoch), &(input->segment->data->epoch),
        sizeof(LIGOTimeGPS) );
    params->rhosqVec->deltaT = input->segment->deltaT;

    for ( j = 0; j < numPoints; ++j )
    {
      REAL4 rhoSq = 0.0;
      
      for( i = 0; i < NACTDTILDEVECS; ++i )
      {
        rhoSq += q[i]->data[j] * q[i]->data[j];
      }

      params->rhosqVec->data->data[j] = rhoSq * normFacSq ;

    }
  }

  /* normal exit */
  DETATCHSTATUSPTR( status );
  RETURN( status );
}

INT4 XLALFindChirpACTDApplyConstraint(
      INT4                        qidx,
      FindChirpFilterInput      * restrict input,
      FindChirpFilterParams     * restrict params )
{
  
  static const char func[] = "XLALFindChirpACTDApplyConstraint";

  REAL4Vector *qSqForConstraint = NULL;
  UINT4        matrixDim;
  INT4         apply = 0;
  UINT4        i, k;

  /* Constraints to apply taking into account SNR */
  REAL4Vector *cons = NULL;

  REAL4 normFacSq;

  /* The SNR-squared as would be calculated from */
  /* non-orthogonal templates */
  REAL4 falseRhoSq = 0.0;

  /* Create the vectors for the constraint */
  qSqForConstraint = XLALCreateREAL4Vector( NACTDTILDEVECS );
  cons             = XLALCreateREAL4Vector( input->fcTmplt->ACTDconstraint->length );

  memset( cons->data, 0, input->fcTmplt->ACTDconstraint->length * sizeof( REAL4 ) );
  memset( qSqForConstraint->data, 0, NACTDTILDEVECS * sizeof( REAL4 ));

  /* Calculate the dimension of the matrix used in the constraint */
  matrixDim = 2 * (input->fcTmplt->stopVecACTD - input->fcTmplt->startVecACTD );

  /* Normalization factor for rho squared */
  normFacSq = 2.0 * params->deltaT / (REAL4)(params->qVecACTD[0]->length);
  normFacSq = normFacSq * normFacSq;

  /* Go through the transformation matrix */
  for ( i = 0; i < NACTDTILDEVECS; i++ )
  {
    for ( k = 0; k < NACTDTILDEVECS; k++ )
    {
      if ( ( i % NACTDVECS - input->fcTmplt->startVecACTD < matrixDim / 2 )
         && ( k % NACTDVECS - input->fcTmplt->startVecACTD < matrixDim / 2 ) )
      {
        UINT4 idx1, idx2, nVec;

        nVec = matrixDim / 2;
        idx1 = (i / NACTDVECS) * nVec + i % NACTDVECS - input->fcTmplt->startVecACTD;
        idx2 = (k / NACTDVECS) * nVec + k % NACTDVECS - input->fcTmplt->startVecACTD;

        qSqForConstraint->data[i] += 
          gsl_matrix_get( input->fcTmplt->ACTDconmatrix, idx1, idx2 ) 
          * params->qVecACTD[k]->data[qidx]
          * sqrt( gsl_vector_get( input->fcTmplt->ACTDconvector, idx2 ) );
      }
    }
    qSqForConstraint->data[i] *= qSqForConstraint->data[i];
    falseRhoSq          += qSqForConstraint->data[i];
  }

  falseRhoSq *= normFacSq;

  /* Having calculated the SNR in this basis, we need to tweak */
  /* the constraint so that it takes into account the SNR.     */
  /* We also have to take into account the variance for noise. */
  for ( i=input->fcTmplt->startVecACTD; i < cons->length; i++ )
  {
    REAL8 noNoiseConstraint = input->fcTmplt->ACTDconstraint->data[i] * falseRhoSq;

    /* Find the SNR to be used in the constraint taking into account noise. */
    /* To give ourselves some headroom, we take the 90% level for the       */
    /* numerator, and the 10% level for the denominator.                    */
    if ( i == 1 )
    {
      cons->data[i] = XLALFindHarmonicRhoSq( noNoiseConstraint, 0.1 );
    }
    else
    {
      cons->data[i] = XLALFindHarmonicRhoSq( noNoiseConstraint, 0.9 );
    }
    if ( XLAL_IS_REAL4_FAIL_NAN( cons->data[i] ) )
    {
      XLAL_ERROR( func, XLAL_EFUNC );
    }
  }

  /* ...and now for something completely different */
  for ( i=input->fcTmplt->startVecACTD; i < cons->length; i++ )
  {
    cons->data[i] /= cons->data[1];
  }

  /* Compare ratio of first harmonic with second in original basis */
  if ( (qSqForConstraint->data[0]
            + qSqForConstraint->data[3])
           / (qSqForConstraint->data[1]
            + qSqForConstraint->data[4])
           > cons->data[0] )
  {
    apply = 1;
  }

  /* Compare ratio of third harmonic with second in original basis */
  if ( ( qSqForConstraint->data[2]
            +  qSqForConstraint->data[5])
           / (qSqForConstraint->data[1]
            + qSqForConstraint->data[4])
           > cons->data[2] )
  {
    apply = 1;
  }


  XLALDestroyREAL4Vector( qSqForConstraint );
  XLALDestroyREAL4Vector( cons );
 
  return apply;
}


static REAL4
XLALFindHarmonicRhoSq( REAL8 expRhoSq, REAL8 cdfLevel )
{

  static const char func[] = "XLALFindHarmonicRhoSq";

  /* For now we hard-code the degrees of freedom to be 2 */
  const INT4 dof = 2;

  const INT4 maxIter = 100;

  const REAL8 epsilon = 1.0e-3;
  /* Non-central chisq CDF and PDF values */
  REAL8 cdf, pdf;

  REAL8 rhoSq, oldRhoSq;

  /* Sometimes the root-finding will overshoot */
  /* If it happens more than once, we will throw an error */
  INT4  overshot = 0;

  INT4 i = 0;

  /* Choose an appropriate starting point */
  rhoSq    = expRhoSq > (REAL8)dof ? expRhoSq : (REAL8)dof;
  oldRhoSq = rhoSq + 1.0;

  while ( fabs( rhoSq - oldRhoSq ) > 1.0e-8 * rhoSq && i <= maxIter )
  {
    oldRhoSq = rhoSq;

    cdf = XLALNonCentralChisqCDF( rhoSq, dof, expRhoSq );
    pdf = XLALNonCentralChisqPDF( rhoSq, dof, expRhoSq );

    if ( XLAL_IS_REAL8_FAIL_NAN( cdf ) || XLAL_IS_REAL8_FAIL_NAN( pdf ) )
    {
      XLAL_ERROR_REAL4( func, XLAL_EFUNC );
    }
    rhoSq = rhoSq - ( cdf - cdfLevel ) / pdf;

    /* If we have overshot, just set rhoSq to a small number */
    if ( rhoSq < 0 )
    {
      if ( !overshot )
      {
        rhoSq = epsilon;
        overshot = 1;
      }
      else
      {
        XLAL_ERROR_REAL4( func, XLAL_ETOL );
      }
    }
    i++;
  }

  if ( i > maxIter )
    XLAL_ERROR_REAL4( func, XLAL_EMAXITER );

  return (REAL4)rhoSq;
}
