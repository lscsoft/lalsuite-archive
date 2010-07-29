#if 0
<lalVerbatim file="FindChirpACTDTemplateCV">
Author: Brown, D. A., Creighton, J. D. E. and Mckechan, D. J. A.
$Id: FindChirpACTDTemplate.c,v 1.2.2.4.2.1 2009/03/31 11:25:18 spxcar Exp $
</lalVerbatim> 

<lalLaTeX>
\subsection{Module \texttt{FindChirpACTDTemplate.c}}
\label{ss:FindChirpACTDTemplate.c}

Provides functions to create time domain inspiral templates in a
form that can be used by the \texttt{FindChirpACTDFilter()} function.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{FindChirpACTDTemplateCP}
\idx{LALFindChirpACTDTemplate()}
\idx{LALFindChirpACTDNormalize()}

The function \texttt{LALFindChirpACTDTemplate()} creates a time
domain template using LALGeneratePPNAmpCorInspiral().

\subsubsection*{Algorithm}

Blah.

\subsubsection*{Uses}
\begin{verbatim}
LALCalloc()
LALFree()
LALCreateVector()
LALDestroyVector()
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{FindChirpACTDTemplateCV}}
</lalLaTeX>
#endif

#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>
#include <lal/LALInspiral.h>
#include <lal/FindChirp.h>
#include <lal/FindChirpTD.h>
#include <lal/FindChirpACTD.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_randist.h>

NRCSID (FINDCHIRPACTDTEMPLATEC, "$Id: FindChirpACTDTemplate.c,v 1.2.2.4.2.1 2009/03/31 11:25:18 spxcar Exp $");


/* <lalVerbatim file="FindChirpACTDTemplateCP"> */
void
LALFindChirpACTDTemplate(
    LALStatus                  *status,
    FindChirpTemplate          *fcTmplt,
    InspiralTemplate           *tmplt,
    FindChirpTmpltParams       *params
    )
/* </lalVerbatim> */
{
  UINT4          i, j;
  UINT4          istart = 0;
  UINT4          istop = NACTDVECS;
  UINT4          shift;
  UINT4          numPoints;
  REAL4Vector    ACTDVecs[NACTDVECS];
  COMPLEX8Vector ACTDtilde[NACTDTILDEVECS];
  REAL8          deltaT;
  REAL8          sampleRate;
  const REAL4    cannonDist = 1.0; /* Mpc */
  CHAR           infomsg[512];
  PPNParamStruc  ppnParams;
  CoherentGW     waveform;

  REAL4Vector  *tmpACTDVec = NULL; /* Used for band-passing */


  INITSTATUS( status, "LALFindChirpACTDTemplate", FINDCHIRPACTDTEMPLATEC );
  ATTATCHSTATUSPTR( status );


  /*
   *
   * check that the arguments are reasonable
   *
   */
 
  ASSERT( NACTDVECS >= 2, status,
      FINDCHIRPACTDH_EACTDV, FINDCHIRPACTDH_MSGEACTDV ); 

  /* check that the output structures exist */
  ASSERT( fcTmplt, status,
      FINDCHIRPTDH_ENULL, FINDCHIRPTDH_MSGENULL );
  ASSERT( fcTmplt->ACTDtilde, status,
      FINDCHIRPTDH_ENULL, FINDCHIRPTDH_MSGENULL );
  ASSERT( fcTmplt->ACTDtilde->length == NACTDTILDEVECS, status, 
      FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT( fcTmplt->ACTDtilde->data, status,
      FINDCHIRPTDH_ENULL, FINDCHIRPTDH_MSGENULL );

  /* check that the parameter structure exists */
  ASSERT( params, status,
      FINDCHIRPTDH_ENULL, FINDCHIRPTDH_MSGENULL );
  ASSERT( params->ACTDVecs, status,
      FINDCHIRPTDH_ENULL, FINDCHIRPTDH_MSGENULL );
  ASSERT( params->ACTDVecs->length == NACTDVECS, status,
      FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT( params->ACTDVecs->data, status,
      FINDCHIRPTDH_ENULL, FINDCHIRPTDH_MSGENULL );


  /* check we have an fft plan for the template */
  ASSERT( params->fwdPlan, status,
      FINDCHIRPTDH_ENULL, FINDCHIRPTDH_MSGENULL );

  /* check that the timestep is positive */
  ASSERT( params->deltaT > 0, status,
      FINDCHIRPTDH_EDELT, FINDCHIRPTDH_MSGEDELT );

  /* check that the input exists */
  ASSERT( tmplt, status, FINDCHIRPTDH_ENULL, FINDCHIRPTDH_MSGENULL );

  /* check that the  approximant is AmpCorPPN */
  if( params->approximant != AmpCorPPN )
  {
    ABORT( status, FINDCHIRPTDH_EMAPX, FINDCHIRPTDH_MSGEMAPX );
  }

  /* store deltaT and zero out the time domain waveform vector */
  deltaT = params->deltaT;
  sampleRate = 1.0 / deltaT;
  numPoints =  params->ACTDVecs->vectorLength;

  ASSERT( numPoints >= (2 * (fcTmplt->ACTDtilde->vectorLength - 1)), status,
      FINDCHIRPTDH_EMISM, FINDCHIRPTDH_MSGEMISM );
  ASSERT( numPoints < (2 * (fcTmplt->ACTDtilde->vectorLength - 1) + 2), status,
      FINDCHIRPTDH_EMISM, FINDCHIRPTDH_MSGEMISM );


  /*
   *
   * generate the waveform using LALGeneratePPNAmpCorInspiral() from inject
   *
   */



  /* input parameters */
  memset( &ppnParams, 0, sizeof(PPNParamStruc) );
  ppnParams.deltaT = deltaT;
  ppnParams.mTot_real8 = tmplt->mass1 + tmplt->mass2;
  ppnParams.eta_real8  = tmplt->mass1 * tmplt->mass2;
  ppnParams.eta_real8 /= ppnParams.mTot_real8;
  ppnParams.eta_real8 /= ppnParams.mTot_real8;

  /* Set distance at 20Mpc for testing, will be normalised anyway */
  ppnParams.d = 1.0e6 * LAL_PC_SI / params->dynRange;
  ppnParams.fStartIn = params->fLow;
  ppnParams.fStopIn = - 1.0 /
              (6.0 * sqrt(6.0) * LAL_PI * ppnParams.mTot_real8 * LAL_MTSUN_SI);

   /* PPN parameter. */
   ppnParams.ppn = NULL;
   LALSCreateVector( status->statusPtr, &(ppnParams.ppn), tmplt->order + 1 );
   ppnParams.ppn->data[0] = 1.0;
   if ( tmplt->order > 0 )
     ppnParams.ppn->data[1] = 0.0;
   for ( i = 2; i <= (UINT4)( tmplt->order ); i++ )
     ppnParams.ppn->data[i] = 1.0;


  /* ACTD specific */
  ppnParams.inc = LAL_PI_2;
  /* 
  ppnParams.ampOrder = ( INT4 )( tmplt->ampOrder );
  ppnParams.ampOrder = 1;
  */
  ppnParams.ampOrder = ( UINT4 )( NACTDVECS - 2 );

  /* XXX Uncomment below for extra testing XXX */  
  /*
  fprintf( stderr, "\n   tmplt->mass1         = %f\n", tmplt->mass1);
  fprintf( stderr, "   tmplt->mass2         = %f\n", tmplt->mass2);
  fprintf( stderr, "   ppnParams.deltaT     = %e\n", ppnParams.deltaT ); 
  fprintf( stderr, "   ppnParams.mTot_real8 = %0.16e\n", ppnParams.mTot_real8 ); 
  fprintf( stderr, "   ppnParams.eta_real8  = %0.16e\n", ppnParams.eta_real8 ); 
  fprintf( stderr, "   ppnParams.fStartIn   = %e\n", ppnParams.fStartIn ); 
  fprintf( stderr, "   ppnParams.fStopIn    = %e\n", ppnParams.fStopIn ); 
  fprintf( stderr, "   ppnParams.inc        = %e\n", ppnParams.inc ); 
  fprintf( stderr, "   ppnParams.amporder   = %d\n", ppnParams.ampOrder ); 
  for( i = 0; i < tmplt->order + 1; ++i )
  {
    fprintf( stderr, "   ppnParams.ppn->data[%d] = %e\n", i, 
                                    ppnParams.ppn->data[i] );
  }
  */  
  /*   XXX Uncomment above for extra testing XXX */ 


  /* generate waveform amplitude and phase */
  memset( &waveform, 0, sizeof(CoherentGW) );
  LALGeneratePPNAmpCorInspiral( status->statusPtr, &waveform, &ppnParams );
  CHECKSTATUSPTR( status );



  /* print termination information and check sampling */
  LALInfo( status, ppnParams.termDescription );
  if ( ppnParams.dfdt > 2.0 )
  {
  /* ABORT( status, FINDCHIRPTDH_ESMPL, FINDCHIRPTDH_MSGESMPL ); */
  }
  if ( waveform.a->data->length > numPoints )
  {
    ABORT( status, FINDCHIRPTDH_ELONG, FINDCHIRPTDH_MSGELONG );
  }

  memset( params->ACTDVecs->data, 0, NACTDVECS * numPoints * sizeof( REAL4 ) );
  memset( fcTmplt->ACTDtilde->data, 0, 
                   NACTDTILDEVECS * (numPoints / 2 + 1) * sizeof( COMPLEX8 ) ); 

  for( i=0; i < NACTDVECS; ++i )
  {
    ACTDVecs[i].length  = numPoints;
    ACTDtilde[i].length = numPoints / 2 + 1;
    ACTDtilde[i+NACTDVECS].length = numPoints / 2 + 1;
    ACTDVecs[i].data  = params->ACTDVecs->data + (i * numPoints);
    ACTDtilde[i].data = fcTmplt->ACTDtilde->data + (i * (numPoints / 2 + 1 ) );
    ACTDtilde[i+NACTDVECS].data = 
            fcTmplt->ACTDtilde->data + ((i+NACTDVECS) * (numPoints / 2 + 1 ) );
  }

  /* compute h(t) */
  /* legacy - length is the length of the vectors rather than vector Length */
  
  /* This will make all 0 if eta = 0.25 and 2nd harmonic not in band! */
  if( ppnParams.eta_real8 > 0.24999999 || fcTmplt->tmplt.ACTDdominantSwitch )
  {
    istart = 1;
    istop  = 2;
  }
  /* If a harmonic is not in band don't loop over it */  
  if( istart == 0 && - ppnParams.fStopIn / 2.0 <= params->fLow )
    istart++;
  if( - ppnParams.fStopIn  <= params->fLow )
    istart++;
  /* Signal generation would have failed if 3rd harmonic not in band */
  
  /* However, we should still check there is anything to filter, */
  /* because we do not filter the 3rd harmonic for equal mass systems */
  if ( istart == istop )
  {
    XLALPrintError( "All applicable signal harmonics out of band!\n" );
    LALSDestroyVectorSequence( status->statusPtr, &(waveform.h->data) );
    LALSDestroyVectorSequence( status->statusPtr, &(waveform.a->data) );
    LALSDestroyVector( status->statusPtr, &(waveform.f->data) );
    LALSDestroyVector( status->statusPtr, &(ppnParams.ppn) );
    LALFree( waveform.h );
    LALFree( waveform.a );
    LALFree( waveform.f );
    LALFree( waveform.phi );
    ABORT( status, FINDCHIRPH_ECHTZ, FINDCHIRPH_MSGECHTZ );
  }

  /* Set tmplt params so that we can track which harmonics are actually used */
  fcTmplt->startVecACTD = istart;
  fcTmplt->stopVecACTD  = istop;


  for ( j = 0; j < waveform.a->data->length; ++j )
  {
    for ( i = istart; i < istop; ++i )
    {
      ACTDVecs[i].data[j] = waveform.a->data->data[3*j + i]
         * cos( ( ( REAL4 )( i ) + 1.0 ) / 2.0 * waveform.phi->data->data[j] );
    }
  }


  /* free the memory allocated by LALGeneratePPNAmpCorInspiral() */
  LALSDestroyVectorSequence( status->statusPtr, &(waveform.h->data) );
  CHECKSTATUSPTR( status );

  LALSDestroyVectorSequence( status->statusPtr, &(waveform.a->data) );
  CHECKSTATUSPTR( status );

  LALSDestroyVector( status->statusPtr, &(waveform.f->data) );


  CHECKSTATUSPTR( status );

  LALDDestroyVector( status->statusPtr, &(waveform.phi->data) );
  CHECKSTATUSPTR( status );

  LALSDestroyVector( status->statusPtr, &(ppnParams.ppn) );
  CHECKSTATUSPTR( status );
  LALFree( waveform.h );
  LALFree( waveform.a );
  LALFree( waveform.f );
  LALFree( waveform.phi );

  /* waveform parameters needed for findchirp filter */
  tmplt->approximant = params->approximant;
  tmplt->tC = ppnParams.tc;
  tmplt->fFinal = ppnParams.fStop;
  tmplt->ACTDdominantSwitch = fcTmplt->tmplt.ACTDdominantSwitch;
  tmplt->ACTDconstraintSwitch = fcTmplt->tmplt.ACTDconstraintSwitch;

  fcTmplt->tmpltNorm = params->dynRange / ( cannonDist * 1.0e6 * LAL_PC_SI );
  fcTmplt->tmpltNorm *= fcTmplt->tmpltNorm;


  /* 
   *
   * Loop over each template and apply tapering/band passing etc
   *
   */

  /* Create a temporary vector to avoid mishaps */
  if( ( tmpACTDVec = XLALCreateREAL4Vector( numPoints ) ) == NULL )
  {
    ABORTXLAL( status );
  }
  
  for( i = istart; i < istop; i++ )
  {
    memcpy( tmpACTDVec->data, ACTDVecs[i].data,
               numPoints * sizeof( *( ACTDVecs[i].data ) ) );

    /* Taper the waveform if required */
    if ( params->taperTmplt != INSPIRAL_TAPER_NONE )
    {
      if ( XLALInspiralWaveTaper( tmpACTDVec, params->taperTmplt )
             == XLAL_FAILURE )
      {
        ABORTXLAL( status );
      }
    }


    /* Find the end of the chirp */
    j = numPoints - 1;
    while ( tmpACTDVec->data[j] == 0 )
    {
      /* search for the end of the chirp but don't fall off the array */
      if ( --j == 0 )
      { 
        ABORT( status, FINDCHIRPTDH_EEMTY, FINDCHIRPTDH_MSGEEMTY );
      }
    }
    ++j;

    if ( params->bandPassTmplt )
    {
      REAL4Vector bpVector; /*Used to save time */
      REAL4 bpfLow = 0.0;
      REAL4 bpfFinal = 0.0;

      /* We want to shift the template to the middle of the vector so */
      /* that band-passing will work properly */
      shift = ( numPoints - j ) / 2;
      memmove( tmpACTDVec->data + shift, tmpACTDVec->data,
                                       j * sizeof( *( tmpACTDVec->data ) ) );

      memset( tmpACTDVec->data, 0, shift * sizeof( *( tmpACTDVec->data ) ) );
      memset( tmpACTDVec->data + ( numPoints + j ) / 2, 0,
    ( numPoints - ( numPoints + j ) / 2 ) * sizeof( *( tmpACTDVec->data ) ) );

      /* Select an appropriate part of the vector to band pass. */
      /* band passing the whole thing takes a lot of time */
      if ( j > 2 * sampleRate && 2 * j <= numPoints )
      {
        bpVector.length = 2 * j;
        bpVector.data   = tmpACTDVec->data + numPoints / 2 - j;
      }
      else if ( j <= 2 * sampleRate && j + 2 * sampleRate <= numPoints )
      {
        bpVector.length = j + 2 * sampleRate;
        bpVector.data   = tmpACTDVec->data
                   + ( numPoints - j ) / 2 - (INT4)sampleRate;
      }
      else
      {
        bpVector.length = numPoints;
        bpVector.data   = tmpACTDVec->data;
      }

      /* Adjust frequencies according to the harmonic */
      bpfLow = 0.98 * tmplt->fLower ;
      bpfFinal = 1.02 * tmplt->fFinal * ( ( REAL4 )( i ) + 1. ) / 2.;


      if ( XLALBandPassInspiralTemplate( &bpVector, bpfLow,
                   tmplt->fFinal, sampleRate ) == XLAL_FAILURE )
      {
        ABORTXLAL( status );
      }

      /* Now we need to do the shift to the end. */
      memcpy( tmpACTDVec->data, tmpACTDVec->data + ( numPoints + j ) / 2,
         ( numPoints - ( numPoints + j ) / 2 )
               * sizeof( *(tmpACTDVec->data) ) );
      memcpy( tmpACTDVec->data + numPoints - ( numPoints + j ) / 2,
         tmpACTDVec->data, ( numPoints + j ) /2
               * sizeof( *( tmpACTDVec->data ) ) );

    }
    else
    {
      /* No need for so much shifting around if not band passing */
      /* shift chirp to end of vector so it is the correct place for filter */
        memmove( tmpACTDVec->data + numPoints - j, tmpACTDVec->data,
                                        j * sizeof( *( tmpACTDVec->data ) ) );
        memset( tmpACTDVec->data, 0,
                        ( numPoints - j ) * sizeof( *( tmpACTDVec->data ) ) );
    }

    memcpy( ACTDVecs[i].data, tmpACTDVec->data,
                           numPoints * sizeof( *( tmpACTDVec->data ) ) );

  }

  XLALDestroyREAL4Vector( tmpACTDVec );
  tmpACTDVec = NULL;
 


  /*
   *
   * create the frequency domain findchirp templates
   *
   */

  /* fft harmonics */
  for( i = istart; i < istop; ++i)
  {
    if ( XLALREAL4ForwardFFT( &ACTDtilde[i], &ACTDVecs[i],
         params->fwdPlan ) == XLAL_FAILURE )
    {
      ABORTXLAL( status );
    }
    ACTDtilde[i].data[0].re = 0.0;
    ACTDtilde[i].data[0].im = 0.0;
    ACTDtilde[i].data[numPoints / 2].re = 0.0;
    ACTDtilde[i].data[numPoints / 2].im = 0.0;
  }

  /* copy the template parameters to the findchirp template structure */
  memcpy( &(fcTmplt->tmplt), tmplt, sizeof(InspiralTemplate) );

  /* print the template normalization constant */
  if ( lalDebugLevel & LALINFO )
  {
    snprintf( infomsg, sizeof(infomsg) / sizeof(*infomsg),
        "tmpltNorm = %e\n", fcTmplt->tmpltNorm );
    LALInfo( status, infomsg );
  }


  /* normal exit */
  DETATCHSTATUSPTR( status );
  RETURN( status );
}


/* <lalVerbatim file="FindChirpACTDTemplate"> */
void
LALFindChirpACTDNormalize(
    LALStatus                  *status,
    FindChirpTemplate          *fcTmplt,
    FindChirpTmpltParams       *tmpltParams,
    FindChirpDataParams        *params
    )
/* </lalVerbatim> */
{
  UINT4           i, j, k;
  UINT4           numPoints;
  UINT4           numTDPoints;
  REAL8           deltaT; 
  COMPLEX8Vector  ACTDtilde[NACTDTILDEVECS];
  COMPLEX8Vector *tmpVec[NACTDTILDEVECS];
  COMPLEX8       *wtilde;
  REAL4           norm;
  REAL4           min12, max12;
  REAL4           min13, max13;
  REAL4           min23, max23;

  REAL4           cons1;
  REAL4           cons2;
  REAL4           cons3;

  /* Variables for GSL */
  UINT4       matrixDim = 0;
  gsl_matrix *innerProd = NULL;
  gsl_matrix *eigenVect = NULL;
  gsl_vector *eigenVal  = NULL;

 


  gsl_eigen_symmv_workspace *workspace = NULL;

  INITSTATUS( status, "LALFindChirpACTDNormalize", FINDCHIRPACTDTEMPLATEC );

  /* check the required input exists */
  ASSERT( fcTmplt, status,
      FINDCHIRPTDH_ENULL, FINDCHIRPTDH_MSGENULL );
  ASSERT( tmpltParams, status,
      FINDCHIRPTDH_ENULL, FINDCHIRPTDH_MSGENULL );

  ASSERT( params, status,
      FINDCHIRPTDH_ENULL, FINDCHIRPTDH_MSGENULL );

  ASSERT( params->wtildeVec, status,
      FINDCHIRPTDH_ENULL, FINDCHIRPTDH_MSGENULL );
  ASSERT( params->wtildeVec->data, status,
      FINDCHIRPTDH_ENULL, FINDCHIRPTDH_MSGENULL );

  ASSERT( params->tmpltPowerVec, status,
      FINDCHIRPTDH_ENULL, FINDCHIRPTDH_MSGENULL );
  ASSERT( params->tmpltPowerVec->data, status,
      FINDCHIRPTDH_ENULL, FINDCHIRPTDH_MSGENULL );

  /* check that we have the correct approximant */
  if( params->approximant != AmpCorPPN )
  {
    ABORT( status, FINDCHIRPTDH_EMAPX, FINDCHIRPTDH_MSGEMAPX );
  }


  wtilde = params->wtildeVec->data;

  numPoints   = fcTmplt->ACTDtilde->vectorLength;
  numTDPoints = 2 * ( numPoints - 1 );

  deltaT = tmpltParams->deltaT;

  matrixDim = 2 * ( fcTmplt->stopVecACTD - fcTmplt->startVecACTD );

  /* Allocate memory for transformation matrices and workspace */
  innerProd = gsl_matrix_alloc( matrixDim, matrixDim );
  eigenVect = gsl_matrix_alloc( matrixDim, matrixDim );
  eigenVal  = gsl_vector_alloc( matrixDim );
  workspace = gsl_eigen_symmv_alloc( matrixDim );


  for( i = 0; i < NACTDTILDEVECS; ++i )
  {
    ACTDtilde[i].length = numPoints;
    ACTDtilde[i].data  = fcTmplt->ACTDtilde->data + (i * numPoints );
  }
  /* Do a little jiggery-pokery so that the dominant 
   * harmonic is the first index.
   * The third harmonic has more power in band so that goes in to 
   * second place. 
   */
/*
  {
    COMPLEX8 *tmp = ACTDtilde[0].data;
    ACTDtilde[0].data = ACTDtilde[1].data;
    ACTDtilde[1].data = tmp;

    tmp = ACTDtilde[1].data;
    ACTDtilde[1].data = ACTDtilde[3].data;
    ACTDtilde[3].data = tmp;

    tmp = ACTDtilde[3].data;
    ACTDtilde[3].data = ACTDtilde[5].data;
    ACTDtilde[5].data = tmp;

    tmp = ACTDtilde[4].data;
    ACTDtilde[4].data = ACTDtilde[5].data;
    ACTDtilde[5].data = tmp;
  }
*/


  /* Norm the templates before we start and use for constraint */
  memset( fcTmplt->ACTDconstraint->data, 0.0, NACTDVECS * sizeof( REAL4 ) );
  for( i = 0; i < NACTDVECS; ++i ) 
  {
    if ( i >= fcTmplt->startVecACTD && i < fcTmplt->stopVecACTD )
    {
      norm =  XLALFindChirpACTDInnerProduct( &ACTDtilde[i], &ACTDtilde[i],
                wtilde, tmpltParams->fLow, deltaT, numTDPoints );
      if( norm != 0.0 )
      {
        norm = sqrt( norm );
      }
      fcTmplt->ACTDconstraint->data[i] = norm;
      for( k = 0; k < numPoints; k++ )
      {
        if( norm != 0.0 )
        {
          ACTDtilde[i].data[k].re /= norm;
          ACTDtilde[i].data[k].im /= norm;
        }          
        /* Set the other quadrature */
        ACTDtilde[i+3].data[k].re = - ACTDtilde[i].data[k].im;
        ACTDtilde[i+3].data[k].im = ACTDtilde[i].data[k].re;   
      }
    }
  }

  /* CALCULATE THE CONSTRAINT */
  if( fcTmplt->tmplt.ACTDconstraintSwitch )
  {
    /* h1 and h2 */
    XLALFindChirpACTDMiniMax( &ACTDtilde[0], &ACTDtilde[3], &ACTDtilde[1], 
            &ACTDtilde[4], wtilde, tmpltParams->fLow, deltaT, &min12, &max12 );
    /* h1 and h3 */
    XLALFindChirpACTDMiniMax( &ACTDtilde[0], &ACTDtilde[3], &ACTDtilde[2], 
            &ACTDtilde[5], wtilde, tmpltParams->fLow, deltaT, &min13, &max13 );
    /* h2 and h3 */
    XLALFindChirpACTDMiniMax( &ACTDtilde[1], &ACTDtilde[4], &ACTDtilde[2], 
            &ACTDtilde[5], wtilde, tmpltParams->fLow, deltaT, &min23, &max23 );

    cons1 = fcTmplt->ACTDconstraint->data[0]
          + max12 * fcTmplt->ACTDconstraint->data[1]
          + max13 * fcTmplt->ACTDconstraint->data[2];

    cons2 = fcTmplt->ACTDconstraint->data[1] 
          + max23 * fcTmplt->ACTDconstraint->data[2]
          + max12 * fcTmplt->ACTDconstraint->data[0];

    cons3 = fcTmplt->ACTDconstraint->data[2] 
          + max13 * fcTmplt->ACTDconstraint->data[0]
          + max23 * fcTmplt->ACTDconstraint->data[1];

    if ( cons2 == 0.0 )
    {
      XLALPrintError( "Second harmonic out of band! Not yet catered for.\n" );
      ABORT( status, FINDCHIRPH_EMASS, FINDCHIRPH_MSGEMASS );
    }

    /* Now repopulate the vector with the appropriate values */
    fcTmplt->ACTDconstraint->data[0] = cons1 / cons2;
    fcTmplt->ACTDconstraint->data[1] = 1.0;
    fcTmplt->ACTDconstraint->data[2] = cons3 / cons2;
  }

  /* Fill the inner product matrix */
  for ( i = 0; i < NACTDTILDEVECS; i++ )
  {
    for ( k = i; k < NACTDTILDEVECS; k++ )
    {
      if ( ( i % NACTDVECS - fcTmplt->startVecACTD < matrixDim / 2 )
         && ( k % NACTDVECS - fcTmplt->startVecACTD < matrixDim / 2 ) )
      {
        UINT4 idx1, idx2, nVec;
         
        nVec = matrixDim / 2;
        idx1 = (i / NACTDVECS) * nVec + (i % NACTDVECS) - fcTmplt->startVecACTD;
        idx2 = (k / NACTDVECS) * nVec + (k % NACTDVECS) - fcTmplt->startVecACTD;

  
        gsl_matrix_set( innerProd, idx1, idx2, 
         (double)XLALFindChirpACTDInnerProduct( &ACTDtilde[i], &ACTDtilde[k],
                             wtilde, tmpltParams->fLow, deltaT, numTDPoints ));

        if ( idx1 != idx2 )
        {
          gsl_matrix_set( innerProd, idx2, idx1, 
               gsl_matrix_get( innerProd, idx1, idx2 ));
        }
      }
    }
  }
/*
  fprintf(stderr,"\n");
  for ( i = 0; i < 2 * NACTDVECS; i++ )
  {
    for ( k = 0; k < 2 * NACTDVECS; k++ )
    {
      if( gsl_matrix_get( innerProd, i, k ) >= 0.0 )
        fprintf(stderr, " ");
      fprintf(stderr," %.2e", gsl_matrix_get( innerProd, i, k ) );
    }
    fprintf(stderr,"\n");
  }
*/

  /* Diagonalize the matrix */
  gsl_eigen_symmv( innerProd, eigenVal, eigenVect, workspace );

  /* Now we perform the co-ordinate transformation */
  for ( i = 0; i < NACTDTILDEVECS; i++ )
  {
    tmpVec[i] = XLALCreateCOMPLEX8Vector( numPoints );
    memset( tmpVec[i]->data, 0.0, numPoints * sizeof( COMPLEX8) );
  }

  for ( i = 0; i < NACTDTILDEVECS; i++ )
  {
    for ( j = 0; j < NACTDTILDEVECS; j++ )
    {
      if ( ( i % NACTDVECS - fcTmplt->startVecACTD < matrixDim / 2 )
         && ( j % NACTDVECS - fcTmplt->startVecACTD < matrixDim / 2 ) )
      {
        UINT4 idx1, idx2, nVec;

        nVec = matrixDim / 2;
        idx1 = (i / NACTDVECS) * nVec + i % NACTDVECS - fcTmplt->startVecACTD;
        idx2 = (j / NACTDVECS) * nVec + j % NACTDVECS - fcTmplt->startVecACTD;

        for ( k = 0; k < numPoints; k++ )
        {
          tmpVec[i]->data[k].re += 
              (REAL4)gsl_matrix_get( eigenVect, idx2, idx1) * ACTDtilde[j].data[k].re;
          tmpVec[i]->data[k].im += 
              (REAL4)gsl_matrix_get( eigenVect, idx2, idx1) * ACTDtilde[j].data[k].im;
        }
      }
    }
  }

  /* Copy in and normalize */
  for ( i = 0; i < NACTDTILDEVECS; i++ )
  {
    memcpy( ACTDtilde[i].data, tmpVec[i]->data, numPoints * sizeof( COMPLEX8 ));
    XLALDestroyCOMPLEX8Vector( tmpVec[i] );
    tmpVec[i] = NULL;

    if ( i % NACTDVECS - fcTmplt->startVecACTD < matrixDim / 2 )
    {
      UINT4 idx1, nVec;
      nVec = matrixDim / 2;
      idx1 = (i / NACTDVECS) * nVec + i % NACTDVECS - fcTmplt->startVecACTD;

      norm = sqrt( gsl_vector_get( eigenVal, idx1 ) );
      for ( k = 0; k < numPoints; k++ )
      {
        ACTDtilde[i].data[k].re /= norm;
        ACTDtilde[i].data[k].im /= norm;
      }
    }
  }


#if 0
  /* Gram-Schmidt works */
  for ( i = 0; i < NACTDTILDEVECS; ++i )
  {
    for ( j = 0; j < i; ++j )
    {
      norm = 0.0;

      REAL4 innerProd = 
                XLALFindChirpACTDInnerProduct( &ACTDtilde[i], &ACTDtilde[j],
                             wtilde, tmpltParams->fLow, deltaT, numTDPoints );

      if( innerProd != 0.0 )
      {
        for ( k = 0; k < numPoints; ++k )
        {
          ACTDtilde[i].data[k].re -= innerProd * ACTDtilde[j].data[k].re;
          ACTDtilde[i].data[k].im -= innerProd * ACTDtilde[j].data[k].im;
        }
      
        /* Now re-norm the vector */
        norm = XLALFindChirpACTDInnerProduct( &ACTDtilde[i], &ACTDtilde[i],
                wtilde, tmpltParams->fLow, deltaT, numTDPoints );
      }
 
      if( norm != 0.0 )
      {
        norm = sqrt(norm);
        for (k = 0; k < numPoints; k++)
        {
          ACTDtilde[i].data[k].re /= norm;
          ACTDtilde[i].data[k].im /= norm;
        }
      }  
    }
  }
#endif
  /* XXX UNCOMMENT BELOW TO TEST ORTHONORMALISATION XXX */
  /*
  fprintf( stderr, "\n\n NORMALIZATION TEST:\n    ");
  for ( i = 0; i < NACTDTILDEVECS; i++ )
  {
    for ( j= 0; j < NACTDTILDEVECS; j++ )
    {
      norm = XLALFindChirpACTDInnerProduct( &ACTDtilde[i], &ACTDtilde[j],
                              wtilde, tmpltParams->fLow, deltaT, numTDPoints );
      fprintf( stderr, "H%dH%d=%.2f ", i, j, fabs(norm) );
    }
    fprintf( stderr, "\n    ");
  }
  fprintf( stderr, "                                    ");
  */
  /* XXX UNCOMMENT ABOVE TO TEST ORTHONORMALIZATION XXX */

  fcTmplt->norm = 2.0 * deltaT / (REAL4)(numTDPoints);
  fcTmplt->norm = fcTmplt->norm * fcTmplt->norm;

  /* Set the transformation matrix to be used to constrain the filter */
  fcTmplt->ACTDconmatrix = eigenVect;

  /* Free memory */
  gsl_matrix_free( innerProd );
  /* gsl_matrix_free( eigenVect ); */
  gsl_vector_free( eigenVal );
  gsl_eigen_symmv_free( workspace );


  /* normal exit */
  RETURN( status );

}

REAL4  XLALFindChirpACTDInnerProduct(
               COMPLEX8Vector  *a,
               COMPLEX8Vector  *b,
               COMPLEX8        *wtilde,
               REAL4            flower,
               REAL4            deltaT,
               UINT4            numPoints
                                    )
{
  INT4  k;
  INT4  cut;
  REAL4 deltaF;
  REAL4 inp = 0.0;
  REAL4 power;

  if ( a->length != b->length )
  {
     XLALPrintError( "Lengths of the vectors do not agree\n" );
     XLAL_ERROR_REAL8( "XLALFindChirpACTDInnerProduct", XLAL_EINVAL );
  }
     

  deltaF = 1.0 / ((REAL4)numPoints * deltaT);
  cut = flower / deltaF > 1 ? flower / deltaF : 1;

  for( k = cut; k < (INT4)a->length-1; ++k )
  {
    power = a->data[k].re * b->data[k].re;
    power += a->data[k].im * b->data[k].im;

    if( isnan( power ) )
    {
      fprintf(stderr, "\nk=%d\ta.re=%.3e\ta.im=%.3e\tb.re=%.3e\tb.im=%.3e\n",
                k, a->data[k].re, b->data[k].re,
                a->data[k].im, b->data[k].im );
      exit( 1 );
    }
    inp += 4.0 * deltaT *  power * wtilde[k].re / (REAL4)(numPoints);
  }

  return inp;
}

INT4 XLALFindChirpACTDCorrelate(
               COMPLEX8Vector  *a,
               COMPLEX8Vector  *b,
               REAL4Vector     *out,
               COMPLEX8        *wtilde,
               REAL4            flower,
               REAL4            deltaT,
               UINT4            numPoints
                               )
{
  INT4  k;
  INT4  cut;
  REAL4 deltaF;
  COMPLEX8Vector *cor = NULL;
  REAL4FFTPlan *revPlan = NULL;

  if ( a->length != b->length )
  {
     XLALPrintError( "Lengths of the vectors do not agree\n" );
     XLAL_ERROR_REAL8( "XLALFindChirpACTDInnerProduct", XLAL_EINVAL );
  }
  if ( a->length != out->length/2+1 )
  {
     XLALPrintError( "Lengths of the vectors do not agree\n" );
     XLAL_ERROR_REAL8( "XLALFindChirpACTDInnerProduct", XLAL_EINVAL );
  }
  
  cor = XLALCreateCOMPLEX8Vector( numPoints/2+1 );
  memset( cor->data, 0.0, (numPoints/2+1) * sizeof( COMPLEX8 ) );
     
  deltaF = 1.0 / ((REAL4)numPoints * deltaT);
  cut = flower / deltaF > 1 ? flower / deltaF : 1;

  for( k = cut; k < (INT4)a->length-1; ++k )
  {
    cor->data[k].re  = a->data[k].re * b->data[k].re;
    cor->data[k].re += a->data[k].im * b->data[k].im;
 
    cor->data[k].im  = a->data[k].re * b->data[k].im;
    cor->data[k].im -= a->data[k].im * b->data[k].re;

    cor->data[k].re *= 2.0 * deltaT * wtilde[k].re / (REAL4)(numPoints);
    cor->data[k].im *= 2.0 * deltaT * wtilde[k].re / (REAL4)(numPoints);
  }
  
  revPlan = XLALCreateReverseREAL4FFTPlan( numPoints, 0 );
  if( !revPlan )
  {
    exit( 1 );
  }
  

  XLALREAL4ReverseFFT( out, cor, revPlan );
  for( k = cut; k < (INT4)(out->length); k++)
  {
    out->data[k] /= out->length;
  }

  XLALDestroyREAL4FFTPlan( revPlan );
  XLALDestroyCOMPLEX8Vector( cor );

  return XLAL_SUCCESS;
}

INT4 XLALFindChirpACTDMiniMax(
               COMPLEX8Vector  *a1,
               COMPLEX8Vector  *a2,
               COMPLEX8Vector  *b1,
               COMPLEX8Vector  *b2,
               COMPLEX8        *wtilde,
               REAL4           fLow,
               REAL4           deltaT,
               REAL4           *min,
               REAL4           *max
                            )
{
  INT4 k, numTDPoints;
  REAL4Vector     *x11 = NULL, *x12 = NULL, *x21 = NULL, *x22 = NULL;
  REAL4 a, b, c, sum, diff, d11, d12, d21, d22;
  REAL4 maxloc, minloc;  

  numTDPoints = 2 * ( a1->length - 1 );

  x11 = XLALCreateREAL4Vector( numTDPoints );
  x12 = XLALCreateREAL4Vector( numTDPoints );
  x21 = XLALCreateREAL4Vector( numTDPoints );
  x22 = XLALCreateREAL4Vector( numTDPoints );
  memset( x11->data, 0.0, numTDPoints * sizeof( REAL4 ) );
  memset( x12->data, 0.0, numTDPoints * sizeof( REAL4 ) );
  memset( x21->data, 0.0, numTDPoints * sizeof( REAL4 ) );
  memset( x22->data, 0.0, numTDPoints * sizeof( REAL4 ) );
  
  XLALFindChirpACTDCorrelate( a1, b1, x11,
          wtilde, fLow, deltaT, numTDPoints );
  XLALFindChirpACTDCorrelate( a1, b2, x12,
          wtilde, fLow, deltaT, numTDPoints );
  XLALFindChirpACTDCorrelate( a2, b1, x21,
          wtilde, fLow, deltaT, numTDPoints );
  XLALFindChirpACTDCorrelate( a2, b2, x22,
          wtilde, fLow, deltaT, numTDPoints );

  d11 = x11->data[0];
  d12 = x12->data[0];
  d21 = x21->data[0];
  d22 = x22->data[0];
  a = d11*d11 + d12*d12;
  b = d21*d21 + d22*d22;
  c = d11*d21 + d12*d22;
  sum  = ( a + b )/2.;
  diff = ( a - b )/2.;
  diff = sqrt( diff*diff + c*c );
  *min = sqrt( sum - diff );
  *max = sqrt( sum + diff );
  for( k=1; k < numTDPoints; k++) 
  {
    d11 = x11->data[k];
    d12 = x12->data[k];
    d21 = x21->data[k];
    d22 = x22->data[k];
    /*
    fprintf( stderr, "\nx11=%.3e x12c=%.3e x21=%.3e z22=%.3e", 
              x11->data[k],
              x12->data[k],
              x21->data[k],
              x22->data[k] );
    */
    a = d11*d11 + d12*d12;
    b = d21*d21 + d22*d22;
    c = d11*d21 + d12*d22;
    sum  = ( a + b )/2.0;
    diff = ( a - b )/2.0;
    diff = sqrt( diff*diff +c*c );
    minloc = sqrt( sum - diff );
    maxloc = sqrt( sum + diff );
    if (maxloc > *max) 
    {
      *max = maxloc;
    }
    if (minloc > *min) 
    {      
      *min = minloc;
    }
  }
  XLALDestroyREAL4Vector( x11 );
  XLALDestroyREAL4Vector( x12 );
  XLALDestroyREAL4Vector( x21 );
  XLALDestroyREAL4Vector( x22 );



 return XLAL_SUCCESS;

}
