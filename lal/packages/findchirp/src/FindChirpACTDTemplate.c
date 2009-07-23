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
#include <lal/DataBuffer.h>
#include <lal/LALInspiral.h>
#include <lal/FindChirp.h>
#include <lal/FindChirpTD.h>
#include <lal/FindChirpACTD.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_eigen.h>

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

  /* check that the template masses are not equal */
  if( tmplt->mass1 == tmplt->mass2 )
  {
    tmplt->mass2 += 1e-05;
    /*
    ABORT( status, FINDCHIRPACTDH_EQMAS, FINDCHIRPACTDH_MSGEQMAS );
    */
  }

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
  ppnParams.inc = LAL_PI_4;
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
    ACTDtilde[i+NACTDVECS].data = fcTmplt->ACTDtilde->data + ((i+NACTDVECS) * (numPoints / 2 + 1 ) );
  }

  /* compute h(t) */
  /* legacy - length is the lentgh of the vectors */
  for ( j = 0; j < waveform.a->data->length; ++j )
  {
    for ( i = 0; i < NACTDVECS; ++i )
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
  
  for( i = 0; i < NACTDVECS; i++ )
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
      

      /* If the harmonic is not in our bandwidth we mega band pass it */
      if( bpfLow >= bpfFinal )
      {
        bpfLow = 0.98 * tmplt->fLower;
        bpfFinal = 1.02* ( tmplt->fLower + 1. );
      }
      else
      {
        /* Adjust frequencies according to the harmonic */
        bpfLow = 0.98 * tmplt->fLower ;
        bpfFinal = 1.02 * tmplt->fFinal * ( ( REAL4 )( i ) + 1. ) / 2.;
      }

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
  for( i = 0; i < NACTDVECS; ++i)
  {
    if ( XLALREAL4ForwardFFT( &ACTDtilde[i], &ACTDVecs[i], 
         params->fwdPlan ) == XLAL_FAILURE )
    {
      ABORTXLAL( status );
    }
  }

  /* Create the tmplt for the other quadrature */
  for( i = 0; i < NACTDVECS; ++i)
  {
    for ( j = 0; j < ACTDtilde[i].length; j++ )
    {
      ACTDtilde[i+NACTDVECS].data[j].re = - ACTDtilde[i].data[j].im;
      ACTDtilde[i+NACTDVECS].data[j].im =  ACTDtilde[i].data[j].re;
    }
  }

  /* copy the template parameters to the findchirp template structure */
  memcpy( &(fcTmplt->tmplt), tmplt, sizeof(InspiralTemplate) );

  /* print the template normalization constant */
  if ( lalDebugLevel & LALINFO )
  {
    LALSnprintf( infomsg, sizeof(infomsg) / sizeof(*infomsg), 
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

  /* gsl variables */
  gsl_matrix *innerProd = NULL;
  gsl_matrix *eigenVect = NULL;
  gsl_vector *eigenVal  = NULL;

  /*
  gsl_matrix *transpose = NULL;
  gsl_matrix *output    = NULL;
  */
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

  innerProd = gsl_matrix_alloc( NACTDTILDEVECS, NACTDTILDEVECS );
  eigenVect = gsl_matrix_alloc( NACTDTILDEVECS, NACTDTILDEVECS );
  /*
  transpose = gsl_matrix_alloc( 2 * NACTDVECS, 2 * NACTDVECS );
  output    = gsl_matrix_alloc( 2 * NACTDVECS, 2 * NACTDVECS );
  */
  eigenVal  = gsl_vector_alloc( NACTDTILDEVECS );
  workspace = gsl_eigen_symmv_alloc( NACTDTILDEVECS );

  for( i = 0; i < NACTDTILDEVECS; ++i )
  {
    ACTDtilde[i].length = numPoints;
    ACTDtilde[i].data  = fcTmplt->ACTDtilde->data + (i * numPoints );
  }

  /* Norm the templates before we start */
  for( i = 0; i < NACTDTILDEVECS; ++i ) 
  {
    norm =  XLALFindChirpACTDInnerProduct( &ACTDtilde[i], &ACTDtilde[i],
                wtilde, tmpltParams->fLow, deltaT, numTDPoints );
    norm = sqrt(norm);

    for (k = 0; k < numPoints; k++)
    {
      ACTDtilde[i].data[k].re /= norm;
      ACTDtilde[i].data[k].im /= norm;
    }
  }

  /* Fill the inner product matrix */
  for ( i = 0; i < NACTDVECS; i++ )
  {
    for ( k = i; k < NACTDVECS; k++ )
    {
   
      gsl_matrix_set( innerProd, i, k, 
         (double)XLALFindChirpACTDInnerProduct( &ACTDtilde[i], &ACTDtilde[k],
                             wtilde, tmpltParams->fLow, deltaT, numTDPoints ));
      gsl_matrix_set( innerProd, i+NACTDVECS, k, 
         (double)XLALFindChirpACTDInnerProduct( &ACTDtilde[i+NACTDVECS], &ACTDtilde[k],
                             wtilde, tmpltParams->fLow, deltaT, numTDPoints ));
      gsl_matrix_set( innerProd, i, k+NACTDVECS, 
         (double)XLALFindChirpACTDInnerProduct( &ACTDtilde[i], &ACTDtilde[k+NACTDVECS],
                             wtilde, tmpltParams->fLow, deltaT, numTDPoints ));
      gsl_matrix_set( innerProd, i+NACTDVECS, k+NACTDVECS,
         (double)XLALFindChirpACTDInnerProduct( &ACTDtilde[i+NACTDVECS], &ACTDtilde[k+NACTDVECS],
                             wtilde, tmpltParams->fLow, deltaT, numTDPoints ));

      gsl_matrix_set( innerProd, k, i, gsl_matrix_get( innerProd, i, k ));
      gsl_matrix_set( innerProd, k, i+NACTDVECS, 
                                  gsl_matrix_get( innerProd, i+NACTDVECS, k ));
      gsl_matrix_set( innerProd, k+NACTDVECS, i, 
                                  gsl_matrix_get( innerProd, i, k+NACTDVECS ));
      gsl_matrix_set( innerProd, k+NACTDVECS, i+NACTDVECS, 
                        gsl_matrix_get( innerProd, i+NACTDVECS, k+NACTDVECS ));
    }
  }

  /* XXX UNCOMMENT BELOW TO SEE INNER PRODUCT MATRIX XXX */
  /*  
  printf("\n");
  for ( i = 0; i < 2 * NACTDVECS; i++ )
  {
    for ( k = 0; k < 2 * NACTDVECS; k++ )
    {
      printf("\t%e", gsl_matrix_get( innerProd, i, k ) );
    }
    printf("\n");
  }
  */
  /* XXX UNCOMMENT ABOVE TO SEE INNER PRODUCT MATRIX XXX */

  /* Diagonalize the matrix */
  gsl_eigen_symmv( innerProd, eigenVal, eigenVect, workspace );

  /* XXX UNCOMMENT BELOW TO SEE EIGENVALUES/VECTORS XXX */
  /*
  printf( "EigenValues:\n");
  for ( i = 0; i < 2 * NACTDVECS; i++ )
  {
    printf("\t%e", gsl_vector_get( eigenVal, i ));
  }
  printf("\n\n");
  printf("Transformation matrix:\n");
  for ( i = 0; i < 2 * NACTDVECS; i++ )
  {
    for ( k = 0; k < 2 * NACTDVECS; k++ )
    {
      printf("\t%e", gsl_matrix_get( eigenVect, i, k ) );
    }
    printf("\n");
  }
  */
  /* XXX UNCOMMENT ABOVE TO SEE EIGENVALUES/VECTORS XXX */

#if 0
  /* Take transpose and check it is the inverse */
  gsl_matrix_transpose_memcpy( transpose, eigenVect );

  gsl_blas_dgemm( CblasNoTrans, CblasNoTrans, 1.0, 
                  transpose, eigenVect, 0.0, output );

  printf("\nA^TA:\n");
  for ( i = 0; i < 2 * NACTDVECS; i++ )
  {
    for ( k = 0; k < 2 * NACTDVECS; k++ )
    {
      printf("\t%e", gsl_matrix_get( output, i, k ) );
    }
    printf("\n");
  }
#endif   

  
  /* Now we perform the co-ordinate transformation */
  for ( i = 0; i < NACTDTILDEVECS; i++ )
  {
    tmpVec[i] = XLALCreateCOMPLEX8Vector( numPoints );
    memset( tmpVec[i]->data, 0, numPoints * sizeof( COMPLEX8) );
  }

  for ( i = 0; i < NACTDTILDEVECS; i++ )
  {
    for ( j = 0; j < NACTDTILDEVECS; j++ )
    {
      for ( k = 0; k < numPoints; k++ )
      {
        tmpVec[i]->data[k].re += 
              (REAL4)gsl_matrix_get( eigenVect, j, i) * ACTDtilde[j].data[k].re;
        tmpVec[i]->data[k].im += 
              (REAL4)gsl_matrix_get( eigenVect, j, i) * ACTDtilde[j].data[k].im;
      }
    }
  }

  /* Copy in and normalize */
  for ( i = 0; i < NACTDTILDEVECS; i++ )
  {
    memcpy( ACTDtilde[i].data, tmpVec[i]->data, numPoints * sizeof( COMPLEX8 ));
    XLALDestroyCOMPLEX8Vector( tmpVec[i] );
    norm = sqrt( gsl_vector_get( eigenVal, i ) );
    for ( k = 0; k < numPoints; k++ )
    {
      ACTDtilde[i].data[k].re /= norm;
      ACTDtilde[i].data[k].im /= norm;
    }
  }

  /* XXX UNCOMMENT BELOW TO TEST ORTHONORMALISATION XXX */
  /*
  printf( " NORMALIZATION TEST:\n\n" );
  for ( i = 0; i < NACTDTILDEVECS; i++ )
  {
    for ( k = 0; k < NACTDTILDEVECS; k++ )
    {
      norm = XLALFindChirpACTDInnerProduct( &ACTDtilde[i], &ACTDtilde[k],
                              wtilde, tmpltParams->fLow, deltaT, numTDPoints );
      printf( "H%dH%d=%.4f ", i, k, fabs(norm) );
    }
    printf("\n");
  }
  */
  /* XXX UNCOMMENT ABOVE TO TEST ORTHONORMALIZATION XXX */

  /* Since the template is now properly normalized, set norm to 1.0 */
  fcTmplt->norm = 2.0 * deltaT / (REAL4)(numPoints);
  fcTmplt->norm = fcTmplt->norm * fcTmplt->norm;

  /* Free memory */
  gsl_matrix_free( innerProd );
  gsl_matrix_free( eigenVect );
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
  REAL4 innerProduct;
  REAL4 deltaF;
  REAL4 sum = 0.0;

  deltaF = 1.0 / ((REAL4)numPoints * deltaT);

  for( k = 1; k < (INT4)a->length; ++k )
  {
    if(  k * deltaF  >= flower )
    {
      REAL4 power;
      power = a->data[k].re * b->data[k].re;
      power += a->data[k].im * b->data[k].im;
      
      sum += 4.0 * deltaT *  power * wtilde[k].re / (REAL4)(numPoints);
    }
  }

  innerProduct = sum;

  return innerProduct;

}
