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

#include <math.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALError.h>
#include <lal/LALConstants.h>
#include <lal/Date.h>
#include <lal/AVFactories.h>
#include <lal/FindChirp.h>
#include <lal/LALInspiral.h>
#include <lal/FindChirp.h>

double rint(double x);

NRCSID (FINDCHIRPPTFFILTERC, "$Id$");

/* LAPACK function for the calculation of the eigenvalues of a NxN matrix */
int
XLALFindChirpFindEigenvalues ( REAL4* WR, REAL4* WI, REAL4* A)
{  
  INT4 errcode = 0;
  INT4 INFO = 0;
  REAL4 vr[25], vl[25], work[25];
  CHAR  n     = 'N';
  INT4  N     = 5;
  INT4  M     = 1;
  INT4  lwork = 25;

  extern void sgeev_ ( CHAR* JOBVLp, CHAR* JOBVRp, INT4* Np, REAL4* A, 
      INT4* LDAp, REAL4* WR, REAL4* WI, REAL4* VL, 
      INT4* LDVLp, REAL4* VR, INT4* LDVRp, REAL4* WORK, 
      INT4* LWORKp, INT4* INFOp);

  sgeev_ ( &n, &n, &N, A, &N, WR, WI, vl, &M, vr, &M, work, 
      &lwork, &INFO);

  if ( INFO != 0)
  {
    XLALPrintError( "Eigenvalue evaluation error: sgeev_ function filed with error %d\n", INFO);
    errcode = XLAL_EFAILED;
  }

  return errcode;
}

/* <lalVerbatim file="FindChirpPTFFilterCP"> */
void
LALFindChirpPTFFilterSegment (
    LALStatus                  *status,
    SnglInspiralTable         **eventList,
    FindChirpFilterInput       *input,
    FindChirpFilterParams      *params
    )
/* </lalVerbatim> */
{
  UINT4                 i, j, k, l, kmax;
  UINT4                 errcode;
  UINT4                 numPoints;
  UINT4                 deltaEventIndex;
  UINT4                 ignoreIndex;
  UINT4                 haveEvent   = 0;
  REAL4                 deltaT, sum, temp, PTFMatrix[25], r, s, x, y;
  REAL4                 deltaF, fFinal;
  REAL4                 snrThresh      = 0;
  REAL4                *snr            = NULL;
  COMPLEX8             *PTFQtilde, *qtilde, *PTFq, *inputData;
  COMPLEX8Vector        qVec;

  /* Variables needed for the eigenvalues finding LAPACK routine */
  REAL4 wr[5], wi[5];

  /*
   *
   * point local pointers to input and output pointers
   *
   */

  /* number of points in a segment */
  numPoints = params->PTFqVec->vectorLength;
  /* workspace vectors */
  snr = params->PTFsnrVec->data;
  qtilde = params->qtildeVec->data;
  PTFq   = params->PTFqVec->data;
  qVec.length = numPoints;

  /* template and data */
  inputData = input->segment->data->data->data;
  PTFQtilde = input->fcTmplt->PTFQtilde->data;

  /* number of points and frequency cutoffs */
  deltaT = (REAL4) params->deltaT;
  deltaF = 1.0 / ( deltaT * (REAL4) numPoints );
  fFinal = (REAL4) input->fcTmplt->tmplt.fFinal;
  kmax =  fFinal / deltaF < numPoints/2 ? fFinal / deltaF : numPoints/2;

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

  /* check that the fft plan exists */
  ASSERT( params->invPlan, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );

  /* check that the workspace vectors exist */

  /* if a rhosqVec vector has been created, check we can store data in it */
  if ( params->rhosqVec ) 
  {
    ASSERT( params->rhosqVec->data->data, status, 
        FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
    ASSERT( params->rhosqVec->data, status, 
        FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  }

  /* make sure that the input structure exists */
  ASSERT( input, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );

  /* make sure that the input structure contains some input */
  ASSERT( input->fcTmplt, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT( input->segment, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );

  /* make sure that the filter has been niitialized for the correct */
  /* approximant                                                    */
  if ( params->approximant != FindChirpPTF )
  {
    ABORT( status, FINDCHIRPH_EAPRX, FINDCHIRPH_MSGEAPRX );
  }

  /* make sure the approximant in the tmplt and segment agree */
  ASSERT( input->fcTmplt->tmplt.approximant == FindChirpPTF, status,
      FINDCHIRPH_EAPRX, FINDCHIRPH_MSGEAPRX );
  ASSERT( input->segment->approximant == FindChirpPTF, status,
      FINDCHIRPH_EAPRX, FINDCHIRPH_MSGEAPRX );

  fprintf(stderr,"LALFindChirpPTFFilterSegment called\n");

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
  ignoreIndex = ( input->segment->invSpecTrunc / 2 ) + deltaEventIndex;

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
  ignoreIndex = numPoints / 4;

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
  memset( params->PTFsnrVec->data, 0, numPoints * sizeof(REAL4) );
  memset( params->PTFqVec->data, 0, 5 * numPoints * sizeof(COMPLEX8) );
  
  for ( i = 0; i < 5; ++i )
  {

    /* compute qtilde using data and Qtilde */

    memset( params->qtildeVec->data, 0, 
        params->qtildeVec->length * sizeof(COMPLEX8) );

    /* qtilde positive frequency, not DC or nyquist */
    for ( k = 1; k < kmax; ++k )
    {
      r = inputData[k].re;
      s = inputData[k].im;
      x = PTFQtilde[i * (numPoints / 2 + 1) + k].re;
      y = 0 - PTFQtilde[i * (numPoints / 2 + 1) + k].im; /* cplx conj */

      qtilde[k].re = r*x - s*y;
      qtilde[k].im = r*y + s*x;
    }

    qVec.data = params->PTFqVec->data + (i * numPoints);

    /* inverse fft to get q */
    LALCOMPLEX8VectorFFT( status->statusPtr, &qVec, params->qtildeVec, 
        params->invPlan );
    CHECKSTATUSPTR( status );
  }

  /* now we have PTFqVec which contains <s|Q^I_0> + i <s|Q^I_\pi/2> */


  for ( k = 0; k < numPoints; ++k ) /* beginning of main loop over time */
  {  
    /* Set PTFMatrxi elements to zero */
    memset(params->PTFA->data, 0 , 25 * sizeof(REAL4));
    memset(params->PTFMatrix->data, 0, 25 * sizeof(REAL4) );
    /* construct A */
    
    for ( i = 0; i < 5; ++i )
    {  
      for ( j = 0; j < i + 1; ++j )
      { 
        params->PTFA->data[5 * i + j] = PTFq[i * numPoints + k].re * 
                                        PTFq[j * numPoints + k].re +
                                        PTFq[i * numPoints + k].im * 
                                        PTFq[j * numPoints + k].im;
        params->PTFA->data[5 * j + i] = params->PTFA->data[ 5 * i + j]; 
      }  
    } 
    
    /* multiply by PTFBinverse to obtain AB^(-1) */

    LALSMatrixMultiply(status->statusPtr, 
        params->PTFMatrix, params->PTFA, 
        input->fcTmplt->PTFBinverse);
    CHECKSTATUSPTR( status );

    /* Transpose PTFMatrix and store it into the corresponding local variable:
     * the input to the LAPACK eigenvalues finding routine must be in 
     * Fortran column-wise format
     */ 

    for ( i = 0; i < 5; ++i ) 
    {
      for ( j = 0; j < 5; ++j )
      {  
        PTFMatrix[i + 5 * j] = params->PTFMatrix->data[j + 5 * i];
      }
    }  

    /* find max eigenvalue and store it in snr vector */
    errcode = XLALFindChirpFindEigenvalues( wr, wi, PTFMatrix);

    if ( errcode != XLAL_SUCCESS )
    {
      ABORT( status, FINDCHIRPH_EIGEN, FINDCHIRPH_MSGEIGEN );
    }

    if (wi[0] == 0) 
      snr[k] = wr[0];
    else
      snr[k] = 0.0;

    for ( i = 1; i < 5; ++i )
    {
      if ( wi[i] == 0) 
      {  
        if ( (temp = wr[i]) > snr[k] ) 
          snr[k] = temp;                
      } 
    }
  } /* End of main loop over time */


  fprintf( stderr, "Ptf filtering data segment\n" );

  /*
   *
   * look for and cluster events in the snr vector
   *
   */


  /* look for an event in the filter output */
  for ( j = ignoreIndex; j < numPoints - ignoreIndex; ++j )
  {
    /* if snrsq exceeds threshold at any point */
    if ( snr[j] > snrThresh )
    {
      haveEvent = 1;        /* mark segment to have events    */
      break;
    }
  }

  /* search the SNR vector for events */
  /* process events in the filter output */
  if ( haveEvent )
  {
#if 0
    LALFindChirpClusterEvents( status->statusPtr, eventList, input,
        params, q, kmax, numPoints, ignoreIndex, 
        norm, modqsqThresh, chisqThreshFac, numChisqBins, searchName );
    CHECKSTATUSPTR( status );
#endif
    fprintf( stderr, "found events!\n" );
  }

  /* normal exit */
  DETATCHSTATUSPTR( status );
  RETURN( status );
}
