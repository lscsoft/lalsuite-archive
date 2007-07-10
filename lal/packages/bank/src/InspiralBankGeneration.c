

#include<lal/LALStdlib.h>
#include<lal/LALStatusMacros.h>
#include<lal/LALInspiral.h>
#include<lal/LALInspiralBank.h>
#include<lal/LIGOMetadataTables.h>


NRCSID(INSPIRALBANKGENERATIONC, "$Id$");

void
LALInspiralBankGeneration(
     LALStatus *status,
     InspiralCoarseBankIn *input,
     SnglInspiralTable **first,
     INT4 *ntiles )
{
  InspiralTemplateList *coarseList = NULL;
  SnglInspiralTable *bank;
  INT4  cnt = 0;
  INT4  chicnt = 0;
  INT4  kappacnt = 0;
  INT4  i;
  REAL4 chi[3], kappa[4];
  
  INITSTATUS(status, "LALInspiralBankGeneration", INSPIRALBANKGENERATIONC);
  ATTATCHSTATUSPTR(status);
    
  ASSERT( input != NULL, status, LALINSPIRALBANKH_ENULL,
          LALINSPIRALBANKH_MSGENULL );
  ASSERT( *first == NULL, status, LALINSPIRALBANKH_ENULL,
          LALINSPIRALBANKH_MSGENULL );

  
  
  /* For nonspinning approximants, call LALInspiralCreateCoarseBank(). */
  switch( input->approximant )
  {
  case BCV:
  case EOB:
  case PadeT1:
  case PadeF1:
  case TaylorF1:
  case TaylorF2:
  case TaylorT1:
  case TaylorT2:
  case TaylorT3:
  case AmpCorPPN:

    /* Use LALInspiralCreateCoarseBank(). */
    TRY( LALInspiralCreateCoarseBank( status->statusPtr, &coarseList, ntiles,
         *input ), status );
 
    /* Convert output data structure. */
    bank = (SnglInspiralTable *) LALCalloc(1, sizeof(SnglInspiralTable));
    if (bank == NULL){
      ABORT( status, LALINSPIRALBANKH_EMEM, LALINSPIRALBANKH_MSGEMEM );
    }
    *first = bank;
    for( cnt = 0; cnt < *ntiles; cnt++ )
    {
      bank = bank->next = (SnglInspiralTable *) LALCalloc( 1, sizeof(
             SnglInspiralTable ) );
      if (bank == NULL)
      {
        ABORT( status, LALINSPIRALBANKH_EMEM, LALINSPIRALBANKH_MSGEMEM );
      }
      bank->mass1 = coarseList[cnt].params.mass1;
      bank->mass2 = coarseList[cnt].params.mass2;
      bank->mchirp = coarseList[cnt].params.chirpMass;
      bank->mtotal = coarseList[cnt].params.totalMass;
      bank->eta = coarseList[cnt].params.eta;
      bank->tau0 = coarseList[cnt].params.t0;
      bank->tau2 = coarseList[cnt].params.t2;
      bank->tau3 = coarseList[cnt].params.t3;
      bank->tau4 = coarseList[cnt].params.t4;
      bank->tau5 = coarseList[cnt].params.t5;
      bank->ttotal = coarseList[cnt].params.tC;
      bank->psi0 = coarseList[cnt].params.psi0;
      bank->psi3 = coarseList[cnt].params.psi3;
      bank->f_final = coarseList[cnt].params.fFinal;
      bank->eta = coarseList[cnt].params.eta;
      bank->beta = coarseList[cnt].params.beta;
      
      /* Copy the 10 metric co-efficients ... */
      memcpy (bank->Gamma, coarseList[cnt].metric.Gamma, 10*sizeof(REAL4));
      
    }
    /* Free first template, which is blank. */
    bank = (*first)->next;
    LALFree( *first );
    *first = bank;
    /* free the coarse list returned by create coarse bank */
    LALFree( coarseList );
    break;
  
  case FindChirpPTF:
    
    for (i=0; i<5; i++)
    {
      if ( i < 3 ) chi[i]     = 0.1 + i * 0.4 ;
      if ( i < 2 ) kappa[i]   = -0.9 + i * 0.4 ;
      if ( i > 2 ) kappa[i-1] = 0.5 + (i-3) * 0.4;
    }
    /* Use LALInspiralCreateCoarseBank(). */
    TRY( LALInspiralCreateCoarseBank( status->statusPtr, &coarseList, ntiles,
         *input ), status );
 
    /* Convert output data structure. */
    bank = (SnglInspiralTable *) LALCalloc(1, sizeof(SnglInspiralTable));
    if (bank == NULL){
      ABORT( status, LALINSPIRALBANKH_EMEM, LALINSPIRALBANKH_MSGEMEM );
    }
    *first = bank;
    for ( chicnt = 0; chicnt < 3; chicnt++ )
    {
      for( kappacnt = 0; kappacnt < 4; kappacnt++ )
      {
        fprintf(stderr,"kappacnt=%d\n",kappacnt);
        for( cnt = 0; cnt < *ntiles; cnt++ )
        {
          bank = bank->next = (SnglInspiralTable *) LALCalloc( 1, sizeof(
                SnglInspiralTable ) );
          if (bank == NULL)
          {
            ABORT( status, LALINSPIRALBANKH_EMEM, LALINSPIRALBANKH_MSGEMEM );
          }
          bank->mass1   = coarseList[cnt].params.mass1;
          bank->mass2   = coarseList[cnt].params.mass2;
          bank->mchirp  = coarseList[cnt].params.chirpMass;
          bank->mtotal  = coarseList[cnt].params.totalMass;
          bank->eta     = coarseList[cnt].params.eta;
          bank->kappa   = kappa[kappacnt];
          fprintf(stderr,"kappa[%d]=%e\n",kappacnt,kappa[kappacnt]);
          bank->chi     = chi[chicnt];
          bank->tau0    = coarseList[cnt].params.t0;
          bank->tau2    = coarseList[cnt].params.t2;
          bank->tau3    = coarseList[cnt].params.t3;
          bank->tau4    = coarseList[cnt].params.t4;
          bank->tau5    = coarseList[cnt].params.t5;
          bank->ttotal  = coarseList[cnt].params.tC;
          bank->psi0    = coarseList[cnt].params.psi0;
          bank->psi3    = coarseList[cnt].params.psi3;
          bank->f_final = coarseList[cnt].params.fFinal;
          bank->eta     = coarseList[cnt].params.eta;
          bank->beta    = coarseList[cnt].params.beta;
          

          /* Copy the 10 metric co-efficients ... */
          memcpy (bank->Gamma, coarseList[cnt].metric.Gamma, 10*sizeof(REAL4));

        }
      }
    }
    /* Free first template, which is blank. */
    bank = (*first)->next;
    LALFree( *first );
    *first = bank;
    /* free the coarse list returned by create coarse bank */
    LALFree( coarseList );
    break;

  case BCVSpin:
    if (input->spinBank==0)
    {
    /* Use LALInspiralSpinBank(); no need to convert output. */
    TRY( LALInspiralSpinBank( status->statusPtr, first, ntiles, input ),
         status );   
    }
    else if (input->spinBank==1)
    {
    /* For extended bank use LALInspiralBCVSpinBank() */
/*
    TRY( LALInspiralBCVSpinBank( status->statusPtr, first, ntiles, input ),
         status );   
*/
    }
    else
    {
      ABORT( status, LALINSPIRALBANKH_ECHOICE, LALINSPIRALBANKH_MSGECHOICE );
    }

    if (*ntiles < 1){       
      ABORT( status, LALINSPIRALBANKH_ENULL, LALINSPIRALBANKH_MSGENULL );
    }
    break;

  default:
    ABORT( status, LALINSPIRALBANKH_ECHOICE, LALINSPIRALBANKH_MSGECHOICE );

  }

  DETATCHSTATUSPTR(status);
  RETURN(status); 
}
