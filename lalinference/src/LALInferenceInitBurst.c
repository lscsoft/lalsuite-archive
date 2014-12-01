/*
 *  LALInferenceCBCInit.c:  Bayesian Followup initialisation routines.
 *
 *  Copyright (C) 2013 James Clark, John Veitch, Salvatore Vitale
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


#include <stdio.h>
#include <lal/Date.h>
#include <lal/GenerateInspiral.h>
#include <lal/LALInference.h>
#include <lal/FrequencySeries.h>
#include <lal/Units.h>
#include <lal/StringInput.h>
#include <lal/LIGOLwXMLInspiralRead.h>
#include <lal/TimeSeries.h>
#include <lal/LALInferencePrior.h>
#include <lal/LALInferenceTemplate.h>
#include <lal/LALInferenceProposal.h>
#include <lal/LALInferenceLikelihood.h>
#include <lal/LALInferenceReadData.h>
#include <lal/LALInferenceReadBurstData.h>
#include <lal/LALInferenceInit.h>
#include <lal/LIGOLwXMLBurstRead.h>
#include <lal/GenerateBurst.h>
#include <lal/LALSimBurst.h>

LALInferenceTemplateFunction LALInferenceInitBurstTemplate(LALInferenceRunState *runState)
{

  char help[]="(--approx [SineGaussian,SineGaussianF,Gaussian,GaussianF,RingdownF]\tSpecify approximant to use (default SineGaussianF)\n";
  ProcessParamsTable *ppt=NULL;
  ProcessParamsTable *commandLine=runState->commandLine;
  /* Print command line arguments if help requested */
  LALInferenceTemplateFunction templt = &LALInferenceTemplateXLALSimInspiralChooseWaveform;
  
  ppt=LALInferenceGetProcParamVal(commandLine,"--approx");
  if(ppt) {
    if(XLALSimBurstImplementedFDApproximants(XLALGetBurstApproximantFromString(ppt->value))){
        templt=&LALInferenceTemplateXLALSimBurstChooseWaveform;
    }
    else if(XLALSimBurstImplementedTDApproximants(XLALGetBurstApproximantFromString(ppt->value))){
        templt=&LALInferenceTemplateXLALSimBurstChooseWaveform;
    }
    else {
      XLALPrintError("Error: unknown template %s\n",ppt->value);
      XLALPrintError(help);
      //XLAL_ERROR_VOID(XLAL_EINVAL);
    }
  }
  return templt;
}

/* Setup the variables to control Burst template generation */
/* Includes specification of prior ranges */

LALInferenceModel * LALInferenceInitBurstModel(LALInferenceRunState *state)
{
  fprintf(stderr,"Using LALInferenceBurstVariables!\n");

  LALStatus status;
  SimBurst *BinjTable=NULL;
  SimInspiralTable *inj_table=NULL;
	state->currentParams=XLALCalloc(1,sizeof(LALInferenceVariables));
	ProcessParamsTable *commandLine=state->commandLine;
	REAL8 endtime=-1;
	REAL8 endtime_from_inj=-1;
  ProcessParamsTable *ppt=NULL;
	memset(&status,0,sizeof(LALStatus));
	INT4 event=0;	
	INT4 i=0;
  BurstApproximant approx = (BurstApproximant) 0;
  char *pinned_params=NULL;
	char help[]="\
Parameter arguments:\n\
(--inj injections.xml)\tSimInspiral or SimBurst Injection XML file to use\n\
(--dt time)\tWidth of time prior, centred around trigger (0.1s)\n\
(--trigtime time)\tTrigger time to use\n\
(--approx Approximant)\tSet approximant (SineGaussianF,SineGaussian,Gaussian,RingdownF)\n\
(--fref fRef)\tSpecify a reference frequency at which parameters are defined (default 0).\n\
(--pinparams [frequency,q,loghrss, etc])\n\tList of parameters to set to injected values\n";
	/* Print command line arguments if help requested */
	ppt=LALInferenceGetProcParamVal(commandLine,"--help");
	if(ppt)
	{
		fprintf(stdout,"%s",help);
		return 0;
	}
  
  LALInferenceModel *model = XLALMalloc(sizeof(LALInferenceModel));
  model->params = XLALCalloc(1, sizeof(LALInferenceVariables));
  memset(model->params, 0, sizeof(LALInferenceVariables));
  LALInferenceVariables *currentParams=model->params;
  
  //state->proposal=&NSWrapMCMCSineGaussProposal;
  
  /* We may have used a CBC injection... test */
  ppt=LALInferenceGetProcParamVal(commandLine,"--trigtime");
  if (ppt)
      endtime=atof(ppt->value);
  ppt=LALInferenceGetProcParamVal(commandLine,"--inj");
  if (ppt) {
    BinjTable=XLALSimBurstTableFromLIGOLw(LALInferenceGetProcParamVal(commandLine,"--inj")->value,0,0);
    if (BinjTable){
      //burst_inj=1;
      ppt=LALInferenceGetProcParamVal(commandLine,"--event");
      if(ppt){
        event = atoi(ppt->value);
        while(i<event) {i++; BinjTable = BinjTable->next;}
      }
      else
        fprintf(stdout,"WARNING: You did not provide an event number with you --inj. Using default event=0 which may not be what you want!!!!\n");
      endtime_from_inj=XLALGPSGetREAL8(&(BinjTable->time_geocent_gps));
    }
    else{
      SimInspiralTableFromLIGOLw(&inj_table,LALInferenceGetProcParamVal(commandLine,"--inj")->value,0,0);
      if (inj_table){
        ppt=LALInferenceGetProcParamVal(commandLine,"--event");
        if(ppt){
          event= atoi(ppt->value);
          fprintf(stderr,"Reading event %d from file\n",event);
          i =0;
          while(i<event) {i++; inj_table=inj_table->next;} /* select event */
          endtime_from_inj=XLALGPSGetREAL8(&(inj_table->geocent_end_time));
        }
        else
            fprintf(stdout,"WARNING: You did not provide an event number with you --inj. Using default event=0 which may not be what you want!!!!\n");
      }
    }
  }
  if(!(BinjTable || inj_table || endtime>=0.0 )){
    printf("Did not provide --trigtime or an xml file and event... Exiting.\n");
    exit(1);
  }
  if (endtime_from_inj!=endtime){
    if(endtime_from_inj>0 && endtime>0)
      fprintf(stderr,"WARNING!!! You set trigtime %lf with --trigtime but event %i seems to trigger at time %lf\n",endtime,event,endtime_from_inj);
      else if(endtime_from_inj>0)
    endtime=endtime_from_inj;
  }
        
  if((ppt=LALInferenceGetProcParamVal(commandLine,"--pinparams"))){
    pinned_params=ppt->value;
    LALInferenceVariables tempParams;
    memset(&tempParams,0,sizeof(tempParams));
    char **strings=NULL;
    UINT4 N;
    LALInferenceParseCharacterOptionString(pinned_params,&strings,&N);
    LALInferenceBurstInjectionToVariables(BinjTable,&tempParams);

    LALInferenceVariableItem *node=NULL;
    while(N>0){
      N--;
      char *name=strings[N];
      node=LALInferenceGetItem(&tempParams,name);
      if(node) {LALInferenceAddVariable(currentParams,node->name,node->value,node->type,node->vary); printf("pinned %s \n",node->name);}
      else {fprintf(stderr,"Error: Cannot pin parameter %s. No such parameter found in injection!\n",node->name);}
    }
  }
    
  /* Over-ride approximant if user specifies */
  ppt=LALInferenceGetProcParamVal(commandLine,"--approximant");
  if(ppt){
    approx = XLALGetBurstApproximantFromString(ppt->value);
  }
  ppt=LALInferenceGetProcParamVal(commandLine,"--approx");
  if(ppt){
    approx = XLALGetBurstApproximantFromString(ppt->value);
    }
   /* Set the model domain appropriately */
  if (XLALSimBurstImplementedFDApproximants(approx)) {
    model->domain = LAL_SIM_DOMAIN_FREQUENCY;
  } else if (XLALSimBurstImplementedTDApproximants(approx)) {
    model->domain = LAL_SIM_DOMAIN_TIME;
  } else {
    fprintf(stderr,"ERROR. Unknown approximant number %i. Unable to choose time or frequency domain model.",approx);
    exit(1);
  }
  
  LALInferenceAddVariable(currentParams, "LAL_APPROXIMANT", &approx,        LALINFERENCE_UINT4_t, LALINFERENCE_PARAM_FIXED);
 
    REAL8 psiMin=0.0,psiMax=LAL_PI; 
    REAL8 raMin=0.0,raMax=LAL_TWOPI; 
    REAL8 decMin=-LAL_PI/2.0,decMax=LAL_PI/2.0; 
    REAL8 qMin=3., qMax=100.0;
    REAL8 ffMin=40., ffMax=1024.0;
    REAL8 durMin=1.0e-4; // min and max value of duration for gaussian templates 
    REAL8 durMax=.5;
    REAL8 hrssMin=1.e-23, hrssMax=1.0e-15;
    REAL8 loghrssMin=log(hrssMin),loghrssMax=log(hrssMax);
    REAL8 dt=0.1;
    REAL8 timeMin=endtime-0.5*dt; 
    REAL8 timeMax=endtime+0.5*dt;
    REAL8 phiMin=0.0,phiMax=LAL_TWOPI;

    REAL8 zero=0.0;

    ppt=LALInferenceGetProcParamVal(commandLine,"--dt");
    if (ppt) dt=atof(ppt->value);
    
    LALInferenceRegisterUniformVariableREAL8(state, model->params, "time", zero, timeMin, timeMax, LALINFERENCE_PARAM_LINEAR);
    
    /* If we are marginalising over the time, remove that variable from the model (having set the prior above) */
    /* Also set the prior in model->params, since Likelihood can't access the state! (ugly hack) */
    if(LALInferenceGetProcParamVal(commandLine,"--margtime") || LALInferenceGetProcParamVal(commandLine, "--margtimephi")){
        LALInferenceVariableItem *p=LALInferenceGetItem(state->priorArgs,"time_min");
        LALInferenceAddVariable(model->params,"time_min",p->value,p->type,p->vary);
        p=LALInferenceGetItem(state->priorArgs,"time_max");
        LALInferenceAddVariable(model->params,"time_max",p->value,p->type,p->vary);
        LALInferenceRemoveVariable(model->params,"time");
        if (LALInferenceGetProcParamVal(commandLine, "--margtimephi")) {
            UINT4 margphi = 1;
            LALInferenceAddVariable(model->params, "margtimephi", &margphi, LALINFERENCE_UINT4_t,LALINFERENCE_PARAM_FIXED);
        }
    }

    LALInferenceRegisterUniformVariableREAL8(state, model->params, "rightascension", zero, raMin, raMax, LALINFERENCE_PARAM_CIRCULAR);
    LALInferenceRegisterUniformVariableREAL8(state, model->params, "declination", zero, decMin, decMax, LALINFERENCE_PARAM_LINEAR);
    LALInferenceRegisterUniformVariableREAL8(state, model->params, "polarisation", zero, psiMin, psiMax, LALINFERENCE_PARAM_LINEAR);

    ppt=LALInferenceGetProcParamVal(commandLine,"--approx");
    if (!strcmp("SineGaussian",ppt->value) || !strcmp("SineGaussianF",ppt->value)|| !strcmp("DampedSinusoid",ppt->value) || !strcmp("DampedSinusoidF",ppt->value)){
     
      LALInferenceRegisterUniformVariableREAL8(state, model->params, "frequency",  zero, ffMin, ffMax,   LALINFERENCE_PARAM_LINEAR);
      LALInferenceRegisterUniformVariableREAL8(state, model->params, "quality",  zero,qMin, qMax,   LALINFERENCE_PARAM_LINEAR);
      if(!LALInferenceGetProcParamVal(commandLine,"--margphi") && !LALInferenceGetProcParamVal(commandLine, "--margtimephi")){
      LALInferenceRegisterUniformVariableREAL8(state, model->params, "phase", zero, phiMin, phiMax, LALINFERENCE_PARAM_CIRCULAR);
      }
    }
    else if (!strcmp("Gaussian",ppt->value) || !strcmp("GaussianF",ppt->value)){
      LALInferenceRegisterUniformVariableREAL8(state, model->params,"duration", zero, durMin,durMax, LALINFERENCE_PARAM_LINEAR);
      
      LALInferenceRegisterUniformVariableREAL8(state, model->params, "phase", zero, phiMin, phiMax, LALINFERENCE_PARAM_FIXED);
    }
    
    if (LALInferenceGetProcParamVal(commandLine,"--use-hrss")){
      LALInferenceRegisterUniformVariableREAL8(state, model->params, "hrss",  zero,hrssMin, hrssMax,   LALINFERENCE_PARAM_LINEAR);
    }
    else
      LALInferenceRegisterUniformVariableREAL8(state, model->params, "loghrss",  zero,loghrssMin, loghrssMax,   LALINFERENCE_PARAM_LINEAR);
    
    LALInferenceRegisterUniformVariableREAL8(state,model->params, "alpha", zero,0.0,2*LAL_PI ,LALINFERENCE_PARAM_CIRCULAR);
    ppt=LALInferenceGetProcParamVal(commandLine,"--cross_only");
    if (ppt){
      printf("Fixing alpha to Pi/2 in template ---> only cross polarization will be used\n");
      LALInferenceRegisterUniformVariableREAL8(state,model->params, "alpha", LAL_PI/2.0,0.0,2*LAL_PI ,LALINFERENCE_PARAM_FIXED);
    }
    ppt=LALInferenceGetProcParamVal(commandLine,"--plus_only");
    if (ppt){
        printf("Fixing alpha to 0 in template ---> only plus polarization will be used\n");
        LALInferenceRegisterUniformVariableREAL8(state,model->params, "alpha", 0.0,0.0,2*LAL_PI, LALINFERENCE_PARAM_FIXED);
    }
    
    /* Needs two condition: must be a burst template and the burst injection must have been provided to do those checks
    if (BinjTable && burst_inj){
        
        if (log(BinjTable->hrss) > loghrssMax || log(BinjTable->hrss) < loghrssMin)
            fprintf(stdout,"WARNING: The injected value of loghrss (%.4e) lies outside the prior range. That may be ok if your template and injection belong to different WF families\n",log(BinjTable->hrss));
        if (BinjTable->q > qMax || BinjTable->q < qMin )
            fprintf(stdout,"WARNING: The injected value of q (%lf) lies outside the prior range\n",BinjTable->q);
         if (BinjTable->frequency > ffMax || BinjTable->frequency < ffMin )
            fprintf(stdout,"WARNING: The injected value of centre_frequency (%lf) lies outside the prior range\n",BinjTable->frequency);
        
        ppt=LALInferenceGetProcParamVal(commandLine,"--approx");
        if (!strcmp("Gaussian",ppt->value) || !strcmp("GaussianF",ppt->value)){
          if (BinjTable->duration >durMax || BinjTable->duration < durMin )
            fprintf(stdout,"WARNING: The injected value of centre_frequency (%lf) lies outside the prior range\n",BinjTable->frequency);
        }
        // Check the max Nyquist frequency for this parameter range
        if ( (ffmax+ 3.0*ffmax/qmin) > state->data->fHigh){
            fprintf(stderr,"WARNING, some of the template in your parameter space will be generated at a frequency higher than Nyquist (%lf). This is bound to produce unwanted consequences. Consider increasing the sampling rate, or reducing (increasing) the max (min) value of frequency (Q). With current setting, srate must be higher than %lf\n",state->data->fHigh,2*(ffmax+ 3.0*ffmax/qmin));
            //exit(1);
        }
            
    }*/
  /* Set model sampling rates to be consistent with data */
  model->deltaT = state->data->timeData->deltaT;
  model->deltaF = state->data->freqData->deltaF;
  UINT4 nifo=0;
  LALInferenceIFOData *dataPtr = state->data;
  while (dataPtr != NULL)
  {
    dataPtr = dataPtr->next;
    nifo++;
  }
  /* Initialize waveform buffers */
  model->timehPlus  = XLALCreateREAL8TimeSeries("timehPlus",
                                                &(state->data->timeData->epoch),
                                                0.0,
                                                model->deltaT,
                                                &lalDimensionlessUnit,
                                                state->data->timeData->data->length);
  model->timehCross = XLALCreateREAL8TimeSeries("timehCross",
                                                &(state->data->timeData->epoch),
                                                0.0,
                                                model->deltaT,
                                                &lalDimensionlessUnit,
                                                state->data->timeData->data->length);
  model->freqhPlus = XLALCreateCOMPLEX16FrequencySeries("freqhPlus",
                                                &(state->data->freqData->epoch),
                                                0.0,
                                                model->deltaF,
                                                &lalDimensionlessUnit,
                                                state->data->freqData->data->length);
  model->freqhCross = XLALCreateCOMPLEX16FrequencySeries("freqhCross",
                                                &(state->data->freqData->epoch),
                                                0.0,
                                                model->deltaF,
                                                &lalDimensionlessUnit,
                                                state->data->freqData->data->length);

  /* Create arrays for holding single-IFO likelihoods, etc. */
  model->ifo_loglikelihoods = XLALCalloc(nifo, sizeof(REAL8));
  model->ifo_SNRs = XLALCalloc(nifo, sizeof(REAL8));

  /* Choose proper template */
  model->templt = LALInferenceInitBurstTemplate(state);

  /* Use same window and FFT plans on model as data */
  model->window = state->data->window;
  model->padding = state->data->padding;
  model->timeToFreqFFTPlan = state->data->timeToFreqFFTPlan;
  model->freqToTimeFFTPlan = state->data->freqToTimeFFTPlan;
  /* Initialize waveform cache */
  model->burstWaveformCache = XLALCreateSimBurstWaveformCache();

  return model;
}

LALInferenceModel *LALInferenceInitModelReviewBurstEvidence_unimod(LALInferenceRunState *state)
{
  ProcessParamsTable *commandLine=state->commandLine;
  ProcessParamsTable *ppt=NULL;
  char **strings=NULL;
  char *pinned_params=NULL;
  UINT4 N=0,i,j;
  if((ppt=LALInferenceGetProcParamVal(commandLine,"--pinparams"))){
    pinned_params=ppt->value;
    LALInferenceVariables tempParams;
    memset(&tempParams,0,sizeof(tempParams));
    LALInferenceParseCharacterOptionString(pinned_params,&strings,&N);
  }
  LALInferenceModel *model = XLALCalloc(1, sizeof(LALInferenceModel));
  model->params = XLALCalloc(1, sizeof(LALInferenceVariables));
  i=0;
  
  struct varSettings {const char *name; REAL8 val, min, max;};
  
  struct varSettings setup[]=
  {
    {.name="time", .val=0.001, .min=-0.006121, .max=0.008121},
    {.name="frequency", .val=210., .min=205.346948, .max=216.653052},
    {.name="quality", .val=6.03626, .min=5.043829, .max=6.956171},
    {.name="loghrss", .val=-46., .min=-46.985195, .max=-45.014805},
    {.name="phase", .val=1.008, .min=0.718919, .max=1.281081},
    {.name="polarisation", .val=0.73, .min=0.427564, .max=0.972436},
    {.name="rightascension", .val=LAL_PI, .min=2.837864, .max=3.445321},
    {.name="declination", .val=0.04, .min= -0.334492, .max=0.334492},
    {.name="alpha", .val=0.58, .min=0.200742, .max=0.799258},
    {.name="END", .val=0., .min=0., .max=0.}
  };
  
  while(strcmp("END",setup[i].name))
  {
    LALInferenceParamVaryType type=LALINFERENCE_PARAM_CIRCULAR;
    /* Check if it is to be fixed */
    for(j=0;j<N;j++) if(!strcmp(setup[i].name,strings[j])) {type=LALINFERENCE_PARAM_FIXED; printf("Fixing parameter %s\n",setup[i].name); break;}
    LALInferenceRegisterUniformVariableREAL8(state, model->params, setup[i].name, setup[i].val, setup[i].min, setup[i].max, type);
    i++;
  }
  return(model);
}

LALInferenceModel *LALInferenceInitModelReviewBurstEvidence_bimod(LALInferenceRunState *state)
{
  ProcessParamsTable *commandLine=state->commandLine;
  ProcessParamsTable *ppt=NULL;
  char **strings=NULL;
  char *pinned_params=NULL;
  UINT4 N=0,i,j;
  if((ppt=LALInferenceGetProcParamVal(commandLine,"--pinparams"))){
    pinned_params=ppt->value;
    LALInferenceVariables tempParams;
    memset(&tempParams,0,sizeof(tempParams));
    LALInferenceParseCharacterOptionString(pinned_params,&strings,&N);
  }
  LALInferenceModel *model = XLALCalloc(1, sizeof(LALInferenceModel));
  model->params = XLALCalloc(1, sizeof(LALInferenceVariables));
  i=0;
  
  struct varSettings {const char *name; REAL8 val, min, max;};
  
  struct varSettings setup[]=
  {
    {.name="time", .val=0.001, .min=-0.006121, .max=0.019514},
    {.name="frequency", .val=211., .min=205.346948, .max=225.697936},
    {.name="quality", .val=6.0, .min=5.043829, .max=8.486044},
    {.name="loghrss", .val=-46., .min=-46.985195, .max=-43.438492},
    {.name="phase", .val=1.0, .min=0.718919, .max=1.730810},
    {.name="polarisation", .val=0.73, .min=0.427564,.max=1.408335},
    {.name="rightascension", .val=LAL_PI, .min=2.837864, .max=3.931287},
    {.name="declination", .val=0.0, .min=-0.334492, .max=0.869678},
    {.name="alpha", .val=0.5, .min=0.200742, .max=1.278070},
    {.name="END", .val=0., .min=0., .max=0.}
  };
  
  while(strcmp("END",setup[i].name))
  {
    LALInferenceParamVaryType type=LALINFERENCE_PARAM_CIRCULAR;
    /* Check if it is to be fixed */
    for(j=0;j<N;j++) if(!strcmp(setup[i].name,strings[j])) {type=LALINFERENCE_PARAM_FIXED; printf("Fixing parameter %s\n",setup[i].name); break;}
    LALInferenceRegisterUniformVariableREAL8(state, model->params, setup[i].name, setup[i].val, setup[i].min, setup[i].max, type);
    i++;
  }
  return(model);
}
