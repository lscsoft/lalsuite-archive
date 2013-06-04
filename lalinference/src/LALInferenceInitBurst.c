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

void LALInferenceInitBurstTemplate(LALInferenceRunState *runState)
{

  char help[]="(--template [SineGauss,SineGaussF,BestIFO,RingdownF,HMNS]\tSpecify template (default LAL)\n";
  ProcessParamsTable *ppt=NULL;
  ProcessParamsTable *commandLine=runState->commandLine;
  /* Print command line arguments if help requested */
  
  runState->templt=&LALInferenceTemplateXLALSimInspiralChooseWaveform;
  ppt=LALInferenceGetProcParamVal(commandLine,"--template");
  if(ppt) {
    if(!strcmp("SineGaussF",ppt->value))
        runState->templt=&LALInferenceTemplateSineGaussianF;
    else if(!strcmp("SineGauss",ppt->value))
        runState->templt=&LALInferenceTemplateSineGaussian;
    else if(!strcmp("BestIFO",ppt->value))
        runState->templt=&LALInferenceTemplateBestIFO;
    else if(!strcmp("RingdownF",ppt->value)){
        printf("Using LALInferenceTemplateXLALSimRingdown: congratulations!\n");
        runState->templt=&LALInferenceTemplateXLALSimRingdown;}
    else if(!strcmp("HMNS",ppt->value)){
        printf("Using LALInferenceTemplateHMNS\n");
        runState->templt=&LALInferenceTemplateHMNS;}
    else {
      XLALPrintError("Error: unknown template %s\n",ppt->value);
      XLALPrintError(help);
      XLAL_ERROR_VOID(XLAL_EINVAL);
    }
  }
  else if(LALInferenceGetProcParamVal(commandLine,"--LALSimulation")){
    fprintf(stderr,"Warning: --LALSimulation is deprecated, the LALSimulation package is now the default. To use LALInspiral specify:\n\
                    --template LALGenerateInspiral (for time-domain templates)\n\
                    --template LAL (for frequency-domain templates)\n");
  }
  else {
    fprintf(stdout,"Template function called is \"LALInferenceTemplateXLALSimInspiralChooseWaveform\"\n");
  }
  return;
}


/* Setup the variables to control Burst template generation */
/* Includes specification of prior ranges */

LALInferenceVariables * LALInferenceInitBurstVariables(LALInferenceRunState *state)
{   printf("---------------------------------Using LALInferenceBurstVariables!\n");

    LALStatus status;
    SimBurst *BinjTable=NULL;
     SimInspiralTable *inj_table=NULL;
	LALInferenceVariables *priorArgs=state->priorArgs;
	state->currentParams=XLALCalloc(1,sizeof(LALInferenceVariables));
	LALInferenceVariables *currentParams=state->currentParams;
	ProcessParamsTable *commandLine=state->commandLine;
	REAL8 endtime=0;
	REAL8 endtime_from_inj=0;
	ProcessParamsTable *ppt=NULL;
    
    REAL8 tmpMax, tmpVal,tmpMin;
	memset(currentParams,0,sizeof(LALInferenceVariables));
	memset(&status,0,sizeof(LALStatus));
	INT4 event=0;	
	INT4 i=0;
    char *pinned_params=NULL;
	char help[]="\
Parameter arguments:\n\
(--inj injections.xml)\tInjection XML file to use\n\
(--dt time)\tWidth of time prior, centred around trigger (0.1s)\n\
(--trigtime time)\tTrigger time to use\n\
(--approx ApproximantphaseOrderPN)\tSet approximant (PhenSpin implicitly enables spin)\n\
(--fref fRef)\tSpecify a reference frequency at which parameters are defined (default 0).\n\
(--pinparams [mchirp,asym_massratio,etc])\n\tList of parameters to set to injected values\n\
(--burst_inj)\t Assume burst signals are injected (for the moment singaussian only)\n";
	/* Print command line arguments if help requested */
	ppt=LALInferenceGetProcParamVal(commandLine,"--help");
	if(ppt)
	{
		fprintf(stdout,"%s",help);
		return 0;
	}
    
    int burst_inj=0;
    state->proposal=&NSWrapMCMCSinGaussProposal;
    
    /* We may have used a CBC injection... test */
    ppt=LALInferenceGetProcParamVal(commandLine,"--trigtime");
    if (ppt)
        endtime=atof(ppt->value);
    ppt=LALInferenceGetProcParamVal(commandLine,"--burst_inj");
    if (ppt) {
        burst_inj=1;
        BinjTable=XLALSimBurstTableFromLIGOLw(LALInferenceGetProcParamVal(commandLine,"--inj")->value,0,0);
         ppt=LALInferenceGetProcParamVal(commandLine,"--event");
        if(ppt){
          event = atoi(ppt->value);
          while(i<event) {i++; BinjTable = BinjTable->next;}
        }
        else
        fprintf(stdout,"WARNING: You did not provide an event number with you --inj. Using default event=0 which may not be what you want!!!!\n");
        endtime_from_inj=XLALGPSGetREAL8(&(BinjTable->time_geocent_gps));
        state->data->modelDomain=LAL_SIM_DOMAIN_TIME; // salvo
    }
    else{
	ppt=LALInferenceGetProcParamVal(commandLine,"--inj");
	if (ppt){
	    SimInspiralTableFromLIGOLw(&inj_table,LALInferenceGetProcParamVal(commandLine,"--inj")->value,0,0);
	    ppt=LALInferenceGetProcParamVal(commandLine,"--event");
	    if(ppt){
	      event= atoi(ppt->value);
	      fprintf(stderr,"Reading event %d from file\n",event);
	      i=0;
	      while(i<event) {i++; inj_table=inj_table->next;} /* select event */
	      endtime_from_inj=XLALGPSGetREAL8(&(inj_table->geocent_end_time));
	      state->data->modelDomain=LAL_SIM_DOMAIN_TIME;
	    }
	    else
	    fprintf(stdout,"WARNING: You did not provide an event number with you --inj. Using default event=0 which may not be what you want!!!!\n");
	}
    }

    if(!(BinjTable || inj_table || endtime )){
            printf("Did not provide --trigtime or an xml file and event... Exiting.\n");
            exit(1);
    }
    if(endtime_from_inj && endtime && (endtime_from_inj!=endtime))
	fprintf(stderr,"WARNING!!! You set trigtime %lf with --trigtime but event %i seems to trigger at time %lf\n",endtime,event,endtime_from_inj);
	
    fprintf(stdout,"Set trigtime %lf for template\n",endtime);
    
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

    gsl_rng *GSLrandom=state->GSLrandom;
    
    REAL8 psimin=0.0,psimax=LAL_PI; 
    REAL8 ramin=0.0,ramax=LAL_TWOPI; 
    REAL8 decmin=-LAL_PI/2.0,decmax=LAL_PI/2.0; 
    REAL8 ffmin=40., ffmax=1300.0;
    REAL8 qmin=3., qmax=100.0;
    REAL8 hrssmin=1.e-23, hrssmax=1.0e-21;
    REAL8 loghrssmin=log(hrssmin),loghrssmax=log(hrssmax);
    REAL8 dt=0.1;
    
    /* With these settings should be ok to run at 1024Hz of srate */

    ppt=LALInferenceGetProcParamVal(commandLine,"--loghrssmin");
    if (ppt){ loghrssmin=atof(ppt->value); fprintf(stdout,"Setting min prior for loghrss to %f\n",atof(ppt->value));}
    ppt=LALInferenceGetProcParamVal(commandLine,"--loghrssmax");
    if (ppt){ loghrssmax=atof(ppt->value); fprintf(stdout,"Setting max prior for loghrss to %f\n",atof(ppt->value));}
    ppt=LALInferenceGetProcParamVal(commandLine,"--hrssmin");
    if (ppt){ loghrssmin=log(atof(ppt->value)); fprintf(stdout,"Setting min prior for loghrss to %f\n",log(atof(ppt->value)));}
    ppt=LALInferenceGetProcParamVal(commandLine,"--hrssmax");
    if (ppt){ loghrssmax=log(atof(ppt->value)); fprintf(stdout,"Setting max prior for loghrss to %f\n",log(atof(ppt->value)));}    
    ppt=LALInferenceGetProcParamVal(commandLine,"--qmin");
    if (ppt){ qmin=atof(ppt->value); fprintf(stdout,"Setting min prior for Q to %f\n",atof(ppt->value));}
    ppt=LALInferenceGetProcParamVal(commandLine,"--qmax");
    if (ppt){ qmax=atof(ppt->value); fprintf(stdout,"Setting max prior for Q to %f\n",atof(ppt->value));}
    ppt=LALInferenceGetProcParamVal(commandLine,"--fmin");
    if (ppt){ ffmin=atof(ppt->value);     fprintf(stdout,"Setting min prior for centre frequency to %f\n",atof(ppt->value));}
    ppt=LALInferenceGetProcParamVal(commandLine,"--fmax");
    if (ppt){ ffmax=atof(ppt->value); fprintf(stdout,"Setting max prior for centre frequency to %f\n",atof(ppt->value));}
    ppt=LALInferenceGetProcParamVal(commandLine,"--dt");
    if (ppt) dt=atof(ppt->value);
	

    REAL8 startpsi=psimin+gsl_rng_uniform(GSLrandom)*(psimax-psimin);
    REAL8 startra=ramin+gsl_rng_uniform(GSLrandom)*(ramax-ramin);
    REAL8 startdec=decmin+gsl_rng_uniform(GSLrandom)*(decmax-decmin);
    REAL8 startf=ffmin+gsl_rng_uniform(GSLrandom)*(ffmax-ffmin);
    REAL8 startq=qmin+gsl_rng_uniform(GSLrandom)*(qmax-qmin);
    REAL8 startloghrss=loghrssmin+gsl_rng_uniform(GSLrandom)*(loghrssmax-loghrssmin);
    
    
    if(!LALInferenceCheckVariable(currentParams,"time")){
        ppt=LALInferenceGetProcParamVal(commandLine,"--t");
        if (ppt) {
            endtime=atof(ppt->value);
            fprintf(stdout,"Fixing time to %lf\n",endtime);
            LALInferenceAddVariable(currentParams, "time",&endtime   ,           LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
        }
        else{
            LALInferenceAddVariable(currentParams, "time",            &endtime   ,           LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR); 
        }
    }
    tmpMin=endtime-0.5*dt; tmpMax=endtime+0.5*dt;
    LALInferenceAddMinMaxPrior(priorArgs, "time",     &tmpMin, &tmpMax,   LALINFERENCE_REAL8_t);	

    if(!LALInferenceCheckVariable(currentParams,"rightascension")){
        ppt=LALInferenceGetProcParamVal(commandLine,"--ra");
        if (ppt) {
            tmpVal=atof(ppt->value);
            fprintf(stdout,"Fixing ra to %lf\n",tmpVal);
            LALInferenceAddVariable(currentParams, "rightascension",  &tmpVal,      LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
            }
        else LALInferenceAddVariable(currentParams, "rightascension",  &startra,      LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_CIRCULAR);
    }
    LALInferenceAddMinMaxPrior(priorArgs, "rightascension",     &ramin, &ramax,   LALINFERENCE_REAL8_t);

    if(!LALInferenceCheckVariable(currentParams,"declination")){
        ppt=LALInferenceGetProcParamVal(commandLine,"--dec");
        if (ppt) {
            tmpVal=atof(ppt->value);
            fprintf(stdout,"Fixing dec to %lf\n",tmpVal);
            LALInferenceAddVariable(currentParams, "declination",     &tmpVal,     LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
        }
        else 
            LALInferenceAddVariable(currentParams, "declination",     &startdec,     LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
    }
    LALInferenceAddMinMaxPrior(priorArgs, "declination",     &decmin, &decmax,   LALINFERENCE_REAL8_t);
    
    
    if(!LALInferenceCheckVariable(currentParams,"polarisation")){
        ppt=LALInferenceGetProcParamVal(commandLine,"--psi");
        if (ppt) {
            tmpVal=atof(ppt->value);
            fprintf(stdout,"Fixing psi to %lf\n",tmpVal);
            LALInferenceAddVariable(currentParams, "polarisation",    &tmpVal,     LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
        }
        else
            LALInferenceAddVariable(currentParams, "polarisation",    &startpsi,     LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_CIRCULAR);
    }
    LALInferenceAddMinMaxPrior(priorArgs, "polarisation",     &psimin, &psimax,   LALINFERENCE_REAL8_t);
       
    if(!LALInferenceCheckVariable(currentParams,"frequency")){
        ppt=LALInferenceGetProcParamVal(commandLine,"--freq");
        if (ppt) {
            tmpVal=atof(ppt->value);
            fprintf(stdout,"Fixing freq to %lf\n",tmpVal);
            if (tmpVal<ffmin || tmpVal>ffmax){
                fprintf(stderr,"ERROR: the value of --freq is outside the prior range for the centre frenquency (%f,%f)! Consider increasing the prior range using --fmin --fmax. Exiting...\n",ffmin,ffmax);
                exit(1);
            }
            LALInferenceAddVariable(currentParams, "frequency",     &tmpVal,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
        }
        else
            LALInferenceAddVariable(currentParams, "frequency",     &startf,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
    }
    LALInferenceAddMinMaxPrior(priorArgs, "frequency",     &ffmin, &ffmax,   LALINFERENCE_REAL8_t);
    
    if(!LALInferenceCheckVariable(currentParams,"loghrss")){
        if (LALInferenceGetProcParamVal(commandLine,"--hrss") || LALInferenceGetProcParamVal(commandLine,"--loghrss")){
            if (LALInferenceGetProcParamVal(commandLine,"--hrss")) 
                tmpVal=log(atof((LALInferenceGetProcParamVal(commandLine,"--hrss"))->value));
            else if (LALInferenceGetProcParamVal(commandLine,"--loghrss"))
                tmpVal=atof((LALInferenceGetProcParamVal(commandLine,"--loghrss"))->value);
            fprintf(stdout,"Fixing loghrss to %lf\n",tmpVal);
            if (tmpVal<loghrssmin || tmpVal>loghrssmax){
                fprintf(stderr,"ERROR: the value of --(log)hrss is outside the prior range for the loghrss (%f,%f)! Consider increasing the prior range using --(log)hrssmin --(log)hrssmax. Exiting...\n",loghrssmin,loghrssmax);
                exit(1);
            }
            LALInferenceAddVariable(currentParams, "loghrss",     &tmpVal,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);   
        }
        else
            LALInferenceAddVariable(currentParams, "loghrss",     &startloghrss,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
    }
    LALInferenceAddMinMaxPrior(priorArgs, "loghrss",     &loghrssmin, &loghrssmax,   LALINFERENCE_REAL8_t);
        
    if(!LALInferenceCheckVariable(currentParams,"Q")){
        ppt=LALInferenceGetProcParamVal(commandLine,"--q");
        if (ppt) {
            tmpVal=atof(ppt->value);
            fprintf(stdout,"Fixing q to %lf\n",tmpVal);
            if (tmpVal<qmin || tmpVal>qmax){
                fprintf(stderr,"ERROR: the value of --q is outside the prior range for the qaulity (%f,%f)! Consider increasing the prior range using --qmin --qmax. Exiting...\n",qmin,qmax);
                exit(1);
            }
            LALInferenceAddVariable(currentParams, "Q",     &tmpVal,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
        }
        else
            LALInferenceAddVariable(currentParams, "Q",     &startq,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
    }
    LALInferenceAddMinMaxPrior(priorArgs, "Q",     &qmin, &qmax,   LALINFERENCE_REAL8_t);
    
    tmpVal=1.0;
    if(!LALInferenceCheckVariable(currentParams,"eccentricity")){
        ppt=LALInferenceGetProcParamVal(commandLine,"--SGonly");
        if (ppt){printf("Fixing eccentricity to 0 in template \n");
            LALInferenceAddVariable(currentParams, "eccentricity",     &tmpVal,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);}
        else
            LALInferenceAddVariable(currentParams, "eccentricity",     &tmpVal,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
    }
    tmpMin=0.0; tmpMax=1.0;//salvo
    LALInferenceAddMinMaxPrior(priorArgs, "eccentricity",     &tmpMin, &tmpMax,   LALINFERENCE_REAL8_t);
	
    tmpVal=LAL_PI_2;
    if(!LALInferenceCheckVariable(currentParams,"polar_angle")) {
        ppt=LALInferenceGetProcParamVal(commandLine,"--SGonly");
	    if (ppt){printf("Fixing polar angle to Pi/2 in template \n");
            LALInferenceAddVariable(currentParams, "polar_angle",    &tmpVal,     LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);}
	    else
            LALInferenceAddVariable(currentParams, "polar_angle",    &tmpVal,     LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_CIRCULAR);
    }
    tmpMin=0.0; tmpMax=LAL_PI;
    LALInferenceAddMinMaxPrior(priorArgs, "polar_angle",     &tmpMin, &tmpMax,   LALINFERENCE_REAL8_t);

    if (BinjTable && burst_inj){
        
        if (log(BinjTable->hrss) > loghrssmax || log(BinjTable->hrss) < loghrssmin)
            fprintf(stdout,"WARNING: The injected value of hrss lies outside the prior range\n");
        if (BinjTable->q > qmax || BinjTable->q < qmin )
            fprintf(stdout,"WARNING: The injected value of q lies outside the prior range\n");
         if (BinjTable->frequency > ffmax || BinjTable->frequency < ffmin )
            fprintf(stdout,"WARNING: The injected value of centre_frequency lies outside the prior range\n");
        // Check the max Nyquist frequency for this parameter range
        
        if ( (ffmax+ 3.0*ffmax/qmin) > state->data->fHigh){
            fprintf(stderr,"WARNING, some of the template in your parameter space will be generated at a frequency higher than Nyquist (%lf). Consider increasing the sampling rate, or reducing (increasing) the max (min) value of frequency (Q). With current setting, srate must be higher than %lf\n",state->data->fHigh,2*(ffmax+ 3.0*ffmax/qmin));
            //exit(1);
        }
            
    }
    return currentParams;
}

LALInferenceVariables * LALInferenceInitBestIFOVariables(LALInferenceRunState *state)
{

    LALStatus status;
    SimBurst *BInjTable=NULL;
    SimInspiralTable *injTable=NULL;
    int CBCTemplate=1;
	LALInferenceVariables *priorArgs=state->priorArgs;
	state->currentParams=XLALCalloc(1,sizeof(LALInferenceVariables));
	LALInferenceVariables *currentParams=state->currentParams;
	ProcessParamsTable *commandLine=state->commandLine;
	REAL8 endtime;
	ProcessParamsTable *ppt=NULL;
	REAL8 dt=0.1;            /* Width of time prior */
    REAL8 tmpMin,tmpMax,tmpVal;
	memset(currentParams,0,sizeof(LALInferenceVariables));
	memset(&status,0,sizeof(LALStatus));
	INT4 event=0;	
	INT4 i=0;
    char *pinned_params=NULL;
	char help[]="\
Parameter arguments:\n\
(--inj injections.xml)\tInjection XML file to use\n\
(--dt time)\tWidth of time prior, centred around trigger (0.1s)\n\
(--trigtime time)\tTrigger time to use\n\
(--approx ApproximantphaseOrderPN)\tSet approximant (PhenSpin implicitly enables spin)\n\
(--fref fRef)\tSpecify a reference frequency at which parameters are defined (default 0).\n\
(--pinparams [mchirp,asym_massratio,etc])\n\tList of parameters to set to injected values\n\
(--burst_inj)\t Assume burst signals are injected (for the moment singaussian only)\n";
	/* Print command line arguments if help requested */
	ppt=LALInferenceGetProcParamVal(commandLine,"--help");
	if(ppt)
	{
		fprintf(stdout,"%s",help);
		return 0;
	}
 
    state->proposal=&NSWrapMCMCLALProposal;
    if(LALInferenceGetProcParamVal(commandLine,"--burst_inj")){
            CBCTemplate=0;
            printf("Injecting from burst Table!\n");
            ppt=LALInferenceGetProcParamVal(commandLine,"--inj");
            BInjTable=XLALSimBurstTableFromLIGOLw(ppt->value,0,0);       
            } 
        else
        {   CBCTemplate=1;
            printf("Injecting from inspiral Table!\n");
            SimInspiralTableFromLIGOLw(&injTable,LALInferenceGetProcParamVal(commandLine,"--inj")->value,0,0);
    }   
		if(!(injTable || BInjTable) ){
			XLALPrintError("Unable to open injection file(LALInferenceReadData) %s\n",LALInferenceGetProcParamVal(commandLine,"--inj")->value);
			//XLAL_ERROR_NULL(XLAL_EFUNC);
		}
    
    if (!CBCTemplate){
     ppt=LALInferenceGetProcParamVal(commandLine,"--event");
        if(ppt){
          event = atoi(ppt->value);
          while(i<event) {i++; BInjTable = BInjTable->next;}
        }
        endtime=XLALGPSGetREAL8(&(BInjTable->time_geocent_gps));
        fprintf(stderr,"Read trig time %lf from injection XML file\n",endtime);
        state->data->modelDomain=LAL_SIM_DOMAIN_TIME; // salvo
    
    if((ppt=LALInferenceGetProcParamVal(commandLine,"--pinparams"))){
            pinned_params=ppt->value;
            LALInferenceVariables tempParams;
            memset(&tempParams,0,sizeof(tempParams));
            char **strings=NULL;
            UINT4 N;
            LALInferenceParseCharacterOptionString(pinned_params,&strings,&N);
            LALInferenceBurstInjectionToVariables(BInjTable,&tempParams);

            LALInferenceVariableItem *node=NULL;
            while(N>0){
                N--;
                char *name=strings[N];
                node=LALInferenceGetItem(&tempParams,name);
                if(node) {LALInferenceAddVariable(currentParams,node->name,node->value,node->type,node->vary); printf("pinned %s \n",node->name);}
                else {fprintf(stderr,"Error: Cannot pin parameter %s. No such parameter found in injection!\n",node->name);}
            }
        }
    }
    else{
        
        ppt=LALInferenceGetProcParamVal(commandLine,"--event");
        if(ppt){
          event = atoi(ppt->value);
          while(i<event) {i++; injTable = injTable->next;}
        }
        endtime=XLALGPSGetREAL8(&injTable->geocent_end_time);
        fprintf(stderr,"Read trig time %lf from injection XML file\n",endtime);
        state->data->modelDomain=LAL_SIM_DOMAIN_FREQUENCY; // salvo
    
    if((ppt=LALInferenceGetProcParamVal(commandLine,"--pinparams"))){
            pinned_params=ppt->value;
            LALInferenceVariables tempParams;
            memset(&tempParams,0,sizeof(tempParams));
            char **strings=NULL;
            UINT4 N;
            LALInferenceParseCharacterOptionString(pinned_params,&strings,&N);
            LALInferenceInjectionToVariables(injTable,&tempParams);

            LALInferenceVariableItem *node=NULL;
            while(N>0){
                N--;
                char *name=strings[N];
                node=LALInferenceGetItem(&tempParams,name);
                if(node) {LALInferenceAddVariable(currentParams,node->name,node->value,node->type,node->vary); printf("pinned %s \n",node->name);}
                else {fprintf(stderr,"Error: Cannot pin parameter %s. No such parameter found in injection!\n",node->name);}
            }
        }
        
        }
    /* Over-ride end time if specified */
    ppt=LALInferenceGetProcParamVal(commandLine,"--trigtime");
    if(ppt){
        endtime=atof(ppt->value);
    }

    /* Over-ride time prior if specified */
    ppt=LALInferenceGetProcParamVal(commandLine,"--dt");
    if(ppt){
        dt=atof(ppt->value);
    }
    
 
        if(!LALInferenceCheckVariable(currentParams,"time")) LALInferenceAddVariable(currentParams, "time",            &endtime   ,           LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR); 
    tmpMin=endtime-0.5*dt; tmpMax=endtime+0.5*dt;
    LALInferenceAddMinMaxPrior(priorArgs, "time",     &tmpMin, &tmpMax,   LALINFERENCE_REAL8_t);	

    tmpVal=1.0;
    if(!LALInferenceCheckVariable(currentParams,"rightascension")) LALInferenceAddVariable(currentParams, "rightascension",  &tmpVal,      LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_CIRCULAR);
    tmpMin=0.0; tmpMax=LAL_TWOPI;
    LALInferenceAddMinMaxPrior(priorArgs, "rightascension",     &tmpMin, &tmpMax,   LALINFERENCE_REAL8_t);

    if(!LALInferenceCheckVariable(currentParams,"declination")) LALInferenceAddVariable(currentParams, "declination",     &tmpVal,     LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
    tmpMin=-LAL_PI/2.0; tmpMax=LAL_PI/2.0;
    LALInferenceAddMinMaxPrior(priorArgs, "declination",     &tmpMin, &tmpMax,   LALINFERENCE_REAL8_t);
     if(!LALInferenceCheckVariable(currentParams,"polarisation")) LALInferenceAddVariable(currentParams, "polarisation",    &tmpVal,     LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_CIRCULAR);
        tmpMin=0.0; tmpMax=LAL_PI;
        LALInferenceAddMinMaxPrior(priorArgs, "polarisation",     &tmpMin, &tmpMax,   LALINFERENCE_REAL8_t);
       
   // Pin the other parameters 
    if (CBCTemplate){
        printf("Pinning CBC intrinsic params for BestIFO template\n");
         REAL8 logDmin=log(1.0);
  REAL8 logDmax=log(100.0);
  REAL8 mcMin=1.0;
  REAL8 mcMax=15.3;
  //REAL8 mMin=1.0,mMax=30.0;
  REAL8 logmcMin=0.0,logmcMax=0.0;
  //REAL8 MTotMax=35.0;
  //REAL8 MTotMin=2.0;
  REAL8 etaMin=0.0312;
  REAL8 etaMax=0.25;
  /*REAL8 qMin=mMin/mMax;
  REAL8 qMax=1.0;
  REAL8 m1min=mMin,m1max=mMax;
  REAL8 m2min=mMin,m2max=mMax;
  REAL8 iotaMin=0.0,iotaMax=LAL_PI;
  REAL8 psiMin=0.0,psiMax=LAL_PI;*/
    tmpVal=injTable->mchirp;
                LALInferenceAddVariable(currentParams,"logmc",&tmpVal, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
                logmcMin=log(mcMin); logmcMax=log(mcMax);
                LALInferenceAddMinMaxPrior(priorArgs,	"logmc",	&logmcMin,	&logmcMax,		LALINFERENCE_REAL8_t);
                tmpVal=injTable->eta;
                LALInferenceAddVariable(currentParams, "massratio",       &tmpVal,             LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
                LALInferenceAddMinMaxPrior(priorArgs,	"massratio",	&etaMin,	&etaMax,	LALINFERENCE_REAL8_t);
                            tmpVal=injTable->coa_phase;
                if(!LALInferenceCheckVariable(currentParams,"phase")) LALInferenceAddVariable(currentParams, "phase",           &tmpVal,             LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED); //SALVO
            tmpMin=0.0; tmpMax=LAL_TWOPI;
            LALInferenceAddMinMaxPrior(priorArgs, "phase",     &tmpMin, &tmpMax,   LALINFERENCE_REAL8_t);
            
            if(LALInferenceGetProcParamVal(commandLine,"--no-logdistance"))
            {
                REAL8 Dmin=exp(logDmin);
                REAL8 Dmax=exp(logDmax);
                tmpVal=injTable->distance;;
                LALInferenceAddVariable(currentParams,"distance", &tmpVal, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
                LALInferenceAddMinMaxPrior(priorArgs, "distance",     &Dmin, &Dmax,   LALINFERENCE_REAL8_t);		
            }
            else 
            {
                tmpVal=log(injTable->distance);
                if(!LALInferenceCheckVariable(currentParams,"logdistance")) LALInferenceAddVariable(currentParams,"logdistance", &tmpVal, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
                LALInferenceAddMinMaxPrior(priorArgs, "logdistance",     &logDmin, &logDmax,   LALINFERENCE_REAL8_t);
            }
            
            tmpVal=injTable->inclination;
          if(!LALInferenceCheckVariable(currentParams,"inclination")) LALInferenceAddVariable(currentParams, "inclination",     &tmpVal,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
        tmpMin=0.0; tmpMax=LAL_PI;
        LALInferenceAddMinMaxPrior(priorArgs, "inclination",     &tmpMin, &tmpMax,   LALINFERENCE_REAL8_t); 
            
   }
       
       
   fprintf(stdout,"WARNING: Using bestIFO template, only the extrinsic parameters are enabled!\n");
            //state->likelihood=&LALInferenceUndecomposedFreqDomainLogLikelihood_BestIFO;

return currentParams;

}

/* Setup variables for HMNS analysis (~same as LALInferenceInitRDVariables) */
LALInferenceVariables * LALInferenceInitHMNSVariables(LALInferenceRunState *state)
{   printf("Using LALInferenceInitHMNSVariables: good luck!\n");
	LALStatus status;
	LALInferenceVariables *priorArgs=state->priorArgs;
	state->currentParams=XLALCalloc(1,sizeof(LALInferenceVariables));
	LALInferenceVariables *currentParams=state->currentParams;
	ProcessParamsTable *commandLine=state->commandLine;
	ProcessParamsTable *ppt=NULL;

    /* Always use phase-marginalised likelihood */
    printf("Using Marginalised Phase Likelihood\n");
    //state->likelihood=&LALInferenceMarginalisedPhaseLogLikelihood_HMNS;

    /* Use sine-Gaussian proposal (should be good for bursts...) */
    state->proposal=&NSWrapMCMCSinGaussProposal;
    /* Prior Ranges */
	REAL8 starttime;
	REAL8 loghrssmin=log(1e-25);  /* amplitude parameter */
	REAL8 loghrssmax=log(1e-15);
	REAL8 f0min=1600.0;   /* dominant frequency */
	REAL8 f0max=4000.0;
	REAL8 qualitymin=10;  /* decay time */
	REAL8 qualitymax=100;
	REAL8 phi0min=0.0;    /* initial phase */
	REAL8 phi0max=LAL_TWOPI;
	REAL8 dt=0.1;            /* Width of time prior */
	REAL8 tmpMin,tmpMax,tmpVal;

	memset(currentParams,0,sizeof(LALInferenceVariables));
	memset(&status,0,sizeof(LALStatus));


	char help[]="\
Parameter arguments:\n\
(--inj injections.xml)\tInjection XML file to use (unsupported)\n\
(--trigtime GPS)\tstart time for ringdown\n\
(--dt time)\tWidth of time prior, centred around trigger (0.1s)\n\
(--trigtime time)\tTrigger time to use\n\
(--hrssmin root-sum-squared amplitude)\tMinimum rss amplitude (1e-24)\n\
(--hrssmax root-sum-squared amplitude)\tMaximum rss amplitude (1e-20)\n\
(--f0min frequency)\tMinimum frequency\n\
(--f0max frequency)\tMaximum frequency\n\
(--qualitymin quality)\tMinimum quality factor\n\
(--qualitymax quality)\tMaximum quality factor\n\
(--phi0max initPhase [degrees])\tMax initial phase\n\
(--phi0min initPhase [degrees])\tMin initial phase\n\
(--pin-RA [radians])\t pin right-ascension to here\n\
(--pin-dec [radians])\t pin declination to here\n";

	/* Print command line arguments if help requested */
	ppt=LALInferenceGetProcParamVal(commandLine,"--help");
	if(ppt)
	{
		fprintf(stdout,"%s",help);
		return 0;
	}

	/* Over-ride end time if specified */
	ppt=LALInferenceGetProcParamVal(commandLine,"--trigtime");
	if(ppt){
		starttime=atof(ppt->value);
	}
	
	/* Over-ride time prior if specified */
	ppt=LALInferenceGetProcParamVal(commandLine,"--dt");
	if(ppt){
		dt=atof(ppt->value);
	}

	/* Over-ride Amplitude min if specified */
	ppt=LALInferenceGetProcParamVal(commandLine,"--hrssmin");
	if(ppt){
		loghrssmin=log(atof(ppt->value));
	}
	
	/* Over-ride Amplitude max if specified */
	ppt=LALInferenceGetProcParamVal(commandLine,"--hrssmax");
	if(ppt){
		loghrssmax=log(atof(ppt->value));
	}

	/* Over-ride freq min if specified */
	ppt=LALInferenceGetProcParamVal(commandLine,"--f0min");
	if(ppt){
		f0min=atof(ppt->value);
	}
	
	/* Over-ride Amplitude max if specified */
	ppt=LALInferenceGetProcParamVal(commandLine,"--f0max");
	if(ppt){
		f0max=atof(ppt->value);
	}

	/* Over-ride quality min if specified */
	ppt=LALInferenceGetProcParamVal(commandLine,"--qualitymin");
	if(ppt){
		qualitymin=atof(ppt->value);
	}
	
	/* Over-ride quality max if specified */
	ppt=LALInferenceGetProcParamVal(commandLine,"--qualitymax");
	if(ppt){
		qualitymax=atof(ppt->value);
	}

	/* Over-ride phi0 min if specified */
	ppt=LALInferenceGetProcParamVal(commandLine,"--phi0min");
	if(ppt){
		phi0min=LAL_PI_180*atof(ppt->value);
	}
	
	/* Over-ride phi0 max if specified */
	ppt=LALInferenceGetProcParamVal(commandLine,"--phi0max");
	if(ppt){
		phi0max=LAL_PI_180*atof(ppt->value);
	}

	/* Pin phi0 */
    ppt=LALInferenceGetProcParamVal(commandLine,"--pin-phi0");
	if(ppt){

        tmpVal=LAL_PI_180*atof(ppt->value);

		LALInferenceAddVariable(currentParams, "phase",  &tmpVal,
				LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
		tmpMin=0.0; tmpMax=LAL_TWOPI;
		LALInferenceAddMinMaxPrior(priorArgs, "phase", &tmpMin, &tmpMax,
				LALINFERENCE_REAL8_t);
    }


    /* Pin right ascension */
    ppt=LALInferenceGetProcParamVal(commandLine,"--pin-RA");
	if(ppt){

        tmpVal=atof(ppt->value);

		LALInferenceAddVariable(currentParams, "rightascension",  &tmpVal,
				LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
		tmpMin=0.0; tmpMax=LAL_TWOPI;
		LALInferenceAddMinMaxPrior(priorArgs, "rightascension", &tmpMin, &tmpMax,
				LALINFERENCE_REAL8_t);
    }

    /* Pin declination */
    ppt=LALInferenceGetProcParamVal(commandLine,"--pin-dec");
	if(ppt){

        tmpVal=atof(ppt->value);

		LALInferenceAddVariable(currentParams, "declination",  &tmpVal,
				LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
		tmpMin=-LAL_PI/2.0; tmpMax=LAL_PI/2.0;
		LALInferenceAddMinMaxPrior(priorArgs, "declination", &tmpMin, &tmpMax,
				LALINFERENCE_REAL8_t);
    }

	if(!LALInferenceCheckVariable(currentParams,"rightascension")) 
	{
		LALInferenceAddVariable(currentParams, "rightascension",  &tmpVal,
				LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_CIRCULAR);
		tmpMin=0.0; tmpMax=LAL_TWOPI;
		LALInferenceAddMinMaxPrior(priorArgs, "rightascension", &tmpMin, &tmpMax,
				LALINFERENCE_REAL8_t);
	}
	
	if(!LALInferenceCheckVariable(currentParams,"declination")) 
	{
		LALInferenceAddVariable(currentParams, "declination",  &tmpVal,
				LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
		tmpMin=-LAL_PI/2.0; tmpMax=LAL_PI/2.0;
		LALInferenceAddMinMaxPrior(priorArgs, "declination", &tmpMin, &tmpMax,
				LALINFERENCE_REAL8_t);
	}
    
	if(!LALInferenceCheckVariable(currentParams,"polarisation")) 
	{
		LALInferenceAddVariable(currentParams, "polarisation",  &tmpVal,
				LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_CIRCULAR);
		tmpMin=0.0; tmpMax=LAL_PI;
		LALInferenceAddMinMaxPrior(priorArgs, "polarisation", &tmpMin, &tmpMax,
				LALINFERENCE_REAL8_t);
	}

 	if(!LALInferenceCheckVariable(currentParams,"inclination")) 
	{
		LALInferenceAddVariable(currentParams, "inclination",  &tmpVal,
				LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
		tmpMin=0.0; tmpMax=LAL_PI;
		LALInferenceAddMinMaxPrior(priorArgs, "inclination", &tmpMin, &tmpMax,
				LALINFERENCE_REAL8_t);
	}

    /* Now Add intrinsic variables to Prior arguments */
    if(!LALInferenceCheckVariable(currentParams,"time"))
        LALInferenceAddVariable(currentParams, "time", &starttime,
                LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR); 
	tmpMin=starttime-0.5*dt; tmpMax=starttime+0.5*dt;
    LALInferenceAddMinMaxPrior(priorArgs, "time", &tmpMin, &tmpMax,
            LALINFERENCE_REAL8_t);	

    tmpVal=loghrssmin+(loghrssmax-loghrssmin)/2.0;
    if(!LALInferenceCheckVariable(currentParams,"loghrss"))
        LALInferenceAddVariable(currentParams,"loghrss",&tmpVal,
                LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
    LALInferenceAddMinMaxPrior(priorArgs, "loghrss", &loghrssmin, &loghrssmax,
            LALINFERENCE_REAL8_t);

    tmpVal=f0min+(f0max-f0min)/2.0;
    if(!LALInferenceCheckVariable(currentParams,"frequency"))
        LALInferenceAddVariable(currentParams, "frequency", &tmpVal,
                LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR); 
    LALInferenceAddMinMaxPrior(priorArgs, "frequency", &f0min, &f0max,
            LALINFERENCE_REAL8_t);	

    tmpVal=qualitymin+(qualitymax-qualitymin)/2.0;
    if(!LALInferenceCheckVariable(currentParams,"Q"))
        LALInferenceAddVariable(currentParams, "Q", &tmpVal,
                LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR); 
    LALInferenceAddMinMaxPrior(priorArgs, "Q", &qualitymin, &qualitymax,
            LALINFERENCE_REAL8_t);	

	tmpVal=1.0;
    if(!LALInferenceCheckVariable(currentParams,"phase")) 
        LALInferenceAddVariable(currentParams, "phase", &tmpVal,
                LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_CIRCULAR);
    LALInferenceAddMinMaxPrior(priorArgs, "phase", &phi0min, &phi0max,
            LALINFERENCE_REAL8_t);


    return currentParams;

}

/* Setup variables for ringdown template generation */
LALInferenceVariables * LALInferenceInitRDVariables(LALInferenceRunState *state)
{   printf("Using LALInferenceInitRDVariables: congratulations!\n");
	LALStatus status;
	LALInferenceVariables *priorArgs=state->priorArgs;
	state->currentParams=XLALCalloc(1,sizeof(LALInferenceVariables));
	LALInferenceVariables *currentParams=state->currentParams;
	ProcessParamsTable *commandLine=state->commandLine;
	ProcessParamsTable *ppt=NULL;
    /* state->likelihood=&LALInferenceUndecomposedFreqDomainLogLikelihood_RD;
    // If --margphi, marginalise over phase //
	if(LALInferenceGetProcParamVal(commandLine,"--margphi")){
	  printf("Using Marginalise Phase Likelihood\n");
	  state->likelihood=&LALInferenceMarginalisedPhaseLogLikelihood_RD;
	}
    */
    state->proposal=&NSWrapMCMCSinGaussProposal;
    /* Prior Ranges */
	REAL8 starttime;
	REAL8 loghrssmin=log(1e-25);  /* amplitude parameter */
	REAL8 loghrssmax=log(1e-15);
	REAL8 f0min=1600.0;   /* dominant frequency */
	REAL8 f0max=4000.0;
	REAL8 qualitymin=10;  /* decay time */
	REAL8 qualitymax=100;
	REAL8 phi0min=0.0;    /* initial phase */
	REAL8 phi0max=LAL_TWOPI;
	REAL8 dt=0.1;            /* Width of time prior */
	REAL8 tmpMin,tmpMax,tmpVal;

	memset(currentParams,0,sizeof(LALInferenceVariables));
	memset(&status,0,sizeof(LALStatus));


	char help[]="\
Parameter arguments:\n\
(--inj injections.xml)\tInjection XML file to use (unsupported)\n\
(--trigtime GPS)\tstart time for ringdown\n\
(--dt time)\tWidth of time prior, centred around trigger (0.1s)\n\
(--trigtime time)\tTrigger time to use\n\
(--hrssmin root-sum-squared amplitude)\tMinimum rss amplitude (1e-24)\n\
(--hrssmax root-sum-squared amplitude)\tMaximum rss amplitude (1e-20)\n\
(--f0min frequency)\tMinimum frequency\n\
(--f0max frequency)\tMaximum frequency\n\
(--qualitymin quality)\tMinimum quality factor\n\
(--qualitymax quality)\tMaximum quality factor\n\
(--phi0max initPhase [degrees])\tMax initial phase\n\
(--phi0min initPhase [degrees])\tMin initial phase\n\
(--pin-RA [radians])\t pin right-ascension to here\n\
(--pin-dec [radians])\t pin declination to here\n";

	/* Print command line arguments if help requested */
	ppt=LALInferenceGetProcParamVal(commandLine,"--help");
	if(ppt)
	{
		fprintf(stdout,"%s",help);
		return 0;
	}

    /*
     *     char *pinned_params=NULL;
    * SALVO: Disabled until it is clear if a BinjTable is defined for these signals
     * if((ppt=LALInferenceGetProcParamVal(commandLine,"--pinparams"))){
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
    */
	/* Over-ride end time if specified */
	ppt=LALInferenceGetProcParamVal(commandLine,"--trigtime");
	if(ppt){
		starttime=atof(ppt->value);
	}
	
	/* Over-ride time prior if specified */
	ppt=LALInferenceGetProcParamVal(commandLine,"--dt");
	if(ppt){
		dt=atof(ppt->value);
	}

	/* Over-ride Amplitude min if specified */
	ppt=LALInferenceGetProcParamVal(commandLine,"--hrssmin");
	if(ppt){
		loghrssmin=log(atof(ppt->value));
	}
	
	/* Over-ride Amplitude max if specified */
	ppt=LALInferenceGetProcParamVal(commandLine,"--hrssmax");
	if(ppt){
		loghrssmax=log(atof(ppt->value));
	}

	/* Over-ride freq min if specified */
	ppt=LALInferenceGetProcParamVal(commandLine,"--f0min");
	if(ppt){
		f0min=atof(ppt->value);
	}
	
	/* Over-ride Amplitude max if specified */
	ppt=LALInferenceGetProcParamVal(commandLine,"--f0max");
	if(ppt){
		f0max=atof(ppt->value);
	}

	/* Over-ride quality min if specified */
	ppt=LALInferenceGetProcParamVal(commandLine,"--qualitymin");
	if(ppt){
		qualitymin=atof(ppt->value);
	}
	
	/* Over-ride quality max if specified */
	ppt=LALInferenceGetProcParamVal(commandLine,"--qualitymax");
	if(ppt){
		qualitymax=atof(ppt->value);
	}

	/* Over-ride phi0 min if specified */
	ppt=LALInferenceGetProcParamVal(commandLine,"--phi0min");
	if(ppt){
		phi0min=LAL_PI_180*atof(ppt->value);
	}
	
	/* Over-ride phi0 max if specified */
	ppt=LALInferenceGetProcParamVal(commandLine,"--phi0max");
	if(ppt){
		phi0max=LAL_PI_180*atof(ppt->value);
	}

	/* Pin phi0 */
    ppt=LALInferenceGetProcParamVal(commandLine,"--pin-phi0");
	if(ppt){

        tmpVal=LAL_PI_180*atof(ppt->value);

		LALInferenceAddVariable(currentParams, "phase",  &tmpVal,
				LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
		tmpMin=0.0; tmpMax=LAL_TWOPI;
		LALInferenceAddMinMaxPrior(priorArgs, "phase", &tmpMin, &tmpMax,
				LALINFERENCE_REAL8_t);
    }


    /* Pin right ascension */
    ppt=LALInferenceGetProcParamVal(commandLine,"--pin-RA");
	if(ppt){

        tmpVal=atof(ppt->value);

		LALInferenceAddVariable(currentParams, "rightascension",  &tmpVal,
				LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
		tmpMin=0.0; tmpMax=LAL_TWOPI;
		LALInferenceAddMinMaxPrior(priorArgs, "rightascension", &tmpMin, &tmpMax,
				LALINFERENCE_REAL8_t);
    }

    /* Pin declination */
    ppt=LALInferenceGetProcParamVal(commandLine,"--pin-dec");
	if(ppt){

        tmpVal=atof(ppt->value);

		LALInferenceAddVariable(currentParams, "declination",  &tmpVal,
				LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
		tmpMin=-LAL_PI/2.0; tmpMax=LAL_PI/2.0;
		LALInferenceAddMinMaxPrior(priorArgs, "declination", &tmpMin, &tmpMax,
				LALINFERENCE_REAL8_t);
    }

	if(!LALInferenceCheckVariable(currentParams,"rightascension")) 
	{
		LALInferenceAddVariable(currentParams, "rightascension",  &tmpVal,
				LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_CIRCULAR);
		tmpMin=0.0; tmpMax=LAL_TWOPI;
		LALInferenceAddMinMaxPrior(priorArgs, "rightascension", &tmpMin, &tmpMax,
				LALINFERENCE_REAL8_t);
	}
	
	if(!LALInferenceCheckVariable(currentParams,"declination")) 
	{
		LALInferenceAddVariable(currentParams, "declination",  &tmpVal,
				LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
		tmpMin=-LAL_PI/2.0; tmpMax=LAL_PI/2.0;
		LALInferenceAddMinMaxPrior(priorArgs, "declination", &tmpMin, &tmpMax,
				LALINFERENCE_REAL8_t);
	}
    
	if(!LALInferenceCheckVariable(currentParams,"polarisation")) 
	{
		LALInferenceAddVariable(currentParams, "polarisation",  &tmpVal,
				LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_CIRCULAR);
		tmpMin=0.0; tmpMax=LAL_PI;
		LALInferenceAddMinMaxPrior(priorArgs, "polarisation", &tmpMin, &tmpMax,
				LALINFERENCE_REAL8_t);
	}

 	if(!LALInferenceCheckVariable(currentParams,"inclination")) 
	{
		LALInferenceAddVariable(currentParams, "inclination",  &tmpVal,
				LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
		tmpMin=0.0; tmpMax=LAL_PI;
		LALInferenceAddMinMaxPrior(priorArgs, "inclination", &tmpMin, &tmpMax,
				LALINFERENCE_REAL8_t);
	}

    /* Now Add intrinsic variables to Prior arguments */
    if(!LALInferenceCheckVariable(currentParams,"time"))
        LALInferenceAddVariable(currentParams, "time", &starttime,
                LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR); 
	tmpMin=starttime-0.5*dt; tmpMax=starttime+0.5*dt;
    LALInferenceAddMinMaxPrior(priorArgs, "time", &tmpMin, &tmpMax,
            LALINFERENCE_REAL8_t);	

    tmpVal=loghrssmin+(loghrssmax-loghrssmin)/2.0;
    if(!LALInferenceCheckVariable(currentParams,"loghrss"))
        LALInferenceAddVariable(currentParams,"loghrss",&tmpVal,
                LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
    LALInferenceAddMinMaxPrior(priorArgs, "loghrss", &loghrssmin, &loghrssmax,
            LALINFERENCE_REAL8_t);

    tmpVal=f0min+(f0max-f0min)/2.0;
    if(!LALInferenceCheckVariable(currentParams,"frequency"))
        LALInferenceAddVariable(currentParams, "frequency", &tmpVal,
                LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR); 
    LALInferenceAddMinMaxPrior(priorArgs, "frequency", &f0min, &f0max,
            LALINFERENCE_REAL8_t);	

    tmpVal=qualitymin+(qualitymax-qualitymin)/2.0;
    if(!LALInferenceCheckVariable(currentParams,"Q"))
        LALInferenceAddVariable(currentParams, "Q", &tmpVal,
                LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR); 
    LALInferenceAddMinMaxPrior(priorArgs, "Q", &qualitymin, &qualitymax,
            LALINFERENCE_REAL8_t);	

	tmpVal=1.0;
    if(!LALInferenceCheckVariable(currentParams,"phase")) 
        LALInferenceAddVariable(currentParams, "phase", &tmpVal,
                LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_CIRCULAR);
    LALInferenceAddMinMaxPrior(priorArgs, "phase", &phi0min, &phi0max,
            LALINFERENCE_REAL8_t);


    return currentParams;

}


LALInferenceVariables * LALInferenceInitPowerBurst(LALInferenceRunState *state)
{
    LALInferenceVariables *currentParams = state->currentParams;
    ProcessParamsTable *ppt=NULL;
    if (LALInferenceGetProcParamVal(state->commandLine,"--powerburst"))
    {
        /* log amplitude at 1 Hz */
        /* N.B. at 1000 Hz 10^-7 smaller */
        REAL8 logampmin=log(1.0e-25),logampmax=log(1.0e-20);
        if(LALInferenceGetProcParamVal(state->commandLine,"--powerburst-amp-min"))
            logampmin=log(atoi(ppt->value));
        if(LALInferenceGetProcParamVal(state->commandLine,"--powerburst-amp-max"))
            logampmax=log(atoi(ppt->value));

        REAL8 logamptmp=0.5*(logampmax+logampmin);
        LALInferenceAddVariable(currentParams,"powerburst_logamp",&logamptmp,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_LINEAR);
        LALInferenceAddMinMaxPrior(state->priorArgs,"powerburst_logamp",&logampmin,&logampmax,LALINFERENCE_REAL8_t);
        REAL8 freqmin=state->data->fLow,freqmax=state->data->fHigh;
        if(LALInferenceGetProcParamVal(state->commandLine,"--powerburst-freq-min"))
            freqmin=atoi(ppt->value);
        if(LALInferenceGetProcParamVal(state->commandLine,"--powerburst-freq-max"))
            freqmax=atoi(ppt->value);

        REAL8 freqtmp=0.5*(freqmin+freqmax);
        LALInferenceAddVariable(currentParams,"powerburst_fstar",&freqtmp,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_LINEAR);
        LALInferenceAddMinMaxPrior(state->priorArgs,"powerburst_fstar",&freqmin,&freqmax,LALINFERENCE_REAL8_t);

       // state->likelihood=&LALInferenceExtraPowerLogLikelihood;
    }

    return currentParams;

}
