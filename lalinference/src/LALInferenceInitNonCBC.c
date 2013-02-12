/*
 *  LALInferenceCBCInit.c:  Bayesian Followup initialisation routines.
 *
 *  Copyright (C) 2012 Vivien Raymond and John Veitch
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
#include <lal/LALInferenceInit.h>

void LALInferenceInitNonCBCTemplate(LALInferenceRunState *runState)
{
  char help[]="(--template [SinGauss,BestIFO]\tSpecify template (default LAL)\n";
  ProcessParamsTable *ppt=NULL;
  ProcessParamsTable *commandLine=runState->commandLine;
  /* Print command line arguments if help requested */
  //Help is taken care of in LALInferenceInitCBCVariables
  //ppt=LALInferenceGetProcParamVal(commandLine,"--help");
  //if(ppt)
  //{
  //	fprintf(stdout,"%s",help);
  //	return;
  //}
  /* This is the LAL template generator for inspiral signals */
  runState->template=&LALInferenceTemplateXLALSimInspiralChooseWaveform;
  ppt=LALInferenceGetProcParamVal(commandLine,"--template");
  if(ppt) {
    if(!strcmp("SinGauss",ppt->value))
        runState->template=&LALInferenceTemplateSineGaussian;
    else if(!strcmp("BestIFO",ppt->value))
        runState->template=&LALInferenceTemplateBestIFO;
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

void LALInferenceInitBurstVariables(LALInferenceRunState *state)
{

    LALStatus status;
    SimBurst *BinjTable=NULL;
	LALInferenceVariables *priorArgs=state->priorArgs;
	state->currentParams=XLALCalloc(1,sizeof(LALInferenceVariables));
	LALInferenceVariables *currentParams=state->currentParams;
	ProcessParamsTable *commandLine=state->commandLine;
	REAL8 endtime;
	ProcessParamsTable *ppt=NULL;
	LALPNOrder PhaseOrder=LAL_PNORDER_THREE_POINT_FIVE;
	INT4 AmpOrder=0;
	Approximant approx=TaylorF2;
	REAL8 logDmin=log(1.0);
	REAL8 logDmax=log(100.0);
    REAL8 phi_spin1_min=0.0;
	REAL8 phi_spin1_max=2.0*LAL_PI;
	REAL8 theta_spin1_min=0.0;
	REAL8 theta_spin1_max=LAL_PI;
	REAL8 fRef=0.; /* freq. at which precessing "initial" cond. specified */
	REAL8 dt=0.1;            /* Width of time prior */
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
		return;
	}
 
    state->likelihood=&LALInferenceUndecomposedFreqDomainLogLikelihood_Burst;
    state->proposal=&NSWrapMCMCLALProposal;
    BinjTable=XLALSimBurstTableFromLIGOLw(LALInferenceGetProcParamVal(commandLine,"--inj")->value,0,0);
    if(!BinjTable){
        fprintf(stderr,"Unable to open injection file %s\n",ppt->value);
        exit(1);
    }
    
     ppt=LALInferenceGetProcParamVal(commandLine,"--event");
        if(ppt){
          event = atoi(ppt->value);
          while(i<event) {i++; BinjTable = BinjTable->next;}
        }
        endtime=XLALGPSGetREAL8(&(BinjTable->time_geocent_gps));
        fprintf(stderr,"Read trig time %lf from injection XML file\n",endtime);
        state->data->modelDomain=LALINFERENCE_DOMAIN_TIME; // salvo
    
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
       
       tmpVal=70.0;
          if(!LALInferenceCheckVariable(currentParams,"frequency")) LALInferenceAddVariable(currentParams, "frequency",     &tmpVal,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
        tmpMin=20.0; tmpMax=280.0;//salvo
        LALInferenceAddMinMaxPrior(priorArgs, "frequency",     &tmpMin, &tmpMax,   LALINFERENCE_REAL8_t);
       tmpVal=-50.0;
          if(!LALInferenceCheckVariable(currentParams,"loghrss")) LALInferenceAddVariable(currentParams, "loghrss",     &tmpVal,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
        tmpMin=-51.0; tmpMax=-46.0;// this is log base e
        LALInferenceAddMinMaxPrior(priorArgs, "loghrss",     &tmpMin, &tmpMax,   LALINFERENCE_REAL8_t);
        
        tmpVal=10.0;
          if(!LALInferenceCheckVariable(currentParams,"Q")) LALInferenceAddVariable(currentParams, "Q",     &tmpVal,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
        tmpMin=1.0; tmpMax=100.0;//salvo
        LALInferenceAddMinMaxPrior(priorArgs, "Q",     &tmpMin, &tmpMax,   LALINFERENCE_REAL8_t);
        tmpVal=0.0;
             tmpVal=0.5;
          if(!LALInferenceCheckVariable(currentParams,"eccentricity")) LALInferenceAddVariable(currentParams, "eccentricity",     &tmpVal,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
        tmpMin=0.0; tmpMax=1.0;//salvo
        LALInferenceAddMinMaxPrior(priorArgs, "eccentricity",     &tmpMin, &tmpMax,   LALINFERENCE_REAL8_t);
         if(!LALInferenceCheckVariable(currentParams,"polar_angle")) LALInferenceAddVariable(currentParams, "polar_angle",    &tmpVal,     LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_CIRCULAR);
        tmpMin=0.0; tmpMax=LAL_PI;
        LALInferenceAddMinMaxPrior(priorArgs, "polar_angle",     &tmpMin, &tmpMax,   LALINFERENCE_REAL8_t);

}

void LALInferenceInitBestIFOVariables(LALInferenceRunState *state)
{

    LALStatus status;
    SimBurst *BinjTable=NULL;
	LALInferenceVariables *priorArgs=state->priorArgs;
	state->currentParams=XLALCalloc(1,sizeof(LALInferenceVariables));
	LALInferenceVariables *currentParams=state->currentParams;
	ProcessParamsTable *commandLine=state->commandLine;
	REAL8 endtime;
	ProcessParamsTable *ppt=NULL;
	LALPNOrder PhaseOrder=LAL_PNORDER_THREE_POINT_FIVE;
	INT4 AmpOrder=0;
	Approximant approx=TaylorF2;
	REAL8 logDmin=log(1.0);
	REAL8 logDmax=log(100.0);
    REAL8 phi_spin1_min=0.0;
	REAL8 phi_spin1_max=2.0*LAL_PI;
	REAL8 theta_spin1_min=0.0;
	REAL8 theta_spin1_max=LAL_PI;
	REAL8 fRef=0.; /* freq. at which precessing "initial" cond. specified */
	REAL8 dt=0.1;            /* Width of time prior */
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
		return;
	}
 
    state->likelihood=&LALInferenceUndecomposedFreqDomainLogLikelihood_Burst;
    state->proposal=&NSWrapMCMCLALProposal;
    BinjTable=XLALSimBurstTableFromLIGOLw(LALInferenceGetProcParamVal(commandLine,"--inj")->value,0,0);
    if(!BinjTable){
        fprintf(stderr,"Unable to open injection file %s\n",ppt->value);
        exit(1);
    }
    
     ppt=LALInferenceGetProcParamVal(commandLine,"--event");
        if(ppt){
          event = atoi(ppt->value);
          while(i<event) {i++; BinjTable = BinjTable->next;}
        }
        endtime=XLALGPSGetREAL8(&(BinjTable->time_geocent_gps));
        fprintf(stderr,"Read trig time %lf from injection XML file\n",endtime);
        state->data->modelDomain=LALINFERENCE_DOMAIN_TIME; // salvo
    
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
       
   fprintf(stdout,"WARNING: Using bestIFO template, only the extrinsic parameters are enabled!\n");
            state->likelihood=&LALInferenceUndecomposedFreqDomainLogLikelihood_BestIFO;

}

