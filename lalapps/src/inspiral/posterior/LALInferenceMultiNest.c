/* 
 *  InferenceNest.c:  MultiNest with LALInference
 *
 *  Copyright (C) 2009 Ilya Mandel, Vivien Raymond, Christian Roever, Marc van der Sluys, John Veitch and Farhan Feroz
 *
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
#include <lalapps.h>
#include <lal/LALInferencePrior.h>
#include <lal/LALInferenceReadData.h>
#include <lal/LALInferenceLikelihood.h>
#include <lal/LALInferenceTemplate.h>
#include <lal/LALInferenceProposal.h>

LALInferenceRunState *initialize(ProcessParamsTable *commandLine);
void initializeMN(LALInferenceRunState *runState);
void initStudentt(LALInferenceRunState *state);
void initVariables(LALInferenceRunState *state);
void initializeTemplate(LALInferenceRunState *runState);
static void mc2masses(double mc, double eta, double *m1, double *m2);
void MultiNestRun(int mmodal, int ceff, int nlive, double tol, double efr, int ndims, int nPar, int nClsPar,  int maxModes,
	int updInt, double Ztol, char root[], int seed, int *pWrap, int fb, int resume, int outfile, int initMPI, double logZero, 
	void (*LogLike)(double *, int *, int *, double *), void (*dumper)(int *, int *, int *, double **, double **, double *, 
	double *, double *, double *), int context);
void getLogLike(double *Cube, int *ndim, int *npars, double *lnew);
void LALInferenceMultiNestAlgorithm(LALInferenceRunState *runState);
void dumper(int *nSamples, int *nlive, int *nPar, double **physLive, double **posterior, double *paramConstr,
	double *maxLogLike, double *logZ, double *logZerr);

LALInferenceRunState *runStateGlobal;


void MultiNestRun(int mmodal, int ceff, int nlive, double tol, double efr, int ndims, int nPar, int nClsPar,  int maxModes,
int updInt, double Ztol, char root[], int seed, int *pWrap, int fb, int resume, int outfile, int initMPI, double logZero, 
void (*getLogLike)(double *, int *, int *, double *), void (*dumper)(int *, int *, int *, double **, double **, double *, 
double *, double *, double *), int context)
{
	int i;
	for (i = strlen(root); i < 100; i++) root[i] = ' ';

        __nested_MOD_nestrun(&mmodal, &ceff, &nlive, &tol, &efr, &ndims, &nPar, &nClsPar, &maxModes, &updInt, &Ztol,
        root, &seed, pWrap, &fb, &resume, &outfile, &initMPI, &logZero, getLogLike, dumper, &context);
}

void getLogLike(double *Cube, int *ndim, int *npars, double *lnew)
{
	// transform the parameter in the unit hypercube to their physical counterparts according to the prior
	LALInferenceVariables *newParams=NULL;
	newParams=calloc(1,sizeof(LALInferenceVariables));
	/* Make a copy of the parameters passed through currentParams */
	LALInferenceCopyVariables(runStateGlobal->currentParams,newParams);
	int i = runStateGlobal->CubeToPrior(runStateGlobal, newParams, Cube);
	//printf("i = %d\n",i);

	// if the parameters violate the prior then set likelihood to log(0);
	if( i == 0 )
	{
		*lnew = -DBL_MAX;
		LALInferenceDestroyVariables(newParams);
		free(newParams);
		return;
	}
	
	// calculate the loglike
	*lnew=runStateGlobal->likelihood(newParams, runStateGlobal->data, runStateGlobal->template);
	*lnew -= (*(REAL8 *)LALInferenceGetVariable(runStateGlobal->algorithmParams, "logZnoise"));
	LALInferenceDestroyVariables(newParams);
	free(newParams);
}

void dumper(int *nSamples, int *nlive, int *nPar, double **physLive, double **posterior, double *paramConstr, double *maxLogLike, double *logZ, double *logZerr)
{
}


/* MultiNestAlgorithm implements the MultiNest algorithm*/
void LALInferenceMultiNestAlgorithm(LALInferenceRunState *runState)
{
	UINT4 Nlive=*(UINT4 *)LALInferenceGetVariable(runState->algorithmParams,"Nlive");
	REAL8 eff=*(REAL8 *)LALInferenceGetVariable(runState->algorithmParams,"eff");
	REAL8 logZnoise;
	UINT4 verbose=0;
	
	if (LALInferenceGetProcParamVal(runState->commandLine, "--correlatedGaussianLikelihood") 
         || LALInferenceGetProcParamVal(runState->commandLine, "--bimodalGaussianLikelihood")
         || LALInferenceGetProcParamVal(runState->commandLine, "--rosenbrockLikelihood")) {
		logZnoise=0.0;
	} else {
		logZnoise=LALInferenceNullLogLikelihood(runState->data);
	}
	LALInferenceAddVariable(runState->algorithmParams,"logZnoise",&logZnoise,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_FIXED);
	//logLikelihoods=(REAL8 *)(*(REAL8Vector **)LALInferenceGetVariable(runState->algorithmParams,"logLikelihoods"))->data;

	verbose=LALInferenceCheckVariable(runState->algorithmParams,"verbose");
	
	/* output file root */
	ProcessParamsTable *ppt=LALInferenceGetProcParamVal(runState->commandLine,"--outfile");
	if(!ppt){
		fprintf(stderr,"Must specify --outfile <filename.dat>\n");
		exit(1);
	}
	char *outfilestr=ppt->value;
	
	runStateGlobal = runState;
	
	// find out the dimensionality of the problem
	int ND = 0;
	LALInferenceVariableItem *item=runState->currentParams->head;
	for(;item;item=item->next)
	{
		if(item->vary==LALINFERENCE_PARAM_LINEAR || item->vary==LALINFERENCE_PARAM_CIRCULAR) ND++;
	}
	
	if( ND==0 )
	{
		double like = runState->likelihood(runState->currentParams,runState->data,runState->template);
		like -= (*(REAL8 *)LALInferenceGetVariable(runState->algorithmParams, "logZnoise"));
		fprintf(stdout,"LOG-LIKELIHOOD VALUE RETURNED = %g\n",like);
		double prior = LALInferenceInspiralSkyLocPrior(runState,runState->currentParams);
		fprintf(stdout,"LOG-PRIOR VALUE RETURNED = %g\n",prior);
		fprintf(stdout,"LOG-POSTERIOR VALUE RETURNED = %g\n",like+prior);
		exit(0);
	}
	
	int mmodal = 0;
	int ceff = 0;
	int nlive = Nlive;
	double efr = eff;
	double tol = 0.5;
	int ndims = ND;
	int nPar = ndims + 3;
	int nClsPar = 2;
	int updInt = 50;
	double Ztol = -1.e90;
	int maxModes = 1;
	int pWrap[ndims];
	item=runState->currentParams->head;
	int k = -1;
	for(;item;item=item->next)
	{
		if(item->vary==LALINFERENCE_PARAM_LINEAR || item->vary==LALINFERENCE_PARAM_CIRCULAR)
		{
			k++;
			if(item->vary==LALINFERENCE_PARAM_CIRCULAR)
				pWrap[k] = 1;
			else
				pWrap[k] = 0;
		}
	}
	char root[100];
	for( int j = 0; j < 100; j++ ) root[j] = outfilestr[j];
	int rseed = -1;
	int fb = verbose;
	int resume = 1;
	int outfile = 1;
	int initMPI = 0;
	double logZero = -1E90;
	int context = 0;


	MultiNestRun(mmodal, ceff, nlive, tol, efr, ndims, nPar, nClsPar, maxModes, updInt, Ztol, root, rseed, pWrap, fb, 
	resume, outfile, initMPI, logZero, getLogLike, dumper, context);
	

	/* Write out the evidence */
	/*fclose(fpout);
	char bayesfile[100];
	sprintf(bayesfile,"%s_B.txt",outfile);
	fpout=fopen(bayesfile,"w");
	fprintf(fpout,"%lf %lf %lf %lf\n",logZ-logZnoise,logZ,logZnoise,logLmax);
	fclose(fpout);*/
}

static void mc2masses(double mc, double eta, double *m1, double *m2)
/*  Compute individual companion masses (m1, m2)   */
/*  for given chirp mass (m_c) & mass ratio (eta)  */
/*  (note: m1 >= m2).                              */
{
  double root = sqrt(0.25-eta);
  double fraction = (0.5+root) / (0.5-root);
  *m2 = mc * (pow(1+fraction,0.2) / pow(fraction,0.6));
  *m1 = mc * (pow(1+1.0/fraction,0.2) / pow(1.0/fraction,0.6));
  return;
}


LALInferenceRunState *initialize(ProcessParamsTable *commandLine)
/* calls the "ReadData()" function to gather data & PSD from files, */
/* and initializes other variables accordingly.                     */
{
  LALInferenceRunState *irs=NULL;
  LALInferenceIFOData *ifoPtr, *ifoListStart;

  irs = calloc(1, sizeof(LALInferenceRunState));
  /* read data from files: */
  fprintf(stdout, " ==== LALInferenceReadData(): started. ====\n");
  irs->commandLine=commandLine;
  irs->data = LALInferenceReadData(commandLine);
  /* (this will already initialise each LALInferenceIFOData's following elements:  */
  /*     fLow, fHigh, detector, timeToFreqFFTPlan, freqToTimeFFTPlan,     */
  /*     window, oneSidedNoisePowerSpectrum, timeDate, freqData         ) */
  fprintf(stdout, " ==== LALInferenceReadData(): finished. ====\n");
  if (irs->data != NULL) {
    fprintf(stdout, " ==== initialize(): successfully read data. ====\n");

    fprintf(stdout, " ==== LALInferenceInjectInspiralSignal(): started. ====\n");
    LALInferenceInjectInspiralSignal(irs->data,commandLine);
    fprintf(stdout, " ==== LALInferenceInjectInspiralSignal(): finished. ====\n");

    ifoPtr = irs->data;
    ifoListStart = irs->data;
    while (ifoPtr != NULL) {
      /*If two IFOs have the same sampling rate, they should have the same timeModelh*,
        freqModelh*, and modelParams variables to avoid excess computation
        in model waveform generation in the future*/
      LALInferenceIFOData * ifoPtrCompare=ifoListStart;
      int foundIFOwithSameSampleRate=0;
      while (ifoPtrCompare != NULL && ifoPtrCompare!=ifoPtr) {
        if(ifoPtrCompare->timeData->deltaT == ifoPtr->timeData->deltaT){
          ifoPtr->timeModelhPlus=ifoPtrCompare->timeModelhPlus;
          ifoPtr->freqModelhPlus=ifoPtrCompare->freqModelhPlus;
          ifoPtr->timeModelhCross=ifoPtrCompare->timeModelhCross;
          ifoPtr->freqModelhCross=ifoPtrCompare->freqModelhCross;
          ifoPtr->modelParams=ifoPtrCompare->modelParams;
          foundIFOwithSameSampleRate=1;
          break;
        }
        ifoPtrCompare = ifoPtrCompare->next;
      }
      if(!foundIFOwithSameSampleRate){
        ifoPtr->timeModelhPlus  = XLALCreateREAL8TimeSeries("timeModelhPlus",
                                                            &(ifoPtr->timeData->epoch),
                                                            0.0,
                                                            ifoPtr->timeData->deltaT,
                                                            &lalDimensionlessUnit,
                                                            ifoPtr->timeData->data->length);
        ifoPtr->timeModelhCross = XLALCreateREAL8TimeSeries("timeModelhCross",
                                                            &(ifoPtr->timeData->epoch),
                                                            0.0,
                                                            ifoPtr->timeData->deltaT,
                                                            &lalDimensionlessUnit,
                                                            ifoPtr->timeData->data->length);
        ifoPtr->freqModelhPlus = XLALCreateCOMPLEX16FrequencySeries("freqModelhPlus",
                                                                    &(ifoPtr->freqData->epoch),
                                                                    0.0,
                                                                    ifoPtr->freqData->deltaF,
                                                                    &lalDimensionlessUnit,
                                                                    ifoPtr->freqData->data->length);
        ifoPtr->freqModelhCross = XLALCreateCOMPLEX16FrequencySeries("freqModelhCross",
                                                                     &(ifoPtr->freqData->epoch),
                                                                     0.0,
                                                                     ifoPtr->freqData->deltaF,
                                                                     &lalDimensionlessUnit,
                                                                     ifoPtr->freqData->data->length);
        ifoPtr->modelParams = calloc(1, sizeof(LALInferenceVariables));
      }
      ifoPtr = ifoPtr->next;
    }
    irs->currentLikelihood=LALInferenceNullLogLikelihood(irs->data);
    printf("Injection Null Log Likelihood: %g\n", irs->currentLikelihood);
  }
  else{
    fprintf(stdout, " initialize(): no data read.\n");
    irs = NULL;
    return(irs);
  }

  return(irs);
}

/***** Initialise MultiNest structures *****/
/************************************************/
void initializeMN(LALInferenceRunState *runState)
{
	char help[]="\
MultiNest arguments:\n\
 --Nlive N\tNumber of live points to use\n\
 --eff e\ttarget efficiency\n\
(--verbose)\tProduce progress information\n";

	ProcessParamsTable *ppt=NULL;
	ProcessParamsTable *commandLine=runState->commandLine;
	
	
	/* Print command line arguments if help requested */
	ppt=LALInferenceGetProcParamVal(commandLine,"--help");
	if(ppt)
	{
		fprintf(stdout,"%s",help);
		return;
	}

	INT4 verbose=0,tmpi=0;
	REAL8 tmpd=0;
	
	
	/* Initialise parameters structure */
	runState->algorithmParams=XLALCalloc(1,sizeof(LALInferenceVariables));
	runState->priorArgs=XLALCalloc(1,sizeof(LALInferenceVariables));
	
	
	/* Set up the appropriate functions for MultiNest */
	runState->algorithm=&LALInferenceMultiNestAlgorithm;
	
	
	/* Choose the template generator for inspiral signals */
	runState->template=&LALInferenceTemplateLAL;
	if(LALInferenceGetProcParamVal(commandLine,"--LALSimulation")){
		runState->template=&LALInferenceTemplateXLALSimInspiralChooseWaveform;
		fprintf(stdout,"Template function called is \"LALInferenceTemplateXLALSimInspiralChooseWaveform\"\n");
	}else{
		ppt=LALInferenceGetProcParamVal(commandLine,"--approximant");
		if(ppt){
			if(strstr(ppt->value,"TaylorF2")) {
				runState->template=&LALInferenceTemplateLAL;
				fprintf(stdout,"Template function called is \"LALInferenceTemplateLAL\"\n");
			}else if(strstr(ppt->value,"35phase_25amp")) {
				runState->template=&LALInferenceTemplate3525TD;
				fprintf(stdout,"Template function called is \"LALInferenceTemplate3525TD\"\n");
			}else{
				runState->template=&LALInferenceTemplateLALGenerateInspiral;
				fprintf(stdout,"Template function called is \"LALInferenceTemplateLALGenerateInspiral\"\n");
			}
		}
	}
	
	
	/* Set up the loglike function */
	/*if (LALInferenceGetProcParamVal(commandLine,"--tdlike")) {
		fprintf(stderr, "Computing likelihood in the time domain.\n");
		runState->likelihood=&LALInferenceTimeDomainLogLikelihood;
	} else */
    if (LALInferenceGetProcParamVal(commandLine, "--zeroLogLike")) {
		/* Use zero log(L) */
		runState->likelihood=&LALInferenceZeroLogLikelihood;
	} else if (LALInferenceGetProcParamVal(commandLine, "--correlatedGaussianLikelihood")) {
		runState->likelihood=&LALInferenceCorrelatedAnalyticLogLikelihood;
	} else if (LALInferenceGetProcParamVal(commandLine, "--bimodalGaussianLikelihood")) {
		runState->likelihood=&LALInferenceBimodalCorrelatedAnalyticLogLikelihood;
	} else if (LALInferenceGetProcParamVal(commandLine, "--rosenbrockLikelihood")) {
		runState->likelihood=&LALInferenceRosenbrockLogLikelihood;
	} else if (LALInferenceGetProcParamVal(commandLine, "--studentTLikelihood")) {
		fprintf(stderr, "Using Student's T Likelihood.\n");
		runState->likelihood=&LALInferenceFreqDomainStudentTLogLikelihood;
	} else {
		runState->likelihood=&LALInferenceUndecomposedFreqDomainLogLikelihood;
	}
	
	
	/* Set up the prior function */
	if(LALInferenceGetProcParamVal(commandLine,"--skyLocPrior")){
		runState->prior=&LALInferenceInspiralSkyLocPrior;
		runState->CubeToPrior = &LALInferenceInspiralSkyLocCubeToPrior;
	} else if (LALInferenceGetProcParamVal(commandLine, "--S6Prior")) {
		runState->prior=&LALInferenceInspiralPriorNormalised;
		runState->CubeToPrior = &LALInferenceInspiralPriorNormalisedCubeToPrior;
	} else if (LALInferenceGetProcParamVal(commandLine, "--AnalyticGaussPrior")) {
		runState->prior = &LALInferenceNullPrior;
		runState->CubeToPrior = &LALInferenceAnalyticGaussianCubeToPrior;
	} else if (LALInferenceGetProcParamVal(commandLine, "--BimodalGaussPrior")) {
		runState->prior = &LALInferenceNullPrior;
		runState->CubeToPrior = &LALInferenceBimodalGaussianCubeToPrior;
	} else if (LALInferenceGetProcParamVal(commandLine, "--RosenbrockPrior")) {
		runState->prior = &LALInferenceNullPrior;
		runState->CubeToPrior = &LALInferenceRosenbrockCubeToPrior;
	} else {
		runState->prior = &LALInferenceInspiralPriorNormalised;
		runState->CubeToPrior = &LALInferenceInspiralPriorNormalisedCubeToPrior;
	}
	
	
	ppt=LALInferenceGetProcParamVal(commandLine,"--verbose");
	if(ppt) {
		verbose=1;
		LALInferenceAddVariable(runState->algorithmParams,"verbose", &verbose , LALINFERENCE_INT4_t,
					LALINFERENCE_PARAM_FIXED);		
	}
	if(verbose) set_debug_level("ERROR|INFO");
	else set_debug_level("NDEBUG");
	
	
	/* Number of live points */
	printf("set number of live points.\n");
	ppt=LALInferenceGetProcParamVal(commandLine,"--Nlive");
	if(ppt)
		tmpi=atoi(ppt->value);
	else {
		fprintf(stderr,"Error, must specify number of live points\n");
		exit(1);
	}
	LALInferenceAddVariable(runState->algorithmParams,"Nlive",&tmpi, LALINFERENCE_INT4_t,LALINFERENCE_PARAM_FIXED);
	
	
	/* Target efficiency */
	ppt=LALInferenceGetProcParamVal(commandLine,"--eff");
	if(ppt)
		tmpd=atof(ppt->value);
	else {
		tmpd=0.5;
	}
	LALInferenceAddVariable(runState->algorithmParams,"eff",&tmpd, LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_FIXED);
	
	return;
	
}

/* Setup the variables to control template generation */
/* Includes specification of prior ranges */

void initVariables(LALInferenceRunState *state)
{

  char help[]="\
               (--inj injections.xml)       Injection XML file to use\n\
               (--tempSkip )                   Number of iterations between proposed temperature swaps (100)\n\
               (--symMassRatio)                Run with symmetric mass ratio eta, instead of q=m2/m1\n\
               (--mc-min mchirp)               Minimum chirp mass\n\
               (--mc-max mchirp)               Maximum chirp mass\n\
               (--eta-min etaMin)              Minimum eta\n\
               (--eta-max etaMax)              Maximum eta\n\
               (--q-min qMin)                  Minimum q\n\
               (--q-max qMax)                  Maximum q\n\
               (--dt time)                     Width of time prior, centred around trigger (0.1s)\n\
               (--trigtime time)               Trigger time to use\n\
               (--mc mchirp)                   Trigger chirpmass to use\n\
               (--fixMc)                       Do not allow chirpmass to vary\n\
               (--eta eta)                     Trigger eta to use\n\
               (--q q)                         Trigger q to use\n\
               (--fixEta)                      Do not allow mass ratio to vary\n\
               (--fixQ)                        Do not allow mass ratio to vary\n\
               (--phi phase)                   Trigger phase to use\n\
               (--fixPhi)                      Do not allow phase to vary\n\
               (--iota inclination)            Trigger inclination to use\n\
               (--fixIota)                     Do not allow inclination to vary\n\
               (--dist dist)                   Trigger distance\n\
               (--fixDist)                     Do not allow distance to vary\n\
               (--ra ra)                       Trigger RA\n\
               (--fixRa)                       Do not allow RA to vary\n\
               (--dec dec)                     Trigger declination\n\
               (--fixDec)                      Do not allow declination to vary\n\
               (--psi psi)                     Trigger psi\n\
               (--fixPsi)                      Do not allow polarization to vary\n\
               (--a1 a1)                       Trigger a1\n\
               (--fixA1)                       Do not allow spin to vary\n\
               (--theta1 theta1)               Trigger theta1\n\
               (--fixTheta1)                   Do not allow spin 1 colatitude to vary\n\
               (--phi1 phi1)                   Trigger phi1\n\
               (--fixPhi1)                     Do not allow spin 1 longitude to vary\n\
               (--a2 a2)                       Trigger a2\n\
               (--fixA2)                       Do not allow spin 2 to vary\n\
               (--theta2 theta2)               Trigger theta2\n\
               (--fixTheta2)                   Do not allow spin 2 colatitude to vary\n\
               (--phi2 phi2)                   Trigger phi2\n\
               (--fixPhi2)                     Do not allow spin 2 longitude to vary\n\
               (--time time)                   Waveform time (overrides random about trigtime)\n\
               (--fixTime)                     Do not allow coalescence time to vary\n\
               (--Dmin dist)                   Minimum distance in Mpc (1)\n\
               (--Dmax dist)                   Maximum distance in Mpc (100)\n\
               (--approximant Approximant)     Specify a template approximant to use, (default TaylorF2)\n\
               (--order PNorder)               Specify a PN order in phase to use, (default threePointFivePN)\n\
               (--ampOrder PNorder)            Specify a PN order in amplitude to use, (default newtonian)\n\
               (--comp-min min)                Minimum component mass (1.0)\n\
               (--comp-max max)                Maximum component mass (30.0)\n\
               (--MTotMax max)                 Maximum total mass (35.0)\n\
               (--tidal)                       Enables tidal corrections, only with LALSimulation\n\
               (--lambda1)                     Trigger lambda1\n\
               (--fixLambda1)                  Do not allow lambda1 to vary\n\
               (--lambda1-min)                 Minimum lambda1 (0)\n\
               (--lambda1-max)                 Maximum lambda1 (80)\n\
               (--lambda2)                     Trigger lambda2\n\
               (--fixLambda2)                  Do not allow lambda2 to vary\n\
               (--lambda2-min)                 Minimum lambda2 (0)\n\
               (--lambda2-max)                 Maximum lambda2 (80)\n\
               (--interactionFlags)            intercation flags, only with LALSimuation (LAL_SIM_INSPIRAL_INTERACTION_ALL)\n";

  /* Print command line arguments if state was not allocated */
  if(state==NULL)
    {
      fprintf(stdout,"%s",help);
      return;
    }

  /* Print command line arguments if help requested */
  if(LALInferenceGetProcParamVal(state->commandLine,"--help"))
    {
      fprintf(stdout,"%s",help);
      return;
    }


  LALStatus status;
  memset(&status,0,sizeof(status));
  SimInspiralTable *injTable=NULL;
  LALInferenceVariables *priorArgs=state->priorArgs;
  state->currentParams=XLALCalloc(1,sizeof(LALInferenceVariables));
  LALInferenceVariables *currentParams=state->currentParams;
  ProcessParamsTable *commandLine=state->commandLine;
  ProcessParamsTable *ppt=NULL;
  //INT4 AmpOrder=0;
  LALPNOrder PhaseOrder=LAL_PNORDER_THREE_POINT_FIVE;
  LALPNOrder AmpOrder=LAL_PNORDER_NEWTONIAN;
  Approximant approx=TaylorF2;
  LALInferenceApplyTaper bookends = LALINFERENCE_TAPER_NONE;
  UINT4 event=0;
  UINT4 i=0;
  REAL8 logDmin=log(1.0);
  REAL8 logDmax=log(100.0);
  REAL8 Dmin=1.0;
  REAL8 Dmax=100.0;
  REAL8 mcMin=1.0;
  REAL8 mcMax=15.3;
  REAL8 mMin=1.0,mMax=30.0;
  REAL8 MTotMax=35.0;
  REAL8 etaMin=0.0312;
  REAL8 etaMax=0.25;
  REAL8 qMin=mMin/mMax;
  REAL8 qMax=1.0;
  REAL8 dt=0.1;            /* Width of time prior */
  REAL8 lambda1Min=0.0;
  REAL8 lambda1Max=80.0;
  REAL8 lambda2Min=0.0;
  REAL8 lambda2Max=80.0;  
  REAL8 tmpMin,tmpMax;//,tmpVal;
  REAL8 endtime=0.0, timeParam=0.0;
  REAL8 start_mc			=4.82;
  REAL8 start_eta			=etaMin;
  REAL8 start_q           =qMin;
  REAL8 start_phase		=0.0;
  REAL8 start_dist		=8.07955;
  REAL8 start_ra			=0.0;
  REAL8 start_dec			=-LAL_PI/2.0;
  REAL8 start_psi			=0.0;
  REAL8 start_iota		=0.0;
  REAL8 start_a_spin1		=0.0;
  REAL8 start_theta_spin1 =0.0;
  REAL8 start_phi_spin1	=0.0;
  REAL8 start_a_spin2		=0.0;
  REAL8 start_theta_spin2 =0.0;
  REAL8 start_phi_spin2	=0.0;
  REAL8 start_lambda1 =lambda1Min;
  REAL8 start_lambda2 =lambda2Min;
  
  memset(currentParams,0,sizeof(LALInferenceVariables));

  if(LALInferenceGetProcParamVal(commandLine,"--skyLocPrior")){
    MTotMax=20.0;
    mMin=1.0;
    mMax=15.0;
    qMin=mMin/mMax;
    Dmin=10.0;
    Dmax=40.0;
    REAL8 densityVNR=1000.0;
    LALInferenceAddVariable(state->priorArgs,"densityVNR", &densityVNR , LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
  }

  /* Read injection XML file for parameters if specified */
  ppt=LALInferenceGetProcParamVal(commandLine,"--inj");
  if(ppt){
    SimInspiralTableFromLIGOLw(&injTable,ppt->value,0,0);
    if(!injTable){
      fprintf(stderr,"Unable to open injection file %s\n",ppt->value);
      exit(1);
    }
    ppt=LALInferenceGetProcParamVal(commandLine,"--event");
    if(ppt){
      event= atoi(ppt->value);
      fprintf(stderr,"Reading event %d from file\n",event);
      i=0;
      while(i<event) {i++; injTable=injTable->next;} /* select event */
      endtime=XLALGPSGetREAL8(&(injTable->geocent_end_time));
      AmpOrder=injTable->amp_order;
      PhaseOrder = XLALGetOrderFromString(injTable->waveform);
	  if( (int) PhaseOrder == XLAL_FAILURE)
	    ABORTXLAL(&status);
	  approx = XLALGetApproximantFromString(injTable->waveform);
	  if( (int) approx == XLAL_FAILURE)
	    ABORTXLAL(&status);
    }
  }

  /* Over-ride approximant if user specifies */
  ppt=LALInferenceGetProcParamVal(commandLine,"--approximant");
  if(ppt){
    if ( ! strcmp( "GeneratePPN", ppt->value ) )
      {
        approx = GeneratePPN;
      }
    else if ( ! strcmp( "TaylorT1", ppt->value ) )
      {
        approx = TaylorT1;
      }
    else if ( ! strcmp( "TaylorT2", ppt->value ) )
      {
        approx = TaylorT2;
      }
    else if ( ! strcmp( "TaylorT3", ppt->value ) )
      {
        approx = TaylorT3;
      }
    else if ( ! strcmp( "TaylorT4", ppt->value ) )
      {
        approx = TaylorT4;
      }
    else if ( ! strcmp( "TaylorF1", ppt->value ) )
      {
        approx = TaylorF1;
      }
    else if ( ! strcmp( "TaylorF2", ppt->value ) )
      {
        approx = TaylorF2;
      }
    else if ( ! strcmp( "EOB", ppt->value ) )
      {
        approx = EOB;
      }
    else if ( ! strcmp( "EOBNR", ppt->value ) )
      {
        approx = EOBNR;
      }
    else if ( ! strcmp( "EOBNRv2", ppt->value ) )
      {
        approx = EOBNRv2;
      }
    else if ( ! strcmp( "EOBNRv2HM", ppt->value ) )
      {
        approx = EOBNRv2HM;
      }
    else if ( ! strcmp( "SpinTaylor", ppt->value ) )
      {
        approx = SpinTaylor;
      }
    else if ( ! strcmp( "SpinTaylorT3", ppt->value ) )
      {
        approx = SpinTaylorT3;
      }
    else if ( ! strcmp( "SpinQuadTaylor", ppt->value ) )
      {
        approx = SpinQuadTaylor;
      }
    else if ( ! strcmp( "SpinTaylorFrameless", ppt->value ) )
      {
        approx = SpinTaylorFrameless;
      }
    else if ( ! strcmp( "PhenSpinTaylorRD", ppt->value ) )
      {
        approx = PhenSpinTaylorRD;
      }
    else if ( ! strcmp( "NumRel", ppt->value ) )
      {
        approx = NumRel;
      }
    else if ( ! strcmp( "IMRPhenomA", ppt->value ) )
      {
        approx = IMRPhenomA;
      }
    else if ( ! strcmp( "IMRPhenomB", ppt->value ) )
      {
        approx = IMRPhenomB;
      }
    else
      {
        fprintf( stderr, "invalid argument to --approximant\n"
                 "unknown approximant %s specified: "
                 "Approximant must be one of: GeneratePPN, TaylorT1, TaylorT2,\n"
                 "TaylorT3, TaylorT4, TaylorF1, TaylorF2,  EOB, EOBNR, EOBNRv2, \n"
                 "EOBNRv2HM, SpinTaylor, SpinTaylorT3, SpinQuadTaylor, SpinTaylorFrameless,\n"
                 "PhenSpinTaylorRD, NumRel, IMRPhenomA, IMRPhenomB \n", ppt->value);
        exit( 1 );
      }
    fprintf(stdout,"Templates will run using Approximant %s (%u)\n",ppt->value,approx);
  }

  /* Over-ride PN order if user specifies */
  ppt=LALInferenceGetProcParamVal(commandLine,"--order");
  if(ppt){
    if ( ! strcmp( "newtonian", ppt->value ) )
      {
        PhaseOrder = LAL_PNORDER_NEWTONIAN;
      }
    else if ( ! strcmp( "oneHalfPN", ppt->value ) )
      {
        PhaseOrder = LAL_PNORDER_HALF;
      }
    else if ( ! strcmp( "onePN", ppt->value ) )
      {
        PhaseOrder = LAL_PNORDER_ONE;
      }
    else if ( ! strcmp( "onePointFivePN", ppt->value ) )
      {
        PhaseOrder = LAL_PNORDER_ONE_POINT_FIVE;
      }
    else if ( ! strcmp( "twoPN", ppt->value ) )
      {
        PhaseOrder = LAL_PNORDER_TWO;
      }
    else if ( ! strcmp( "twoPointFivePN", ppt->value ) )
      {
        PhaseOrder = LAL_PNORDER_TWO_POINT_FIVE;
      }
    else if ( ! strcmp( "threePN", ppt->value ) )
      {
        PhaseOrder = LAL_PNORDER_THREE;
      }
    else if ( ! strcmp( "threePointFivePN", ppt->value ) )
      {
        PhaseOrder = LAL_PNORDER_THREE_POINT_FIVE;
      }
    else if ( ! strcmp( "pseudoFourPN", ppt->value ) )
      {
        PhaseOrder = LAL_PNORDER_PSEUDO_FOUR;
      }
    else
      {
        fprintf( stderr, "invalid argument to --order:\n"
                 "unknown order specified: "
                 "PN order must be one of: newtonian, oneHalfPN, onePN,\n"
                 "onePointFivePN, twoPN, twoPointFivePN, threePN or\n"
                 "threePointFivePN\n");
        exit( 1 );
      }
    fprintf(stdout,"Templates will be generated at %.1f PN order in phase\n",((float)(PhaseOrder))/2.0);
  }
  
  /* Over-ride amplitude PN order if user specifies */
  ppt=LALInferenceGetProcParamVal(commandLine,"--ampOrder");
  if(ppt){
    if ( ! strcmp( "newtonian", ppt->value ) )
      {
        AmpOrder = LAL_PNORDER_NEWTONIAN;
      }
    else if ( ! strcmp( "oneHalfPN", ppt->value ) )
      {
        AmpOrder = LAL_PNORDER_HALF;
      }
    else if ( ! strcmp( "onePN", ppt->value ) )
      {
        AmpOrder = LAL_PNORDER_ONE;
      }
    else if ( ! strcmp( "onePointFivePN", ppt->value ) )
      {
        AmpOrder = LAL_PNORDER_ONE_POINT_FIVE;
      }
    else if ( ! strcmp( "twoPN", ppt->value ) )
      {
        AmpOrder = LAL_PNORDER_TWO;
      }
    else if ( ! strcmp( "twoPointFivePN", ppt->value ) )
      {
        AmpOrder = LAL_PNORDER_TWO_POINT_FIVE;
      }
    else if ( ! strcmp( "threePN", ppt->value ) )
      {
        AmpOrder = LAL_PNORDER_THREE;
      }
    else if ( ! strcmp( "threePointFivePN", ppt->value ) )
      {
        AmpOrder = LAL_PNORDER_THREE_POINT_FIVE;
      }
    else if ( ! strcmp( "pseudoFourPN", ppt->value ) )
      {
        AmpOrder = LAL_PNORDER_PSEUDO_FOUR;
      }
    else
      {
        fprintf( stderr, "invalid argument to --ampOrder:\n"
                 "unknown order specified: "
                 "PN order must be one of: newtonian, oneHalfPN, onePN,\n"
                 "onePointFivePN, twoPN, twoPointFivePN, threePN or\n"
                 "threePointFivePN\n");
        exit( 1 );
      }
    fprintf(stdout,"Templates will be generated at %.1f PN order in amplitude\n",((float)(AmpOrder))/2.0);
  }
  
  /* This flag was added to account for the broken Big Dog
     injection, which had the opposite sign in H and L compared
     to Virgo. */
  if (LALInferenceGetProcParamVal(commandLine, "--crazyInjectionHLSign")) {
    INT4 flag = 1;
    LALInferenceAddVariable(currentParams, "crazyInjectionHLSign", &flag, LALINFERENCE_UINT4_t, LALINFERENCE_PARAM_FIXED);
  } else {
    INT4 flag = 0;
    LALInferenceAddVariable(currentParams, "crazyInjectionHLSign", &flag, LALINFERENCE_UINT4_t, LALINFERENCE_PARAM_FIXED);
  }

  /* Over-ride taper if specified */
  ppt=LALInferenceGetProcParamVal(commandLine,"--taper");
  if(ppt){
    if(strstr(ppt->value,"STARTEND")) bookends=LALINFERENCE_TAPER_STARTEND;
    if(strstr(ppt->value,"STARTONLY")) bookends=LALINFERENCE_TAPER_START;
    if(strstr(ppt->value,"ENDONLY")) bookends=LALINFERENCE_TAPER_END;
    if(strstr(ppt->value,"RING")) bookends=LALINFERENCE_RING;
    if(strstr(ppt->value,"SMOOTH")) bookends=LALINFERENCE_SMOOTH;
  }

  /* Over-ride end time if specified */
  ppt=LALInferenceGetProcParamVal(commandLine,"--trigtime");
  if(ppt){
    endtime=atof(ppt->value);
  }

  /* Over-ride chirp mass if specified */
  ppt=LALInferenceGetProcParamVal(commandLine,"--mc");
  if(ppt){
    start_mc=atof(ppt->value);
  }

  /* Over-ride eta if specified */
  ppt=LALInferenceGetProcParamVal(commandLine,"--eta");
  if(ppt){
    start_eta=atof(ppt->value);
  }

  /* Over-ride q if specified */
  ppt=LALInferenceGetProcParamVal(commandLine,"--q");
  if(ppt){
    start_q=atof(ppt->value);
  }

  /* Over-ride phase if specified */
  ppt=LALInferenceGetProcParamVal(commandLine,"--phi");
  if(ppt){
    start_phase=atof(ppt->value);
  }

  /* Over-ride inclination if specified */
  ppt=LALInferenceGetProcParamVal(commandLine,"--iota");
  if(ppt){
    start_iota=atof(ppt->value);
  }

  /* Over-ride distance if specified */
  ppt=LALInferenceGetProcParamVal(commandLine,"--dist");
  if (ppt) {
    start_dist = atof(ppt->value);
  }

  ppt=LALInferenceGetProcParamVal(commandLine,"--ra");
  if (ppt) {
    start_ra = atof(ppt->value);
  }

  ppt=LALInferenceGetProcParamVal(commandLine,"--dec");
  if (ppt) {
    start_dec = atof(ppt->value);
  }

  ppt=LALInferenceGetProcParamVal(commandLine,"--psi");
  if (ppt) {
    start_psi = atof(ppt->value);
  }

  ppt=LALInferenceGetProcParamVal(commandLine,"--a1");
  if (ppt) {
    start_a_spin1 = atof(ppt->value);
  }

  ppt=LALInferenceGetProcParamVal(commandLine,"--theta1");
  if (ppt) {
    start_theta_spin1 = atof(ppt->value);
  }

  ppt=LALInferenceGetProcParamVal(commandLine,"--phi1");
  if (ppt) {
    start_phi_spin1 = atof(ppt->value);
  }

  ppt=LALInferenceGetProcParamVal(commandLine,"--a2");
  if (ppt) {
    start_a_spin2 = atof(ppt->value);
  }

  ppt=LALInferenceGetProcParamVal(commandLine,"--theta2");
  if (ppt) {
    start_theta_spin2 = atof(ppt->value);
  }

  ppt=LALInferenceGetProcParamVal(commandLine,"--phi2");
  if (ppt) {
    start_phi_spin2 = atof(ppt->value);
  }

  ppt=LALInferenceGetProcParamVal(commandLine,"--lambda1");
  if (ppt) {
    start_lambda1 = atof(ppt->value);
  }

  ppt=LALInferenceGetProcParamVal(commandLine,"--lambda2");
  if (ppt) {
    start_lambda2 = atof(ppt->value);
  }
  
  /* Over-ride time prior if specified */
  ppt=LALInferenceGetProcParamVal(commandLine,"--dt");
  if(ppt){
    dt=atof(ppt->value);
  }

  /* Over-ride Distance min if specified */
  ppt=LALInferenceGetProcParamVal(commandLine,"--Dmin");
  if(ppt){
    logDmin=log(atof(ppt->value));
    Dmin=atof(ppt->value);
  }

  /* Over-ride Distance max if specified */
  ppt=LALInferenceGetProcParamVal(commandLine,"--Dmax");
  if(ppt){
    logDmax=log(atof(ppt->value));
    Dmax=atof(ppt->value);
  }

  /* Over-ride lambda1 min if specified */
  ppt=LALInferenceGetProcParamVal(commandLine,"--lambda1-min");
  if(ppt){
    lambda1Min=atof(ppt->value);
  }
  
  /* Over-ride lambda1 max if specified */
  ppt=LALInferenceGetProcParamVal(commandLine,"--lambda1-max");
  if(ppt){
    lambda1Max=atof(ppt->value);
  }
  
  /* Over-ride lambda2 min if specified */
  ppt=LALInferenceGetProcParamVal(commandLine,"--lambda2-min");
  if(ppt){
    lambda2Min=atof(ppt->value);
  }
  
  /* Over-ride lambda2 max if specified */
  ppt=LALInferenceGetProcParamVal(commandLine,"--lambda2-max");
  if(ppt){
    lambda2Max=atof(ppt->value);
  }
  
  
  /* Over-ride Mass priors if specified */
  ppt=LALInferenceGetProcParamVal(commandLine,"--mc-min");
  if(ppt)
    {
      mcMin=atof(ppt->value);
      if (mcMin < 0)
        {
          fprintf(stderr,"ERROR: Minimum value of mchirp must be > 0");
          exit(1);
        }
    }

  ppt=LALInferenceGetProcParamVal(commandLine,"--mc-max");
  if(ppt)
    {
      mcMax=atof(ppt->value);
      if (mcMax <= 0)
        {
          fprintf(stderr,"ERROR: Maximum value of mchirp must be > 0");
          exit(1);
        }
    }

  ppt=LALInferenceGetProcParamVal(commandLine,"--eta-min");
  if(ppt)
    {
      etaMin=atof(ppt->value);
      if (etaMin < 0.0)
        {
          fprintf(stderr,"ERROR: Minimum value of eta must be > 0");
          exit(1);
        }
    }

  ppt=LALInferenceGetProcParamVal(commandLine,"--eta-max");
  if(ppt)
    {
      etaMax=atof(ppt->value);
      if (etaMax > 0.25 || etaMax <= 0.0)
        {
          fprintf(stderr,"ERROR: Maximum value of eta must be between 0 and 0.25\n");
          exit(1);
        }
    }

  /* Over-ride component masses */
  ppt=LALInferenceGetProcParamVal(commandLine,"--comp-min");
  if(ppt)	mMin=atof(ppt->value);
  LALInferenceAddVariable(priorArgs,"component_min",&mMin,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_FIXED);

  ppt=LALInferenceGetProcParamVal(commandLine,"--comp-max");
  if(ppt)	mMax=atof(ppt->value);
  LALInferenceAddVariable(priorArgs,"component_max",&mMax,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_FIXED);

  ppt=LALInferenceGetProcParamVal(commandLine,"--MTotMax");
  if(ppt)	MTotMax=atof(ppt->value);
  LALInferenceAddVariable(priorArgs,"MTotMax",&MTotMax,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_FIXED);

  ppt=LALInferenceGetProcParamVal(commandLine,"--q-min");
  if(ppt)
    {
      qMin=atof(ppt->value);
      if (qMin <= 0.0 || qMin < mMin/mMax || qMin < mMin/(MTotMax-mMin) || qMin > 1.0)
        {
          fprintf(stderr,"ERROR: invalid qMin ( max{0,mMin/mMax,mMin/(MTotMax-mMin) < q < 1.0} )");
          exit(1);
        }
    }

  ppt=LALInferenceGetProcParamVal(commandLine,"--q-max");
  if(ppt)
    {
      qMax=atof(ppt->value);
      if (qMax > 1.0 || qMax <= 0.0 || qMax < mMin/mMax || qMax < mMin/(MTotMax-mMin))
        {
          fprintf(stderr,"ERROR: invalid qMax ( max{0,mMin/mMax,mMin/(MTotMax-mMin) < q < 1.0} )");
          exit(1);
        }
    }

  printf("Read end time %f\n",endtime);

  LALInferenceAddVariable(currentParams, "LAL_APPROXIMANT", &approx,        LALINFERENCE_UINT4_t, LALINFERENCE_PARAM_FIXED);
  LALInferenceAddVariable(currentParams, "LAL_PNORDER",     &PhaseOrder,        LALINFERENCE_UINT4_t, LALINFERENCE_PARAM_FIXED);
  if(LALInferenceGetProcParamVal(commandLine,"--ampOrder")) 
    LALInferenceAddVariable(currentParams, "LAL_AMPORDER",     &AmpOrder,        LALINFERENCE_UINT4_t, LALINFERENCE_PARAM_FIXED);

  ppt=LALInferenceGetProcParamVal(commandLine,"--taper");
  if(ppt){
    LALInferenceAddVariable(currentParams, "LALINFERENCE_TAPER",     &bookends,        LALINFERENCE_UINT4_t, LALINFERENCE_PARAM_FIXED);
  }
  ppt=LALInferenceGetProcParamVal(commandLine,"--newswitch");
  int newswitch=0;
  if(ppt){
    newswitch=1;
    LALInferenceAddVariable(currentParams, "newswitch", &newswitch, LALINFERENCE_UINT4_t, LALINFERENCE_PARAM_FIXED);
  }
  /* Set up the variable parameters */
  ppt=LALInferenceGetProcParamVal(commandLine,"--fixMc");
  if(ppt){
    LALInferenceAddVariable(currentParams, "chirpmass",    &start_mc,    LALINFERENCE_REAL8_t,	LALINFERENCE_PARAM_FIXED);
    fprintf(stdout,"chirpmass fixed and set to %f\n",start_mc);
  }else{
    LALInferenceAddVariable(currentParams, "chirpmass",    &start_mc,    LALINFERENCE_REAL8_t,	LALINFERENCE_PARAM_LINEAR);
  }
  LALInferenceAddMinMaxPrior(priorArgs,	"chirpmass",	&mcMin,	&mcMax,		LALINFERENCE_REAL8_t);

  /* Check if running with symmetric (eta) or asymmetric (q) mass ratio.*/
  ppt=LALInferenceGetProcParamVal(commandLine,"--symMassRatio");
  if(ppt){
    ppt=LALInferenceGetProcParamVal(commandLine,"--fixEta");
    if(ppt){
      LALInferenceAddVariable(currentParams, "massratio",       &start_eta,             LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
      fprintf(stdout,"eta fixed and set to %f\n",start_eta);
    }else{
      LALInferenceAddVariable(currentParams, "massratio",       &start_eta,             LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
    }
    LALInferenceAddMinMaxPrior(priorArgs,	"massratio",	&etaMin,	&etaMax,	LALINFERENCE_REAL8_t);
  }else{
    ppt=LALInferenceGetProcParamVal(commandLine,"--fixQ");
    if(ppt){
      LALInferenceAddVariable(currentParams, "asym_massratio",       &start_q,             LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
      fprintf(stdout,"q fixed and set to %f\n",start_q);
    }else{
      LALInferenceAddVariable(currentParams, "asym_massratio",       &start_q,             LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
    }
    LALInferenceAddMinMaxPrior(priorArgs,	"asym_massratio",	&qMin,	&qMax,	LALINFERENCE_REAL8_t);
  }

  tmpMin=endtime-dt; tmpMax=endtime+dt;

  /* Set up start time. */
  ppt=LALInferenceGetProcParamVal(commandLine, "--time");
  if (ppt) {
    /* User has specified start time. */
    timeParam = atof(ppt->value);
  } else {
    timeParam = endtime;
  }

  ppt=LALInferenceGetProcParamVal(commandLine,"--fixTime");
  if(ppt){
    LALInferenceAddVariable(currentParams, "time",            &timeParam   ,           LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
    fprintf(stdout,"time fixed and set to %f\n",timeParam);
  }else{
    LALInferenceAddVariable(currentParams, "time",            &timeParam   ,           LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
  }
  LALInferenceAddMinMaxPrior(priorArgs, "time",     &tmpMin, &tmpMax,   LALINFERENCE_REAL8_t);

  tmpMin=0.0; tmpMax=LAL_TWOPI;
  ppt=LALInferenceGetProcParamVal(commandLine,"--fixPhi");
  if(ppt){
    LALInferenceAddVariable(currentParams, "phase",           &start_phase,        LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
    fprintf(stdout,"phase fixed and set to %f\n",start_phase);
  }else{
    LALInferenceAddVariable(currentParams, "phase",           &start_phase,        LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_CIRCULAR);
  }
  LALInferenceAddMinMaxPrior(priorArgs, "phase",     &tmpMin, &tmpMax,   LALINFERENCE_REAL8_t);

  ppt=LALInferenceGetProcParamVal(commandLine,"--fixDist");
  if(ppt){
    LALInferenceAddVariable(currentParams,"distance", &start_dist, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
    fprintf(stdout,"distance fixed and set to %f\n",start_dist);
  }else{
    LALInferenceAddVariable(currentParams,"distance", &start_dist, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
  }
  LALInferenceAddMinMaxPrior(priorArgs, "distance",     &Dmin, &Dmax,   LALINFERENCE_REAL8_t);


  tmpMin=0.0; tmpMax=LAL_TWOPI;
  ppt=LALInferenceGetProcParamVal(commandLine,"--fixRa");
  if(ppt){
    LALInferenceAddVariable(currentParams, "rightascension",  &start_ra,      LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
    fprintf(stdout,"R.A. fixed and set to %f\n",start_ra);
  }else{
    LALInferenceAddVariable(currentParams, "rightascension",  &start_ra,      LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_CIRCULAR);
  }
  LALInferenceAddMinMaxPrior(priorArgs, "rightascension",     &tmpMin, &tmpMax,   LALINFERENCE_REAL8_t);

  tmpMin=-LAL_PI/2.0; tmpMax=LAL_PI/2.0;
  ppt=LALInferenceGetProcParamVal(commandLine,"--fixDec");
  if(ppt){
    LALInferenceAddVariable(currentParams, "declination",     &start_dec,     LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
    fprintf(stdout,"declination fixed and set to %f\n",start_dec);
  }else{
    LALInferenceAddVariable(currentParams, "declination",     &start_dec,     LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
  }
  LALInferenceAddMinMaxPrior(priorArgs, "declination",     &tmpMin, &tmpMax,   LALINFERENCE_REAL8_t);

  tmpMin=0.0; tmpMax=LAL_PI;
  ppt=LALInferenceGetProcParamVal(commandLine,"--fixPsi");
  if(ppt){
    LALInferenceAddVariable(currentParams, "polarisation",    &start_psi,     LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
    fprintf(stdout,"polarisation fixed and set to %f\n",start_psi);
  }else{
    LALInferenceAddVariable(currentParams, "polarisation",    &start_psi,     LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_CIRCULAR);
  }
  LALInferenceAddMinMaxPrior(priorArgs, "polarisation",     &tmpMin, &tmpMax,   LALINFERENCE_REAL8_t);

  tmpMin=0.0; tmpMax=LAL_PI;

  ppt=LALInferenceGetProcParamVal(commandLine,"--max-iota");
  if (ppt) {
    tmpMax = atof(ppt->value);
  }

  ppt=LALInferenceGetProcParamVal(commandLine,"--fixIota");
  if(ppt){
    LALInferenceAddVariable(currentParams, "inclination",     &start_iota,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
    fprintf(stdout,"iota fixed and set to %f\n",start_iota);
  }else{
    LALInferenceAddVariable(currentParams, "inclination",     &start_iota,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
  }
  LALInferenceAddMinMaxPrior(priorArgs, "inclination",     &tmpMin, &tmpMax,   LALINFERENCE_REAL8_t);

  ppt=LALInferenceGetProcParamVal(commandLine, "--noSpin");
  if((approx==SpinTaylor || approx==SpinTaylorFrameless || approx==PhenSpinTaylorRD) && !ppt){


    ppt=LALInferenceGetProcParamVal(commandLine, "--spinAligned");
    if(ppt) tmpMin=-1.0;
    else tmpMin=0.0;
    tmpMax=1.0;
    ppt=LALInferenceGetProcParamVal(commandLine,"--fixA1");
    if(ppt){
      LALInferenceAddVariable(currentParams, "a_spin1",     &start_a_spin1,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
      fprintf(stdout,"spin 1 fixed and set to %f\n",start_a_spin1);
    }else{
      LALInferenceAddVariable(currentParams, "a_spin1",     &start_a_spin1,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
    }
    LALInferenceAddMinMaxPrior(priorArgs, "a_spin1",     &tmpMin, &tmpMax,   LALINFERENCE_REAL8_t);

    ppt=LALInferenceGetProcParamVal(commandLine, "--spinAligned");
    if(ppt) fprintf(stdout,"Running with spin1 aligned to the orbital angular momentum.\n");
    else {
      tmpMin=0.0; tmpMax=LAL_PI;
      ppt=LALInferenceGetProcParamVal(commandLine,"--fixTheta1");
      if(ppt){
        LALInferenceAddVariable(currentParams, "theta_spin1",     &start_theta_spin1,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
        fprintf(stdout,"theta 1 fixed and set to %f\n",start_theta_spin1);
      }else{
        LALInferenceAddVariable(currentParams, "theta_spin1",     &start_theta_spin1,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
      }
      LALInferenceAddMinMaxPrior(priorArgs, "theta_spin1",     &tmpMin, &tmpMax,   LALINFERENCE_REAL8_t);

      tmpMin=0.0; tmpMax=LAL_TWOPI;
      ppt=LALInferenceGetProcParamVal(commandLine,"--fixPhi1");
      if(ppt){
        LALInferenceAddVariable(currentParams, "phi_spin1",     &start_phi_spin1,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
        fprintf(stdout,"phi 1 fixed and set to %f\n",start_phi_spin1);
      }else{
        LALInferenceAddVariable(currentParams, "phi_spin1",     &start_phi_spin1,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_CIRCULAR);
      }
      LALInferenceAddMinMaxPrior(priorArgs, "phi_spin1",     &tmpMin, &tmpMax,   LALINFERENCE_REAL8_t);
    }
    ppt=LALInferenceGetProcParamVal(commandLine, "--singleSpin");
    if(ppt) fprintf(stdout,"Running with first spin set to 0\n");
    else {
      ppt=LALInferenceGetProcParamVal(commandLine, "--spinAligned");
      if(ppt) tmpMin=-1.0;
      else tmpMin=0.0;
      tmpMax=1.0;
      ppt=LALInferenceGetProcParamVal(commandLine,"--fixA2");
      if(ppt){
        LALInferenceAddVariable(currentParams, "a_spin2",     &start_a_spin2,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
        fprintf(stdout,"spin 2 fixed and set to %f\n",start_a_spin2);
      }else{
        LALInferenceAddVariable(currentParams, "a_spin2",     &start_a_spin2,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
      }
      LALInferenceAddMinMaxPrior(priorArgs, "a_spin2",     &tmpMin, &tmpMax,   LALINFERENCE_REAL8_t);

      ppt=LALInferenceGetProcParamVal(commandLine, "--spinAligned");
      if(ppt) fprintf(stdout,"Running with spin2 aligned to the orbital angular momentum.\n");
      else {
        tmpMin=0.0; tmpMax=LAL_PI;
        ppt=LALInferenceGetProcParamVal(commandLine,"--fixTheta2");
        if(ppt){
          LALInferenceAddVariable(currentParams, "theta_spin2",     &start_theta_spin2,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
          fprintf(stdout,"theta spin 2 fixed and set to %f\n",start_theta_spin2);
        }else{
          LALInferenceAddVariable(currentParams, "theta_spin2",     &start_theta_spin2,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
        }
        LALInferenceAddMinMaxPrior(priorArgs, "theta_spin2",     &tmpMin, &tmpMax,   LALINFERENCE_REAL8_t);

        tmpMin=0.0; tmpMax=LAL_TWOPI;
        ppt=LALInferenceGetProcParamVal(commandLine,"--fixPhi2");
        if(ppt){
          LALInferenceAddVariable(currentParams, "phi_spin2",     &start_phi_spin2,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
          fprintf(stdout,"phi 2 fixed and set to %f\n",start_phi_spin2);
        }else{
          LALInferenceAddVariable(currentParams, "phi_spin2",     &start_phi_spin2,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_CIRCULAR);
        }
        LALInferenceAddMinMaxPrior(priorArgs, "phi_spin2",     &tmpMin, &tmpMax,   LALINFERENCE_REAL8_t);
      }
    }
  }
  ppt=LALInferenceGetProcParamVal(commandLine, "--spinAligned");
  if(approx==TaylorF2 && ppt){

    tmpMin=-1.0; tmpMax=1.0;
    ppt=LALInferenceGetProcParamVal(commandLine,"--fixA1");
    if(ppt){
      LALInferenceAddVariable(currentParams, "spin1",     &start_a_spin1,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
      fprintf(stdout,"spin 1 fixed and set to %f\n",start_a_spin1);
    }else{
      LALInferenceAddVariable(currentParams, "spin1",     &start_a_spin1,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
    }
    LALInferenceAddMinMaxPrior(priorArgs, "spin1",     &tmpMin, &tmpMax,   LALINFERENCE_REAL8_t);

    tmpMin=-1.0; tmpMax=1.0;
    ppt=LALInferenceGetProcParamVal(commandLine,"--fixA2");
    if(ppt){
      LALInferenceAddVariable(currentParams, "spin2",     &start_a_spin2,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
      fprintf(stdout,"spin 2 fixed and set to %f\n",start_a_spin2);
    }else{
      LALInferenceAddVariable(currentParams, "spin2",     &start_a_spin2,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
    }
    LALInferenceAddMinMaxPrior(priorArgs, "spin2",     &tmpMin, &tmpMax,   LALINFERENCE_REAL8_t);

  }

  if (LALInferenceGetProcParamVal(commandLine, "--studentTLikelihood")) {
    
  }

  ppt=LALInferenceGetProcParamVal(commandLine, "--TaylorF2ppE");
  if(approx==TaylorF2 && ppt){

    REAL8 start_alpha, start_A, start_a, start_beta, start_B, start_b;

    tmpMin = -1000;
    tmpMax = 1000;
    start_alpha = tmpMin;
    ppt=LALInferenceGetProcParamVal(commandLine,"--ppealpha");
    if (ppt) {
      start_alpha = atof(ppt->value);
    }
    ppt=LALInferenceGetProcParamVal(commandLine,"--fixppealpha");
    if(ppt){
      LALInferenceAddVariable(currentParams, "ppealpha",     &start_alpha,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
      fprintf(stdout,"ppE alpha fixed and set to %f\n",start_alpha);
    }else{
      LALInferenceAddVariable(currentParams, "ppealpha",     &start_alpha,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
    }
    LALInferenceAddMinMaxPrior(priorArgs, "ppealpha",     &tmpMin, &tmpMax,   LALINFERENCE_REAL8_t);

    start_beta = tmpMin;
    ppt=LALInferenceGetProcParamVal(commandLine,"--ppebeta");
    if (ppt) {
      start_beta = atof(ppt->value);
    }
    ppt=LALInferenceGetProcParamVal(commandLine,"--fixppebeta");
    if(ppt){
      LALInferenceAddVariable(currentParams, "ppebeta",     &start_beta,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
      fprintf(stdout,"ppE beta fixed and set to %f\n",start_beta);
    }else{
      LALInferenceAddVariable(currentParams, "ppebeta",     &start_beta,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
    }
    LALInferenceAddMinMaxPrior(priorArgs, "ppebeta",     &tmpMin, &tmpMax,   LALINFERENCE_REAL8_t);

    tmpMin = -3;
    tmpMax = 3;
    start_A = tmpMin;
    ppt=LALInferenceGetProcParamVal(commandLine,"--ppeuppera");
    if (ppt) {
      start_A = atof(ppt->value);
    }
    ppt=LALInferenceGetProcParamVal(commandLine,"--fixppeuppera");
    if(ppt){
      LALInferenceAddVariable(currentParams, "ppeuppera",     &start_A,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
      fprintf(stdout,"ppE A fixed and set to %f\n",start_A);
    }else{
      LALInferenceAddVariable(currentParams, "ppeuppera",     &start_A,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
    }
    LALInferenceAddMinMaxPrior(priorArgs, "ppeuppera",     &tmpMin, &tmpMax,   LALINFERENCE_REAL8_t);

    start_B = tmpMin;
    ppt=LALInferenceGetProcParamVal(commandLine,"--ppeupperb");
    if (ppt) {
      start_B = atof(ppt->value);
    }
    ppt=LALInferenceGetProcParamVal(commandLine,"--fixppeupperb");
    if(ppt){
      LALInferenceAddVariable(currentParams, "ppeupperb",     &start_B,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
      fprintf(stdout,"ppE B fixed and set to %f\n",start_B);
    }else{
      LALInferenceAddVariable(currentParams, "ppeupperb",     &start_B,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
    }
    LALInferenceAddMinMaxPrior(priorArgs, "ppeupperb",     &tmpMin, &tmpMax,   LALINFERENCE_REAL8_t);

    tmpMin = -3.0;
    tmpMax = 2.0/3.0;
    start_a = tmpMin;
    ppt=LALInferenceGetProcParamVal(commandLine,"--ppelowera");
    if (ppt) {
      start_a = atof(ppt->value);
    }
    ppt=LALInferenceGetProcParamVal(commandLine,"--fixppelowera");
    if(ppt){
      LALInferenceAddVariable(currentParams, "ppelowera",     &start_a,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
      fprintf(stdout,"ppE a fixed and set to %f\n",start_a);
    }else{
      LALInferenceAddVariable(currentParams, "ppelowera",     &start_a,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
    }
    LALInferenceAddMinMaxPrior(priorArgs, "ppelowera",     &tmpMin, &tmpMax,   LALINFERENCE_REAL8_t);

    tmpMin = -4.5;
    tmpMax = 1.0;
    start_b = tmpMin;
    ppt=LALInferenceGetProcParamVal(commandLine,"--ppelowerb");
    if (ppt) {
      start_b = atof(ppt->value);
    }
    ppt=LALInferenceGetProcParamVal(commandLine,"--fixppelowerb");
    if(ppt){
      LALInferenceAddVariable(currentParams, "ppelowerb",     &start_b,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
      fprintf(stdout,"ppE b fixed and set to %f\n",start_b);
    }else{
      LALInferenceAddVariable(currentParams, "ppelowerb",     &start_b,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
    }
    LALInferenceAddMinMaxPrior(priorArgs, "ppelowerb",     &tmpMin, &tmpMax,   LALINFERENCE_REAL8_t);

  }
  
  ppt=LALInferenceGetProcParamVal(commandLine,"--tidal");
  if(ppt){
    ppt=LALInferenceGetProcParamVal(commandLine,"--fixLambda1");
    if(ppt){
      LALInferenceAddVariable(currentParams, "lambda1",           &start_lambda1,        LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
      fprintf(stdout,"phase fixed and set to %f\n",start_lambda1);
    }else{
      LALInferenceAddVariable(currentParams, "lambda1",           &start_lambda1,        LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
    }
    LALInferenceAddMinMaxPrior(priorArgs, "lambda1",     &lambda1Min, &lambda1Max,   LALINFERENCE_REAL8_t);
  
    ppt=LALInferenceGetProcParamVal(commandLine,"--fixLambda2");
    if(ppt){
    LALInferenceAddVariable(currentParams, "lambda2",           &start_lambda2,        LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
      fprintf(stdout,"phase fixed and set to %f\n",start_lambda2);
    }else{
      LALInferenceAddVariable(currentParams, "lambda2",           &start_lambda2,        LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
    }
    LALInferenceAddMinMaxPrior(priorArgs, "lambda2",     &lambda2Min, &lambda2Max,   LALINFERENCE_REAL8_t);
  }
  
  LALSimInspiralInteraction interactionFlags=LAL_SIM_INSPIRAL_INTERACTION_ALL;
  ppt=LALInferenceGetProcParamVal(commandLine,"--interactionFlags");
  if(ppt){
    if(strstr(ppt->value,"LAL_SIM_INSPIRAL_INTERACTION_NONE")) interactionFlags=LAL_SIM_INSPIRAL_INTERACTION_NONE;
    if(strstr(ppt->value,"LAL_SIM_INSPIRAL_INTERACTION_SPIN_ORBIT_15PN")) interactionFlags=LAL_SIM_INSPIRAL_INTERACTION_SPIN_ORBIT_15PN;
    if(strstr(ppt->value,"LAL_SIM_INSPIRAL_INTERACTION_SPIN_SPIN_2PN")) interactionFlags=LAL_SIM_INSPIRAL_INTERACTION_SPIN_SPIN_2PN;
    if(strstr(ppt->value,"LAL_SIM_INSPIRAL_INTERACTION_SPIN_SPIN_SELF_2PN")) interactionFlags=LAL_SIM_INSPIRAL_INTERACTION_SPIN_SPIN_SELF_2PN;
    if(strstr(ppt->value,"LAL_SIM_INSPIRAL_INTERACTION_QUAD_MONO_2PN")) interactionFlags=LAL_SIM_INSPIRAL_INTERACTION_QUAD_MONO_2PN;
    if(strstr(ppt->value,"LAL_SIM_INSPIRAL_INTERACTION_SPIN_ORBIT_25PN")) interactionFlags=LAL_SIM_INSPIRAL_INTERACTION_SPIN_ORBIT_25PN;
    if(strstr(ppt->value,"LAL_SIM_INSPIRAL_INTERACTION_TIDAL_5PN")) interactionFlags=LAL_SIM_INSPIRAL_INTERACTION_TIDAL_5PN;
    if(strstr(ppt->value,"LAL_SIM_INSPIRAL_INTERACTION_TIDAL_6PN")) interactionFlags=LAL_SIM_INSPIRAL_INTERACTION_TIDAL_6PN;
    if(strstr(ppt->value,"LAL_SIM_INSPIRAL_INTERACTION_ALL_SPIN")) interactionFlags=LAL_SIM_INSPIRAL_INTERACTION_ALL_SPIN;
    if(strstr(ppt->value,"LAL_SIM_INSPIRAL_INTERACTION_ALL")) interactionFlags=LAL_SIM_INSPIRAL_INTERACTION_ALL;
    LALInferenceAddVariable(currentParams, "interactionFlags", &interactionFlags,        LALINFERENCE_UINT4_t, LALINFERENCE_PARAM_FIXED); 
 }

  return;
}

/** Initialise student-t extra variables, set likelihood */
void initStudentt(LALInferenceRunState *state)
{	
        char help[]="\
Student T Likelihood Arguments:\n\
(--studentTLikelihood)\tUse student-t likelihood function\n";
	
	ProcessParamsTable *ppt=NULL;
	LALInferenceIFOData *ifo=state->data;

	/* Print command line arguments if help requested */
        if(LALInferenceGetProcParamVal(state->commandLine,"--help"))
        {
                fprintf(stdout,"%s",help);
		while(ifo) {
			fprintf(stdout,"(--dof-%s DoF)\tDegrees of freedom for %s\n",ifo->name,ifo->name);
			ifo=ifo->next;
		}
		return;
        }
	/* Don't do anything unless asked */
	if(!LALInferenceGetProcParamVal(state->commandLine,"--studentTLikelihood")) return;

	/* initialise degrees of freedom parameters for each IFO */
	while(ifo){
		CHAR df_argument_name[128];
		CHAR df_variable_name[64];
		REAL8 dof=10.0; /* Degrees of freedom parameter */
		
		sprintf(df_argument_name,"--dof-%s",ifo->name);
		if((ppt=LALInferenceGetProcParamVal(state->commandLine,df_argument_name)))
			dof=atof(ppt->value);
    		sprintf(df_variable_name,"df_%s",ifo->name);
    		LALInferenceAddVariable(state->currentParams,df_variable_name,&dof,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_FIXED);
		fprintf(stdout,"Setting %lf degrees of freedom for %s\n",dof,ifo->name);
		ifo=ifo->next;
	}

	/* Set likelihood to student-t */
	state->likelihood = &LALInferenceFreqDomainStudentTLogLikelihood;
	
	/* Set the noise model evidence to the student t model value */
	LALInferenceTemplateNullFreqdomain(state->data);
	REAL8 noiseZ=LALInferenceFreqDomainStudentTLogLikelihood(state->currentParams,state->data,&LALInferenceTemplateNullFreqdomain);
	LALInferenceAddVariable(state->algorithmParams,"logZnoise",&noiseZ,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_FIXED);
	fprintf(stdout,"Student-t Noise evidence %lf\n",noiseZ);

	return;
}

/*************** MAIN **********************/


int main(int argc, char *argv[]){
        char help[]="LALInferenceNest:\n\
Bayesian analysis tool using MultiNest algorithm\n\
for CBC analysis. Uses LALInference library for back-end.\n\n\
Arguments for each section follow:\n\n";

	LALInferenceRunState *state;
	ProcessParamsTable *procParams=NULL;

	/* Read command line and parse */
	procParams=LALInferenceParseCommandLine(argc,argv);
	
	/* initialise runstate based on command line */
	/* This includes reading in the data */
	/* And performing any injections specified */
	/* And allocating memory */
	state = initialize(procParams);
	
	/* Set up structures for MultiNest */
	initializeMN(state);
	
	/* Set up currentParams with variables to be used */
	initVariables(state);
	
	/* Check for student-t and apply */
	initStudentt(state);

       /* Print command line arguments if help requested */
        if(LALInferenceGetProcParamVal(state->commandLine,"--help"))
        {
                fprintf(stdout,"%s",help);
		exit(0);
        }

	/* Call MultiNest algorithm */
	state->algorithm(state);

	/* end */
	return 0;
}
