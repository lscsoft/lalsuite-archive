/* Nested sampling algorithm */
/* And support functions */
/* (C) John Veitch 2009 */

#define LAL_USE_OLD_COMPLEX_STRUCTS
#include <lal/Units.h>
#include <lal/LALStdlib.h>
#include "LALInspiralMCMC.h"
#include "LALInspiralMCMCUser.h"
#include <lal/LALError.h>
#include <lal/TimeDelay.h>
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_linalg.h>
#include "nest_calc.h"
#include <float.h>

//#define TOLERANCE 0.1
#define TOLERANCE 0.01

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

gsl_matrix *cov_mat;

CHAR outfile[FILENAME_MAX];
double etawindow;

INT4 seed;

double logadd(double a,double b){
	if(a>b) return(a+log(1.0+exp(b-a)));
	else return(b+log(1.0+exp(a-b)));
}

void NestInit2PN(LALMCMCParameter *parameter, void *iT){
	REAL8 trg_time;
	SnglInspiralTable *inspiralTable = (SnglInspiralTable *)iT;
	REAL4 UNUSED mtot, eta, UNUSED mwindow;
	trg_time = (REAL8)inspiralTable->end_time.gpsSeconds + (REAL8)inspiralTable->end_time.gpsNanoSeconds *1.0e-9;
	parameter->param=NULL;
	parameter->dimension=0;
	mwindow=0.5;
	REAL8 m1=inspiralTable->mass1;
	REAL8 m2=inspiralTable->mass2;
	mtot=inspiralTable->mass1 + inspiralTable->mass2;
	eta=inspiralTable->eta;
	if (eta==0.0) eta=m1*m2/((m1+m2)*(m1+m2));
	double etamin = eta-0.5*etawindow;
	etamin = etamin<0.01?0.01:etamin;
	double etamax = eta+0.5*etawindow;
	etamax = etamax>0.25?0.25:etamax;
	etamin=0.05; etamax=0.25;
	REAL4 localetawin=etamax-etamin;
	fprintf(stderr,"clipped eta prior %lf\n",localetawin);
	
	XLALMCMCAddParam(parameter,"mtotal",2.0+33.0*gsl_rng_uniform(RNG),2.0,35.0,0); /* Low mass */
	/*XLALMCMCAddParam(parameter, "mtotal",	mtot*(1.0 + mwindow*(gsl_rng_uniform(RNG)-0.5)),mtot*(1.0-0.5*mwindow),mtot*(1.0+0.5*mwindow),0);*/
	/*XLALMCMCAddParam(parameter, "eta", gsl_rng_uniform(RNG)*0.25 , 0, 0.25, 0);*/
	XLALMCMCAddParam(parameter, "eta", gsl_rng_uniform(RNG)*localetawin+etamin , etamin, etamax, 0);
	/*XLALMCMCAddParam(parameter, "eta",	eta*(1.0 + etawindow*(gsl_rng_uniform(RNG)-0.5)),eta*(1.0-0.5*etawindow),eta*(1.0+0.5*etawindow),0);*/
	XLALMCMCAddParam(parameter, "time",		(gsl_rng_uniform(RNG)-0.5)*timewindow + trg_time ,trg_time-0.5*timewindow,trg_time+0.5*timewindow,0);
	XLALMCMCAddParam(parameter, "phi",		LAL_TWOPI*gsl_rng_uniform(RNG),0.0,LAL_TWOPI,1);
	XLALMCMCAddParam(parameter, "distMpc", 99.0*gsl_rng_uniform(RNG)+1.0, 1.0, 100.0, 0);
	XLALMCMCAddParam(parameter,"ra",LAL_TWOPI*gsl_rng_uniform(RNG),0,LAL_TWOPI,1);
	XLALMCMCAddParam(parameter,"dec",LAL_PI*(gsl_rng_uniform(RNG)-0.5),-LAL_PI*0.5,LAL_PI*0.5,0);
	XLALMCMCAddParam(parameter,"psi",0.5*LAL_PI*gsl_rng_uniform(RNG),0,LAL_PI*0.5,0);
	XLALMCMCAddParam(parameter,"iota",LAL_PI*gsl_rng_uniform(RNG),0,LAL_PI,0);
	
	return;
}

void Inject2PN(LALMCMCParameter *parameter, LALMCMCInput *inputMCMC, double SNR){
	static LALStatus status;
	InspiralTemplate template;
	REAL8 real,imag,chisq;
	UINT4 Nmodel;
	UINT4 i;
	double SNR1,mul_factor=1.0;
	
	memset(&template,0,sizeof(InspiralTemplate));
	template.totalMass = XLALMCMCGetParameter(parameter,"mtotal");
	template.eta = XLALMCMCGetParameter(parameter,"eta");
	template.massChoice = totalMassAndEta;
	template.fLower = inputMCMC->fLow;
	template.distance = XLALMCMCGetParameter(parameter,"distMpc");
	template.order = LAL_PNORDER_TWO;
	template.approximant=inputMCMC->approximant;
	template.tSampling = 1.0/inputMCMC->deltaT;
	template.fCutoff = 0.5/inputMCMC->deltaT -1.0;
	template.nStartPad = 0;
	template.nEndPad =0;
	template.startPhase = XLALMCMCGetParameter(parameter,"phi");
	template.startTime = XLALMCMCGetParameter(parameter,"time");
	template.startTime -= inputMCMC->stilde[0]->epoch.gpsSeconds + 1e-9*inputMCMC->stilde[0]->epoch.gpsNanoSeconds;
	template.ieta = 1;
	template.next = NULL;
	template.fine = NULL;
	
	LALInspiralParameterCalc(&status,&template);
	
	template.startTime-=template.tC;
	
	LALInspiralRestrictedAmplitude(&status,&template);
	Nmodel=inputMCMC->stilde[0]->data->length*2; /* *2 for real/imag packing format */
	
	LALCreateVector(&status,&model,Nmodel);
	
	LALInspiralWave(&status,model,&template); /* Create the model */
	
	int lowBin = (int)(inputMCMC->fLow / inputMCMC->stilde[0]->deltaF);
	/* Find the SNR of this wave */
	for(chisq=0.0,i=lowBin;i<Nmodel/2;i++){
		real=model->data[i]; imag=model->data[Nmodel-i];
		chisq+=(real*real + imag*imag)*inputMCMC->invspec[0]->data->data[i];
	}
	SNR1=sqrt(chisq);
	if(SNR>=0.0) mul_factor = SNR/SNR1; /* Multiplicative factor to achieve desired SNR */
	/* Inject the wave */
	for(chisq=0.0,i=lowBin;i<Nmodel/2;i++){
		inputMCMC->stilde[0]->data->data[i].re+=mul_factor*(REAL8) model->data[i];
		inputMCMC->stilde[0]->data->data[i].im+=mul_factor*(REAL8) model->data[Nmodel-i];
		real=model->data[i]; imag=model->data[Nmodel-i];
		chisq+=mul_factor*mul_factor*(real*real + imag*imag)*inputMCMC->invspec[0]->data->data[i];
	}
	fprintf(stderr,"Injected wave, SNR = %e\n",sqrt(chisq));
	model=NULL;
	/*	free(model); */
}


REAL8 computeZ(LALMCMCInput *MCMCinput)
{
  UINT4 i=0;
  UINT4 j;

  REAL8 logZnoise=0.0;

  topdown_sum=calloc((size_t)MCMCinput->numberDataStreams,sizeof(REAL8Vector *));
  for (i=0;i<MCMCinput->numberDataStreams;i++) {
    topdown_sum[i]=XLALCreateREAL8Vector(MCMCinput->stilde[i]->data->length);
    topdown_sum[i]->data[topdown_sum[i]->length-1] = (pow(MCMCinput->stilde[i]->data->data[topdown_sum[i]->length-1].re,2.0)+pow(MCMCinput->stilde[i]->data->data[topdown_sum[i]->length-1].im,2.0))*MCMCinput->invspec[i]->data->data[topdown_sum[i]->length-1];
    for(j=topdown_sum[i]->length-2;j>0;j--) {
      topdown_sum[i]->data[j]=topdown_sum[i]->data[j+1]+(pow(MCMCinput->stilde[i]->data->data[j].re,2.0)+pow(MCMCinput->stilde[i]->data->data[j].im,2.0))*MCMCinput->invspec[i]->data->data[j];
    }
  }

  /* Likelihood of the noise model */
  logZnoise=0.0;
  for (j=0;j<MCMCinput->numberDataStreams;j++){
    int lowBin=(int)MCMCinput->fLow/MCMCinput->deltaF;
    logZnoise+=topdown_sum[j]->data[lowBin];
  }
  logZnoise*=-2.0*MCMCinput->deltaF;

  return logZnoise;
}

REAL8 logLmax=-DBL_MAX;


REAL8 nestZ(UINT4 Nruns, UINT4 Nlive, LALMCMCParameter **Live, LALMCMCInput *MCMCinput)
{
	UINT4 i=0;
	UINT4 j,minpos;
	static LALStatus status;
	REAL4 accept;
	REAL8 *logZarray,*logwarray,*Harray,*oldZarray,*Wtarray;
	REAL8 logw, H=0.0, logLmin, UNUSED logWt, logZ=-DBL_MAX, logZnew, UNUSED deltaZ;
	REAL8 MCMCfail=0;
	REAL8 logZnoise=0.0;
//	REAL8 logLmax=-DBL_MAX;
	FILE *fpout=NULL;
	CHAR outEnd[FILENAME_MAX];
	LALMCMCParameter *temp=(LALMCMCParameter *)malloc(sizeof(LALMCMCParameter));
	LALMCMCParam *param_ptr;

	if(!(MCMCinput->randParams)) LALCreateRandomParams(&status,&(MCMCinput->randParams),seed);
	
	MCMCinput->Live=Live;
	MCMCinput->Nlive=(UINT4)Nlive;
	/* Initialise the RNG */
	gsl_rng_env_setup();
	RNG=gsl_rng_alloc(gsl_rng_default);
	gsl_rng_set(RNG,seed==0 ? (unsigned long int)time(NULL) : (unsigned int long) seed);
	
	/* Initialise the optimised tables declared in LALInspiralMCMCUser.h*/
	
	normalisations = calloc((size_t)MCMCinput->numberDataStreams,sizeof(REAL8));
	
	for(i=0;i<MCMCinput->numberDataStreams;i++) {
		normalisations[i]=0.0;
		if(MCMCinput->invspec[i]!=NULL){
			for(j=(int)(MCMCinput->fLow/MCMCinput->invspec[i]->deltaF);j<MCMCinput->invspec[i]->data->length-1;j++) normalisations[i]+=0.5*log(MCMCinput->invspec[i]->data->data[j]);
		}
	}
	
	topdown_sum=calloc((size_t)MCMCinput->numberDataStreams,sizeof(REAL8Vector *));
	for (i=0;i<MCMCinput->numberDataStreams;i++){
		topdown_sum[i]=XLALCreateREAL8Vector(MCMCinput->stilde[i]->data->length);
		topdown_sum[i]->data[topdown_sum[i]->length-1]=
		(pow(MCMCinput->stilde[i]->data->data[topdown_sum[i]->length-1].re,2.0)+pow(MCMCinput->stilde[i]->data->data[topdown_sum[i]->length-1].im,2.0))*MCMCinput->invspec[i]->data->data[topdown_sum[i]->length-1];
		for(j=topdown_sum[i]->length-2;j>0;j--) topdown_sum[i]->data[j]=topdown_sum[i]->data[j+1]+(pow(MCMCinput->stilde[i]->data->data[j].re,2.0)+pow(MCMCinput->stilde[i]->data->data[j].im,2.0))*MCMCinput->invspec[i]->data->data[j];
	}
	
    if (MCMCinput->injectionTable!=NULL){		
        LALMCMCParameter *injected  =(LALMCMCParameter *)malloc(sizeof(LALMCMCParameter));    
        NestInitInjectedParam(injected,(void *)MCMCinput->injectionTable, MCMCinput);
        REAL8 injlogL=MCMCinput->funcLikelihood(MCMCinput,injected);
        fprintf(stdout,"Injected logL  %10.15e \n",injlogL);  
        XLALMCMCDestroyPara(&injected);     
     }
     //exit(1);
	if(MCMCinput->injectionTable!=NULL) MCMCinput->funcInit(temp,(void *)MCMCinput->injectionTable);
	else MCMCinput->funcInit(temp,(void *)MCMCinput->inspiralTable);
	if(!PriorIsSane(temp))
		{fprintf(stderr,"ERROR: Prior is not sane, check ranges specified\n"); exit(1);}
   
	/* Likelihood of the noise model */
	logZnoise=0.0;
	for (j=0;j<MCMCinput->numberDataStreams;j++){
		int lowBin=(int)MCMCinput->fLow/MCMCinput->deltaF;
		logZnoise+=topdown_sum[j]->data[lowBin];
	}
	logZnoise*=-2.0*MCMCinput->deltaF;
	
	fprintf(stdout,"Noise evidence: %lf\n",logZnoise);
	fprintf(stderr,"Sprinkling initial points, may take a while");
	/* Set up the parameters for the live points */
	for(i=0;i<Nlive;i++) {
		do{
			if(MCMCinput->injectionTable!=NULL) MCMCinput->funcInit(Live[i],(void *)MCMCinput->injectionTable);
			else MCMCinput->funcInit(Live[i],(void *)MCMCinput->inspiralTable);
			MCMCinput->dim=Live[i]->dimension;
			MCMCinput->funcPrior(MCMCinput,Live[i]);
			if(!PriorIsSane(Live[i]))
				{fprintf(stderr,"ERROR: Prior is not sane, check ranges specified\n"); exit(1);}
			if(Live[i]->logPrior==-DBL_MAX) XLALMCMCFreePara(Live[i]);
		} while(Live[i]->logPrior==-DBL_MAX);
		MCMCinput->funcLikelihood(MCMCinput,Live[i]);
	}
	/* Set up covariance matrix */
	cov_mat = gsl_matrix_alloc(MCMCinput->dim,MCMCinput->dim);
	calcCVM(cov_mat,Live,Nlive);
	
	for(i=0;i<Nlive;i++) {
		accept=MCMCSampleLimitedPrior(Live[i],temp,MCMCinput,-DBL_MAX,cov_mat,MCMCinput->numberDraw,10.0);
		if(i%50==0)fprintf(stderr,".");
        Live[i]->dZ=1000.0;//logadd(logZ,logLmax+logw+log((REAL8)Nlive))-logZ;
        Live[i]->accept=accept;
	}
	if(MCMCinput->verbose) fprintf(stderr,"Set up %i live points\n",Nlive);
	
	/* Find max likelihood currently */
	for(i=1,logLmax=Live[0]->logLikelihood;i<Nlive;i++) logLmax=logLmax>Live[i]->logLikelihood ? logLmax : Live[i]->logLikelihood;
	
	/* Set up arrays for parallel runs */
	logZarray = calloc(Nruns,sizeof(REAL8));
	oldZarray = calloc(Nruns,sizeof(REAL8));
	Harray = calloc(Nruns,sizeof(REAL8));
	logwarray = calloc(Nruns,sizeof(REAL8));
	Wtarray = calloc(Nruns,sizeof(REAL8));
	if(logZarray==NULL || Harray==NULL || oldZarray==NULL || logwarray==NULL) {fprintf(stderr,"Unable to allocate RAM\n"); exit(-1);}
	
	logw=log(1.0-exp(-1.0/Nlive));
	for(i=0;i<Nruns;i++)  {logwarray[i]=logw; logZarray[i]=-DBL_MAX; oldZarray[i]=-DBL_MAX; Harray[i]=0.0;}
	i=0;
	
    /* Write list of parameter names */
    sprintf(outEnd,"%s_params.txt",outfile);
    fpout=fopen(outEnd,"w");
    for(param_ptr=temp->param;param_ptr;param_ptr=param_ptr->next)
    {
        fprintf(fpout,"%s\t",param_ptr->core->name);
    }
    fprintf(fpout,"logl");
    fclose(fpout);

	/* open outfile */
	fpout=fopen(outfile,"w");
	if(fpout==NULL) {
		fprintf(stderr,"Unable to open output file %s\n",outfile);
		exit(1);
	}
	if(setvbuf(fpout,NULL,_IOFBF, 0x10000)) /* set output buffer to avoid NFS thrashing */
		fprintf(stderr,"Warning: unable to set output file buffer size\n");
	
	/* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- Nested sampling loop -=-=-=-=--=-=-=-=-==-=-=-=-=-=-= */
	/*	while(((REAL8)i)<=((REAL8)Nlive)*infosafe*H || i<3*Nlive) */
	deltaZ=1.0;
	/*	while((REAL8)i<=((REAL8)Nlive)*infosafe*H ? 1 : Nlive*fabs(deltaZ/logZ)>1e-6)*/
	/*while(((REAL8)i)<=((REAL8)Nlive) || logLmax+logw > logZ-5)*/  /* This termination condition: when remaining prior can't
	 account for more than exp(-5) of the evidence, even
	 if entire support is at Lmax */
        while(((REAL8)i)<=((REAL8)Nlive) || logadd(logZ,logLmax-((double)i/(double)Nlive))-logZ > TOLERANCE )
	{
		minpos=0;
		/* Find minimum likelihood sample to replace */
		for(j=0;j<Nlive;j++) {if(Live[j]->logLikelihood <= Live[minpos]->logLikelihood) minpos=j;}
		logLmin = Live[minpos]->logLikelihood;
		logWt = logw + Live[minpos]->logLikelihood;
		fprintSample(fpout,Live[minpos]);
		/* update evidence array */
		for(j=0;j<Nruns;j++){
			logZarray[j]=logadd(logZarray[j],Live[minpos]->logLikelihood + logwarray[j]);
			Wtarray[j]=logwarray[j]+Live[minpos]->logLikelihood;
			Harray[j]= exp(Wtarray[j]-logZarray[j])*Live[minpos]->logLikelihood
			+ exp(oldZarray[j]-logZarray[j])*(Harray[j]+oldZarray[j])-logZarray[j];
		}
		logZnew=mean(logZarray,Nruns);
		deltaZ=logZnew-logZ;
		H=mean(Harray,Nruns);
		logZ=logZnew;
		for(j=0;j<Nruns;j++) oldZarray[j]=logZarray[j];
		MCMCfail=0.0;
		
		/* Replace the previous sample by a new one */
		
		/* Update covariance matrix every so often */
		if(!(i%(Nlive/4)))	calcCVM(cov_mat,Live,Nlive);
		/* generate new live point */
/*UINT4 kk;       
 for(kk=0;kk<Nlive;kk++) {
		printf("before LP i=%i logL=%10.10e dist=%10.10e \n",kk,Live[kk]->logLikelihood,XLALMCMCGetParameter(Live[kk],"distMpc"));
	}
    int once=0;*/
		do{
			while((j=gsl_rng_uniform_int(RNG,Nlive))==minpos){};
           // printf("Changing logMin=%10.20e address %p j=%d minpos=%d\n",logLmin,&(Live[minpos]),j,minpos);
			XLALMCMCCopyPara(&(Live[minpos]),Live[j]);
           	accept = MCMCSampleLimitedPrior(Live[minpos],temp,MCMCinput,logLmin,cov_mat,MCMCinput->numberDraw,logadd(logZ,logLmax+logw+log((REAL8)Nlive))-logZ);
			MCMCinput->funcLikelihood(MCMCinput,Live[minpos]);
			MCMCinput->funcPrior(MCMCinput,Live[minpos]);
			MCMCfail+=1.0;
            /*if (accept==0 && once==0) 
            {fprintf(stdout,"accept end NS =%lf \n",accept);
            once=1;}
            else fprintf(stdout,"accept end NS =%lf but not zero\n",accept);*/
		}while((Live[minpos]->logLikelihood<=logLmin)||(accept==0.0));
     
    /*    for(kk=0;kk<Nlive;kk++) {
		printf("after LP i=%i logL=%10.10e dist=%10.10e \n",kk,Live[kk]->logLikelihood,XLALMCMCGetParameter(Live[kk],"distMpc"));
	}
    once=0 ;
    for(kk=0;kk<Nlive;kk++) {
       
	 if (fabs((REAL8) ((REAL8)Live[kk]->logLikelihood+ (REAL8)35011.8127710000 ))<5.E-5) {once++;printf("logL FOUND %d ##########################################\n",kk);}
       
	}
    
    if (once==1) 
    printf("1\n");
    else if (once==2)
        printf("2\n");
        else if(once==3)
            printf("3\n");
            else if(once==4)
                exit(1);
      */          
		if(Live[minpos]->logLikelihood > logLmax) logLmax = Live[minpos]->logLikelihood;
		for(j=0;j<Nruns;j++) logwarray[j]+=sample_logt(Nlive);
		logw=mean(logwarray,Nruns);
        Live[minpos]->dZ=logadd(logZ,logLmax+logw+log((REAL8)Nlive))-logZ;
        Live[minpos]->accept=accept/MCMCfail;
		if(MCMCinput->verbose) fprintf(stderr,"%i: (%2.1lf%%) accpt: %1.3f H: %3.3lf nats (%3.3lf b) logL:%lf ->%lf dZ: %lf logZ: %lf Zratio: %lf db\n",
									   i,100.0*((REAL8)i)/(((REAL8) Nlive)*H),accept/MCMCfail,H,H/log(2.0),logLmin,Live[minpos]->logLikelihood,logadd(logZ,logLmax+logw+log((REAL8)Nlive))-logZ,logZ,10.0*log10(exp(1.0))*(logZ-logZnoise));
		if(fpout && !(i%50)) fflush(fpout);
		i++;
	}
	
	/* Sort the remaining points (not essential, just nice)*/
	for(i=0;i<Nlive-1;i++){
		minpos=i;
		logLmin=Live[i]->logLikelihood;
		for(j=i+1;j<Nlive;j++){
			if(Live[j]->logLikelihood<logLmin) {minpos=j; logLmin=Live[j]->logLikelihood;}
		}
		temp=Live[minpos]; /* Put the minimum remaining point in the current position */
		Live[minpos]=Live[i];
		Live[i]=temp;
	}
	
	/* final corrections */
	for(i=0;i<Nlive;i++){
		logZ=logadd(logZ,Live[i]->logLikelihood+logw);
		for(j=0;j<Nruns;j++){
			logwarray[j]+=sample_logt(Nlive);
			logZarray[j]=logadd(logZarray[j],Live[i]->logLikelihood+logwarray[j]);
		}
		fprintSample(fpout,Live[i]);
	}
	
	/* Output the maximum template, data, etc */
	sprintf(outEnd,"%s_maxLdata.dat",outfile);
	MCMCinput->dumpfile=outEnd;
	MCMCinput->funcLikelihood(MCMCinput,Live[Nlive-1]);

	/* Output some statistics */
	double Npoints = MCMCinput->numberDataStreams*MCMCinput->stilde[0]->data->length-(int)(MCMCinput->fLow/MCMCinput->deltaF);
	fprintf(stdout,"MaxL = %lf\nReduced chi squared = %lf\n",Live[Nlive-1]->logLikelihood,-Live[Nlive-1]->logLikelihood/Npoints);
	logZ=mean(logZarray,Nruns);
	fprintf(stdout,"deltaLmax = %lf\n",Live[Nlive-1]->logLikelihood-logZnoise);
	double zscore =( -2.0*Live[Nlive-1]->logLikelihood - Npoints) / sqrt(2.0*Npoints);
	fprintf(stdout,"Z-score = %lf\n",zscore);
	fclose(fpout);
	sprintf(outEnd,"%s_B.txt",outfile);
	fpout=fopen(outEnd,"w");
	fprintf(fpout,"%lf %lf %lf %lf %lf\n",logZ-logZnoise,logZ,logZnoise,Live[Nlive-1]->logLikelihood-logZnoise,zscore);
	fclose(fpout);
	free(Harray); free(logwarray); free(Wtarray); free(oldZarray); free(logZarray);
	fprintf(stdout,"lodds ratio %lf\n",logZ-logZnoise);
	return logZ;
	
}

REAL8 mean(REAL8 *array,int N){
	REAL8 sum=0.0;
	int i;
	for(i=0;i<N;i++) sum+=array[i];
	return sum/((REAL8) N);
}

REAL4 MCMCSampleLimitedPrior(LALMCMCParameter *sample, LALMCMCParameter *temp, LALMCMCInput *MCMCInput,REAL8 minL,gsl_matrix *covM,INT4 N,REAL8 dZ)
{
	/* Sample from prior using MCMC to evolve the existing value of sample subject to the new likelihood being >minL*/
	/* Returns acceptance ratio */
#define ROTATEFRAC 0.1
#define REFLECTFRAC 0.1
	int i=0;
	int a_cnt=0;
	int accept=0;
	int nreflect=0;
	REAL8 jump_select=0;
	int ret=0;

int sky_rotate_i=0;
int normal_jump_i=0;
int diffusive_i=0;
int spin_rotate_i=0;
int spin_jump_i=0;

int sky_rotate_t=0;
int normal_jump_t=0;
int diffusive_t=0;
int spin_rotate_t=0;
int spin_jump_t=0;

int sky_rotate_s=0;
int normal_jump_s=0;
int diffusive_s=0;
int spin_rotate_s=0;
int spin_jump_s=0;

int failed_logP=0;
    REAL8 spin_only_ratio=0.01;
    REAL8 sky_jump=0.05;
    REAL8 differential_ratio=0.2;
    REAL8 spin_jump_ratio=0.01;
    int dont_skip=1;
if (dZ> 100.0){differential_ratio=0.2;spin_only_ratio=0.15;sky_jump=0.15;}    
else if(dZ>80.0){dont_skip=0;differential_ratio=0.0;spin_only_ratio=.5;spin_jump_ratio=.5;sky_jump=0.5;}
else if (dZ<70.0 && dZ>60.0 ){dont_skip=0;differential_ratio=0.0;spin_only_ratio=.5;spin_jump_ratio=.5;sky_jump=0.5;}
else if (dZ>30.0){differential_ratio=0.2;spin_only_ratio=0.;sky_jump=0.;}
else if (dZ>20.0){dont_skip=0;differential_ratio=0.0;spin_only_ratio=.5;spin_jump_ratio=.5;sky_jump=0.5;}
else if (dZ<10.0 && dZ>8.0 ){dont_skip=0;differential_ratio=0.0;spin_only_ratio=.5;spin_jump_ratio=.5;sky_jump=0.5;}
else if (dZ>0.0) {spin_only_ratio=0.00;sky_jump=0.00;differential_ratio=0.4;}
	

MCMCInput->funcPrior(MCMCInput,sample);
	
	XLALMCMCCopyPara(&temp,sample);
	
	i=0;
	while (i<N || (nreflect==a_cnt && nreflect>0 && nreflect%2==0)){
        sky_rotate_i=0;
        normal_jump_i=0;
        diffusive_i=0;
        spin_rotate_i=0;
        spin_jump_i=0;
		i++;
		jump_select = gsl_rng_uniform(RNG);
		if(jump_select<sky_jump && MCMCInput->numberDataStreams>1 && XLALMCMCCheckWrapping(sample,"ra")!=-1 && XLALMCMCCheckWrapping(sample,"dec")!=-1){
			if(MCMCInput->numberDataStreams>1) jump_select = gsl_rng_uniform(RNG);
			else jump_select=0;
			if(jump_select>0.5) {
			    ret=0;
			    /* Custom proposal for symmetry of the sky sphere */
			    if(MCMCInput->numberDataStreams>=3 && jump_select>0.9){ 
                    sky_rotate:
                    ret=XLALMCMCReflectDetPlane(MCMCInput,temp);sky_rotate_t+=1;sky_rotate_i=1; 
                    //printf("Reflect Sky proposal iteration %d\n",i);
                }
                else { //printf("Im in the pitch iteration %d\n",i);
                    if (dont_skip==1) goto diff_jump;
                    else {goto spin_jump;}}
			    if(ret==0) nreflect++;
			    if(ret==-1 || jump_select <=0.9 || MCMCInput->numberDataStreams==2) {XLALMCMCRotateSky(MCMCInput,temp);
                //printf("RotateSky 2 proposal iteration %d\n",i);
                }
			    /* Custom proposal for mass/eta1 surface of constant 1PN term */
			    /* if(gsl_rng_uniform(RNG)<0.1) XLALMCMC1PNMasseta(MCMCInput,temp); */
			}
            else if (dont_skip!=1) goto spin_jump;
            else { diff_jump: {XLALMCMCDifferentialEvolution(MCMCInput,temp);diffusive_t+=1;diffusive_i=1;
                   //printf("Differential proposal iteration %d\n",i);
                   }
            }
		}
		else    
		{
		    if( (jump_select=gsl_rng_uniform(RNG))<differential_ratio/*0.2*/){ XLALMCMCDifferentialEvolution(MCMCInput,temp);diffusive_t+=1;diffusive_i=1;
                    //printf("Differential2 proposal iteration %d\n",i);
            }
		    else {
		    /* Check for higher harmonics present */
		        if((jump_select=gsl_rng_uniform(RNG))<0.1 && MCMCInput->ampOrder!=0)
                    {XLALMCMCJumpHarmonic(MCMCInput,temp);printf("Higher Harmonic  proposal iteration %d\n",i);}
                    /* Spin jumps */
                    // else if ((jump_select=gsl_rng_uniform(RNG))<spin_only_ratio && XLALMCMCCheckParameter(sample,"a1") && XLALMCMCCheckParameter(sample,"a1")){
                else if ((jump_select=gsl_rng_uniform(RNG))<spin_only_ratio && XLALMCMCCheckWrapping(sample,"theta1")!=-1 && XLALMCMCCheckWrapping(sample,"theta2")!=-1 && XLALMCMCCheckWrapping(sample,"phi1")!=-1 && XLALMCMCCheckWrapping(sample,"phi2")!=-1 && XLALMCMCCheckWrapping(sample,"a1")!=-1 && XLALMCMCCheckWrapping(sample,"a2")!=-1)
                    {spin_jump:
                    XLALMCMCRotateSpins(MCMCInput,temp);spin_rotate_t+=1;spin_rotate_i=1; 
                     //printf("Rotate Spins iteration %d ============= SPINS \n",i);
                }
                else if ((jump_select=gsl_rng_uniform(RNG))<spin_jump_ratio && XLALMCMCCheckWrapping(sample,"theta1")!=-1 && XLALMCMCCheckWrapping(sample,"theta2")!=-1 && XLALMCMCCheckWrapping(sample,"phi1")!=-1 && XLALMCMCCheckWrapping(sample,"phi2")!=-1 && XLALMCMCCheckWrapping(sample,"a1")!=-1 && XLALMCMCCheckWrapping(sample,"a2")!=-1)
                {
                    XLALMCMCJumpSpins(MCMCInput,temp,covM);spin_jump_t+=1;spin_jump_i=1; 
                     //printf("Rotate Spins iteration %d ============= SPINS \n",i);
                }
                else if (dont_skip!=1) {goto sky_rotate;}
                else /* Otherwise just perform a regular jump */
                    {XLALMCMCJump(MCMCInput,temp,covM);normal_jump_t+=1;normal_jump_i=1;
                    //printf("Normal Jumps proposal iteration %d\n",i);
                }
            }
		}
		/* Evoluate the MH ratio */	
        //printf("Entering calculation prior temp \n");	
		MCMCInput->funcPrior(MCMCInput,temp);
		if(temp->logPrior!=-DBL_MAX && ( (temp->logPrior - sample->logPrior) > log(gsl_rng_uniform(RNG)) )) {
			/* this would be accepted based on the priors, we can now confirm that its likelihood is above the limit
			 with the more expensive calculation */
			MCMCInput->funcLikelihood(MCMCInput,temp);
			if(temp->logLikelihood>minL) accept = 1;
		}
        //else if (temp->logPrior==-DBL_MAX) {failed_logP+=1;}//printf("logPrior was -inf. Jump tried %d %d %d %d %d %d \n",sky_rotate_i,normal_jump_i,diffusive_i,spin_jump_i,spin_diffusive_i,spin_rotate_i);}

sky_rotate_s+=accept*sky_rotate_i;
normal_jump_s+=accept*normal_jump_i;
diffusive_s+=accept*diffusive_i;
spin_jump_s+=accept*spin_jump_i;
//spin_diffusive_s+=accept*spin_diffusive_i;
spin_rotate_s+=accept*spin_rotate_i;

		if(accept==1) {
        // printf("LogL GOOD PRIOR AND LOGL %d (old=%lf proposed %lf) ACCEPT:1\n",i,sample->logLikelihood,temp->logLikelihood);
            XLALMCMCCopyPara(&sample,temp); a_cnt++; accept=0;
            if(temp->logLikelihood > logLmax) logLmax = temp->logLikelihood;
                     
        }
		else {XLALMCMCCopyPara(&temp,sample); 
        //printf("LogL was REFUSED at the iteration %d (old=%lf proposed %lf) ACCEPT:%d\n",i,sample->logLikelihood,temp->logLikelihood,accept);
        if (i==N && a_cnt==0) printf("Could not find a better point in %d MCMC jumps\n",N);
        }
    }   
sample->sky_rotate_s=(sky_rotate_t==0? 41.: (REAL8) sky_rotate_s/ (REAL8)sky_rotate_t);
sample->normal_jump_s=(normal_jump_t==0?41.: (REAL8) normal_jump_s/ (REAL8)normal_jump_t);
sample->diffusive_s=(diffusive_t==0?41.:(REAL8)diffusive_s/ (REAL8)diffusive_t);
sample->spin_jump_s=(spin_jump_t==0?41.:(REAL8) spin_jump_s/ (REAL8)spin_jump_t);
sample->spin_diffusive_s=41.;//(spin_diffusive_t==0?41.:(REAL8) spin_diffusive_s/ (REAL8)spin_diffusive_t);
sample->spin_rotate_s=(spin_rotate_s==0?41.:(REAL8) spin_rotate_s/ (REAL8)spin_rotate_t);
//sample->failed_logP=failed_logP;

	return(((REAL4) a_cnt)/((REAL4) i));
}

void calcCVM(gsl_matrix *cvm, LALMCMCParameter **samples,UINT4 N)
{ UINT4 i,j,k;
	UINT4 ND=samples[0]->dimension;
	REAL8 *means;
	LALMCMCParam *p;
	LALMCMCParam *jp;
	LALMCMCParam *kp;
	gsl_matrix *oldcvm;
	
	oldcvm = gsl_matrix_alloc(ND,ND);
	gsl_matrix_memcpy (oldcvm, cvm);
	
	/* clear the matrix */
	for(i=0;i<cvm->size1;i++) for(j=0;j<cvm->size2;j++) gsl_matrix_set(cvm,i,j,0.0);
	
	/* Find the means */
	if(NULL==(means = malloc((size_t)ND*sizeof(REAL8)))){fprintf(stderr,"Can't allocate RAM"); exit(-1);}
	for(i=0;i<ND;i++) means[i]=0.0;
	for(i=0;i<N;i++){
		p=samples[i]->param;
		for(j=0;j<ND;j++) {if(p->core->wrapping==0) {means[j]+=p->value;} p=p->next;}
	}
	for(j=0;j<ND;j++) means[j]/=(REAL8)N;
	
	/* Find the (co)-variances */
	for(i=0;i<N;i++){
		kp=jp=p=samples[i]->param;
		for(j=0;j<ND;j++){
			for(k=0,kp=p;k<=j;k++){
				gsl_matrix_set(cvm,j,k,gsl_matrix_get(cvm,j,k) + (kp->value - means[k])*(jp->value - means[j]));
				kp=kp->next;
			}
			jp=jp->next;
		}
	}
	
	/* Normalise */
	for(i=0;i<ND;i++) for(j=0;j<ND;j++) gsl_matrix_set(cvm,i,j,gsl_matrix_get(cvm,i,j)/((REAL8) N));
	free(means);
	/* Fill in variances for angle parameters */
	for(p=samples[0]->param,j=0;j<ND;j++,p=p->next) {
		if(p->core->wrapping==1) {
			for(k=0;k<j;k++) gsl_matrix_set(cvm,j,k,0.0);
			gsl_matrix_set(cvm,j,j,ang_var(samples,p->core->name,N));
			for(k=j+1;k<ND;k++) gsl_matrix_set(cvm,k,j,0.0);
		}
		/* Set variance to default 1 for non-varying params */
		else if (p->core->wrapping==-1) {
			for(k=0;k<j;k++) gsl_matrix_set(cvm,j,k,0.0);
			gsl_matrix_set(cvm,j,j,1.0);
			for(k=j+1;k<ND;k++) gsl_matrix_set(cvm,k,j,0.0);
		}
	}
	
	/* the other half */
	for(i=0;i<ND;i++) for(j=0;j<i;j++) gsl_matrix_set(cvm,j,i,gsl_matrix_get(cvm,i,j));
	
	/*fprintf(stderr,"shrinkage: ");
	 for(i=0;i<ND;i++) fprintf(stderr,"%lf ",sqrt(gsl_matrix_get(oldcvm,i,i)/gsl_matrix_get(cvm,i,i)));
	 fprintf(stderr,"\n");
	 */
	/* Check the 2nd order first moment - indicates non-gaussian structure */
	/* double *twopt=calloc(ND,sizeof(double));
	 double *max=calloc(ND,sizeof(double));
	 double *min=calloc(ND,sizeof(double));
	 for(k=0;k<ND;k++) {min[k]=DBL_MAX; max[k]=-DBL_MAX;}
	 for(i=0;i<N;i++){
	 p=samples[i]->param;
	 for(k=0;k<ND;k++){
	 min[k]=min[k]<p->value?min[k]:p->value;
	 max[k]=max[k]>p->value?max[k]:p->value;
	 p=p->next;
	 }
	 }
	 for(i=0;i<N;i++) for(j=i+1;j<N;j++) {
	 p=samples[i]->param;
	 jp=samples[j]->param;
	 for(k=0;k<ND;k++) {
	 twopt[k]+=fabs(p->value-jp->value);
	 p=p->next; jp=jp->next;}
	 }
	 
	 
	 for(k=0;k<ND;k++) twopt[k]/= (1.0/(2.0*sqrt(2.0)))*((double)(N*(N-1)/2))*(max[k]-min[k]);
	 fprintf(stderr,"structure indicator: ");
	 for(k=0;k<ND;k++) fprintf(stderr,"%lf ",twopt[k]);
	 fprintf(stderr,"\n");
	 free(twopt);
	 free(max); free(min);
	 */
	/* Debug output the matrix
	 for(i=0;i<ND;i++){
	 for(j=0;j<ND;j++){ fprintf(stderr,"%e ",gsl_matrix_get(cvm,i,j));}
	 fprintf(stderr,"\n");
	 }
	 */
	
	return;
}

/* Calculate shortest angular distance between a1 and a2 */
REAL8 ang_dist(REAL8 a1, REAL8 a2){
	double raw = (a2>a1 ? a2-a1 : a1-a2);
	return(raw>LAL_PI ? 2.0*LAL_PI - raw : raw);
}

/* Calculate the variance of a modulo-2pi distribution */
REAL8 ang_var(LALMCMCParameter **list,const char *pname, int N){
	int i=0;
	REAL8 ang_mean=0.0;
	REAL8 var=0.0;
	REAL8 ms,mc;
	/* Calc mean */
	for(i=0,ms=0.0,mc=0.0;i<N;i++) {
		ms+=sin(XLALMCMCGetParameter(list[i],pname));
		mc+=cos(XLALMCMCGetParameter(list[i],pname));
	}
	ms/=N; mc/=N;
	ang_mean=atan2(ms,mc);
	ang_mean = ang_mean<0? 2.0*LAL_PI + ang_mean : ang_mean;
	/* calc variance */
	for(i=0;i<N;i++) var+=ang_dist(XLALMCMCGetParameter(list[i],pname),ang_mean)*ang_dist(XLALMCMCGetParameter(list[i],pname),ang_mean);
	return(var/(REAL8)N);
}

void fprintSample(FILE *fp,LALMCMCParameter *sample){
	LALMCMCParam *p=sample->param;
	if(fp==NULL) return;
	while(p!=NULL) {fprintf(fp,"%15.15lf\t",p->value); p=p->next;}
    fprintf(fp,"%lf\t",sample->dZ);
    fprintf(fp,"%lf\t",sample->accept);
    fprintf(fp,"%d\t",sample->failed_logP);    
    fprintf(fp,"%lf\t",sample->sky_rotate_s);
    fprintf(fp,"%lf\t",sample->normal_jump_s);
    fprintf(fp,"%lf\t",sample-> diffusive_s);
    fprintf(fp,"%lf\t",sample->spin_jump_s);
    fprintf(fp,"%lf\t",sample->spin_diffusive_s);
    fprintf(fp,"%lf\t",sample->spin_rotate_s);
	fprintf(fp,"%lf\n",sample->logLikelihood);
	return;
}

REAL8 sample_logt(int Nlive){
	REAL8 t=0.0;
	REAL8 a=0.0;
	while((Nlive--)>1) {a=gsl_rng_uniform(RNG); t = t>a ? t : a;}
	return(log(t));
}
