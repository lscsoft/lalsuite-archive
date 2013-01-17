/* Wrapper for MultiNest */
/* And support functions */
/* (C) Farhan Feroz 2010 */

#include "MultiNest_calc.h"

LALMCMCInput *MultiNestInput;
LALMCMCParameter *MultiNestParam;
double nullZ;

void MultiNestRun(int mmodal, int ceff, int nlive, double tol, double efr, int ndims, int nPar, int nClsPar,  int maxModes,
int updInt, double Ztol, char root[], int seed, int *pWrap, int fb, int resume, int outfile, int initMPI, double logZero, 
void (*LogLike)(double *, int *, int *, double *), void (*dumper)(int *, int *, int *, double **, double **, double *, 
double *, double *, double *), int context)
{
	int i;
	for (i = strlen(root); i < 100; i++) root[i] = ' ';

        __nested_MOD_nestrun(&mmodal, &ceff, &nlive, &tol, &efr, &ndims, 
&nPar, &nClsPar, &maxModes, &updInt, &Ztol,
        root, &seed, pWrap, &fb, &resume, &outfile, &initMPI, &logZero, LogLike, dumper, &context);
}

void LogLike(double *Cube, int *ndim, int *npars, double *lnew)
{
	// transform the parameter in the unit hypercube to their physical counterparts according to the prior
	int i = MultiNestInput->funcMultiNestPrior(Cube, MultiNestInput, MultiNestParam);

	// if the parameters violate the prior then set likelihood to log(0);
	if( i == 0 )
	{
		*lnew = -DBL_MAX;
		return;
	}
	
	// calculate the loglike
	MultiNestInput->funcLikelihood(MultiNestInput, MultiNestParam);
	*lnew = MultiNestParam->logLikelihood - nullZ;
	if( isnan(*lnew) )
	{
		printf("likelihood is NaN. Aborting!\n");
		abort();
	}
}

void dumper(int *nSamples, int *nlive, int *nPar, double **physLive, double **posterior, double *paramConstr, double *maxLogLike, double *logZ, double *logZerr)
{
}

void MultiNestZ(UINT4 Nlive, LALMCMCInput *MCMCinput)
{
	/*if( multinest_seg < 1 || multinest_seg > 4 )
	{
		fprintf(stdout,"multinest_seg can not be less then 1 or greater than 4\n");
		abort();
	}*/

	UINT4 i=0;
	UINT4 j;
	static LALStatus status;
	REAL8 logZnoise=0.0;
	LALMCMCParameter *MCMCParam=(LALMCMCParameter *)malloc(sizeof(LALMCMCParameter));
	
	if(!(MCMCinput->randParams)) LALCreateRandomParams(&status,&(MCMCinput->randParams),seed);
	
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
		//topdown_sum[i]->data[topdown_sum[i]->length-1]=
		//(pow(MCMCinput->stilde[i]->data->data[topdown_sum[i]->length-1].re,2.0)+pow(MCMCinput->stilde[i]->data->data[topdown_sum[i]->length-1].im,2.0))*MCMCinput->invspec[i]->data->data[topdown_sum[i]->length-1];
		//for(j=topdown_sum[i]->length-2;j>0;j--) topdown_sum[i]->data[j]=topdown_sum[i]->data[j+1]+(pow(MCMCinput->stilde[i]->data->data[j].re,2.0)+pow(MCMCinput->stilde[i]->data->data[j].im,2.0))*MCMCinput->invspec[i]->data->data[j];
	}
	
	if(MCMCinput->injectionTable!=NULL) MCMCinput->funcInit(MCMCParam,(void *)MCMCinput->injectionTable);
	else MCMCinput->funcInit(MCMCParam,(void *)MCMCinput->inspiralTable);
	MCMCinput->dim=MCMCParam->dimension;
	
	MultiNestInput = MCMCinput;
	MultiNestParam = MCMCParam;
	
	/* Likelihood of the noise model */
	logZnoise=0.0;
	for (j=0;j<MCMCinput->numberDataStreams;j++){
		int lowBin=(int)MCMCinput->fLow/MCMCinput->deltaF;
		logZnoise+=topdown_sum[j]->data[lowBin];
	}
	logZnoise*=-2.0*MCMCinput->deltaF;
	nullZ = logZnoise;
	
	fprintf(stdout,"Noise evidence: %lf\n",logZnoise);
	
	// call MultiNest
	
	int mmodal = 0;
	int ceff = 0;
	int nlive = Nlive;
	double efr = 0.01;
	double tol = 0.5;
	int ndims = MCMCinput->dim;
	int nPar = ndims + 2;
	int nClsPar = 2;
	int updInt = 100;
	double Ztol = -1.e90;
	int maxModes = 100;
	int pWrap[ndims];
	for( int k = 0; k < ndims; k++ ) pWrap[k] = 0;
	pWrap[1] = pWrap[3] = pWrap[4] = 1;
	char root[100];
	for( int k = 0; k < 100; k++ ) root[k] = outfile[k];
	int rseed = -1;
	int fb = 1;
	int resume = 1;
	int outfile = 1;		// write output files?
	int initMPI = 1;		// initialize MPI routines?, relevant only if compiling with MPI
	double logZero = -1E10;		// points with loglike < logZero will be ignored by MultiNest
	int context = 0;


	MultiNestRun(mmodal, ceff, nlive, tol, efr, ndims, nPar, nClsPar, maxModes, updInt, Ztol, root, rseed, pWrap, fb, resume, outfile, initMPI, logZero,
	LogLike, dumper, context);
	
}
