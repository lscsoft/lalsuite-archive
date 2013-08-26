/* Setup the variables to control template generation for the Ringdown model */
/* Includes specification of prior ranges. Sets runState->currentParams and
 returns address of new LALInferenceVariables */

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



LALInferenceVariables *LALInferenceInitRingdownVariables(LALInferenceRunState *state)
{

  char help[]="\
               \n\
               ------------------------------------------------------------------------------------------------------------------\n\
               --- Injection Arguments ------------------------------------------------------------------------------------------\n\
               ------------------------------------------------------------------------------------------------------------------\n\
               (--inj injections.xml)          Injection XML file to use.\n\
               (--event N)                     Event number from Injection XML file to use.\n\
               \n\
               ------------------------------------------------------------------------------------------------------------------\n\
               --- Template Arguments -------------------------------------------------------------------------------------------\n\
               ------------------------------------------------------------------------------------------------------------------\n\
               (--symMassRatio)                Jump in symmetric mass ratio eta, instead of q=m2/m1.\n\
               (--logdistance)                 Jump in log(distance) instead of distance.\n\
               (--system-frame                 Jump in spin parameters defined in the system coordinates, relative to total angular momentum\n\
               (--approx)                      Specify a template approximant to use.\n\
                                               (default RingdownTD). Available approximants:\n\
                                               default modeldomain=\"time\": RingdownTD, RingdownLeaver, RingdownNumRel.\n\
                                               default modeldomain=\"frequency\": RingdownFD.\n\
               (--fref fRef)                   Specify a reference frequency at which parameters are defined (default 0).\n\
               (--modeldomain)                 domain the waveform template will be computed in (\"time\" or \"frequency\").\n\
               \n\
               ------------------------------------------------------------------------------------------------------------------\n\
               --- Starting Parameters ------------------------------------------------------------------------------------------\n\
               ------------------------------------------------------------------------------------------------------------------\n\
               (--trigtime time)               Trigger time to use.\n\
               (--time time)                   Waveform time (overrides random about trigtime).\n\
               (--mass mass)                   Trigger mass (total mass of the final black hole).\n\
               (--a a_spin)                    Trigger a (dimensionless spin) to use.\n\
               (--eta eta)                     Trigger eta (symmetric mass ratio) to use.\n\
               (--q q)                         Trigger q (asymmetric mass ratio) to use.\n\
               (--phi phase)                   Trigger phase to use.\n\
               (--iota inclination)            Trigger inclination to use.\n\
               (--dist dist)                   Trigger distance.\n\
               (--ra ra)                       Trigger RA.\n\
               (--dec dec)                     Trigger declination.\n\
               (--psi psi)                     Trigger psi.\n\
               \n\
               ------------------------------------------------------------------------------------------------------------------\n\
               --- Prior Arguments ----------------------------------------------------------------------------------------------\n\
               ------------------------------------------------------------------------------------------------------------------\n\
               (--mass-min massMin)            Minimum mass.\n\
               (--mass-max massMax)            Maximum mass.\n\
               (--eta-min etaMin)              Minimum eta.\n\
               (--eta-max etaMax)              Maximum eta.\n\
               (--q-min qMin)                  Minimum q.\n\
               (--q-max qMax)                  Maximum q.\n\
               (--a-min max)                   Minimum component spin (-1.0).\n\
               (--a-max max)                   Maximum component spin (1.0).\n\
               (--iota-max max)                Maximum inclination angle (pi).\n\
               (--Dmin dist)                   Minimum distance in Mpc (1).\n\
               (--Dmax dist)                   Maximum distance in Mpc (100).\n\
               (--dt time)                     Width of time prior, centred around trigger (0.1s).\n\
               \n\
               ------------------------------------------------------------------------------------------------------------------\n\
               --- Fix Parameters -----------------------------------------------------------------------------------------------\n\
               ------------------------------------------------------------------------------------------------------------------\n\
               (--fixMass)                     Do not allow final BH to vary.\n\
               (--fixEta)                      Do not allow mass ratio to vary.\n\
               (--fixQ)                        Do not allow mass ratio to vary.\n\
               (--fixPhi)                      Do not allow phase to vary.\n\
               (--fixIota)                     Do not allow inclination to vary.\n\
               (--fixDist)                     Do not allow distance to vary.\n\
               (--fixRa)                       Do not allow RA to vary.\n\
               (--fixDec)                      Do not allow declination to vary.\n\
               (--fixPsi)                      Do not allow polarization to vary.\n\
               (--fixTime)                     Do not allow coalescence time to vary.\n\
               (--fixSpin)                     Do not allow final BH spin (a) to vary.\n\
               (--pinparams)                   List of parameters to set to injected values [mchirp,asym_massratio,etc].\n";


  /* Print command line arguments if state was not allocated */
  if(state==NULL)
    {
      fprintf(stdout,"%s",help);
      return(NULL);
    }

  /* Print command line arguments if help requested */
  if(LALInferenceGetProcParamVal(state->commandLine,"--help"))
    {
      fprintf(stdout,"%s",help);
      return(NULL);
    }


  LALStatus status;
  memset(&status,0,sizeof(status));
  // int errnum;
  SimInspiralTable *injTable=NULL;
  LALInferenceVariables *priorArgs=state->priorArgs;
  state->currentParams=XLALCalloc(1,sizeof(LALInferenceVariables));
  LALInferenceVariables *currentParams=state->currentParams;
  ProcessParamsTable *commandLine=state->commandLine;
  ProcessParamsTable *ppt=NULL;
  // ProcessParamsTable *ppt_order=NULL;
  Approximant approx=NumApproximants;
  REAL8 fRef = 0.0;
  LALInferenceApplyTaper bookends = LALINFERENCE_TAPER_NONE;
  UINT4 analytic=0;
  LALInferenceIFOData *dataPtr;
  LALSimulationDomain modelDomain;
  UINT4 event=0;
  UINT4 i=0;
  REAL8 Dmin=1.0;
  REAL8 Dmax=100.0;
  REAL8 mMin=2.0,mMax=30.0;
  REAL8 etaMin=0.0312;
  REAL8 etaMax=0.25;
  REAL8 qMin=1.0;
  REAL8 qMax=30.0;
  REAL8 iotaMin=0.0,iotaMax=LAL_PI;
  REAL8 psiMin=0.0,psiMax=LAL_PI;
  REAL8 decMin=-LAL_PI/2.0,decMax=LAL_PI/2.0;
  REAL8 raMin=0.0,raMax=LAL_TWOPI;
  REAL8 phiMin=0.0,phiMax=LAL_TWOPI;
  REAL8 amin=-1.0,amax=1.0;
  REAL8 dt=0.1;            /* Width of time prior */
  gsl_rng *GSLrandom=state->GSLrandom;
  REAL8 starttime=0.0, timeParam=0.0;
  REAL8 timeMin=starttime-dt,timeMax=starttime+dt;

  memset(currentParams,0,sizeof(LALInferenceVariables));

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
    }
  }

  /* See if there are any parameters pinned to injection values */
  if((ppt=LALInferenceGetProcParamVal(commandLine,"--pinparams"))){
    char *pinned_params=ppt->value;
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
      fprintf(stdout,"Pinning parameter %s\n",node->name);
      node=LALInferenceGetItem(&tempParams,name);
      if(node) LALInferenceAddVariable(currentParams,node->name,node->value,node->type,node->vary);
      else {fprintf(stderr,"Error: Cannot pin parameter %s. No such parameter found in injection!\n",node->name);}
    }
  }
  
  
  /* Over-ride approximant if user specifies */
  
  ppt=LALInferenceGetProcParamVal(commandLine,"--approximant");
  if(ppt){
    approx = XLALGetApproximantFromString(ppt->value);
  }
  ppt=LALInferenceGetProcParamVal(commandLine,"--approx");
  if(ppt){
    approx = XLALGetApproximantFromString(ppt->value);
  }


  
  
  
  if(approx==NumApproximants && injTable){ /* Read aproximant from injection file */
    approx=XLALGetApproximantFromString(injTable->waveform);
  }
  if(approx==NumApproximants){
      
       approx=TaylorF2; /* Defaults to TF2 */
       XLALPrintWarning("You did not provide an approximant for the templates. Using default %s, which might now be what you want!\n",XLALGetStringFromApproximant(approx));
  }

  ppt=LALInferenceGetProcParamVal(commandLine, "--fref");
  if (ppt) fRef = atof(ppt->value);
  
  /* Set the modeldomain appropriately */
  ppt=LALInferenceGetProcParamVal(commandLine,"--modeldomain");
  if(ppt){
    if ( ! strcmp( "time", ppt->value ) )
    {
      modelDomain = LAL_SIM_DOMAIN_TIME;
    }
    else if ( ! strcmp( "frequency", ppt->value ) )
    {
      modelDomain = LAL_SIM_DOMAIN_FREQUENCY;
    }
    else
    {
      fprintf( stderr, "invalid argument to --modeldomain:\n"
      "unknown domain specified: "
      "domain must be one of: time, frequency\n");
      exit( 1 );
    }
  } else {
    /* default domain */
    modelDomain = LAL_SIM_DOMAIN_FREQUENCY;
    XLALPrintWarning("You did not provide a modeldomain for the templates. Using default frequency, which might now be what you want!\n");
  }
  
  dataPtr = state->data;
  while (dataPtr != NULL) {
    dataPtr->modelDomain = modelDomain;
    dataPtr = dataPtr->next;
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



  /************ Prior Related Settings START ******************************/
  
  /* Over-ride time prior if specified */
  ppt=LALInferenceGetProcParamVal(commandLine,"--dt");
  if(ppt){
    dt=atof(ppt->value);
  }

  /* Over-ride Distance min if specified */
  ppt=LALInferenceGetProcParamVal(commandLine,"--Dmin");
  if(ppt){
    Dmin=atof(ppt->value);
  }

  /* Over-ride Distance max if specified */
  ppt=LALInferenceGetProcParamVal(commandLine,"--Dmax");
  if(ppt){
    Dmax=atof(ppt->value);
  }
  
  ppt=LALInferenceGetProcParamVal(commandLine,"--mass-min");
  if(!ppt) ppt=LALInferenceGetProcParamVal(commandLine,"--massmin");
  if(ppt)
  {
    mMin=atof(ppt->value);
    if (mMin < 0.0)
    {
      fprintf(stderr,"ERROR: Minimum value of mass must be > 0");
      exit(1);
    }
  }

  ppt=LALInferenceGetProcParamVal(commandLine,"--mass-max");
  if(!ppt) ppt=LALInferenceGetProcParamVal(commandLine,"--massmax");
  if(ppt)
  {
    mMax=atof(ppt->value);
  } 

  ppt=LALInferenceGetProcParamVal(commandLine,"--eta-min");
  if(!ppt) ppt=LALInferenceGetProcParamVal(commandLine,"--etamin");
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
  if(!ppt) ppt=LALInferenceGetProcParamVal(commandLine,"--etamax");
  if(ppt)
    {
      etaMax=atof(ppt->value);
      if (etaMax > 0.25 || etaMax <= 0.0)
        {
          fprintf(stderr,"ERROR: Maximum value of eta must be between 0 and 0.25\n");
          exit(1);
        }
    }
  
  
  ppt=LALInferenceGetProcParamVal(commandLine,"--q-min");
  if(ppt)
    {
      qMin=atof(ppt->value);
      if (qMin <= 0.0 || qMin < 1.0/mMax || qMin < 1.0/(mMax-1.0) || qMin > 1.0)
        {
          fprintf(stderr,"ERROR: invalid qMin ( max{0,1.0/mMax,1.0/(mMax-1.0) < q < 1.0} )");
          exit(1);
        }
    }

  ppt=LALInferenceGetProcParamVal(commandLine,"--q-max");
  if(ppt)
    {
      qMax=atof(ppt->value);
      if (qMax > 1.0 || qMax <= 0.0 || qMax < 1.0/mMax || qMax < 1.0/(mMax-1.0))
        {
          fprintf(stderr,"ERROR: invalid qMax ( max{0,1.0/mMax,1.0/(mMax-1.0) < q < 1.0} )");
          exit(1);
        }
    }

  ppt=LALInferenceGetProcParamVal(commandLine,"--a-min");
  if (ppt) {
    amin = atof(ppt->value);
  }

  ppt=LALInferenceGetProcParamVal(commandLine,"--a-max");
  if (ppt) {
    amax = atof(ppt->value);
  }

  ppt=LALInferenceGetProcParamVal(commandLine,"--iota-max");
  if (ppt) {
    iotaMax = atof(ppt->value);
  }

  /* Prior related arguments END */

  

  /* Initial values set after setting up prior */
  /************ Initial Value Related Argument START *************/

  REAL8 start_mass    = mMin+gsl_rng_uniform(GSLrandom)*(mMax-mMin) ;
  REAL8 start_eta     =etaMin+gsl_rng_uniform(GSLrandom)*(etaMax-etaMin);
  REAL8 start_q       =qMin+gsl_rng_uniform(GSLrandom)*(qMax-qMin);
  REAL8 start_phase   =0.0+gsl_rng_uniform(GSLrandom)*(LAL_TWOPI-0.0);
  REAL8 start_dist    =Dmin+gsl_rng_uniform(GSLrandom)*(Dmax-Dmin);
  REAL8 start_ra      =0.0+gsl_rng_uniform(GSLrandom)*(LAL_TWOPI-0.0);
  REAL8 start_dec     =-LAL_PI/2.0+gsl_rng_uniform(GSLrandom)*(LAL_PI_2-(-LAL_PI_2));
  REAL8 start_psi     =0.0+gsl_rng_uniform(GSLrandom)*(LAL_PI-0.0);
  REAL8 start_iota    =0.0+gsl_rng_uniform(GSLrandom)*(LAL_PI-0.0);
  REAL8 start_a_spin  =0.0+gsl_rng_uniform(GSLrandom)*(amax-amin);
  
  
  /* Read time parameter from injection file */
  if(injTable)
  {
    starttime=XLALGPSGetREAL8(&(injTable->geocent_end_time));
    // Ringdown start_time is inspiral end_time
    fprintf(stdout,"Using start time from injection file: %lf\n", starttime);
  }
  /* Over-ride end time if specified */
  ppt=LALInferenceGetProcParamVal(commandLine,"--trigtime");
  if(ppt && !analytic){
    starttime=atof(ppt->value);
    printf("Read start time %f\n",starttime);
  }
  /* Adjust prior accordingly */
  if (!analytic) {
    timeMin=starttime-dt; timeMax=starttime+dt;
  }
  
  /* Over-ride mass if specified */
  ppt=LALInferenceGetProcParamVal(commandLine,"--mass");
  if(ppt){
    start_mass=atof(ppt->value);
  }
  
  /* Over-ride eta if specified */
  ppt=LALInferenceGetProcParamVal(commandLine,"--eta");
  if(ppt){
    start_eta=atof(ppt->value);
    REAL8 m1 = 0.0, m2 = 0.0 ;
    REAL8 mc = start_mass * pow(start_eta,3.0/5.0);
    LALInferenceMcEta2Masses(mc, start_eta, &m1, &m2);
    start_q=m2/m1;
  }
  
  /* Over-ride q if specified */
  ppt=LALInferenceGetProcParamVal(commandLine,"--q");
  if(ppt){
    start_q=atof(ppt->value);
    LALInferenceQ2Eta(start_q, &start_eta);
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
  
  ppt=LALInferenceGetProcParamVal(commandLine,"--a");
  if (ppt) {
    start_a_spin = atof(ppt->value);
  }

  /* Set up start time. */
  ppt=LALInferenceGetProcParamVal(commandLine, "--time");
  if (ppt) {
    /* User has specified start time. */
    timeParam = atof(ppt->value);
  } else {
    timeParam = starttime+dt*(1.0-2.0*gsl_rng_uniform(GSLrandom));
    timeParam = starttime+gsl_ran_gaussian(GSLrandom,0.01);
  }

  /* Initial Value Related END */
  

  LALInferenceAddVariable(currentParams, "LAL_APPROXIMANT", &approx,        LALINFERENCE_UINT4_t, LALINFERENCE_PARAM_FIXED);
  
  LALInferenceAddVariable(currentParams, "fRef", &fRef, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);

  ppt=LALInferenceGetProcParamVal(commandLine,"--taper");
  if(ppt){
    LALInferenceAddVariable(currentParams, "LALINFERENCE_TAPER",     &bookends,        LALINFERENCE_UINT4_t, LALINFERENCE_PARAM_FIXED);
  }  


  
  /* Set up the variable parameters */

  /********************* TBL: Adding noise-fitting parameters  *********************/
  UINT4 nscale_block; //number of noise parameters per IFO (1 per frequency block)
  UINT4 nscale_bin;   //number of Fourier bins in each noise block
  REAL8 nscale_dflog; //logarithmic spacing for noise parameters

  UINT4 nscale_dim;   //total dimension of noise model (params X detectors)
  UINT4 nscale_flag;  //flag to tell likelihood function if psd fitting is in use

  REAL8Vector *nscale_prior = NULL; //std. dev. of prior distribution
  REAL8Vector *nscale_sigma = NULL; //std. dev. of prior distribution

  //assume no noise fitting
  nscale_flag=0;

  //set Nblock to default unless specified at command line
  ppt = LALInferenceGetProcParamVal(commandLine, "--psdNblock");
  if(ppt) nscale_block = atoi(ppt->value);
  else nscale_block = 8;

  //First, figure out sizes of dataset to set up noise blocks
  UINT4 nifo; //number of data channels
  UINT4 imin; //minimum Fourier bin for integration in IFO
  UINT4 imax; //maximum Fourier bin for integration in IFO
  UINT4 f_min = 1; //minimum Fourier bin for integration over network
  UINT4 f_max = 1; //maximum Fourier bin for integration over network
  REAL8 df = 1.0; //frequency resolution

  //compute imin,imax for each IFO -- may be different
  nifo=0;
  dataPtr = state->data;
  while (dataPtr != NULL)
  {
    dt      = dataPtr->timeData->deltaT;
    df      = 1.0 / (((double)dataPtr->timeData->data->length) * dt);
    imin    = (UINT4)ceil( dataPtr->fLow  / df);
    imax    = (UINT4)floor(dataPtr->fHigh / df);

    if(nifo==0)
    {
      f_min=imin;
      f_max=imax;
    }
    else
    {
      if(imin<f_min)
      {
        fprintf(stderr,"Warning: Different IFO's have different minimum frequencies -- bad for noise fitting\n");
        f_min=imin;
      }
      if(imax>f_max)
      {
        fprintf(stderr,"Warning: Different IFO's have different minimum frequencies -- bad for noise fitting\n");
        f_max=imax;
      }
    }

    dataPtr = dataPtr->next;
    nifo++;
  }

  UINT4 j = 0;

  ppt = LALInferenceGetProcParamVal(commandLine, "--psdFit");
  if(ppt)//MARK: Here is where noise PSD parameters are being added to the model
  {
 
    printf("Setting up PSD fitting for %i ifos...\n",nifo);

    dataPtr = state->data;
    UINT4 nscale_block_temp   = 10000;
    gsl_matrix *bands_min_temp      = gsl_matrix_alloc(nifo,nscale_block_temp);
    gsl_matrix *bands_max_temp      = gsl_matrix_alloc(nifo,nscale_block_temp);

    i=0;
    while (dataPtr != NULL)
    {

      for (j = 0; j < nscale_block_temp; j++)
      {
        gsl_matrix_set(bands_min_temp,i,j,-1.0);
        gsl_matrix_set(bands_max_temp,i,j,-1.0);
      }

      printf("ifo=%i  %s\n",i,dataPtr->name);fflush(stdout);

      char ifoPSDFitBands[500];
      snprintf (ifoPSDFitBands,500, "--%s-psdFit",dataPtr->name);

      ppt = LALInferenceGetProcParamVal(commandLine,ifoPSDFitBands);
      if(ppt || LALInferenceGetProcParamVal(commandLine, "--xcorrbands"))
      /* Load in values from file if requested */
      {
        char line [ 128 ];
        char bands_tempfile[500];

        if (LALInferenceGetProcParamVal(commandLine, "--xcorrbands")) {
           snprintf (bands_tempfile,500, "%s-XCorrBands.dat",dataPtr->name);
        }
        else {
           char *bands_tempfile_temp = ppt->value;
           strcpy( bands_tempfile, bands_tempfile_temp );
        }
        printf("Reading bands_temp from %s\n",bands_tempfile);

        UINT4 band_min = 0, band_max = 0;
 
        nscale_block = 0;
        char * pch;
        j = 0;

        FILE *file = fopen ( bands_tempfile, "r" );
        if ( file != NULL )
        {
          while ( fgets ( line, sizeof line, file ) != NULL )
          {
              pch = strtok (line," ");
              int count = 0;
              while (pch != NULL)
              {
                  if (count==0) {band_min = atof(pch);}
                  if (count==1) {band_max = atof(pch);}
                  pch = strtok (NULL, " ");
                  count++;
              }

              gsl_matrix_set(bands_min_temp,i,j,band_min/df);
              gsl_matrix_set(bands_max_temp,i,j,band_max/df);
 
              nscale_block++;
              j++;

          }
        fclose ( file );
        }

        else
        {
          perror ( bands_tempfile ); /* why didn't the file open? */
        }
  

      }
      else // Otherwise use defaults
      {

        nscale_bin   = (f_max+1-f_min)/nscale_block;
        nscale_dflog = log( (double)(f_max+1)/(double)f_min )/(double)nscale_block;

        int freq_min, freq_max;

        for (j = 0; j < nscale_block; j++)
        {

            freq_min = (int) exp(log((double)f_min ) + nscale_dflog*j);
            freq_max = (int) exp(log((double)f_min ) + nscale_dflog*(j+1));

            gsl_matrix_set(bands_min_temp,i,j,freq_min);
            gsl_matrix_set(bands_max_temp,i,j,freq_max);
        }

      }  

      dataPtr = dataPtr->next;
      i++;

    }


    gsl_matrix *bands_min      = gsl_matrix_alloc(nifo,nscale_block);
    gsl_matrix *bands_max      = gsl_matrix_alloc(nifo,nscale_block);

    for (i = 0; i < nifo; i++)
    {
      for (j = 0; j < nscale_block; j++)
      {
        gsl_matrix_set(bands_min,i,j,gsl_matrix_get(bands_min_temp,i,j));
        gsl_matrix_set(bands_max,i,j,gsl_matrix_get(bands_max_temp,i,j));

      }
    }

    printf("Running PSD fitting with bands (Hz)...\n");
    dataPtr = state->data;
    i=0;
    while (dataPtr != NULL)
    {
      printf("%s:",dataPtr->name);
      for (j = 0; j < nscale_block; j++)
      {
        printf(" %f-%f ",gsl_matrix_get(bands_min,i,j)*df,gsl_matrix_get(bands_max,i,j)*df);
      }

      dataPtr = dataPtr->next;
      i++;
    }

    nscale_bin   = (f_max+1-f_min)/nscale_block;
    nscale_dflog = log( (double)(f_max+1)/(double)f_min )/(double)nscale_block;

    //nscale_min   = 0.00;
    //nscale_max   = 10.0;
    nscale_dim   = nscale_block*nifo;
    nscale_flag  = 1;

    // Set noise parameter arrays.
    nscale_prior = XLALCreateREAL8Vector(nscale_block);
    nscale_sigma = XLALCreateREAL8Vector(nscale_block);
    for(i=0; i<nscale_block; i++)
    {
      nscale_prior->data[i] = 1.0/sqrt( f_min*exp( (double)(i+1)*nscale_dflog ) );
      nscale_sigma->data[i] = nscale_prior->data[i]/sqrt((double)(nifo*nscale_block));

    }

    gsl_matrix *nscale = gsl_matrix_alloc(nifo,nscale_block);
    gsl_matrix *nstore = gsl_matrix_alloc(nifo,nscale_block);

    gsl_matrix_set_all(nscale, 1.0);
    gsl_matrix_set_all(nstore, 1.0);


    LALInferenceAddVariable(currentParams, "psdscale", &nscale, LALINFERENCE_gslMatrix_t, LALINFERENCE_PARAM_LINEAR);
    LALInferenceAddVariable(currentParams, "psdstore", &nstore, LALINFERENCE_gslMatrix_t, LALINFERENCE_PARAM_FIXED);
    LALInferenceAddVariable(currentParams, "logdeltaf", &nscale_dflog, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);

    LALInferenceAddVariable(currentParams, "psdBandsMin", &bands_min, LALINFERENCE_gslMatrix_t, LALINFERENCE_PARAM_FIXED);
    LALInferenceAddVariable(currentParams, "psdBandsMax", &bands_max, LALINFERENCE_gslMatrix_t, LALINFERENCE_PARAM_FIXED);

    //Set up noise priors
    LALInferenceAddVariable(priorArgs,      "psddim",   &nscale_dim,  LALINFERENCE_INT4_t,  LALINFERENCE_PARAM_FIXED);
    //LALInferenceAddMinMaxPrior(priorArgs,   "psdscale", &nscale_min,  &nscale_max,   LALINFERENCE_REAL8_t);
    LALInferenceAddVariable(priorArgs,      "psdsigma", &nscale_prior, LALINFERENCE_REAL8Vector_t, LALINFERENCE_PARAM_FIXED);

    //Store meta data for noise model in proposal
    LALInferenceAddVariable(state->proposalArgs, "psdblock", &nscale_block, LALINFERENCE_INT4_t,  LALINFERENCE_PARAM_FIXED);
    LALInferenceAddVariable(state->proposalArgs, "psdbin",   &nscale_bin,   LALINFERENCE_INT4_t,  LALINFERENCE_PARAM_FIXED);
    LALInferenceAddVariable(state->proposalArgs, "psdsigma", &nscale_sigma, LALINFERENCE_REAL8Vector_t, LALINFERENCE_PARAM_FIXED);


  }//End of noise model initialization
  LALInferenceAddVariable(currentParams, "psdScaleFlag", &nscale_flag, LALINFERENCE_UINT4_t, LALINFERENCE_PARAM_FIXED);

  UINT4 psdGaussianPrior=1;
  ppt = LALInferenceGetProcParamVal(commandLine, "--psdFlatPrior");
  if(ppt)psdGaussianPrior=0;
  LALInferenceAddVariable(priorArgs, "psdGaussianPrior", &psdGaussianPrior,  LALINFERENCE_INT4_t,  LALINFERENCE_PARAM_FIXED);

  /********************* TBL: Adding line-removal parameters  *********************/
  UINT4 lines_flag  = 0;   //flag tells likelihood if line-removal is turned on

  #define max(x, y) (((x) > (y)) ? (x) : (y))
  ppt = LALInferenceGetProcParamVal(commandLine, "--removeLines");

  if(ppt)//MARK: Here is where noise line removal parameters are being added to the model
  {
    lines_flag = 1;
    dataPtr = state->data;
    UINT4 lines_num_temp   = 10000;
    UINT4 lines_num_ifo;
    UINT4 lines_num = 0;
    gsl_matrix *lines_temp      = gsl_matrix_alloc(nifo,lines_num_temp);
    gsl_matrix *linewidth_temp  = gsl_matrix_alloc(nifo,lines_num_temp);

    i=0;
    while (dataPtr != NULL)
    {

      for (j = 0; j < lines_num_temp; j++)
      {
        gsl_matrix_set(lines_temp,i,j,-1.0);
        gsl_matrix_set(linewidth_temp,i,j,1.0);
      }

      printf("ifo=%i  %s\n",i,dataPtr->name);fflush(stdout);

      char ifoRemoveLines[500];
      snprintf (ifoRemoveLines,500, "--%s-removeLines",dataPtr->name);

      ppt = LALInferenceGetProcParamVal(commandLine,ifoRemoveLines);
      if(ppt || LALInferenceGetProcParamVal(commandLine, "--chisquaredlines") || LALInferenceGetProcParamVal(commandLine, "--KSlines") || LALInferenceGetProcParamVal(commandLine, "--powerlawlines"))
      /* Load in values from file if requested */
      {
        char line [ 128 ];
        char lines_tempfile[500];

        if (LALInferenceGetProcParamVal(commandLine, "--chisquaredlines")) {
             snprintf (lines_tempfile,500, "%s-ChiSquaredLines.dat",dataPtr->name);
        }
        else if (LALInferenceGetProcParamVal(commandLine, "--KSlines")) {
             snprintf (lines_tempfile,500, "%s-KSLines.dat",dataPtr->name);
        }
        else if (LALInferenceGetProcParamVal(commandLine, "--powerlawlines")) {
             snprintf (lines_tempfile,500, "%s-PowerLawLines.dat",dataPtr->name);
        }
        else {
             char *lines_tempfile_temp = ppt->value;
             strcpy( lines_tempfile, lines_tempfile_temp );
        }
        printf("Reading lines_temp from %s\n",lines_tempfile);

        char * pch;
        j = 0;
        double freqline = 0, freqlinewidth = 0;
        lines_num_ifo = 0;
        FILE *file = fopen ( lines_tempfile, "r" );
        if ( file != NULL )
        {
          while ( fgets ( line, sizeof line, file ) != NULL )
          {

            pch = strtok (line," ");
            int count = 0;
            while (pch != NULL)
            {
                if (count==0) {freqline = atof(pch);}
                if (count==1) {freqlinewidth = atof(pch);}
                pch = strtok (NULL, " ");
                count++;
            }

            gsl_matrix_set(lines_temp,i,j,freqline/df);
            gsl_matrix_set(linewidth_temp,i,j,freqlinewidth/df);
            j++;
            lines_num_ifo++;
          }
          fclose ( file );
        }

      }


      else // Otherwise use defaults
      {
        lines_num_ifo = 10;
        /* top lines_temp_num lines_temp for each interferometer, and widths */
        if(!strcmp(dataPtr->name,"H1"))
        {
          gsl_matrix_set(lines_temp,i,0,35.0/df);   gsl_matrix_set(linewidth_temp,i,0,3.0/df);
          gsl_matrix_set(lines_temp,i,1,45.0/df);   gsl_matrix_set(linewidth_temp,i,1,1.5/df);
          gsl_matrix_set(lines_temp,i,2,51.0/df);   gsl_matrix_set(linewidth_temp,i,2,2.5/df);
          gsl_matrix_set(lines_temp,i,3,60.0/df);   gsl_matrix_set(linewidth_temp,i,3,3.0/df);
          gsl_matrix_set(lines_temp,i,4,72.0/df);   gsl_matrix_set(linewidth_temp,i,4,3.0/df);
          gsl_matrix_set(lines_temp,i,5,87.0/df);   gsl_matrix_set(linewidth_temp,i,5,0.5/df);
          gsl_matrix_set(lines_temp,i,6,108.0/df);  gsl_matrix_set(linewidth_temp,i,6,0.5/df);
          gsl_matrix_set(lines_temp,i,7,117.0/df);  gsl_matrix_set(linewidth_temp,i,7,0.5/df);
          gsl_matrix_set(lines_temp,i,8,122.0/df);  gsl_matrix_set(linewidth_temp,i,8,5.0/df);
          gsl_matrix_set(lines_temp,i,9,180.0/df);  gsl_matrix_set(linewidth_temp,i,9,2.0/df);
        }

        if(!strcmp(dataPtr->name,"L1"))
        {
          gsl_matrix_set(lines_temp,i,0,35.0/df);   gsl_matrix_set(linewidth_temp,i,0,3.0/df);
          gsl_matrix_set(lines_temp,i,1,60.0/df);   gsl_matrix_set(linewidth_temp,i,1,4.0/df);
          gsl_matrix_set(lines_temp,i,2,69.0/df);   gsl_matrix_set(linewidth_temp,i,2,2.5/df);
          gsl_matrix_set(lines_temp,i,3,106.4/df);  gsl_matrix_set(linewidth_temp,i,3,0.8/df);
          gsl_matrix_set(lines_temp,i,4,113.0/df);  gsl_matrix_set(linewidth_temp,i,4,1.5/df);
          gsl_matrix_set(lines_temp,i,5,120.0/df);  gsl_matrix_set(linewidth_temp,i,5,2.5/df);
          gsl_matrix_set(lines_temp,i,6,128.0/df);  gsl_matrix_set(linewidth_temp,i,6,3.5/df);
          gsl_matrix_set(lines_temp,i,7,143.0/df);  gsl_matrix_set(linewidth_temp,i,7,1.0/df);
          gsl_matrix_set(lines_temp,i,8,180.0/df);  gsl_matrix_set(linewidth_temp,i,8,2.5/df);
          gsl_matrix_set(lines_temp,i,9,191.5/df);  gsl_matrix_set(linewidth_temp,i,9,4.0/df);
        }

        if(!strcmp(dataPtr->name,"V1"))
        {
          gsl_matrix_set(lines_temp,i,0,35.0/df);   gsl_matrix_set(linewidth_temp,i,0,3.0/df);
          gsl_matrix_set(lines_temp,i,1,60.0/df);   gsl_matrix_set(linewidth_temp,i,1,4.0/df);
          gsl_matrix_set(lines_temp,i,2,69.0/df);   gsl_matrix_set(linewidth_temp,i,2,2.5/df);
          gsl_matrix_set(lines_temp,i,3,106.4/df);  gsl_matrix_set(linewidth_temp,i,3,0.8/df);
          gsl_matrix_set(lines_temp,i,4,113.0/df);  gsl_matrix_set(linewidth_temp,i,4,1.5/df);
          gsl_matrix_set(lines_temp,i,5,120.0/df);  gsl_matrix_set(linewidth_temp,i,5,2.5/df);
          gsl_matrix_set(lines_temp,i,6,128.0/df);  gsl_matrix_set(linewidth_temp,i,6,3.5/df);
          gsl_matrix_set(lines_temp,i,7,143.0/df);  gsl_matrix_set(linewidth_temp,i,7,1.0/df);
          gsl_matrix_set(lines_temp,i,8,180.0/df);  gsl_matrix_set(linewidth_temp,i,8,2.5/df);
          gsl_matrix_set(lines_temp,i,9,191.5/df);  gsl_matrix_set(linewidth_temp,i,9,4.0/df);
        }
      }
      dataPtr = dataPtr->next;
      i++;

      lines_num = max(lines_num,lines_num_ifo);

    }

    gsl_matrix *lines      = gsl_matrix_alloc(nifo,lines_num);
    gsl_matrix *linewidth  = gsl_matrix_alloc(nifo,lines_num);

    for (i = 0; i < nifo; i++)
    {
      for (j = 0; j < lines_num; j++)
      {
        gsl_matrix_set(lines,i,j,gsl_matrix_get(lines_temp,i,j));
        gsl_matrix_set(linewidth,i,j,gsl_matrix_get(linewidth_temp,i,j));
      }
    }

    /* Add line matrices to variable lists */
    LALInferenceAddVariable(currentParams, "line_center", &lines,     LALINFERENCE_gslMatrix_t, LALINFERENCE_PARAM_FIXED);
    LALInferenceAddVariable(currentParams, "line_width",  &linewidth, LALINFERENCE_gslMatrix_t, LALINFERENCE_PARAM_FIXED);


  }//End of line-removal initialization
   
  LALInferenceAddVariable(currentParams, "removeLinesFlag", &lines_flag, LALINFERENCE_UINT4_t, LALINFERENCE_PARAM_FIXED);
  /*********************************************************************************/



  /* Check if running with symmetric (eta) or asymmetric (q) mass ratio.*/
  ppt=LALInferenceGetProcParamVal(commandLine,"--fixQ");
  if(ppt){
    LALInferenceRegisterUniformVariableREAL8(state, currentParams, "asym_massratio", start_q, qMin, qMax, LALINFERENCE_PARAM_FIXED);
    if(lalDebugLevel>0) fprintf(stdout,"q fixed and set to %f\n",start_q);
  }else{
    ppt=LALInferenceGetProcParamVal(commandLine,"--fixEta");
    if(ppt){
      LALInferenceRegisterUniformVariableREAL8(state, currentParams, "massratio", start_eta, etaMin, etaMax, LALINFERENCE_PARAM_FIXED);
      if(lalDebugLevel>0) fprintf(stdout,"eta fixed and set to %f\n",start_eta);
    }else{
      ppt=LALInferenceGetProcParamVal(commandLine,"--symMassRatio");
      if(ppt){
        LALInferenceRegisterUniformVariableREAL8(state, currentParams, "massratio", start_eta, etaMin, etaMax, LALINFERENCE_PARAM_LINEAR);
      }else{
        LALInferenceRegisterUniformVariableREAL8(state, currentParams, "asym_massratio", start_q, qMin, qMax, LALINFERENCE_PARAM_LINEAR);
      }
    }
  }

  ppt = LALInferenceGetProcParamVal(commandLine, "--fixMass");
  if (ppt) {
	  LALInferenceRegisterUniformVariableREAL8(state, currentParams, "rdMass", start_mass, mMin, mMax, LALINFERENCE_PARAM_FIXED);
	  if (lalDebugLevel>0) fprintf(stdout, "final black hole mass fixed and set to %f\n", start_mass);
  } 
  else {
	  LALInferenceRegisterUniformVariableREAL8(state, currentParams, "rdMass", start_mass, mMin, mMax, LALINFERENCE_PARAM_LINEAR);
  }


  ppt=LALInferenceGetProcParamVal(commandLine,"--fixTime");
  if(ppt){
    LALInferenceRegisterUniformVariableREAL8(state, currentParams, "time", timeParam, timeMin, timeMax, LALINFERENCE_PARAM_FIXED);
    if(lalDebugLevel>0) fprintf(stdout,"time fixed and set to %f\n",timeParam);
  }else{
    LALInferenceRegisterUniformVariableREAL8(state, currentParams, "time", timeParam, timeMin, timeMax, LALINFERENCE_PARAM_LINEAR);
  }

  // if(!LALInferenceGetProcParamVal(commandLine,"--margphi")){
    ppt=LALInferenceGetProcParamVal(commandLine,"--fixPhi");
    if(ppt){
      LALInferenceRegisterUniformVariableREAL8(state, currentParams, "phase", start_phase, phiMin, phiMax, LALINFERENCE_PARAM_FIXED);
      if(lalDebugLevel>0) fprintf(stdout,"phase fixed and set to %f\n",start_phase);
    }else{
      LALInferenceRegisterUniformVariableREAL8(state, currentParams, "phase", start_phase, phiMin, phiMax, LALINFERENCE_PARAM_CIRCULAR);
    }
  // }

  /* Jump in log distance if requested, otherwise use distance */
  if(LALInferenceGetProcParamVal(commandLine,"--use-logdistance")){
    LALInferenceRegisterUniformVariableREAL8(state, currentParams, "distance", log(start_dist), log(Dmin), log(Dmax), LALInferenceGetProcParamVal(commandLine,"--fixDist")?LALINFERENCE_PARAM_FIXED:LALINFERENCE_PARAM_LINEAR);
  } else {
    LALInferenceRegisterUniformVariableREAL8(state, currentParams, "distance", start_dist, Dmin, Dmax, LALInferenceGetProcParamVal(commandLine,"--fixDist")?LALINFERENCE_PARAM_FIXED:LALINFERENCE_PARAM_LINEAR);
  }

  LALInferenceRegisterUniformVariableREAL8(state, currentParams, "rdSpin", start_a_spin, amin, amax, LALInferenceGetProcParamVal(commandLine,"--fixSpin")?LALINFERENCE_PARAM_FIXED:LALINFERENCE_PARAM_LINEAR);
  
  LALInferenceRegisterUniformVariableREAL8(state, currentParams, "rightascension", start_ra, raMin, raMax, LALInferenceGetProcParamVal(commandLine,"--fixRa")?LALINFERENCE_PARAM_FIXED:LALINFERENCE_PARAM_CIRCULAR);
  
  LALInferenceRegisterUniformVariableREAL8(state, currentParams, "declination", start_dec, decMin, decMax, LALInferenceGetProcParamVal(commandLine,"--fixDec")?LALINFERENCE_PARAM_FIXED:LALINFERENCE_PARAM_LINEAR);
  
  LALInferenceRegisterUniformVariableREAL8(state, currentParams, "polarisation", start_psi, psiMin, psiMax, LALInferenceGetProcParamVal(commandLine,"--fixPsi")?LALINFERENCE_PARAM_FIXED:LALINFERENCE_PARAM_LINEAR);
  

  LALInferenceRegisterUniformVariableREAL8(state, currentParams, "inclination", start_iota, iotaMin, iotaMax, LALInferenceGetProcParamVal(commandLine,"--fixIota")?LALINFERENCE_PARAM_FIXED:LALINFERENCE_PARAM_LINEAR);

     /* Print info about orders and waveflags used for templates */
     fprintf(stdout,"\n\n---\t\t ---\n");
     fprintf(stdout,"Templates will run using ringdown Approximant %i (%s) in the %s domain.\n",approx,XLALGetStringFromApproximant(approx),modelDomain==LAL_SIM_DOMAIN_TIME?"time":"frequency");
     fprintf(stdout,"---\t\t ---\n\n");
  return(currentParams);
}
