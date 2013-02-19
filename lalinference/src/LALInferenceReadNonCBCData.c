/* 
 *  LALInferenceReadData.c:  Bayesian Followup functions
 *
 *  Copyright (C) 2009,2012 Ilya Mandel, Vivien Raymond, Christian
 *  Roever, Marc van der Sluys, John Veitch, Salvatore Vitale, and
 *  Will M. Farr
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
#include <stdlib.h>

#define LAL_USE_OLD_COMPLEX_STRUCTS
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALInspiral.h>
#include <lal/FrameCache.h>
#include <lal/FrameStream.h>
#include <lal/TimeFreqFFT.h>
#include <lal/LALDetectors.h>
#include <lal/AVFactories.h>
#include <lal/ResampleTimeSeries.h>
#include <lal/TimeSeries.h>
#include <lal/FrequencySeries.h>
#include <lal/Units.h>
#include <lal/Date.h>
#include <lal/StringInput.h>
#include <lal/VectorOps.h>
#include <lal/Random.h>
#include <lal/LALNoiseModels.h>
#include <lal/XLALError.h>
#include <lal/GenerateInspiral.h>
#include <lal/LIGOLwXMLRead.h>
#include <lal/LIGOLwXMLInspiralRead.h>
#include <lal/SeqFactories.h>
#include <lal/DetectorSite.h>
#include <lal/GenerateInspiral.h>
#include <lal/GeneratePPNInspiral.h>
#include <lal/SimulateCoherentGW.h>
#include <lal/Inject.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/LIGOMetadataUtils.h>
#include <lal/LIGOMetadataInspiralUtils.h>
#include <lal/LIGOMetadataRingdownUtils.h>
#include <lal/LALInspiralBank.h>
#include <lal/FindChirp.h>
#include <lal/LALInspiralBank.h>
#include <lal/GenerateInspiral.h>
#include <lal/NRWaveInject.h>
#include <lal/GenerateInspRing.h>
#include <lal/LALErrno.h>
#include <math.h>
#include <lal/LALInspiral.h>
#include <lal/LALSimulation.h>
#include <lal/LALInference.h>
#include <lal/LALInferenceLikelihood.h>
#include <lal/LALInferenceTemplate.h>
#include <lal/LIGOLwXMLBurstRead.h>
#include <lal/GenerateBurst.h>
#include <lal/LALSimBurst.h>
#include <lal/LALInferenceReadNonCBCData.h>
#define LALINFERENCE_DEFAULT_FLOW "40.0"
//typedef void (NoiseFunc)(LALStatus *statusPtr,REAL8 *psd,REAL8 f);

static const LALUnit strainPerCount={0,{0,0,0,0,0,1,-1},{0,0,0,0,0,0,0}};
 struct fvec {
	REAL8 f;
	REAL8 x;
};

typedef void (NoiseFunc)(LALStatus *statusPtr,REAL8 *psd,REAL8 f);
static struct fvec *interpFromFile(char *filename){
	UINT4 fileLength=0;
	UINT4 i=0;
	UINT4 minLength=100; /* size of initial file buffer, and also size of increment */
	FILE *interpfile=NULL;
	struct fvec *interp=NULL;
	interp=calloc(minLength,sizeof(struct fvec)); /* Initialise array */
	if(!interp) {printf("Unable to allocate memory buffer for reading interpolation file\n");}
	fileLength=minLength;
	REAL8 f=0.0,x=0.0;
	interpfile = fopen(filename,"r");
	if (interpfile==NULL){
		printf("Unable to open file %s\n",filename);
		exit(1);
	}
	while(2==fscanf(interpfile," %lf %lf ", &f, &x )){
		interp[i].f=f; interp[i].x=x*x;
		i++;
		if(i>fileLength-1){ /* Grow the array */
			interp=realloc(interp,(fileLength+minLength)*sizeof(struct fvec));
			fileLength+=minLength;
		}
	}
	interp[i].f=0; interp[i].x=0;
	fileLength=i+1;
	interp=realloc(interp,fileLength*sizeof(struct fvec)); /* Resize array */
	fclose(interpfile);
	printf("Read %i records from %s\n",fileLength-1,filename);
	return interp;
}

static REAL8 interpolate(struct fvec *fvec, REAL8 f);
static REAL8 interpolate(struct fvec *fvec, REAL8 f){
	int i=0;
	REAL8 a=0.0; /* fractional distance between bins */
	REAL8 delta=0.0;
	if(f<fvec[0].f) return(0.0);
	while(fvec[i].f<f && (fvec[i].x!=0.0 )){i++;}; //&& fvec[i].f!=0.0)){i++;};
  //printf("%d\t%lg\t%lg\t%lg\n",i,fvec[i].f,f,fvec[i].x);
	if (fvec[i].f==0.0 && fvec[i].x==0.0) /* Frequency above moximum */
	{
		return (fvec[i-1].x);
	}
  //if(i==0){return (fvec[0].x);}
	a=(fvec[i].f-f)/(fvec[i].f-fvec[i-1].f);
	delta=fvec[i].x-fvec[i-1].x;
	return (fvec[i-1].x + delta*a);
}

static void MetaNoiseFunc(LALStatus *status, REAL8 *psd, REAL8 f, struct fvec *interp, NoiseFunc *noisefunc){
	if(interp==NULL&&noisefunc==NULL){
		printf("ERROR: Trying to calculate PSD with NULL inputs\n");
		exit(1);
	}
	if(interp!=NULL && noisefunc!=NULL){
		printf("ERROR: You have specified both an interpolation vector and a function to calculate the PSD\n");
		exit(1);
	}
	if(noisefunc!=NULL){
		noisefunc(status,psd,f);
		return;
	}
	if(interp!=NULL){ /* Use linear interpolation of the interp vector */
		*psd=interpolate(interp,f);
		return;
	}
}

static REAL8TimeSeries *readTseries (CHAR *cachefile, CHAR *channel, LIGOTimeGPS start, REAL8 length);
static void makeWhiteData(LALInferenceIFOData *IFOdata);
static void makeWhiteData(LALInferenceIFOData *IFOdata) {
  REAL8 deltaF = IFOdata->freqData->deltaF;
  REAL8 deltaT = IFOdata->timeData->deltaT;

  IFOdata->whiteFreqData = 
    XLALCreateCOMPLEX16FrequencySeries("whitened frequency data", 
                                       &(IFOdata->freqData->epoch),
                                       0.0,
                                       deltaF,
                                       &lalDimensionlessUnit,
                                       IFOdata->freqData->data->length);
	if(!IFOdata->whiteFreqData) XLAL_ERROR_VOID(XLAL_EFUNC);
  IFOdata->whiteTimeData = 
    XLALCreateREAL8TimeSeries("whitened time data",
                              &(IFOdata->timeData->epoch),
                              0.0,
                              deltaT,
                              &lalDimensionlessUnit,
                              IFOdata->timeData->data->length);
	if(!IFOdata->whiteTimeData) XLAL_ERROR_VOID(XLAL_EFUNC);

  REAL8 iLow = IFOdata->fLow / deltaF;
  REAL8 iHighDefaultCut = 0.95 * IFOdata->freqData->data->length;
  REAL8 iHighFromFHigh = IFOdata->fHigh / deltaF;
  REAL8 iHigh = (iHighDefaultCut < iHighFromFHigh ? iHighDefaultCut : iHighFromFHigh);
  REAL8 windowSquareSum = 0.0;

  UINT4 i;

  for (i = 0; i < IFOdata->freqData->data->length; i++) {
    IFOdata->whiteFreqData->data->data[i].re = IFOdata->freqData->data->data[i].re / IFOdata->oneSidedNoisePowerSpectrum->data->data[i];
    IFOdata->whiteFreqData->data->data[i].im = IFOdata->freqData->data->data[i].im / IFOdata->oneSidedNoisePowerSpectrum->data->data[i];
		
    if (i == 0) {
      /* Cut off the average trend in the data. */
      IFOdata->whiteFreqData->data->data[i].re = 0.0;
      IFOdata->whiteFreqData->data->data[i].im = 0.0;
    }
    if (i <= iLow) {
      /* Need to taper to implement the fLow cutoff.  Tukey window
			 that starts at zero, and reaches 100% at fLow. */
      REAL8 weight = 0.5*(1.0 + cos(M_PI*(i-iLow)/iLow)); /* Starts at -Pi, runs to zero at iLow. */
			
      IFOdata->whiteFreqData->data->data[i].re *= weight;
      IFOdata->whiteFreqData->data->data[i].im *= weight;
			
      windowSquareSum += weight*weight;
    } else if (i >= iHigh) {
      /* Also taper at high freq end, Tukey window that starts at 100%
			 at fHigh, then drops to zero at Nyquist.  Except that we
			 always taper at least 5% of the data at high freq to avoid a
			 sharp edge in freq space there. */
      REAL8 NWind = IFOdata->whiteFreqData->data->length - iHigh;
      REAL8 weight = 0.5*(1.0 + cos(M_PI*(i-iHigh)/NWind)); /* Starts at 0, runs to Pi at i = length */
			
      IFOdata->whiteFreqData->data->data[i].re *= weight;
      IFOdata->whiteFreqData->data->data[i].im *= weight;
			
      windowSquareSum += weight*weight;
    } else {
      windowSquareSum += 1.0;
    }
  }
	
  REAL8 norm = sqrt(IFOdata->whiteFreqData->data->length / windowSquareSum);
  for (i = 0; i < IFOdata->whiteFreqData->data->length; i++) {
    IFOdata->whiteFreqData->data->data[i].re *= norm;
    IFOdata->whiteFreqData->data->data[i].im *= norm;
  }
	
  XLALREAL8FreqTimeFFT(IFOdata->whiteTimeData, IFOdata->whiteFreqData, IFOdata->freqToTimeFFTPlan);
}



static REAL8TimeSeries *readTseries(CHAR *cachefile, CHAR *channel, LIGOTimeGPS start, REAL8 length)
{
	LALStatus status;
	memset(&status,0,sizeof(status));
	FrCache *cache = NULL;
	FrStream *stream = NULL;
	REAL8TimeSeries *out = NULL;
	
	cache  = XLALFrImportCache( cachefile );
        int err;
        err = *XLALGetErrnoPtr();
	if(cache==NULL) {fprintf(stderr,"ERROR: Unable to import cache file \"%s\",\n       XLALError: \"%s\".\n",cachefile, XLALErrorString(err)); exit(-1);}
	stream = XLALFrCacheOpen( cache );
	if(stream==NULL) {fprintf(stderr,"ERROR: Unable to open stream from frame cache file\n"); exit(-1);}
	out = XLALFrInputREAL8TimeSeries( stream, channel, &start, length , 0 );
	if(out==NULL) fprintf(stderr,"ERROR: unable to read channel %s from %s at time %i\nCheck the specified data duration is not too long\n",channel,cachefile,start.gpsSeconds);
	LALDestroyFrCache(&status,&cache);
	LALFrClose(&status,&stream);
	return out;
}

/** Parse the command line looking for options of the kind --ifo H1 --H1-channel H1:LDAS_STRAIN --H1-cache H1.cache --H1-flow 40.0 --H1-fhigh 4096.0 --H1-timeslide 100.0 ...
    It is necessary to use this method instead of the old method for the pipeline to work in DAX mode. Warning: do not mix options between
    the old and new style.
*/
static INT4 getDataOptionsByDetectors(ProcessParamsTable *commandLine, char ***ifos, char ***caches, char ***channels, char ***flows , char ***fhighs, char ***timeslides, UINT4 *N)
{
    /* Check that the input has no lists with [ifo,ifo] */
    ProcessParamsTable *this=commandLine;
    UINT4 i=0;
    *caches=*ifos=*channels=*flows=*fhighs=*timeslides=NULL;
    *N=0;
    char tmp[128];
    if(!this) {fprintf(stderr,"No command line arguments given!\n"); exit(1);}
    while(this)
    {
        if(!strcmp(this->param,"--ifo") || !strcmp(this->param,"--IFO"))
        for(i=0;this->value[i]!='\0';i++)
            if(this->value[i]=='[' || this->value[i]==']')
            {
                fprintf(stderr,"Found old-style input arguments for %s\n",this->param);
                return(0);
            }
        this=this->next;
    }
    /* Construct a list of IFOs */
    for(this=commandLine;this;this=this->next)
    {
        if(!strcmp(this->param,"--ifo")||!strcmp(this->param,"--IFO"))
        {
            (*N)++;
            *ifos=realloc(*ifos,*N*sizeof(char *));
            (*ifos)[*N-1]=strdup(this->value);
        }
    }
    *caches=calloc(*N,sizeof(char *));
    *channels=calloc(*N,sizeof(char *));
    *flows=calloc(*N,sizeof(REAL8));
    *fhighs=calloc(*N,sizeof(REAL8));
    *timeslides=calloc(*N,sizeof(REAL8));
    /* For each IFO, fetch the other options if available */
    for(i=0;i<*N;i++)
    {
        /* Cache */
        sprintf(tmp,"--%s-cache",(*ifos)[i]);
        this=LALInferenceGetProcParamVal(commandLine,tmp);
        if(!this){fprintf(stderr,"ERROR: Must specify a cache file for %s with --%s-cache\n",(*ifos)[i],(*ifos)[i]); exit(1);}
        (*caches)[i]=strdup(this->value);
        
        /* Channel */
        sprintf(tmp,"--%s-channel",(*ifos)[i]);
        this=LALInferenceGetProcParamVal(commandLine,tmp);
        (*channels)[i]=strdup(this?this->value:"Unknown channel");

        /* flow */
        sprintf(tmp,"--%s-flow",(*ifos)[i]);
        this=LALInferenceGetProcParamVal(commandLine,tmp);
        (*flows)[i]=strdup(this?this->value:LALINFERENCE_DEFAULT_FLOW);
        
        /* fhigh */
        sprintf(tmp,"--%s-fhigh",(*ifos)[i]);
        this=LALInferenceGetProcParamVal(commandLine,tmp);
        (*fhighs)[i]=this?strdup(this->value):NULL;

        /* timeslides */
        sprintf(tmp,"--%s-timeslide",(*ifos)[i]);
        this=LALInferenceGetProcParamVal(commandLine,tmp);
        (*timeslides)[i]=strdup(this?this->value:"0.0");
    }
    return(1);
}




#define USAGE "\
 --ifo IFO1 [--ifo IFO2 ...]    IFOs can be H1,L1,V1\n\
 --IFO1-cache cache1 [--IFO2-cache2 cache2 ...]    LAL cache files (LALLIGO, LALAdLIGO, LALVirgo to simulate these detectors)\n\
 --psdstart GPStime             GPS start time of PSD estimation data\n\
 --psdlength length             length of PSD estimation data in seconds\n\
 --seglen length                length of segments for PSD estimation and analysis in seconds\n\
 --trigtime GPStime             GPS time of the trigger to analyse\n\
(--srate rate)                  Downsample data to rate in Hz (4096.0,)\n\
(--injectionsrate rate)         Downsample injection signal to rate in Hz (--srate)\n\
(--IFO1-flow freq1 [--IFO2-flow freq2 ...])      Specify lower frequency cutoff for overlap integral (40.0)\n\
(--IFO1-fhigh freq1 [--IFO2-fhigh freq2 ...])     Specify higher frequency cutoff for overlap integral (2048.0)\n\
(--IFO1-channel chan1 [--IFO2-channel chan2 ...])   Specify channel names when reading cache files\n\
(--dataseed number)             Specify random seed to use when generating data\n\
(--lalinspiralinjection)      Enables injections via the LALInspiral package\n\
(--inj-lambda1)                 value of lambda1 to be injected, LALSimulation only (0)\n\
(--inj-lambda2)                 value of lambda1 to be injected, LALSimulation only (0)\n\
(--inj-spinOrder PNorder)           Specify twice the PN order (e.g. 5 <==> 2.5PN) of spin effects to use, only for LALSimulation (default: -1 <==> Use all spin effects).\n\
(--inj-tidalOrder PNorder)          Specify twice the PN order (e.g. 10 <==> 5PN) of tidal effects to use, only for LALSimulation (default: -1 <==> Use all tidal effects).\n\
(--snrpath) 			Set a folder where to write a file with the SNRs being injected\n\
(--0noise)                      Sets the noise realisation to be identically zero (for the fake caches above only)\n"

LALInferenceIFOData *LALInferenceReadBurstData(ProcessParamsTable *commandLine)
/* Read in the data and store it in a LALInferenceIFOData structure */
{

    LALStatus status;
    INT4 dataseed=0;
    memset(&status,0,sizeof(status));
    ProcessParamsTable *procparam=NULL,*ppt=NULL;
    LALInferenceIFOData *headIFO=NULL,*IFOdata=NULL;
    REAL8 SampleRate=4096.0,SegmentLength=0;
    if(LALInferenceGetProcParamVal(commandLine,"--srate")) SampleRate=atof(LALInferenceGetProcParamVal(commandLine,"--srate")->value);
    REAL8 defaultFLow = atof(LALINFERENCE_DEFAULT_FLOW);
    int nSegs=0;
    size_t seglen=0;
    REAL8TimeSeries *PSDtimeSeries=NULL;
    REAL8 padding=0.4;//Default was 1.0 second. However for The Event the Common Inputs specify a Tukey parameter of 0.1, so 0.4 second of padding for 8 seconds of data.
    UINT4 Ncache=0,Npsd=0,Nifo=0,Nchannel=0,NfLow=0,NfHigh=0;
    UINT4 i,j;
    //int FakeFlag=0; - set but not used
    char strainname[]="LSC-STRAIN";
    UINT4 q=0;	
    //typedef void (NoiseFunc)(LALStatus *statusPtr,REAL8 *psd,REAL8 f);
    NoiseFunc *PSD=NULL;
    REAL8 scalefactor=1;
    SimBurst *BInjTable=NULL;
    RandomParams *datarandparam;
    UINT4 event=0;
    char *chartmp=NULL;
    char **channels=NULL;
    char **caches=NULL;
    char **psds=NULL;
    char **IFOnames=NULL;
    char **fLows=NULL,**fHighs=NULL;
    char **timeslides=NULL;
    UINT4 Ntimeslides=0;
    LIGOTimeGPS GPSstart,GPStrig,segStart;
    REAL8 PSDdatalength=0;
    REAL8 AIGOang=0.0; //orientation angle for the proposed Australian detector.
    procparam=LALInferenceGetProcParamVal(commandLine,"--aigoang");
    if(!procparam) procparam=LALInferenceGetProcParamVal(commandLine,"--AIGOang");
    if(procparam)
        AIGOang=atof(procparam->value)*LAL_PI/180.0;

    struct fvec *interp;
    int interpFlag=0;

    /* Check if the new style command line arguments are used */
    INT4 dataOpts=getDataOptionsByDetectors(commandLine, &IFOnames, &caches, &channels, &fLows , &fHighs, &timeslides, &Nifo);
    /* Check for options if not given in the new style */
    if(!dataOpts){
        if(!LALInferenceGetProcParamVal(commandLine,"--cache")||!(LALInferenceGetProcParamVal(commandLine,"--IFO")||LALInferenceGetProcParamVal(commandLine,"--ifo")))
            {fprintf(stderr,USAGE); return(NULL);}
        if(LALInferenceGetProcParamVal(commandLine,"--channel")){
            LALInferenceParseCharacterOptionString(LALInferenceGetProcParamVal(commandLine,"--channel")->value,&channels,&Nchannel);
        }
        LALInferenceParseCharacterOptionString(LALInferenceGetProcParamVal(commandLine,"--cache")->value,&caches,&Ncache);
        ppt=LALInferenceGetProcParamVal(commandLine,"--ifo");
        if(!ppt) ppt=LALInferenceGetProcParamVal(commandLine,"--IFO");
        LALInferenceParseCharacterOptionString(ppt->value,&IFOnames,&Nifo);

        ppt=LALInferenceGetProcParamVal(commandLine,"--flow");
        if(!ppt) ppt=LALInferenceGetProcParamVal(commandLine,"--fLow");
        if(ppt){
            LALInferenceParseCharacterOptionString(ppt->value,&fLows,&NfLow);
        }
        ppt=LALInferenceGetProcParamVal(commandLine,"--fhigh");
        if(!ppt) ppt=LALInferenceGetProcParamVal(commandLine,"--fHigh");
        if(ppt){
            LALInferenceParseCharacterOptionString(ppt->value,&fHighs,&NfHigh);
        }
        if((ppt=LALInferenceGetProcParamVal(commandLine,"--timeslide"))) LALInferenceParseCharacterOptionString(ppt->value,&timeslides,&Ntimeslides);
        if(Nifo!=Ncache) {fprintf(stderr,"ERROR: Must specify equal number of IFOs and Cache files\n"); exit(1);}
        if(Nchannel!=0 && Nchannel!=Nifo) {fprintf(stderr,"ERROR: Please specify a channel for all caches, or omit to use the defaults\n"); exit(1);}
    }
    else
    {
        NfHigh=Ntimeslides=Ncache=Nchannel=NfLow=Nifo;

    }
    /* Check for remaining required options */
	if(!(LALInferenceGetProcParamVal(commandLine,"--PSDstart")||LALInferenceGetProcParamVal(commandLine,"--psdstart")) ||
            !(LALInferenceGetProcParamVal(commandLine,"--PSDlength")||LALInferenceGetProcParamVal(commandLine,"--psdlength")) ||!LALInferenceGetProcParamVal(commandLine,"--seglen"))
    {fprintf(stderr,USAGE); return(NULL);}

    /* ET detectors */
    LALDetector dE1,dE2,dE3;
    /* response of the detectors */
    dE1.type = dE2.type = dE3.type = LALDETECTORTYPE_IFODIFF;
    dE1.location[0] = dE2.location[0] = dE3.location[0] = 4.5464e6;
    dE1.location[1] = dE2.location[1] = dE3.location[1] = 8.4299e5;
    dE1.location[2] = dE2.location[2] = dE3.location[2] = 4.3786e6;
    sprintf(dE1.frDetector.name,"ET-1");
    sprintf(dE1.frDetector.prefix,"E1");
    dE1.response[0][0] = 0.1666;
    dE1.response[1][1] = -0.2484;
    dE1.response[2][2] = 0.0818;
    dE1.response[0][1] = dE1.response[1][0] = -0.2188;
    dE1.response[0][2] = dE1.response[2][0] = -0.1300;
    dE1.response[1][2] = dE1.response[2][1] = 0.2732;
    sprintf(dE2.frDetector.name,"ET-2");
    sprintf(dE2.frDetector.prefix,"E2");
    dE2.response[0][0] = -0.1992;
    dE2.response[1][1] = 0.4234;
    dE2.response[2][2] = 0.0818;
    dE2.response[0][1] = dE2.response[1][0] = -0.0702;
    dE2.response[0][2] = dE2.response[2][0] = 0.2189;
    dE2.response[1][2] = dE2.response[2][1] = -0.0085;
    sprintf(dE3.frDetector.name,"ET-3");
    sprintf(dE3.frDetector.prefix,"E3");
    dE3.response[0][0] = 0.0326;
    dE3.response[1][1] = -0.1750;
    dE3.response[2][2] = 0.1423;
    dE3.response[0][1] = dE3.response[1][0] = 0.2891;
    dE3.response[0][2] = dE3.response[2][0] = -0.0889;
    dE3.response[1][2] = dE3.response[2][1] = -0.2647;  

    //TEMPORARY. JUST FOR CHECKING USING SPINSPIRAL PSD
    char **spinspiralPSD=NULL;
    UINT4 NspinspiralPSD = 0;
    if (LALInferenceGetProcParamVal(commandLine, "--spinspiralPSD")) {
        LALInferenceParseCharacterOptionString(LALInferenceGetProcParamVal(commandLine,"--spinspiralPSD")->value,&spinspiralPSD,&NspinspiralPSD);
    }    

    if(LALInferenceGetProcParamVal(commandLine,"--psd")){
        LALInferenceParseCharacterOptionString(LALInferenceGetProcParamVal(commandLine,"--psd")->value,&psds,&Npsd);
    }

    if(LALInferenceGetProcParamVal(commandLine,"--dataseed")){
        procparam=LALInferenceGetProcParamVal(commandLine,"--dataseed");
        dataseed=atoi(procparam->value);
    }

    IFOdata=headIFO=calloc(sizeof(LALInferenceIFOData),Nifo);
    if(!IFOdata) XLAL_ERROR_NULL(XLAL_ENOMEM);

    if(LALInferenceGetProcParamVal(commandLine,"--injXML"))
    {
        XLALPrintError("ERROR: --injXML option is deprecated. Use --inj and update your scripts\n");
        exit(1);
    }
    procparam=LALInferenceGetProcParamVal(commandLine,"--inj");
	if(procparam){
        if(LALInferenceGetProcParamVal(commandLine,"--burst_inj")){
        printf("Injecting from burst Table!\n");

        BInjTable=XLALSimBurstTableFromLIGOLw(LALInferenceGetProcParamVal(commandLine,"--inj")->value,0,0);       
        } 
        if(!BInjTable){
                XLALPrintError("Unable to open injection file(LALInferenceReadData) %s\n",procparam->value);
                XLAL_ERROR_NULL(XLAL_EFUNC);
            }
        if(LALInferenceGetProcParamVal(commandLine,"--burst_inj")){
            procparam=LALInferenceGetProcParamVal(commandLine,"--event");
            if(procparam) event=atoi(procparam->value);
            while(q<event) {q++; BInjTable=BInjTable->next;}
        }
    }
    procparam=LALInferenceGetProcParamVal(commandLine,"--psdstart");
    if (!procparam) procparam=LALInferenceGetProcParamVal(commandLine,"--PSDstart");
    LALStringToGPS(&status,&GPSstart,procparam->value,&chartmp);
    if(status.statusCode) REPORTSTATUS(&status);

    if(LALInferenceGetProcParamVal(commandLine,"--trigtime")){
        procparam=LALInferenceGetProcParamVal(commandLine,"--trigtime");
        LALStringToGPS(&status,&GPStrig,procparam->value,&chartmp);
    }
    else{
        if(BInjTable) memcpy(&GPStrig,&(BInjTable->time_geocent_gps),sizeof(GPStrig));
        else {
            XLALPrintError("Error: No trigger time specifed and no injection given \n");
            XLAL_ERROR_NULL(XLAL_EINVAL);
        }
    }
    if(status.statusCode) REPORTSTATUS(&status);
    ppt=LALInferenceGetProcParamVal(commandLine,"--psdlength");
    if(!ppt) ppt=LALInferenceGetProcParamVal(commandLine,"--PSDlength");
    PSDdatalength=atof(ppt->value);
    SegmentLength=atof(LALInferenceGetProcParamVal(commandLine,"--seglen")->value);
    seglen=(size_t)(SegmentLength*SampleRate);
    nSegs=(int)floor(PSDdatalength/SegmentLength);

    for(i=0;i<Nifo;i++) {
        IFOdata[i].fLow=fLows?atof(fLows[i]):defaultFLow; 
        if(fHighs) IFOdata[i].fHigh=fHighs[i]?atof(fHighs[i]):(SampleRate/2.0-(1.0/SegmentLength));
        else IFOdata[i].fHigh=(SampleRate/2.0-(1.0/SegmentLength));
        strncpy(IFOdata[i].name, IFOnames[i], DETNAMELEN);
        IFOdata[i].STDOF = 4.0 / M_PI * nSegs;
        fprintf(stderr, "Detector %s will run with %g DOF if Student's T likelihood used.\n",
                IFOdata[i].name, IFOdata[i].STDOF);
    }

    /* Only allocate this array if there weren't channels read in from the command line */
    if(!dataOpts && !Nchannel) channels=calloc(Nifo,sizeof(char *));
    for(i=0;i<Nifo;i++) {
        if(!dataOpts && !Nchannel) channels[i]=malloc(VARNAME_MAX);
        IFOdata[i].detector=calloc(1,sizeof(LALDetector));

        if(!strcmp(IFOnames[i],"H1")) {			
            memcpy(IFOdata[i].detector,&lalCachedDetectors[LALDetectorIndexLHODIFF],sizeof(LALDetector));
            if(!Nchannel) sprintf((channels[i]),"H1:%s",strainname); continue;}
        if(!strcmp(IFOnames[i],"H2")) {
            memcpy(IFOdata[i].detector,&lalCachedDetectors[LALDetectorIndexLHODIFF],sizeof(LALDetector));
            if(!Nchannel) sprintf((channels[i]),"H2:%s",strainname); continue;}
        if(!strcmp(IFOnames[i],"LLO")||!strcmp(IFOnames[i],"L1")) {
            memcpy(IFOdata[i].detector,&lalCachedDetectors[LALDetectorIndexLLODIFF],sizeof(LALDetector));
            if(!Nchannel) sprintf((channels[i]),"L1:%s",strainname); continue;}
        if(!strcmp(IFOnames[i],"V1")||!strcmp(IFOnames[i],"VIRGO")) {
            memcpy(IFOdata[i].detector,&lalCachedDetectors[LALDetectorIndexVIRGODIFF],sizeof(LALDetector));
            if(!Nchannel) sprintf((channels[i]),"V1:h_16384Hz"); continue;}
        if(!strcmp(IFOnames[i],"GEO")||!strcmp(IFOnames[i],"G1")) {
            memcpy(IFOdata[i].detector,&lalCachedDetectors[LALDetectorIndexGEO600DIFF],sizeof(LALDetector));
            if(!Nchannel) sprintf((channels[i]),"G1:DER_DATA_H"); continue;}
        /*		if(!strcmp(IFOnames[i],"TAMA")||!strcmp(IFOnames[i],"T1")) {memcpy(IFOdata[i].detector,&lalCachedDetectors[LALDetectorIndexTAMA300DIFF]); continue;}*/

        if(!strcmp(IFOnames[i],"E1")){
            memcpy(IFOdata[i].detector,&dE1,sizeof(LALDetector));
            if(!Nchannel) sprintf((channels[i]),"E1:STRAIN"); continue;}
        if(!strcmp(IFOnames[i],"E2")){
            memcpy(IFOdata[i].detector,&dE2,sizeof(LALDetector));
            if(!Nchannel) sprintf((channels[i]),"E2:STRAIN"); continue;}
        if(!strcmp(IFOnames[i],"E3")){
            memcpy(IFOdata[i].detector,&dE3,sizeof(LALDetector));
            if(!Nchannel) sprintf((channels[i]),"E3:STRAIN"); continue;}
        if(!strcmp(IFOnames[i],"HM1")){
            /* Note, this is a sqrt(2)*7.5-km 3rd gen detector */
            LALFrDetector ETHomestakeFr;
            sprintf(ETHomestakeFr.name,"ET-HomeStake1");
            sprintf(ETHomestakeFr.prefix,"M1");
            /* Location of Homestake Mine vertex is */
            /* 44d21'23.11" N, 103d45'54.71" W */
            ETHomestakeFr.vertexLatitudeRadians = (44.+ 21./60  + 23.11/3600)*LAL_PI/180.0;
            ETHomestakeFr.vertexLongitudeRadians = - (103. +45./60 + 54.71/3600)*LAL_PI/180.0;
            ETHomestakeFr.vertexElevation=0.0;
            ETHomestakeFr.xArmAltitudeRadians=0.0;
            ETHomestakeFr.xArmAzimuthRadians=LAL_PI/2.0;
            ETHomestakeFr.yArmAltitudeRadians=0.0;
            ETHomestakeFr.yArmAzimuthRadians=0.0;
            ETHomestakeFr.xArmMidpoint = ETHomestakeFr.yArmMidpoint = sqrt(2.0)*7.5/2.0;
            IFOdata[i].detector=calloc(1,sizeof(LALDetector));
            XLALCreateDetector(IFOdata[i].detector,&ETHomestakeFr,LALDETECTORTYPE_IFODIFF);
            printf("Created Homestake Mine ET detector, location: %lf, %lf, %lf\n",IFOdata[i].detector->location[0],IFOdata[i].detector->location[1],IFOdata[i].detector->location[2]);
            printf("detector tensor:\n");
            for(int jdx=0;jdx<3;jdx++){
                for(j=0;j<3;j++) printf("%f ",IFOdata[i].detector->response[jdx][j]);
                printf("\n");
            }
            continue;
        }
        if(!strcmp(IFOnames[i],"HM2")){
            /* Note, this is a sqrt(2)*7.5-km 3rd gen detector */
            LALFrDetector ETHomestakeFr;
            sprintf(ETHomestakeFr.name,"ET-HomeStake2");
            sprintf(ETHomestakeFr.prefix,"M2");
            /* Location of Homestake Mine vertex is */
            /* 44d21'23.11" N, 103d45'54.71" W */
            ETHomestakeFr.vertexLatitudeRadians = (44.+ 21./60  + 23.11/3600)*LAL_PI/180.0;
            ETHomestakeFr.vertexLongitudeRadians = - (103. +45./60 + 54.71/3600)*LAL_PI/180.0;
            ETHomestakeFr.vertexElevation=0.0;
            ETHomestakeFr.xArmAltitudeRadians=0.0;
            ETHomestakeFr.xArmAzimuthRadians=3.0*LAL_PI/4.0;
            ETHomestakeFr.yArmAltitudeRadians=0.0;
            ETHomestakeFr.yArmAzimuthRadians=LAL_PI/4.0;
            ETHomestakeFr.xArmMidpoint = ETHomestakeFr.yArmMidpoint = sqrt(2.0)*7500./2.0;
            IFOdata[i].detector=calloc(1,sizeof(LALDetector));
            XLALCreateDetector(IFOdata[i].detector,&ETHomestakeFr,LALDETECTORTYPE_IFODIFF);
            printf("Created Homestake Mine ET detector, location: %lf, %lf, %lf\n",IFOdata[i].detector->location[0],IFOdata[i].detector->location[1],IFOdata[i].detector->location[2]);
            printf("detector tensor:\n");
            for(int jdx=0;jdx<3;jdx++){
                for(j=0;j<3;j++) printf("%f ",IFOdata[i].detector->response[jdx][j]);
                printf("\n");
            }
            continue;
        }
        if(!strcmp(IFOnames[i],"EM1")){
            LALFrDetector ETmic1;
            sprintf(ETmic1.name,"ET_Michelson_1");
            sprintf(ETmic1.prefix,"F1");
            ETmic1.vertexLatitudeRadians = (43. + 37./60. + 53.0921/3600)*LAL_PI/180.0;
            ETmic1.vertexLongitudeRadians = (10. + 30./60. + 16.1878/3600.)*LAL_PI/180.0;
            ETmic1.vertexElevation = 0.0;
            ETmic1.xArmAltitudeRadians = ETmic1.yArmAltitudeRadians = 0.0;
            ETmic1.xArmAzimuthRadians = LAL_PI/2.0;
            ETmic1.yArmAzimuthRadians = 0.0;
            ETmic1.xArmMidpoint = ETmic1.yArmMidpoint = sqrt(2.0)*7500./2.;
            IFOdata[i].detector=calloc(1,sizeof(LALDetector));
            XLALCreateDetector(IFOdata[i].detector,&ETmic1,LALDETECTORTYPE_IFODIFF);
            printf("Created ET L-detector 1 (N/E) arms, location: %lf, %lf, %lf\n",IFOdata[i].detector->location[0],IFOdata[i].detector->location[1],IFOdata[i].detector->location[2]);
            printf("detector tensor:\n");
            for(int jdx=0;jdx<3;jdx++){
                for(j=0;j<3;j++) printf("%f ",IFOdata[i].detector->response[jdx][j]);
                printf("\n");
            }
            continue;
        }
        if(!strcmp(IFOnames[i],"EM2")){
            LALFrDetector ETmic2;
            sprintf(ETmic2.name,"ET_Michelson_2");
            sprintf(ETmic2.prefix,"F2");
            ETmic2.vertexLatitudeRadians = (43. + 37./60. + 53.0921/3600)*LAL_PI/180.0;
            ETmic2.vertexLongitudeRadians = (10. + 30./60. + 16.1878/3600.)*LAL_PI/180.0;
            ETmic2.vertexElevation = 0.0;
            ETmic2.xArmAltitudeRadians = ETmic2.yArmAltitudeRadians = 0.0;
            ETmic2.xArmAzimuthRadians = 3.0*LAL_PI/4.0;
            ETmic2.yArmAzimuthRadians = LAL_PI/4.0;
            ETmic2.xArmMidpoint = ETmic2.yArmMidpoint = sqrt(2.0)*7500./2.;
            IFOdata[i].detector=calloc(1,sizeof(LALDetector));
            XLALCreateDetector(IFOdata[i].detector,&ETmic2,LALDETECTORTYPE_IFODIFF);
            printf("Created ET L-detector 2 (NE/SE) arms, location: %lf, %lf, %lf\n",IFOdata[i].detector->location[0],IFOdata[i].detector->location[1],IFOdata[i].detector->location[2]);
            printf("detector tensor:\n");
            for(int jdx=0;jdx<3;jdx++){
                for(j=0;j<3;j++) printf("%f ",IFOdata[i].detector->response[jdx][j]);
                printf("\n");
            }
            continue;
        }
        if(!strcmp(IFOnames[i],"I1")||!strcmp(IFOnames[i],"LIGOIndia")){
            /* Detector in India with 4k arms */
            LALFrDetector LIGOIndiaFr;
            sprintf(LIGOIndiaFr.name,"LIGO_India");
            sprintf(LIGOIndiaFr.prefix,"I1");
            /* Location of India site is */
            /* 14d14' N 76d26' E */
            LIGOIndiaFr.vertexLatitudeRadians = (14. + 14./60.)*LAL_PI/180.0;
            LIGOIndiaFr.vertexLongitudeRadians = (76. + 26./60.)*LAL_PI/180.0;
            LIGOIndiaFr.vertexElevation = 0.0;
            LIGOIndiaFr.xArmAltitudeRadians = 0.0;
            LIGOIndiaFr.yArmAltitudeRadians = 0.0;
            LIGOIndiaFr.yArmMidpoint = 2000.;
            LIGOIndiaFr.xArmMidpoint = 2000.;
            LIGOIndiaFr.xArmAzimuthRadians = LAL_PI/2.;
            LIGOIndiaFr.yArmAzimuthRadians = 0.;
            IFOdata[i].detector=malloc(sizeof(LALDetector));
            memset(IFOdata[i].detector,0,sizeof(LALDetector));
            XLALCreateDetector(IFOdata[i].detector,&LIGOIndiaFr,LALDETECTORTYPE_IFODIFF);
            printf("Created LIGO India Detector, location %lf, %lf, %lf\n",IFOdata[i].detector->location[0],IFOdata[i].detector->location[1],IFOdata[i].detector->location[2]);
            printf("Detector tensor:\n");
            for(int jdx=0;jdx<3;jdx++){
                for(j=0;j<3;j++) printf("%f ",IFOdata[i].detector->response[jdx][j]);
                printf("\n");
            }
            continue;
        }
        if(!strcmp(IFOnames[i],"A1")||!strcmp(IFOnames[i],"LIGOSouth")){
            /* Construct a detector at AIGO with 4k arms */
            LALFrDetector LIGOSouthFr;
            sprintf(LIGOSouthFr.name,"LIGO-South");
            sprintf(LIGOSouthFr.prefix,"A1");
            /* Location of the AIGO detector vertex is */
            /* 31d21'27.56" S, 115d42'50.34"E */
            LIGOSouthFr.vertexLatitudeRadians = - (31. + 21./60. + 27.56/3600.)*LAL_PI/180.0;
            LIGOSouthFr.vertexLongitudeRadians = (115. + 42./60. + 50.34/3600.)*LAL_PI/180.0;
            LIGOSouthFr.vertexElevation=0.0;
            LIGOSouthFr.xArmAltitudeRadians=0.0;
            LIGOSouthFr.xArmAzimuthRadians=AIGOang+LAL_PI/2.;
            LIGOSouthFr.yArmAltitudeRadians=0.0;
            LIGOSouthFr.yArmAzimuthRadians=AIGOang;
            LIGOSouthFr.xArmMidpoint=2000.;
            LIGOSouthFr.yArmMidpoint=2000.;
            IFOdata[i].detector=malloc(sizeof(LALDetector));
            memset(IFOdata[i].detector,0,sizeof(LALDetector));
            XLALCreateDetector(IFOdata[i].detector,&LIGOSouthFr,LALDETECTORTYPE_IFODIFF);
            printf("Created LIGO South detector, location: %lf, %lf, %lf\n",IFOdata[i].detector->location[0],IFOdata[i].detector->location[1],IFOdata[i].detector->location[2]);
            printf("Detector tensor:\n");
            for(int jdx=0;jdx<3;jdx++){
                for(j=0;j<3;j++) printf("%f ",IFOdata[i].detector->response[jdx][j]);
                printf("\n");
            }
            continue;
        }
        if(!strcmp(IFOnames[i],"J1")||!strcmp(IFOnames[i],"LCGT")){
            /* Construct the LCGT telescope */
            REAL8 LCGTangle=19.0*(LAL_PI/180.0);
            LALFrDetector LCGTFr;
            sprintf(LCGTFr.name,"LCGT");
            sprintf(LCGTFr.prefix,"J1");
            LCGTFr.vertexLatitudeRadians  = 36.25 * LAL_PI/180.0;
            LCGTFr.vertexLongitudeRadians = (137.18 * LAL_PI/180.0);
            LCGTFr.vertexElevation=0.0;
            LCGTFr.xArmAltitudeRadians=0.0;
            LCGTFr.xArmAzimuthRadians=LCGTangle+LAL_PI/2.;
            LCGTFr.yArmAltitudeRadians=0.0;
            LCGTFr.yArmAzimuthRadians=LCGTangle;
            LCGTFr.xArmMidpoint=1500.;
            LCGTFr.yArmMidpoint=1500.;
            IFOdata[i].detector=malloc(sizeof(LALDetector));
            memset(IFOdata[i].detector,0,sizeof(LALDetector));
            XLALCreateDetector(IFOdata[i].detector,&LCGTFr,LALDETECTORTYPE_IFODIFF);
            printf("Created LCGT telescope, location: %lf, %lf, %lf\n",IFOdata[i].detector->location[0],IFOdata[i].detector->location[1],IFOdata[i].detector->location[2]);
            printf("Detector tensor:\n");
            for(int jdx=0;jdx<3;jdx++){
                for(j=0;j<3;j++) printf("%f ",IFOdata[i].detector->response[jdx][j]);
                printf("\n");
            }
            continue;
        }
        fprintf(stderr,"Unknown interferometer %s. Valid codes: H1 H2 L1 V1 GEO A1 J1 I1 E1 E2 E3 HM1 HM2 EM1 EM2\n",IFOnames[i]); exit(-1);
    }

    /* Set up FFT structures and window */
    for (i=0;i<Nifo;i++){
        /* Create FFT plans */
        IFOdata[i].timeToFreqFFTPlan = XLALCreateForwardREAL8FFTPlan((UINT4) seglen, 0 );
        if(!IFOdata[i].timeToFreqFFTPlan) XLAL_ERROR_NULL(XLAL_EFUNC);
        IFOdata[i].freqToTimeFFTPlan = XLALCreateReverseREAL8FFTPlan((UINT4) seglen,0);
        if(!IFOdata[i].freqToTimeFFTPlan) XLAL_ERROR_NULL(XLAL_EFUNC);		
        /* Setup windows */
        IFOdata[i].window=XLALCreateTukeyREAL8Window(seglen,(REAL8)2.0*padding*SampleRate/(REAL8)seglen);
        if(!IFOdata[i].window) XLAL_ERROR_NULL(XLAL_EFUNC);
    }

    /* Trigger time = 2 seconds before end of segment (was 1 second, but Common Inputs for The Events are -6 +2*/
    memcpy(&segStart,&GPStrig,sizeof(LIGOTimeGPS));
    XLALGPSAdd(&segStart,-SegmentLength+2);


    /* Read the PSD data */
    for(i=0;i<Nifo;i++) {
        memcpy(&(IFOdata[i].epoch),&segStart,sizeof(LIGOTimeGPS));
        /* Check to see if an interpolation file is specified */
        interpFlag=0;
        interp=NULL;
        if(strstr(caches[i],"interp:")==caches[i]){
            /* Extract the file name */
            char *interpfilename=&(caches[i][7]);
            printf("Looking for interpolation file %s\n",interpfilename);
            interpFlag=1;
            interp=interpFromFile(interpfilename);
        }    
        /* Check if fake data is requested */
        if(interpFlag || (!(strcmp(caches[i],"LALLIGO") && strcmp(caches[i],"LALVirgo") && strcmp(caches[i],"LALGEO") && strcmp(caches[i],"LALEGO")
                        && strcmp(caches[i],"LALAdLIGO"))))
        {
            //FakeFlag=1; - set but not used
            if (!LALInferenceGetProcParamVal(commandLine,"--dataseed")){
                fprintf(stderr,"Error: You need to specify a dataseed when generating data with --dataseed <number>.\n\
                        (--dataseed 0 uses a non-reproducible number from the system clock, and no parallel run is then possible.)\n" );
                exit(-1);
            }
            datarandparam=XLALCreateRandomParams(dataseed?dataseed+(int)i:dataseed);
            if(!datarandparam) XLAL_ERROR_NULL(XLAL_EFUNC);
            /* Selection of the noise curve */
            if(!strcmp(caches[i],"LALLIGO")) {PSD = &LALLIGOIPsd; scalefactor=9E-46;}
            if(!strcmp(caches[i],"LALVirgo")) {PSD = &LALVIRGOPsd; scalefactor=1.0;}
            if(!strcmp(caches[i],"LALGEO")) {PSD = &LALGEOPsd; scalefactor=1E-46;}
            if(!strcmp(caches[i],"LALEGO")) {PSD = &LALEGOPsd; scalefactor=1.0;}
            if(!strcmp(caches[i],"LALAdLIGO")) {PSD = &LALAdvLIGOPsd; scalefactor = 1E-49;}
            if(interpFlag) {PSD=NULL; scalefactor=1.0;}
            //if(!strcmp(caches[i],"LAL2kLIGO")) {PSD = &LALAdvLIGOPsd; scalefactor = 36E-46;}
            if(PSD==NULL && !interpFlag) {fprintf(stderr,"Error: unknown simulated PSD: %s\n",caches[i]); exit(-1);}


            IFOdata[i].oneSidedNoisePowerSpectrum=(REAL8FrequencySeries *)
                XLALCreateREAL8FrequencySeries("spectrum",&GPSstart,0.0,
                        (REAL8)(SampleRate)/seglen,&lalDimensionlessUnit,seglen/2 +1);
            if(!IFOdata[i].oneSidedNoisePowerSpectrum) XLAL_ERROR_NULL(XLAL_EFUNC);
            for(j=0;j<IFOdata[i].oneSidedNoisePowerSpectrum->data->length;j++)
            {
                MetaNoiseFunc(&status,&(IFOdata[i].oneSidedNoisePowerSpectrum->data->data[j]),j*IFOdata[i].oneSidedNoisePowerSpectrum->deltaF,interp,PSD);
                //PSD(&status,&(IFOdata[i].oneSidedNoisePowerSpectrum->data->data[j]),j*IFOdata[i].oneSidedNoisePowerSpectrum->deltaF);
                IFOdata[i].oneSidedNoisePowerSpectrum->data->data[j]*=scalefactor;
            }
            IFOdata[i].freqData = (COMPLEX16FrequencySeries *)XLALCreateCOMPLEX16FrequencySeries("stilde",&segStart,0.0,IFOdata[i].oneSidedNoisePowerSpectrum->deltaF,&lalDimensionlessUnit,seglen/2 +1);
            if(!IFOdata[i].freqData) XLAL_ERROR_NULL(XLAL_EFUNC);

            /* Create the fake data */
            int j_Lo = (int) IFOdata[i].fLow/IFOdata[i].freqData->deltaF;
            if(LALInferenceGetProcParamVal(commandLine,"--0noise")){
                for(j=j_Lo;j<IFOdata[i].freqData->data->length;j++){
                    IFOdata[i].freqData->data->data[j].re=IFOdata[i].freqData->data->data[j].im=0.0;
                }
            } else {
                for(j=j_Lo;j<IFOdata[i].freqData->data->length;j++){
                    IFOdata[i].freqData->data->data[j].re=XLALNormalDeviate(datarandparam)*(0.5*sqrt(IFOdata[i].oneSidedNoisePowerSpectrum->data->data[j]/IFOdata[i].freqData->deltaF));
                    IFOdata[i].freqData->data->data[j].im=XLALNormalDeviate(datarandparam)*(0.5*sqrt(IFOdata[i].oneSidedNoisePowerSpectrum->data->data[j]/IFOdata[i].freqData->deltaF));
                }
            }
            IFOdata[i].freqData->data->data[0].re=0; 			IFOdata[i].freqData->data->data[0].im=0;
            const char timename[]="timeData";
            IFOdata[i].timeData=(REAL8TimeSeries *)XLALCreateREAL8TimeSeries(timename,&segStart,0.0,(REAL8)1.0/SampleRate,&lalDimensionlessUnit,(size_t)seglen);
            if(!IFOdata[i].timeData) XLAL_ERROR_NULL(XLAL_EFUNC);
            XLALREAL8FreqTimeFFT(IFOdata[i].timeData,IFOdata[i].freqData,IFOdata[i].freqToTimeFFTPlan);
            if(*XLALGetErrnoPtr()) printf("XLErr: %s\n",XLALErrorString(*XLALGetErrnoPtr()));
            XLALDestroyRandomParams(datarandparam);
        }
        else{ /* Not using fake data, load the data from a cache file */
            if (LALInferenceGetProcParamVal(commandLine, "--psd")){
                interp=NULL;
                char *interpfilename=&(psds[i][0]);
                fprintf(stderr,"Reading PSD for %s using %s\n",IFOnames[i],interpfilename);
                printf("Looking for psd interpolation file %s\n",interpfilename);
                interp=interpFromFile(interpfilename);
                IFOdata[i].oneSidedNoisePowerSpectrum=(REAL8FrequencySeries *)
                    XLALCreateREAL8FrequencySeries("spectrum",&GPSstart,0.0,
                            (REAL8)(SampleRate)/seglen,&lalDimensionlessUnit,seglen/2 +1);
                if(!IFOdata[i].oneSidedNoisePowerSpectrum) XLAL_ERROR_NULL(XLAL_EFUNC);
                for(j=0;j<IFOdata[i].oneSidedNoisePowerSpectrum->data->length;j++)
                {
                    MetaNoiseFunc(&status,&(IFOdata[i].oneSidedNoisePowerSpectrum->data->data[j]),j*IFOdata[i].oneSidedNoisePowerSpectrum->deltaF,interp,NULL);
                }
            }else{
                fprintf(stderr,"Estimating PSD for %s using %i segments of %i samples (%lfs)\n",IFOnames[i],nSegs,(int)seglen,SegmentLength);
                PSDtimeSeries=readTseries(caches[i],channels[i],GPSstart,PSDdatalength);
                if(!PSDtimeSeries) {XLALPrintError("Error reading PSD data for %s\n",IFOnames[i]); XLAL_ERROR_NULL(XLAL_EFUNC);}
                XLALResampleREAL8TimeSeries(PSDtimeSeries,1.0/SampleRate);
                PSDtimeSeries=(REAL8TimeSeries *)XLALShrinkREAL8TimeSeries(PSDtimeSeries,(size_t) 0, (size_t) seglen*nSegs);
                if(!PSDtimeSeries) {
                    fprintf(stderr,"ERROR while estimating PSD for %s\n",IFOnames[i]);
                    XLAL_ERROR_NULL(XLAL_EFUNC);
                }
                IFOdata[i].oneSidedNoisePowerSpectrum=(REAL8FrequencySeries *)XLALCreateREAL8FrequencySeries("spectrum",&PSDtimeSeries->epoch,0.0,(REAL8)(SampleRate)/seglen,&lalDimensionlessUnit,seglen/2 +1);
                if(!IFOdata[i].oneSidedNoisePowerSpectrum) XLAL_ERROR_NULL(XLAL_EFUNC);
                if (LALInferenceGetProcParamVal(commandLine, "--PSDwelch"))
                    XLALREAL8AverageSpectrumWelch(IFOdata[i].oneSidedNoisePowerSpectrum ,PSDtimeSeries, seglen, (UINT4)seglen, IFOdata[i].window, IFOdata[i].timeToFreqFFTPlan);
                else
                    XLALREAL8AverageSpectrumMedian(IFOdata[i].oneSidedNoisePowerSpectrum ,PSDtimeSeries, seglen, (UINT4)seglen, IFOdata[i].window, IFOdata[i].timeToFreqFFTPlan);	

                XLALDestroyREAL8TimeSeries(PSDtimeSeries);
            }

            /* Read the data segment */
            LIGOTimeGPS truesegstart=segStart;
            if(Ntimeslides) {
                REAL4 deltaT=-atof(timeslides[i]);
                XLALGPSAdd(&segStart, deltaT);
                fprintf(stderr,"Slid %s by %f s from %10.10lf to %10.10lf\n",IFOnames[i],deltaT,truesegstart.gpsSeconds+1e-9*truesegstart.gpsNanoSeconds,segStart.gpsSeconds+1e-9*segStart.gpsNanoSeconds);
            }
            IFOdata[i].timeData=readTseries(caches[i],channels[i],segStart,SegmentLength);
            segStart=truesegstart;
            if(Ntimeslides) IFOdata[i].timeData->epoch=truesegstart;
            /* FILE *out; */
            /* char fileName[256]; */
            /* snprintf(fileName, 256, "readTimeData-%d.dat", i); */
            /* out = fopen(fileName, "w"); */
            /* for (j = 0; j < IFOdata[i].timeData->data->length; j++) { */
            /*   fprintf(out, "%g %g\n", j*IFOdata[i].timeData->deltaT, IFOdata[i].timeData->data->data[j]); */
            /* } */
            /* fclose(out); */

            if(!IFOdata[i].timeData) {
                XLALPrintError("Error reading segment data for %s at %i\n",IFOnames[i],segStart.gpsSeconds);
                XLAL_ERROR_NULL(XLAL_EFUNC);
            }
            XLALResampleREAL8TimeSeries(IFOdata[i].timeData,1.0/SampleRate);	 
            if(!IFOdata[i].timeData) {XLALPrintError("Error reading segment data for %s\n",IFOnames[i]); XLAL_ERROR_NULL(XLAL_EFUNC);}
            IFOdata[i].freqData=(COMPLEX16FrequencySeries *)XLALCreateCOMPLEX16FrequencySeries("freqData",&(IFOdata[i].timeData->epoch),0.0,1.0/SegmentLength,&lalDimensionlessUnit,seglen/2+1);
            if(!IFOdata[i].freqData) XLAL_ERROR_NULL(XLAL_EFUNC);
            IFOdata[i].windowedTimeData=(REAL8TimeSeries *)XLALCreateREAL8TimeSeries("windowed time data",&(IFOdata[i].timeData->epoch),0.0,1.0/SampleRate,&lalDimensionlessUnit,seglen);
            if(!IFOdata[i].windowedTimeData) XLAL_ERROR_NULL(XLAL_EFUNC);
            XLALDDVectorMultiply(IFOdata[i].windowedTimeData->data,IFOdata[i].timeData->data,IFOdata[i].window->data);
            XLALREAL8TimeFreqFFT(IFOdata[i].freqData,IFOdata[i].windowedTimeData,IFOdata[i].timeToFreqFFTPlan);

            for(j=0;j<IFOdata[i].freqData->data->length;j++){
                IFOdata[i].freqData->data->data[j].re/=sqrt(IFOdata[i].window->sumofsquares / IFOdata[i].window->data->length);
                IFOdata[i].freqData->data->data[j].im/=sqrt(IFOdata[i].window->sumofsquares / IFOdata[i].window->data->length);
                IFOdata[i].windowedTimeData->data->data[j] /= sqrt(IFOdata[i].window->sumofsquares / IFOdata[i].window->data->length);
            }
        } /* End of data reading process */

        //		/* Now that the PSD is set up, make the TDW. */
        //    IFOdata[i].timeDomainNoiseWeights = 
        //                  (REAL8TimeSeries *)XLALCreateREAL8TimeSeries("time domain weights", 
        //                                                               &(IFOdata[i].oneSidedNoisePowerSpectrum->epoch),
        //                                                               0.0,
        //                                                               1.0/SampleRate,
        //                                                               &lalDimensionlessUnit,
        //                                                               seglen);
        //		if(!IFOdata[i].timeDomainNoiseWeights) XLAL_ERROR_NULL(XLAL_EFUNC);
        //		LALInferencePSDToTDW(IFOdata[i].timeDomainNoiseWeights, IFOdata[i].oneSidedNoisePowerSpectrum, IFOdata[i].freqToTimeFFTPlan,
        //                         IFOdata[i].fLow, IFOdata[i].fHigh);

        makeWhiteData(&(IFOdata[i]));

        if (LALInferenceGetProcParamVal(commandLine, "--spinspiralPSD")) {
            FILE *in;
            //char fileNameIn[256];
            //snprintf(fileNameIn, 256, spinspiralPSD);
            double freq_temp, psd_temp, temp;
            int n=0;
            int k=0;
            int templen=0;
            char buffer[256];
            char * line=buffer;

            //in = fopen(fileNameIn, "r");
            in = fopen(spinspiralPSD[i], "r");
            while(fgets(buffer, 256, in)){
                templen++;
            }

            // REAL8 *tempPSD = NULL;
            // REAL8 *tempfreq = NULL;
            // tempPSD=calloc(sizeof(REAL8),templen+1);
            // tempfreq=calloc(sizeof(REAL8),templen+1);

            rewind(in);
            IFOdata[i].oneSidedNoisePowerSpectrum->data->data[0] = 1.0;
            while(fgets(buffer, 256, in)){
                line=buffer;

                sscanf(line, "%lg%n", &freq_temp,&n);
                line+=n;
                sscanf(line, "%lg%n", &psd_temp,&n);
                line+=n;
                sscanf(line, "%lg%n", &temp,&n);
                line+=n;

                // tempfreq[k]=freq_temp;
                // tempPSD[k]=psd_temp*psd_temp;

                IFOdata[i].oneSidedNoisePowerSpectrum->data->data[k+1]=psd_temp*psd_temp;

                k++;
                //fprintf(stdout, "%g %g \n",freq_temp, psd_temp); fflush(stdout);
            }
            fclose(in);
        }

        if (LALInferenceGetProcParamVal(commandLine, "--data-dump")) {
            const UINT4 nameLength=256;
            char filename[nameLength];
            FILE *out;

            snprintf(filename, nameLength, "%s-PSD.dat", IFOdata[i].name);
            out = fopen(filename, "w");
            for (j = 0; j < IFOdata[i].oneSidedNoisePowerSpectrum->data->length; j++) {
                REAL8 f = IFOdata[i].oneSidedNoisePowerSpectrum->deltaF*j;
                REAL8 psd = IFOdata[i].oneSidedNoisePowerSpectrum->data->data[j];

                fprintf(out, "%g %g\n", f, psd);
            }
            fclose(out);

            snprintf(filename, nameLength, "%s-timeData.dat", IFOdata[i].name);
            out = fopen(filename, "w");
            for (j = 0; j < IFOdata[i].timeData->data->length; j++) {
                REAL8 t = XLALGPSGetREAL8(&(IFOdata[i].timeData->epoch)) + 
                    j * IFOdata[i].timeData->deltaT;
                REAL8 d = IFOdata[i].timeData->data->data[j];

                fprintf(out, "%.6f %g\n", t, d);
            }
            fclose(out);

            snprintf(filename, nameLength, "%s-freqData.dat", IFOdata[i].name);
            out = fopen(filename, "w");
            for (j = 0; j < IFOdata[i].freqData->data->length; j++) {
                REAL8 f = IFOdata[i].freqData->deltaF * j;
                REAL8 dre = IFOdata[i].freqData->data->data[j].re;
                REAL8 dim = IFOdata[i].freqData->data->data[j].im;

                fprintf(out, "%g %g %g\n", f, dre, dim);
            }
            fclose(out);

        }

    }

    for (i=0;i<Nifo;i++) IFOdata[i].SNR=0.0; //SNR of the injection ONLY IF INJECTION. Set to 0.0 by default.

    for (i=0;i<Nifo-1;i++) IFOdata[i].next=&(IFOdata[i+1]);

    for(i=0;i<Nifo;i++) {
        if(channels) if(channels[i]) free(channels[i]);
        if(caches) if(caches[i]) free(caches[i]);
        if(IFOnames) if(IFOnames[i]) free(IFOnames[i]);
        if(fLows) if(fLows[i]) free(fLows[i]);
        if(fHighs) if(fHighs[i]) free(fHighs[i]);
    }
    if(channels) free(channels);
    if(caches) free(caches);
    if(IFOnames) free(IFOnames);
    if(fLows) free(fLows);
    if(fHighs) free(fHighs);

    return headIFO;
}

void LALInferenceInjectBurstSignal(LALInferenceRunState *irs, ProcessParamsTable *commandLine)
{
	LALStatus status;
	memset(&status,0,sizeof(status));
	SimBurst *injTable=NULL;
    SimBurst *injEvent=NULL;
	UINT4 Ninj=0;
	UINT4 event=0;
	UINT4 i=0,j=0;
    char *SNRpath = NULL;
    //REAL8 responseScale=1.0;
	//CoherentGW InjectGW;
	//PPNParamStruc InjParams;
	LIGOTimeGPS injstart;
	REAL8 SNR=0.0,NetworkSNR=0.0, previous_snr=0.0; 
	UINT4 best_ifo_snr=0;//best_ifo_snr+=1;
	UINT4 highest_snr_index=0;
	DetectorResponse det;
	memset(&injstart,0,sizeof(LIGOTimeGPS));
	//memset(&InjParams,0,sizeof(PPNParamStruc));
	COMPLEX16FrequencySeries *injF=NULL;
	//FILE *rawWaveform=NULL;
	ProcessParamsTable *ppt=NULL;
	REAL8 bufferLength = 2048.0; /* Default length of buffer for injections (seconds) */
	UINT4 bufferN=0;
	LIGOTimeGPS bufferStart;

	LALInferenceIFOData *IFOdata=irs->data;
	LALInferenceIFOData *thisData=IFOdata->next;
	REAL8 minFlow=IFOdata->fLow;
	REAL8 MindeltaT=IFOdata->timeData->deltaT;
    REAL8 InjSampleRate=1.0/MindeltaT;
	REAL4TimeSeries *injectionBuffer=NULL;
    REAL8 padding=0.4; //default, set in LALInferenceReadData()
	
  
	while(thisData){
          minFlow   = minFlow>thisData->fLow ? thisData->fLow : minFlow;
          MindeltaT = MindeltaT>thisData->timeData->deltaT ? thisData->timeData->deltaT : MindeltaT;
          thisData  = thisData->next;
	}
	thisData=IFOdata;
	//InjParams.deltaT = MindeltaT;
	//InjParams.fStartIn=(REAL4)minFlow;

	if(!LALInferenceGetProcParamVal(commandLine,"--inj")) {fprintf(stdout,"No injection file specified, not injecting\n"); return;}
	if(LALInferenceGetProcParamVal(commandLine,"--event")){
    event= atoi(LALInferenceGetProcParamVal(commandLine,"--event")->value);
    fprintf(stdout,"Injecting event %d\n",event);
	}
        if(LALInferenceGetProcParamVal(commandLine,"--snrpath")){
                ppt = LALInferenceGetProcParamVal(commandLine,"--snrpath");
		SNRpath = calloc(strlen(ppt->value)+1,sizeof(char));
		memcpy(SNRpath,ppt->value,strlen(ppt->value)+1);
                fprintf(stdout,"Writing SNRs in %s\n",SNRpath)     ;

	}
	injTable=XLALSimBurstTableFromLIGOLw(LALInferenceGetProcParamVal(commandLine,"--inj")->value,0,0);
	REPORTSTATUS(&status);
    Ninj=-1;
    while(injTable){Ninj++;injTable=injTable->next;}
	printf("Ninj %d\n", Ninj);
	if(Ninj<event) fprintf(stderr,"Error reading event %d from %s\n",event,LALInferenceGetProcParamVal(commandLine,"--inj")->value);
	injTable=XLALSimBurstTableFromLIGOLw(LALInferenceGetProcParamVal(commandLine,"--inj")->value,0,0);
	while(i<event) {i++; injTable = injTable->next;} /* Select event */
	injEvent = injTable;
	injEvent->next = NULL;
	
	//memset(&InjectGW,0,sizeof(InjectGW));
	
	REPORTSTATUS(&status);
	//LALGenerateInspiral(&status,&InjectGW,injTable,&InjParams);
	//if(status.statusCode!=0) {fprintf(stderr,"Error generating injection!\n"); REPORTSTATUS(&status); }
	
	/* Begin loop over interferometers */
	while(thisData){
		
		InjSampleRate=1.0/thisData->timeData->deltaT;
		if(LALInferenceGetProcParamVal(commandLine,"--injectionsrate")) InjSampleRate=atof(LALInferenceGetProcParamVal(commandLine,"--injectionsrate")->value);
		
		memset(&det,0,sizeof(det));
		det.site=thisData->detector;
		COMPLEX8FrequencySeries *resp = XLALCreateCOMPLEX8FrequencySeries("response",&thisData->timeData->epoch,
																		  0.0,
																		  thisData->freqData->deltaF,
																		  &strainPerCount,
																		  thisData->freqData->data->length);
		
		for(i=0;i<resp->data->length;i++) {resp->data->data[i].re=(REAL4)1.0; resp->data->data[i].im=0.0;}
		/* Originally created for injecting into DARM-ERR, so transfer function was needed.  
		But since we are injecting into h(t), the transfer function from h(t) to h(t) is 1.*/

		/* We need a long buffer to inject into so that FindChirpInjectSignals() works properly
		 for low mass systems. Use 100 seconds here */
		bufferN = (UINT4) (bufferLength*InjSampleRate);// /thisData->timeData->deltaT);
		memcpy(&bufferStart,&thisData->timeData->epoch,sizeof(LIGOTimeGPS));
		XLALGPSAdd(&bufferStart,(REAL8) thisData->timeData->data->length * thisData->timeData->deltaT);
		XLALGPSAdd(&bufferStart,-bufferLength);
		injectionBuffer=(REAL4TimeSeries *)XLALCreateREAL4TimeSeries(thisData->detector->frDetector.prefix,
																	 &bufferStart, 0.0, 1.0/InjSampleRate,//thisData->timeData->deltaT,
																	 &lalADCCountUnit, bufferN);
		REAL8TimeSeries *inj8Wave=(REAL8TimeSeries *)XLALCreateREAL8TimeSeries("injection8",
                                                                           &thisData->timeData->epoch,
                                                                           0.0,
                                                                           thisData->timeData->deltaT,
                                                                           //&lalDimensionlessUnit,
                                                                           &lalStrainUnit,
                                                                           thisData->timeData->data->length);
		if(!inj8Wave) XLAL_ERROR_VOID(XLAL_EFUNC);
		/* This marks the sample in which the real segment starts, within the buffer */
		for(i=0;i<injectionBuffer->data->length;i++) injectionBuffer->data->data[i]=0.0;
		for(i=0;i<inj8Wave->data->length;i++) inj8Wave->data->data[i]=0.0;
		INT4 realStartSample=(INT4)((thisData->timeData->epoch.gpsSeconds - injectionBuffer->epoch.gpsSeconds)/thisData->timeData->deltaT);
		realStartSample+=(INT4)((thisData->timeData->epoch.gpsNanoSeconds - injectionBuffer->epoch.gpsNanoSeconds)*1e-9/thisData->timeData->deltaT);

		/*LALSimulateCoherentGW(&status,injWave,&InjectGW,&det);*/
    //LALFindChirpInjectSignals(&status,injectionBuffer,injEvent,resp);

    if(LALInferenceGetProcParamVal(commandLine,"--lalsimulationinjection")){
      printf("Using LALSimulation for injection\n");
      REAL8TimeSeries *hplus=NULL;  /**< +-polarization waveform */
      REAL8TimeSeries *hcross=NULL; /**< x-polarization waveform */
      REAL8TimeSeries       *signalvecREAL8=NULL;
      REAL8 Q, centre_frequency,hrss,eccentricity,polar_angle;
      
      Q=injEvent->q;
      centre_frequency=injEvent->frequency;
      hrss=injEvent->hrss;
     // psi=injEvent->psi;
      polar_angle=injEvent->pol_ellipse_angle;
      eccentricity=injEvent->pol_ellipse_e; // salvo
      //delta_t= LAL_SQRT2*Q/(LAL_TWOPI*centre_frequency); // salvo
      printf("READ polar angle %lf\n",polar_angle); 
      printf("READ ecc %lf\n",injEvent->pol_ellipse_e); 

      XLALSimBurstSineGaussian(&hplus,&hcross, Q, centre_frequency,hrss,eccentricity,polar_angle,thisData->timeData->deltaT);

 
      XLALResampleREAL8TimeSeries(hplus,thisData->timeData->deltaT);
      XLALResampleREAL8TimeSeries(hcross,thisData->timeData->deltaT);
      if(!hplus || !hcross) {
        fprintf(stderr,"Error: XLALSimInspiralChooseWaveform() failed to produce waveform.\n");
        exit(-1);
        //XLALPrintError("XLALSimInspiralChooseWaveform() failed to produce waveform.\n");
        //XLAL_ERROR_VOID(XLAL_EFUNC);
      }
      /* XLALSimInspiralChooseTDWaveform always ends the waveform at t=0 */
      /* So we can adjust the epoch so that the end time is as desired */
      XLALGPSAddGPS(&(hplus->epoch), &(injEvent->time_geocent_gps));       // SALVO check this
      XLALGPSAddGPS(&(hcross->epoch), &(injEvent->time_geocent_gps));
      
      signalvecREAL8=XLALSimDetectorStrainREAL8TimeSeries(hplus, hcross, injEvent->ra, injEvent->dec, injEvent->psi, det.site);
      if (!signalvecREAL8) XLAL_ERROR_VOID(XLAL_EFUNC);
      
      for(i=0;i<signalvecREAL8->data->length;i++){
        if(isnan(signalvecREAL8->data->data[i])) signalvecREAL8->data->data[i]=0.0;
      }
      
      if(signalvecREAL8->data->length > thisData->timeData->data->length-(UINT4)ceil((2.0*padding+2.0)/thisData->timeData->deltaT)){
        fprintf(stderr, "WARNING: waveform length = %u is longer than thisData->timeData->data->length = %d minus the window width = %d and the 2.0 seconds after tc (total of %d points available).\n", signalvecREAL8->data->length, thisData->timeData->data->length, (INT4)ceil((2.0*padding)/thisData->timeData->deltaT) , thisData->timeData->data->length-(INT4)ceil((2.0*padding+2.0)/thisData->timeData->deltaT));
        fprintf(stderr, "The waveform injected is %f seconds long. Consider increasing the %f seconds segment length (--seglen) to be greater than %f. (in %s, line %d)\n",signalvecREAL8->data->length * thisData->timeData->deltaT , thisData->timeData->data->length * thisData->timeData->deltaT, signalvecREAL8->data->length * thisData->timeData->deltaT + 2.0*padding + 2.0, __FILE__, __LINE__);
      }
      XLALSimAddInjectionREAL8TimeSeries(inj8Wave, signalvecREAL8, NULL);
      
      if ( hplus ) XLALDestroyREAL8TimeSeries(hplus);
      if ( hcross ) XLALDestroyREAL8TimeSeries(hcross);
      
    }
    XLALDestroyREAL4TimeSeries(injectionBuffer);
    injF=(COMPLEX16FrequencySeries *)XLALCreateCOMPLEX16FrequencySeries("injF",
										&thisData->timeData->epoch,
										0.0,
										thisData->freqData->deltaF,
										&lalDimensionlessUnit,
										thisData->freqData->data->length);
    if(!injF) {
      XLALPrintError("Unable to allocate memory for injection buffer\n");
      XLAL_ERROR_VOID(XLAL_EFUNC);
    }
    /* Window the data */
    REAL4 WinNorm = sqrt(thisData->window->sumofsquares/thisData->window->data->length);
        for(j=0;j<inj8Wave->data->length;j++) inj8Wave->data->data[j]*=thisData->window->data->data[j]; /* /WinNorm; */ /* Window normalisation applied only in freq domain */
    XLALREAL8TimeFreqFFT(injF,inj8Wave,thisData->timeToFreqFFTPlan);
    /*for(j=0;j<injF->data->length;j++) printf("%lf\n",injF->data->data[j].re);*/
    if(thisData->oneSidedNoisePowerSpectrum){
	for(SNR=0.0,j=thisData->fLow/injF->deltaF;j<injF->data->length;j++){
	  SNR+=pow(injF->data->data[j].re,2.0)/thisData->oneSidedNoisePowerSpectrum->data->data[j];
	  SNR+=pow(injF->data->data[j].im,2.0)/thisData->oneSidedNoisePowerSpectrum->data->data[j];
	}
        SNR*=4.0*injF->deltaF;
    }
    thisData->SNR=sqrt(SNR);
    NetworkSNR+=SNR;
    
    if (thisData->SNR > previous_snr) {best_ifo_snr=highest_snr_index;    previous_snr=thisData->SNR;}
    highest_snr_index++;

    if (!(SNRpath==NULL)){ /* If the user provided a path with --snrpath store a file with injected SNRs */
      //PrintSNRsToFile(IFOdata , &injTable);
    }
    /* Actually inject the waveform */
    for(j=0;j<inj8Wave->data->length;j++) thisData->timeData->data->data[j]+=inj8Wave->data->data[j];
      fprintf(stdout,"Injected SNR in detector %s = %10.10e\n",thisData->name,thisData->SNR);
      char filename[256];
      sprintf(filename,"%s_timeInjection.dat",thisData->name);
      FILE* file=fopen(filename, "w");
      for(j=0;j<inj8Wave->data->length;j++){   
	fprintf(file, "%.6f\t%lg\n", XLALGPSGetREAL8(&thisData->timeData->epoch) + thisData->timeData->deltaT*j, inj8Wave->data->data[j]);
      }
      fclose(file);
      sprintf(filename,"%s_freqInjection.dat",thisData->name);
      file=fopen(filename, "w");
      for(j=0;j<injF->data->length;j++){   
	thisData->freqData->data->data[j].re+=injF->data->data[j].re/WinNorm;
	thisData->freqData->data->data[j].im+=injF->data->data[j].im/WinNorm;
	fprintf(file, "%lg %lg \t %lg\n", thisData->freqData->deltaF*j, injF->data->data[j].re, injF->data->data[j].im);
      }
      fclose(file);
    
      XLALDestroyREAL8TimeSeries(inj8Wave);
      XLALDestroyCOMPLEX16FrequencySeries(injF);
      thisData=thisData->next;
    }
    NetworkSNR=sqrt(NetworkSNR);
    fprintf(stdout,"Network SNR of event %d = %g\n",event,NetworkSNR);
    
    thisData=IFOdata;
    //LALInferenceIFOData *IFOdataRed=NULL;
    UINT4 Nifo=3;
    //IFOdataRed=calloc(sizeof(LALInferenceIFOData),Nifo-1);
    highest_snr_index=0;
    i=0;
    
    while(thisData){
    thisData->BestIFO=LALMalloc (sizeof(LALInferenceBestIFO));
    thisData->BestIFO->detector= LALMalloc (sizeof(LALDetector));
    thisData->BestIFO->TemplateFromInjection=(COMPLEX16FrequencySeries *)XLALCreateCOMPLEX16FrequencySeries("injastempl",&thisData->timeData->epoch,
										0.0,
										thisData->freqData->deltaF,
										&lalDimensionlessUnit,
										thisData->freqData->data->length);
    thisData=thisData->next;
    }
    thisData=IFOdata;
    LALInferenceIFOData *thisData2=NULL;
    thisData2=IFOdata;
    
   while(thisData){
	    if (best_ifo_snr==highest_snr_index){
            while(thisData2){
                for(j=0;j<injF->data->length;j++){   
                    thisData2->BestIFO->TemplateFromInjection->data->data[j].re=thisData->freqData->data->data[j].re;
                    thisData2->BestIFO->TemplateFromInjection->data->data[j].im=thisData->freqData->data->data[j].im;
                  }
                memcpy(thisData2->BestIFO->detector,thisData->detector,sizeof(LALDetector));
                thisData2=thisData2->next;
                printf("doing \n");
            }
            break;
            }
      else
		highest_snr_index++;
		thisData=thisData->next;
		//Salvo: I have to remove the IFO which gave the best snr from the poll
		
		}
    
    printf("HSNR %d, nifo %d \n",highest_snr_index,Nifo);
    
    if (highest_snr_index==0)
    irs->data=&(IFOdata[1]);
    else if (highest_snr_index==2)
    IFOdata[highest_snr_index-1].next=NULL;
    else{
    IFOdata[highest_snr_index-1].next=&(IFOdata[highest_snr_index+1]);
    }
    return;
}

/** Fill the variables passed in vars with the parameters of the injection passed in event
    will over-write and destroy any existing parameters. Param vary type will be fixed */
void LALInferenceBurstInjectionToVariables(SimBurst *theEventTable, LALInferenceVariables *vars)
{
    if(!vars) {
	XLALPrintError("Encountered NULL variables pointer");
   	XLAL_ERROR_VOID(XLAL_EINVAL);
	}
    /* Destroy existing parameters */
    if(vars->head!=NULL) LALInferenceDestroyVariables(vars);
    REAL8 q = theEventTable->q;
    REAL8 psi = theEventTable->psi;
    REAL8 injGPSTime = XLALGPSGetREAL8(&(theEventTable->time_geocent_gps));
    REAL8 hrss = theEventTable->hrss;
    REAL8 loghrss=log(hrss);
    REAL8 f0 = theEventTable->frequency;
    REAL8 pol_angle = theEventTable->pol_ellipse_angle;
    REAL8 eccentricity = theEventTable->pol_ellipse_e;

    REAL8 dec = theEventTable->dec;
    REAL8 ra = theEventTable->ra;
    
    LALInferenceAddVariable(vars, "Q", &q, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
    LALInferenceAddVariable(vars, "frequency", &f0, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
    LALInferenceAddVariable(vars, "time", &injGPSTime, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
   // LALInferenceAddVariable(vars, "hrss", &hrss, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
    LALInferenceAddVariable(vars, "polar_angle", &pol_angle, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
    LALInferenceAddVariable(vars, "eccentricity", &eccentricity, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
    LALInferenceAddVariable(vars, "polarisation", &(psi), LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
    LALInferenceAddVariable(vars, "declination", &dec, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
    LALInferenceAddVariable(vars, "rightascension", &ra, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
    LALInferenceAddVariable(vars, "loghrss", &loghrss, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);


}
