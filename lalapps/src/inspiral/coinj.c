/**********************************************************************
*        Copyright (C) 2009 John Veitch, Stephen Fairhurst
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
**********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <config.h>
#include <math.h>
#include <getopt.h>
#include <string.h>

//#define LAL_USE_OLD_COMPLEX_STRUCTS
#include <lalapps.h>
#include <processtable.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LIGOLwXML.h>
#include <lal/LIGOLwXMLInspiralRead.h>
#include <lal/Date.h>
#include <lal/FindChirp.h>
#include <lal/Units.h>
#include <lal/FrequencySeries.h>
#include <lal/TimeSeries.h>
#include <lal/TimeFreqFFT.h>
#include <lal/InspiralInjectionParams.h>
#include <lal/VectorOps.h>
#include <lal/FrameStream.h>
#include <lal/LALDetectors.h>
#include <lal/LALFrameIO.h>
#include <lal/LALNoiseModels.h>
#include <gsl/gsl_rng.h>

#include <lal/LALInspiralStationaryPhaseApprox2Test.h>
#include <lal/LALInspiralMassiveGraviton.h>
#include <lal/LALInspiralBransDicke.h>
#include <lal/LALInspiralPPE.h>
#include <lal/GenerateInspiral.h>


#define PROGRAM_NAME "coinj"

#define USAGE \
"lalpps_coinj [options]\n\
--help                       display this message \n\
--input <injection.xml>      Specify input SimInspiralTable xml file\n\
--output-path OUTPUTPATH directory path where frames are going to be written to\n\
--skip-ascii-output           skip generation of  ASCII txt output file\n\
--response-type TYPE         TYPE of injection, [ strain | etmx | etmy ]\n\
--frames                     Create h(t) frame files\n\n\
[--AdvNoise\t\t Use AdvLIGO noise to estimate the SNR]\n\
[--maxSNR snrhigh --minSNR snrlow                     Adjust injections to have combined SNR between snrlow and snrhigh in the H1H2L1V1 network]\n\
[--SNR snr      adjust distance to get precisely this snr]\n\
[--GPSstart A --GPSend B     Only generate waveforms for injection between GPS seconds A and B (int)]\n\
[--max-chirp-dist DIST	     Set the maximum chirp distance in H, L or V to DIST]\n\
lalapps_coinj: create coherent injection files for LIGO and VIRGO\n"

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif
RCSID("$Id");

extern int vrbflg;
extern int lalDebugLevel;
int AdvNoises=0;
REAL8 srate=4096.0;
gsl_rng *RNG;
typedef enum
{
  noResponse,
  unityResponse,
  design,
  actuationX,
  actuationY
} ResponseType;

typedef struct actuationparameters
{
  REAL4        ETMXcal;
  REAL4        ETMYcal;
  REAL4        pendFX;
  REAL4        pendFY;
  REAL4        pendQX;
  REAL4        pendQY;
  REAL4        length;
} ActuationParametersType;

typedef void (NoiseFunc)(LALStatus *status,REAL8 *psd,REAL8 f);
void InjectFD_singleIFO(LALStatus status, COMPLEX16FrequencySeries *injFD_out,LALDetector *detector, SimInspiralTable *inj_table);

double chirpDist(SimInspiralTable *inj, char ifo);

double chirpDist(SimInspiralTable *inj, char ifo){
  /* eff_dist(IFO)*(2.8*0.25^(3/5)/mchirp)^(5/6) */
  double eff_dist=0.0;
  if(!inj) {printf("Null pointer passed to chirpDist!\n"); exit(1);}
  switch (ifo)
  {
	case 'H': {eff_dist=inj->eff_dist_h; break;}
	case 'L': {eff_dist=inj->eff_dist_l; break;}
	case 'G': {eff_dist=inj->eff_dist_g; break;}
	case 'T': {eff_dist=inj->eff_dist_t; break;}
	case 'V': {eff_dist=inj->eff_dist_v; break;}
        default: {printf("ERROR: Unknown IFO code %c\n",ifo);}
  }
  return( eff_dist*pow(1.218/inj->mchirp , (5./6.)));
}

int main(int argc, char *argv[])
{
  LALStatus        status=blank_status;
  CHAR                inputfile[FILENAME_MAX];
  CHAR                outputpath[1000];
  CHAR            adjustedfile[FILENAME_MAX+10];
  CHAR                injtype[30];
  CHAR                det_name[10];
  INT4                detectorFlags;
  LIGOTimeGPS inj_epoch;
  REAL8                deltaT= 1.0/srate;
  REAL8                injLength=210.0; /* Ten seconds at end */
    REAL8 deltaF=1.0/injLength;

  REAL8                LeadupTime=95.0;
  REAL8                dynRange=1.0/3.0e-23;
  INT4 seed=0;
  //UINT4 det_id=0;
  UINT4                Nsamples,det_idx,i,inj_num=0;
  /* Initialise the RNG */
  gsl_rng_env_setup();
  RNG=gsl_rng_alloc(gsl_rng_default);
  gsl_rng_set(RNG,seed==0 ? (unsigned long int)time(NULL) : (unsigned int long) seed);

  ActuationParametersType actuationParams[LAL_NUM_IFO];
  ActuationParametersType actData = {0,0,0,0,0,0,0};
  ResponseType injectionResponse=noResponse;
  FILE *                outfile;
  LIGOLwXMLStream                *xmlfp=NULL;
  CHAR                outfilename[FILENAME_MAX];
  CHAR                fname[FILENAME_MAX];

  const LALUnit strainPerCount={0,{0,0,0,0,0,1,-1},{0,0,0,0,0,0,0}};
  const LALUnit countPerStrain={0,{0,0,0,0,0,-1,1},{0,0,0,0,0,0,0}};
  NoiseFunc *PSD = NULL;
  REAL8 PSDscale=1.0;
  int c;
  int repeatLoop=0;
  int hitTarget=0;
  SimInspiralTable *injTable=NULL;
  SimInspiralTable this_injection;
  SimInspiralTable *headTable=NULL;
  MetadataTable MDT;
  REAL4TimeSeries *TimeSeries;
  double max_chirp_dist=0.0;

  REAL4TimeSeries *actuationTimeSeries;
  COMPLEX8FrequencySeries *resp;
  COMPLEX8FrequencySeries *actuationResp;
  COMPLEX8FrequencySeries *transfer;
  COMPLEX8Vector *unity;

  FrOutPar        UNUSED VirgoOutPars;
  CHAR                VirgoParsSource[100];
  CHAR                VirgoParsInfo[100];
  REAL8                NetworkSNR=0.0;
  INT4  makeFrames=0;
  INT4 outputRaw=0;
  INT4 skipASCIIoutput=1;
  COMPLEX8FrequencySeries *fftData;
  REAL8 mySNRsq,mySNR;
  REAL4FFTPlan *fwd_plan;
  REAL8 minSNR=0.0,maxSNR=0.0;
  REAL8 UNUSED maxRatio=1.0, UNUSED minRatio=1.0;
  REAL8 targetSNR=0.0;
  INT4 GPSstart=0,GPSend=2147483647;
  int UNUSED SNROK=1;
  int rewriteXML=0;
  FrameH *frame;
  LALDetector*	 detector;
  LALDetector *detector2=calloc(LAL_NUM_DETECTORS,sizeof(LALDetector));
	for(i=0;i<LAL_NUM_DETECTORS;i++) memcpy(&(detector2[i]),&lalCachedDetectors[i],sizeof(LALDetector));

  /*vrbflg=6;
    lalDebugLevel=6; */

  struct option long_options[]=
    {
      {"help", no_argument, 0, 'h'},
      {"input",required_argument,0, 'i'},
      {"output-path",required_argument,0,'o'},
      {"skip-ascii-output",no_argument,&skipASCIIoutput,'A'},
      {"response-type",required_argument,0,'r'},
      {"frames",no_argument,&makeFrames,'F'},
      {"rawstrain",no_argument,&outputRaw,'s'},
      {"verbose",no_argument,&vrbflg,1},
      {"minSNR",required_argument,0,2},
      {"maxSNR",required_argument,0,3},
      {"SNR",required_argument,0,6},
      {"GPSstart",required_argument,0,4},
      {"GPSend",required_argument,0,5},
      {"max-chirp-dist",required_argument,0,'d'},
      {"AdvNoise",no_argument,0,11},
      {0,0,0,0}
    };

  /*taken from Calibration CVS file: 
   * calibration/frequencydomain/runs/S5/H1/model/V3/H1DARMparams_849677446.m */
  actuationParams[LAL_IFO_H1].ETMXcal = -0.795e-9;
  actuationParams[LAL_IFO_H1].pendFX  = 0.767;
  actuationParams[LAL_IFO_H1].pendQX  = 10.0;
  actuationParams[LAL_IFO_H1].ETMYcal = -0.827e-9;
  actuationParams[LAL_IFO_H1].pendFY  = 0.761;
  actuationParams[LAL_IFO_H1].pendQY  = 10.0;
  actuationParams[LAL_IFO_H1].length  = 4000.0;

  /*taken from Calibration CVS file: 
   * calibration/frequencydomain/runs/S5/H2/model/V3/H2DARMparams_849678155.m */
  actuationParams[LAL_IFO_H2].ETMXcal = -0.876e-9;
  actuationParams[LAL_IFO_H2].pendFX  = 0.749;
  actuationParams[LAL_IFO_H2].pendQX  = 10.0;
  actuationParams[LAL_IFO_H2].ETMYcal = -0.912e-9;
  actuationParams[LAL_IFO_H2].pendFY  = 0.764;
  actuationParams[LAL_IFO_H2].pendQY  = 10.0;
  actuationParams[LAL_IFO_H2].length  = 2000.0;

  /*taken from Calibration CVS file: 
   * calibration/frequencydomain/runs/S5/L1/model/V3/L1DARMparams_841930071.m */
  actuationParams[LAL_IFO_L1].ETMXcal = -0.447e-9;
  actuationParams[LAL_IFO_L1].pendFX  = 0.766;
  actuationParams[LAL_IFO_L1].pendQX  = 100.0;
  actuationParams[LAL_IFO_L1].ETMYcal = -0.438e-9;
  actuationParams[LAL_IFO_L1].pendFY  = 0.756;
  actuationParams[LAL_IFO_L1].pendQY  = 100.0;
  actuationParams[LAL_IFO_L1].length  = 4000.0;


/*******************************************************************************/

/* Process input arguments */
  while(1)
    {
      int option_idx=0;
      c=getopt_long_only(argc,argv,"hFi:d:",long_options,&option_idx);
      if(c==-1) break;
      switch(c)
        {
        case 'h':
          fprintf(stderr,USAGE);
          exit(0);
          break;
        case 'i':
          strncpy(inputfile,optarg,FILENAME_MAX-1);
          break;
        case 'o':
          strncpy(outputpath,optarg,999);
          break;          
        case 'r':
          if(!strcmp("strain",optarg))          injectionResponse = unityResponse;
          else if(!strcmp("etmx",optarg))   injectionResponse = actuationX;
          else if(!strcmp("etmy",optarg))   injectionResponse = actuationY;
          else {fprintf(stderr,"Invalid argument to response-type: %s\nResponse type must be strain, etmy or etmx\n", \
                        optarg); exit(1);}
          break;
        case 2:
          minSNR=atof(optarg);
          fprintf(stderr,"Using minimum SNR of %f\n",minSNR);
          break;
        case 3:
          maxSNR=atof(optarg);
          fprintf(stderr,"Using maximum SNR of %f\n",maxSNR);
          break;
        case 4:
          GPSstart=atoi(optarg);
          break;
        case 5:
          GPSend=atoi(optarg);
          break;
        case 6:
          targetSNR=atof(optarg);
          fprintf(stderr,"Target SNR = %lf\n",targetSNR);
          break;
        case 'd':
	  max_chirp_dist=atof(optarg);
          fprintf(stderr,"Using maximum chirp distance of %lf\n",max_chirp_dist);
	  break;
        case 11:
          AdvNoises=1;
          fprintf(stdout,"Using Advanced LIGO noise\n");
          break;
        default:
	     fprintf(stdout,USAGE); exit(0);
	     break;	
       }
    }

  if(minSNR!=0 && maxSNR!=0 && (maxSNR<minSNR)){
    fprintf(stderr,"Error: minSNR must be less than maxSNR\n");
    exit(1);
    if(targetSNR!=0.0 && (targetSNR<minSNR || targetSNR>maxSNR)){
      fprintf(stderr,"Target SNR %lf is not in range %lf to %lf, ignoring min and max\n",targetSNR,minSNR,maxSNR);
    }
  }

  memset(&status,0,sizeof(status));

  /* Read in the input XML */
  SimInspiralTableFromLIGOLw(&injTable,inputfile,0,0);
  headTable=injTable;
  Nsamples = (UINT4)injLength/deltaT;
    REAL8 singleIFO_threshold=5.5;
    REAL8 * SNRs=NULL;
    SNRs=calloc(LAL_NUM_IFO+1 ,sizeof(REAL8));

  do{
    memcpy(&this_injection,injTable,sizeof(SimInspiralTable));
    this_injection.next=NULL;
    NetworkSNR=0.0;

    /* Set epoch */
    memcpy(&inj_epoch,&(this_injection.geocent_end_time),sizeof(LIGOTimeGPS));
    inj_epoch = this_injection.geocent_end_time;
    XLALGPSAdd(&inj_epoch, -LeadupTime);
    inj_epoch.gpsNanoSeconds=0;
    SNROK=0; /* Reset this to 0 = OK */
    minRatio=2.0;
    maxRatio=0.0;
    repeatLoop=0;
    /* Loop over detectors */
    for(det_idx=0;det_idx<LAL_NUM_IFO;det_idx++){
      /* Only generate within chosen bounds, if specified */
      if((this_injection.geocent_end_time.gpsSeconds-(int)LeadupTime )<GPSstart || (this_injection.geocent_end_time.gpsSeconds-(int)LeadupTime)>GPSend) continue;

      if(det_idx==LAL_IFO_T1||det_idx==LAL_IFO_G1||det_idx==LAL_IFO_H2) continue; /* Don't generate for GEO or TAMA */
      if (AdvNoises==0){
      switch(det_idx)
        {
        case LAL_IFO_H1: sprintf(det_name,"H1"); PSD=&LALLIGOIPsd; PSDscale=9E-46; detectorFlags = LAL_LHO_4K_DETECTOR_BIT; break;
        case LAL_IFO_H2: sprintf(det_name,"H2"); PSD=&LALLIGOIPsd; PSDscale=9E-46; detectorFlags = LAL_LHO_2K_DETECTOR_BIT; break;
        case LAL_IFO_L1: sprintf(det_name,"L1"); PSD=&LALLIGOIPsd; PSDscale=9E-46; detectorFlags = LAL_LLO_4K_DETECTOR_BIT; break; 
        case LAL_IFO_V1: sprintf(det_name,"V1"); PSD=&LALVIRGOPsd; PSDscale=1.0; detectorFlags = LAL_VIRGO_DETECTOR_BIT;  break;
        case LAL_IFO_G1: sprintf(det_name,"G1"); PSD=&LALGEOPsd; PSDscale=1E-46; detectorFlags = LAL_GEO_600_DETECTOR_BIT;  break;
        case LAL_IFO_T1: sprintf(det_name,"T1"); PSD=&LALTAMAPsd; PSDscale=75E-46; detectorFlags = LAL_TAMA_300_DETECTOR_BIT; break;
        default: fprintf(stderr,"Unknown IFO\n"); exit(1); break;
        }
      }
      else{
      /*Uses ALIGO noise */
      switch(det_idx)
        {
        case LAL_IFO_H1: sprintf(det_name,"H1"); PSD=&LALAdvLIGOPsd; PSDscale=1.35e-50; detectorFlags = LAL_LHO_4K_DETECTOR_BIT; break;
        case LAL_IFO_L1: sprintf(det_name,"L1"); PSD=&LALAdvLIGOPsd; PSDscale=1.35e-50; detectorFlags = LAL_LLO_4K_DETECTOR_BIT; break; 
        case LAL_IFO_V1: sprintf(det_name,"V1"); PSD=&LALAdvVIRGOPsd; PSDscale=1E-47; detectorFlags = LAL_VIRGO_DETECTOR_BIT;  break;
        default: fprintf(stderr,"Unknown IFO\n"); exit(1); break;
        }
      }
      TimeSeries=XLALCreateREAL4TimeSeries(det_name,&inj_epoch,0.0,deltaT,&lalADCCountUnit,(size_t)Nsamples);

      if (strstr(this_injection.waveform,"TaylorT2")){
        printf("I'm in TaylorT2 WF\n");
      for(i=0;i<Nsamples;i++) TimeSeries->data->data[i]=0.0;
      resp = XLALCreateCOMPLEX8FrequencySeries("response",&inj_epoch,0.0,1.0/injLength,&strainPerCount,(size_t)Nsamples/2+1);
      for(i=0;i<resp->data->length;i++) {resp->data->data[i].re=(REAL4)1.0/dynRange; resp->data->data[i].im=0.0;}

      /* Create h(t) time series for this detector */
      LAL_CALL( LALFindChirpInjectSignals(&status,TimeSeries,&this_injection,resp) , &status);

      XLALDestroyCOMPLEX8FrequencySeries(resp);

      if(det_idx==LAL_IFO_V1 && injectionResponse!=unityResponse){
        actuationTimeSeries=NULL;
        goto calcSNRandwriteFrames;
      }

      /* -=-=-=-=-=-=- Prepare actuations -=-=-=-=-=-=- */

      if(injectionResponse==actuationX || injectionResponse==actuationY) actData=actuationParams[det_idx];
      actuationResp = XLALCreateCOMPLEX8FrequencySeries("actuationResponse",&inj_epoch,0.0,1.0/(2.0*injLength),&strainPerCount,(size_t)Nsamples/2+1);
      /* Create actuation response */
      switch(injectionResponse){
      case unityResponse:
        sprintf(injtype,"STRAIN");
        for(i=0;i<actuationResp->data->length;i++){actuationResp->data->data[i].re=1.0; actuationResp->data->data[i].im=0.0;}
        break;
      case actuationX:
        sprintf(injtype,"ETMX");
        actuationResp=generateActuation(actuationResp,actData.ETMXcal/actData.length,actData.pendFX,actData.pendQX);
        break;
      case actuationY:
        sprintf(injtype,"ETMY");
        actuationResp=generateActuation(actuationResp,actData.ETMYcal/actData.length,actData.pendFY,actData.pendQY);
        break;
      default:
        fprintf(stderr,"Must specify response function: strain, etmy or etmx\n"); exit(1);
        break;
      }



      if(injectionResponse!=unityResponse) {
        actuationTimeSeries=XLALCreateREAL4TimeSeries(det_name,&inj_epoch,0.0,deltaT,&lalADCCountUnit,(size_t)Nsamples);
        unity = XLALCreateCOMPLEX8Vector(actuationResp->data->length);
        transfer = XLALCreateCOMPLEX8FrequencySeries("transfer",&inj_epoch,0.0,1.0/(2.0*injLength),&countPerStrain,(size_t)Nsamples/2+1);
        for(i=0;i<unity->length;i++) {unity->data[i].re=1.0; unity->data[i].im=0.0;}
        XLALCCVectorDivide(transfer->data,unity,actuationResp->data);
        for(i=0;i<Nsamples;i++) actuationTimeSeries->data->data[i]=TimeSeries->data->data[i];
        actuationTimeSeries = XLALRespFilt(actuationTimeSeries,transfer);
        XLALDestroyCOMPLEX8FrequencySeries(transfer);
        XLALDestroyCOMPLEX8Vector(unity);
        for(i=0;i<actuationTimeSeries->data->length;i++) actuationTimeSeries->data->data[i]/=dynRange;
      }
      else actuationTimeSeries=TimeSeries;

      XLALDestroyCOMPLEX8FrequencySeries(actuationResp);

      /* Output the actuation time series */
      /*sprintf(outfilename,"HWINJ_%i_%s_%i_%s.out",inj_num,injtype,inj_epoch.gpsSeconds,det_name);*/
      char massBin[5];
      if(this_injection.mass1<2.0 && this_injection.mass2<2.0) sprintf(massBin,"BNS");
      else if(this_injection.mass1<2.0 || this_injection.mass2<2.0) sprintf(massBin,"NSBH");
      else sprintf(massBin,"BBH");
      
      if (!(skipASCIIoutput)){
        sprintf(outfilename,"%s%s%i_CBC_%s_%i_%s_%s.txt",outputpath,"/",inj_epoch.gpsSeconds,massBin,inj_num,injtype,det_name);
        outfile=fopen(outfilename,"w");
        fprintf(stdout,"Injected signal %i for %s into file %s\n",inj_num,det_name,outfilename);
        for(i=0;i<actuationTimeSeries->data->length;i++) fprintf(outfile,"%10.10e\n",actuationTimeSeries->data->data[i]);
        fclose(outfile);
      };

    calcSNRandwriteFrames:

      /* Calculate SNR for this injection */
      fwd_plan = XLALCreateForwardREAL4FFTPlan( TimeSeries->data->length, 0 );
      fftData = XLALCreateCOMPLEX8FrequencySeries(TimeSeries->name,&(TimeSeries->epoch),0,1.0/TimeSeries->deltaT,&lalDimensionlessUnit,TimeSeries->data->length/2 +1);
      XLALREAL4TimeFreqFFT(fftData,TimeSeries,fwd_plan);
      XLALDestroyREAL4FFTPlan(fwd_plan);
      
        mySNRsq = 0.0;
      mySNR=0.0;
      for(i=1;i<fftData->data->length;i++){
        REAL8 freq;
        REAL8 sim_psd_value=0;
        freq = fftData->deltaF * i;
        PSD( &status, &sim_psd_value, freq );
        mySNRsq += fftData->data->data[i].re * fftData->data->data[i].re /
          (sim_psd_value*PSDscale);
        mySNRsq += fftData->data->data[i].im * fftData->data->data[i].im /
          (sim_psd_value*PSDscale);
      }
      mySNRsq *= 4.0*fftData->deltaF;
      XLALDestroyCOMPLEX8FrequencySeries( fftData );
      if(det_idx==LAL_IFO_H2) mySNRsq/=4.0;
      mySNR = sqrt(mySNRsq)/dynRange;
      SNRs[det_idx]=mySNR;
      fprintf(stdout,"SNR in design %s of injection %i = %lf\n",det_name,inj_num,mySNR);

      for(i=0;i<TimeSeries->data->length;i++) {
        TimeSeries->data->data[i]=TimeSeries->data->data[i]/dynRange +0.0;
      }
      NetworkSNR+=mySNR*mySNR;

      if(makeFrames){ /* Also output frames for Virgo */
        sprintf(VirgoParsSource,"%s-INSP%i",det_name,inj_num);
        VirgoOutPars.source=VirgoParsSource;
        sprintf(VirgoParsInfo,"HWINJ-STRAIN");
        VirgoOutPars.description=VirgoParsInfo;
        VirgoOutPars.type=ProcDataChannel;
        VirgoOutPars.nframes=(UINT4)injLength;
        VirgoOutPars.frame=0;
        VirgoOutPars.run=2;
        fprintf(stdout,"Generating frame file for %s-%s-%i\n",VirgoParsSource,VirgoParsInfo,TimeSeries->epoch.gpsSeconds);
       /* LALFrWriteREAL4TimeSeries(&status,TimeSeries,&VirgoOutPars); */
        /* define frame */
        frame = XLALFrameNew( &TimeSeries->epoch, injLength, "HWINJ-STRAIN", 0, 1, detectorFlags );
        /* add time series as a channel to the frame */
        XLALFrameAddREAL4TimeSeriesSimData( frame, TimeSeries );
        /* write frame */
        sprintf(fname,"%s%s%s-INSP%i_HWINJ_STRAIN-%i-%i.gwf",outputpath,"/",det_name, inj_num, inj_epoch.gpsSeconds, (UINT4)injLength);
        /*sprintf(fname, "%s%s%s",outputpath, "/", fname);*/
        if (XLALFrameWrite( frame, fname, 8) != 0)
        {
          fprintf( stderr, "ERROR: Cannot save frame file: '%s'\n", fname );
          exit( 1 );
        }

        /* clear frame */
        FrameFree( frame );

        
      }
        
      if(TimeSeries==actuationTimeSeries) XLALDestroyREAL4TimeSeries(TimeSeries);
      else {
        if(injectionResponse) XLALDestroyREAL4TimeSeries(actuationTimeSeries);
        XLALDestroyREAL4TimeSeries(TimeSeries);
      }
} //end if TaylorT2
else if (strstr(this_injection.waveform,"TaylorF2")){
          printf("I'm in TaylorF2 WF deltaF %lf \n",(REAL8) 1.0/injLength);

      COMPLEX16FrequencySeries *FData=(COMPLEX16FrequencySeries *)XLALCreateCOMPLEX16FrequencySeries("InjFD",  &inj_epoch,0.0,(REAL8)1.0/injLength,&lalDimensionlessUnit,(Nsamples-1)*2);
      
       for (UINT4 Midx=0;Midx<FData->data->length;Midx++){
         FData->data->data[Midx].re=0.0;
            FData->data->data[Midx].im=0.0;
         }
         
if(!strcmp(det_name,"H1")) detector=&detector2[LALDetectorIndexLHODIFF];
else if (!strcmp(det_name,"L1")) detector=&detector2[LALDetectorIndexLLODIFF];
else if (!strcmp(det_name,"V1")) detector=&detector2[LALDetectorIndexVIRGODIFF];

InjectFD_singleIFO(status, FData,detector,injTable);
 mySNRsq = 0.0;
      mySNR=0.0;
      UINT4 lowBin = (UINT4)(injTable->f_lower/deltaF);
      for(i=lowBin;i<FData->data->length-1;i++){
        REAL8 freq;
        REAL8 sim_psd_value=0;
        freq = FData->deltaF * i;
        PSD( &status, &sim_psd_value, freq );
        mySNRsq += FData->data->data[i].re * FData->data->data[i].re /
          (sim_psd_value*PSDscale);
        mySNRsq += FData->data->data[i].im * FData->data->data[i].im /
          (sim_psd_value*PSDscale);
          //printf("F %lf Re %10.10e IM %10.10e\n",freq,FData->data->data[i].re,FData->data->data[i].im);
      }
      mySNRsq *= 4.0*FData->deltaF;
      XLALDestroyCOMPLEX16FrequencySeries( FData );
      
      mySNR = sqrt(mySNRsq);
      SNRs[det_idx]=mySNR;
      fprintf(stdout,"SNR in design %s of injection %i = %lf\n",det_name,inj_num,mySNR);
      NetworkSNR+=mySNR*mySNR;


} //end if TF2
    
} /* End loop over detectors */

    /*fprintf(stdout,"Finished injecting signal %i, network SNR %f\n",inj_num,sqrt(NetworkSNR));*/
    NetworkSNR=sqrt(NetworkSNR);
    if(NetworkSNR!=0.0){ /* Check the SNR if we did this injection */
      if(targetSNR!=0.0 && repeatLoop==0 && hitTarget==0) {
        injTable->distance*=(REAL4)(NetworkSNR/targetSNR);
        injTable->eff_dist_h*=(REAL4)(NetworkSNR/targetSNR);
        injTable->eff_dist_l*=(REAL4)(NetworkSNR/targetSNR);
        injTable->eff_dist_v*=(REAL4)(NetworkSNR/targetSNR);
        injTable->eff_dist_g*=(REAL4)(NetworkSNR/targetSNR);
        injTable->eff_dist_t*=(REAL4)(NetworkSNR/targetSNR);
        rewriteXML=1; repeatLoop=1; hitTarget=1;}
      else {repeatLoop=0; hitTarget=0;}
      if(targetSNR==0.0 && (minSNR>NetworkSNR || maxSNR<NetworkSNR)) {
      REAL8 randomSNR;
      INT4 above_threshold;
      REAL8 local_min=minSNR;
      do{
      // We want the distances to be uniform in the volume i.e. p(D)\propto D^2. Changing variable, we need to have p(SNR)\propto 1/SNR^4.
      // This implies that p(1/SNR^3) is uniform. We can extract from this uniform distribution and then find SNR taking the inverse of the cube root.

      randomSNR=1.0/(maxSNR*maxSNR*maxSNR)+(1.0/(local_min*local_min*local_min)- 1.0/(maxSNR*maxSNR*maxSNR))*gsl_rng_uniform(RNG);
      randomSNR=1.0/cbrt(randomSNR);
      printf("Setting random SNR to %lf\n",randomSNR); 
      local_min=randomSNR;
      above_threshold=0;
      for(i=0;i<LAL_NUM_IFO;i++)
      {
      if(SNRs[i]*randomSNR/NetworkSNR>singleIFO_threshold) above_threshold+=1;
      }
      if ((maxSNR-local_min)<0.001) goto here;
      }while(above_threshold<2);       
      printf("DONE!\n");
      here:
 injTable->distance*=(REAL4)(1.*NetworkSNR/randomSNR);
        injTable->eff_dist_h*=(REAL4)(1.*NetworkSNR/randomSNR);
        injTable->eff_dist_l*=(REAL4)(1.*NetworkSNR/randomSNR);
        injTable->eff_dist_v*=(REAL4)(1.*NetworkSNR/randomSNR);
        injTable->eff_dist_t*=(REAL4)(1.*NetworkSNR/randomSNR);
        injTable->eff_dist_g*=(REAL4)(1.*NetworkSNR/randomSNR);
        rewriteXML=1; repeatLoop=0; //SALVO
fprintf(stderr,"Multiplying distance by %lf to get from Distance %lf to target %lf \n",0.99*(NetworkSNR/randomSNR),injTable->distance*(randomSNR/NetworkSNR),injTable->distance);
fprintf(stderr,"Multiplying SNR by %lf to get from SNR %lf to %lf \n",1/(0.99*(NetworkSNR/randomSNR)),NetworkSNR,randomSNR);

}
     
    }
    if(max_chirp_dist!=0.0 && (maxSNR==0 || (maxSNR!=0 && (maxSNR>NetworkSNR) ) )  ){
	double this_max=0.0;
	double second_max=0.0;
	char ifostr[]="HLV";
	for(int cidx=0;cidx<3;cidx++)
	{
		if(this_max<chirpDist(injTable,ifostr[cidx])) this_max=chirpDist(injTable,ifostr[cidx]);
	}
	for(int cidx=0;cidx<3;cidx++)
		if(second_max<chirpDist(injTable,ifostr[cidx]) && chirpDist(injTable,ifostr[cidx])<this_max) second_max=chirpDist(injTable,ifostr[cidx]);
	if(second_max>max_chirp_dist){
	injTable->distance*=0.95*(max_chirp_dist/second_max);
	injTable->eff_dist_h*=0.95*max_chirp_dist/second_max;
	injTable->eff_dist_l*=0.95*max_chirp_dist/second_max;
	injTable->eff_dist_v*=0.95*max_chirp_dist/second_max;
	injTable->eff_dist_t*=0.95*max_chirp_dist/second_max;
	injTable->eff_dist_g*=0.95*max_chirp_dist/second_max;
	rewriteXML=1; repeatLoop=1;
	fprintf(stderr,"MCD: Multiplying distance by %lf to get from %lf to target\n",0.95*max_chirp_dist/second_max,second_max);
    	}
    }
    if(repeatLoop==1) fprintf(stderr,"Reinjecting with new distance %f for desired SNR\n\n",injTable->distance);

    if(repeatLoop==0){
      fprintf(stderr,"\nNetwork SNR of %i = %lf\n",inj_num,NetworkSNR);
      injTable=injTable->next;
      inj_num++;
    }
 

  }while(injTable!=NULL);

  /* If the distances were adjusted, re-write the SimInspiral table */
  if(rewriteXML){
    memset(&MDT,0,sizeof(MDT));
    MDT.simInspiralTable = headTable;
        
    fprintf(stdout,"Overwriting %s with adjusted distances\n",inputfile);
    strncat(adjustedfile,inputfile,strlen(inputfile)-4); /* Cut off the .xml */
    sprintf(inputfile,"%s_adj.xml",adjustedfile);
    xmlfp=XLALOpenLIGOLwXMLFile((const char *)inputfile);
    if(xmlfp==NULL) fprintf(stderr,"Error! Cannot open %s for writing\n",inputfile);
    LAL_CALL( LALBeginLIGOLwXMLTable( &status, xmlfp, sim_inspiral_table ), &status );
    LAL_CALL( LALWriteLIGOLwXMLTable( &status, xmlfp, MDT,sim_inspiral_table ), &status );
    LAL_CALL( LALEndLIGOLwXMLTable ( &status, xmlfp ), &status );
    LAL_CALL( LALCloseLIGOLwXMLFile ( &status, xmlfp ), &status );
  }

  return(0);
}

///*-----------------------------------------------------------*/
void InjectFD_singleIFO(LALStatus status, COMPLEX16FrequencySeries *injFD_out,LALDetector *detector, SimInspiralTable *inj_table)
///*-------------- Inject in Frequency domain -----------------*/
{
	/* Inject a gravitational wave into the data in the frequency domain */
	REAL4Vector *injWaveFD=NULL;
    InspiralTemplate template;
	UINT4 idx;
	REAL8 end_time = 0.0;
	//REAL8 deltaF = inputMCMC->deltaF;
	REAL8 TimeFromGC,resp_r,resp_i;
	UINT4 Nmodel; /* Length of the model */
	LALDetAMResponse det_resp;
    memset(&template,0,sizeof(InspiralTemplate));
    REAL8 dphis[10]={0.0};
    /* Populate the template */
	  REAL8 ChirpISCOLength;
  	expnFunc expnFunction;
	  expnCoeffs ak;
    TofVIn TofVparams;
    //REAL8 * SNRs=NULL;
    //SNRs=calloc(nIFO+1 ,sizeof(REAL8));
    REAL8 fLow=(REAL8) inj_table->f_lower;
    REAL8 injLength=210.0;
    REAL8 deltaT= 1.0/srate;
    REAL8 deltaF=1.0/injLength; 
    UINT4 Nsamples = (UINT4)injLength/deltaT;
    LIGOTimeGPS epoch;
    Nmodel=(UINT4) (Nsamples-1)*2;
    REAL8  LeadupTime=95.0;
    memcpy(&epoch,&(inj_table->geocent_end_time),sizeof(LIGOTimeGPS));
    epoch = inj_table->geocent_end_time;
    XLALGPSAdd(&epoch, -LeadupTime);
    epoch.gpsNanoSeconds=0;
printf("Nmodel %d Nsamples %d \n",Nmodel,Nsamples);
    /* read in the injection approximant and determine whether is TaylorF2 or something else*/
    Approximant injapprox;
    LALPNOrder phase_order;
    LALGetApproximantFromString(&status,inj_table->waveform,&injapprox);
    LALGetOrderFromString(&status,inj_table->waveform,&phase_order);
	template.totalMass = inj_table->mass1+inj_table->mass2;
	template.eta = inj_table->eta;
	template.massChoice = totalMassAndEta;
	template.fLower = inj_table->f_lower;
    template.distance = inj_table->distance; /* This must be in Mpc, contrary to the docs */
	template.order=phase_order;
	template.approximant=injapprox;
	template.tSampling = 1.0/deltaT;
	template.fCutoff = 0.5/deltaT -1.0;
    //fprintf(stdout,"%f \n", template.fCutoff);
	template.nStartPad = 0;
	template.nEndPad =0;
    template.startPhase = inj_table->coa_phase;
	template.startTime = 0.0;
	template.ieta = 1;
	template.next = NULL;
	template.fine = NULL;
	//Nmodel = (inputMCMC->stilde[0]->data->length-1)*2; /* *2 for real/imag packing format */

    COMPLEX16FrequencySeries *injF=(COMPLEX16FrequencySeries *)XLALCreateCOMPLEX16FrequencySeries("InjF",  &epoch,0.0,deltaF,&lalDimensionlessUnit,Nmodel);
   
	if(injWaveFD==NULL)	LALCreateVector(&status,&injWaveFD,Nmodel); /* Allocate storage for the waveform */
        
	/* Create the wave */
	LALInspiralParameterCalc(&status,&template);
	LALInspiralRestrictedAmplitude(&status,&template);
    
    printf("Injection Approx: %i\n",template.approximant);
    if (template.approximant==IMRPhenomFBTest) {
        dphis[0]=inj_table->dphi0;
        dphis[1]=inj_table->dphi1;
        dphis[2]=inj_table->dphi2;
        dphis[3]=inj_table->dphi3;
        dphis[4]=inj_table->dphi4;
        dphis[5]=inj_table->dphi5;
        dphis[6]=inj_table->dphi6;
        dphis[7]=inj_table->dphi7;
        dphis[8]=inj_table->dphi8;
        dphis[9]=inj_table->dphi9;
        printf("Using approximant IMRPhenomFBTest\n");
        for (int k=0;k<10;k++) {
			fprintf(stderr,"Injecting dphi%i = %e\n",k,dphis[k]);
		}
		if (dphis[9]!=0.) {
			fprintf(stderr,"Coefficient psi_9 is not available in IMRPhenomB. Value is set to 0.");
			dphis[9]=0.;
		}
        LALBBHPhenWaveFreqDomTest(&status, injWaveFD, &template, dphis, 0.0);
    }
    else if (template.approximant==TaylorF2Test){
        dphis[0]=inj_table->dphi0;
        dphis[1]=inj_table->dphi1;
        dphis[2]=inj_table->dphi2;
        dphis[3]=inj_table->dphi3;
        dphis[4]=inj_table->dphi4;
        dphis[5]=inj_table->dphi5;
        dphis[6]=inj_table->dphi5l;
        dphis[7]=inj_table->dphi6;
        dphis[8]=inj_table->dphi6l;
        dphis[9]=inj_table->dphi7;
        template.spin1[0]=inj_table->spin1x;
        template.spin1[1]=inj_table->spin1y;
        template.spin1[2]=inj_table->spin1z;
        template.spin2[0]=inj_table->spin2x;
        template.spin2[1]=inj_table->spin2y;
        template.spin2[2]=inj_table->spin2z;
        for (int k=0;k<10;k++) fprintf(stderr,"Injecting dphi%i = %e\n",k,dphis[k]);
        fprintf(stderr, "Injecting spin1 : (%e, %e, %e)\n", template.spin1[0], template.spin1[1], template.spin1[2]);
        fprintf(stderr, "Injecting spin2 : (%e, %e, %e)\n", template.spin2[0], template.spin2[1], template.spin2[2]);
        LALInspiralStationaryPhaseApprox2Test(&status, injWaveFD, &template, dphis, 0.0);
    }
    else if (template.approximant==MassiveGraviton) {
		fprintf(stderr,"Injecting logLambdaG = %e\n",inj_table->loglambdaG);
        template.loglambdaG=inj_table->loglambdaG;
		LALInspiralMassiveGraviton(&status, injWaveFD, &template);
	} 
    else if (template.approximant==PPE) {
		fprintf(stderr,"Injecting aPPE = %e\n",inj_table->aPPE);
        fprintf(stderr,"Injecting alphaPPE = %e\n",inj_table->alphaPPE);
        fprintf(stderr,"Injecting bPPE = %e\n",inj_table->bPPE);
        fprintf(stderr,"Injecting betaPPE = %e\n",inj_table->betaPPE);
        template.aPPE=inj_table->aPPE;
        template.alphaPPE=inj_table->alphaPPE;
        template.bPPE=inj_table->bPPE;
        template.betaPPE=inj_table->betaPPE;
		LALInspiralPPE(&status, injWaveFD, &template, 0.0);    
    }
    else if (template.approximant==BransDicke) {
		fprintf(stderr,"Injecting Scalar Charge 1 = %e\n",inj_table->ScalarCharge1);
        fprintf(stderr,"Injecting Scalar Charge 2 = %e\n",inj_table->ScalarCharge2);
        fprintf(stderr,"Injecting OmegaBD = %e\n",inj_table->omegaBD);
        template.ScalarCharge1=inj_table->ScalarCharge1;
        template.ScalarCharge2=inj_table->ScalarCharge2;
        template.omegaBD=inj_table->omegaBD;
		LALInspiralBransDicke(&status, injWaveFD, &template);    
    }
 /*   else if (template.approximant==IMRPhenomFB) {
		fprintf(stderr,"GR injection");
        LALBBHPhenWaveFreqDom(&status, injWaveFD, &template);
    } */
    else {
		fprintf(stderr,"GR injection\n");
        LALInspiralWave(&status,injWaveFD,&template);
    }
    
    memset(&ak,0,sizeof(expnCoeffs));
	memset(&TofVparams,0,sizeof(TofVparams));

    LALInspiralSetup(&status,&ak,&template);
	LALInspiralChooseModel(&status,&expnFunction,&ak,&template);
	TofVparams.coeffs=&ak;
	TofVparams.dEnergy=expnFunction.dEnergy;
	TofVparams.flux=expnFunction.flux;
	TofVparams.v0= ak.v0;
	TofVparams.t0= ak.t0;
	TofVparams.vlso= ak.vlso;
	TofVparams.totalmass=ak.totalmass;
/*	LALInspiralTofV(&status,&ChirpISCOLength,pow(6.0,-0.5),(void *)&TofVparams);*/
	ChirpISCOLength=ak.tn;
   
    if(template.approximant == IMRPhenomB || template.approximant==IMRPhenomFB || template.approximant==IMRPhenomFBTest){
		ChirpISCOLength = template.tC;
	}

    end_time = (REAL8) inj_table->geocent_end_time.gpsSeconds + (REAL8) inj_table->geocent_end_time.gpsNanoSeconds*1.0e-9;
    end_time-=(REAL8) epoch.gpsSeconds + 1.0e-9*epoch.gpsNanoSeconds;
    
    if(!(template.approximant == IMRPhenomB || template.approximant==IMRPhenomFB || template.approximant==IMRPhenomFBTest)){
       /* IMR FD is created with the end of the waveform at the time of epoch. So we don't need to shift by the length of the WF. */
        end_time-=ChirpISCOLength;
	}
    
	/* Calculate response of the detectors */
	LALSource source;
	memset(&source,0,sizeof(LALSource));
	source.equatorialCoords.longitude = (REAL8) inj_table->longitude;
	source.equatorialCoords.latitude = (REAL8) inj_table->latitude;
	source.equatorialCoords.system = COORDINATESYSTEM_EQUATORIAL;
	source.orientation = (REAL8) inj_table->polarization;
	strncpy(source.name,"blah",sizeof(source.name));

	LALDetAndSource det_source;
	det_source.pSource = &source;

	REAL8 ci = cos((REAL8) inj_table->inclination);
//	REAL8 SNRinj=0;

	REAL8 time_sin,time_cos;
    //inputMCMC->numberDataStreams=nIFO;
   
	//for (det_i=0;det_i<nIFO;det_i++){ //nIFO
        UINT4 lowBin = (UINT4)(fLow /deltaF);
        UINT4 highBin = (UINT4)(template.fFinal / deltaF);
        
        if(highBin==0 || highBin>Nmodel-1) highBin=Nmodel-1;
		
        if(template.approximant==IMRPhenomFB || template.approximant==IMRPhenomB || template.approximant==IMRPhenomFBTest || template.approximant==EOBNR) highBin=Nmodel-1;
		//char InjFileName[50];
		//sprintf(InjFileName,"injection_%i.dat",det_i);
		//FILE *outInj=fopen(InjFileName,"w");

		/* Compute detector amplitude response */
		det_source.pDetector = (detector); /* select detector */
		LALComputeDetAMResponse(&status,&det_resp,&det_source,&epoch); /* Compute det_resp */
        /* Time delay from geocentre */
        TimeFromGC = (REAL8) XLALTimeDelayFromEarthCenter(detector->location, source.equatorialCoords.longitude, source.equatorialCoords.latitude, &(epoch));
		det_resp.plus*=0.5*(1.0+ci*ci);
		det_resp.cross*=-ci;
		//REAL8 chisq=0.0;
        
        
        for(idx=lowBin;idx<=highBin;idx++){
        /* Calculate the WF to be injected for each frequency bin and fill injF, injFnoError and injFwithError with it. 
         * Nothing is yet added to the data */
			time_sin = sin(LAL_TWOPI*(end_time+TimeFromGC)*((REAL8) idx)*(deltaF));
			time_cos = cos(LAL_TWOPI*(end_time+TimeFromGC)*((REAL8) idx)*(deltaF));
			REAL8 hc = (REAL8)injWaveFD->data[idx]*time_cos + (REAL8)injWaveFD->data[Nmodel-idx]*time_sin;
			REAL8 hs = (REAL8)injWaveFD->data[Nmodel-idx]*time_cos - (REAL8)injWaveFD->data[idx]*time_sin;
      //printf("idx %d injWaveRE %10.10e injWaveIM %10.10e time sin %lf time cos %lf\n",idx,injWaveFD->data[idx],injWaveFD->data[Nmodel-idx],time_sin,time_cos);
			resp_r = det_resp.plus * hc - det_resp.cross * hs;
			resp_i = det_resp.cross * hc + det_resp.plus * hs;
			resp_r/=deltaF; resp_i/=deltaF;
            injF->data->data[idx].re=resp_r;
            injF->data->data[idx].im=resp_i;
        }
        
                
        for(idx=lowBin;idx<=highBin;idx++){
		 /* Copy the WF. */  
            injFD_out->data->data[idx].re=injF->data->data[idx].re;
            injFD_out->data->data[idx].im=injF->data->data[idx].im;

		}
   
    
	return;
}
