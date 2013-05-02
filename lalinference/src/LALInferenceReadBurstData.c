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
#include <lal/LALInferenceReadBurstData.h>
#include <lal/LALSimNoise.h>

#define LALINFERENCE_DEFAULT_FLOW "40.0"
//typedef void (NoiseFunc)(LALStatus *statusPtr,REAL8 *psd,REAL8 f);
static void PrintBurstSNRsToFile(LALInferenceIFOData *IFOdata ,SimBurst *inj_table);
void InjectSineGaussianFD(LALInferenceIFOData *IFOdata, SimBurst *inj_table, ProcessParamsTable *commandLine);
char *BurstSNRpath = NULL;

static const LALUnit strainPerCount={0,{0,0,0,0,0,1,-1},{0,0,0,0,0,0,0}};
 struct fvec {
	REAL8 f;
	REAL8 x;
};

typedef void (NoiseFunc)(LALStatus *statusPtr,REAL8 *psd,REAL8 f);

void LALInferenceInjectBurstSignal(LALInferenceRunState *irs, ProcessParamsTable *commandLine)
{
    lalDebugLevel=1;
	LALStatus status;
	memset(&status,0,sizeof(status));
	SimBurst *injTable=NULL;
    SimBurst *injEvent=NULL;
	INT4 Ninj=0;
	INT4 event=0;
	UINT4 i=0,j=0;
	int si=0;
    //nchar *BurstSNRpath = NULL;
    //REAL8 responseScale=1.0;
	//CoherentGW InjectGW;
	//PPNParamStruc InjParams;
	LIGOTimeGPS injstart;
	REAL8 SNR=0.0,NetworkSNR=0.0;// previous_snr=0.0; 
	//UINT4 best_ifo_snr=0;//best_ifo_snr+=1;
	//UINT4 highest_snr_index=0;
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
	else
	fprintf(stdout,"WARNING: you did not give --event. Injecting event 0 of the xml table, which may not be what you want!\n");
        if(LALInferenceGetProcParamVal(commandLine,"--snrpath")){
                ppt = LALInferenceGetProcParamVal(commandLine,"--snrpath");
		BurstSNRpath = calloc(strlen(ppt->value)+1,sizeof(char));
		memcpy(BurstSNRpath,ppt->value,strlen(ppt->value)+1);
                fprintf(stdout,"Writing SNRs in %s\n",BurstSNRpath)     ;

	}
	injTable=XLALSimBurstTableFromLIGOLw(LALInferenceGetProcParamVal(commandLine,"--inj")->value,0,0);
	REPORTSTATUS(&status);
    Ninj=-1;
    while(injTable){Ninj++;injTable=injTable->next;}
	if(Ninj < event){ 
	    fprintf(stderr,"Error reading event %d from %s\n",event,LALInferenceGetProcParamVal(commandLine,"--inj")->value);
	    exit(1);
	    }
	injTable=XLALSimBurstTableFromLIGOLw(LALInferenceGetProcParamVal(commandLine,"--inj")->value,0,0);
	while(si<event) {si++; injTable = injTable->next;} /* Select event */
	injEvent = injTable;
	injEvent->next = NULL;
	
	//memset(&InjectGW,0,sizeof(InjectGW));
	
	REPORTSTATUS(&status);
	//LALGenerateInspiral(&status,&InjectGW,injTable,&InjParams);
	//if(status.statusCode!=0) {fprintf(stderr,"Error generating injection!\n"); REPORTSTATUS(&status); }
    
    /* Inject burst in the FreqDomain */
    int FDinj=0;
    if (injTable){
        if(!strcmp("SineGaussianF",injEvent->waveform)) FDinj=1;
        }
    
	if (LALInferenceGetProcParamVal(commandLine,"--FDinjections") || FDinj==1)
    {
         InjectSineGaussianFD(thisData, injEvent, commandLine);
         return;
    }
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
      printf("--------------------------Using LALSimulation for injection\n");
      REAL8TimeSeries *hplus=NULL;  /**< +-polarization waveform */
      REAL8TimeSeries *hcross=NULL; /**< x-polarization waveform */
      REAL8TimeSeries       *signalvecREAL8=NULL;
      REAL8 Q, centre_frequency,hrss,eccentricity,polar_angle;
      
      Q=injEvent->q;
      centre_frequency=injEvent->frequency;
      hrss=injEvent->hrss;
     // psi=injEvent->psi;
      polar_angle=injEvent->pol_ellipse_angle;
      eccentricity=injEvent->pol_ellipse_e; 
    /* Check that 2*width_gauss_envelope is inside frequency range */
    if ((centre_frequency + 3.0*centre_frequency/Q)>=  1.0/(2.0*thisData->timeData->deltaT))
{
    fprintf(stdout, "WARNING: Your sample rate is too low to ensure a good analysis for a SG centered at f0=%lf and with Q=%lf. Consider increasing it to more than %lf. Exiting...\n",centre_frequency,Q,2.0*(centre_frequency + 3.0*centre_frequency/Q));
//exit(1);
}
    if ((centre_frequency -3.0*centre_frequency/Q)<=  thisData->fLow)
{
    fprintf(stdout, "WARNING: The low frenquency tail of your SG centered at f0=%lf and with Q=%lf will lie below the low frequency cutoff. Whit your current settings and parameters the minimum f0 you can analyze without cuts is %lf.\n Continuing... \n",centre_frequency,Q,centre_frequency -3.0*centre_frequency/Q);
//exit(1);
}
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
      XLALGPSAddGPS(&(hplus->epoch), &(injEvent->time_geocent_gps));  
      XLALGPSAddGPS(&(hcross->epoch), &(injEvent->time_geocent_gps));
//for(i=0;i<hcross->data->length;i++){      
    
    //printf("%10.10e \n",hcross->data->data[i]);
    //}
      signalvecREAL8=XLALSimDetectorStrainREAL8TimeSeries(hplus, hcross, injEvent->ra, injEvent->dec, injEvent->psi, det.site);
      if (!signalvecREAL8) XLAL_ERROR_VOID(XLAL_EFUNC);
      
      for(i=0;i<signalvecREAL8->data->length;i++){
         // printf("%10.10e \n",signalvecREAL8->data->data[i]);
        if(isnan(signalvecREAL8->data->data[i])) {signalvecREAL8->data->data[i]=0.0;printf("isnan %d\n",i);}
      }
      
      if(signalvecREAL8->data->length > thisData->timeData->data->length-(UINT4)ceil((2.0*padding)/thisData->timeData->deltaT)){
        fprintf(stderr, "WARNING: waveform length = %u is longer than thisData->timeData->data->length = %d minus the window width = %d (total of %d points available).\n", signalvecREAL8->data->length, thisData->timeData->data->length, (INT4)ceil((2.0*padding)/thisData->timeData->deltaT) , thisData->timeData->data->length-(INT4)ceil((2.0*padding)/thisData->timeData->deltaT));
        fprintf(stderr, "The waveform injected is %f seconds long. Consider increasing the %f seconds segment length (--seglen) to be greater than %f. (in %s, line %d)\n",signalvecREAL8->data->length * thisData->timeData->deltaT , thisData->timeData->data->length * thisData->timeData->deltaT, signalvecREAL8->data->length * thisData->timeData->deltaT + 2.0*padding , __FILE__, __LINE__);
      }
      //for(i=0;i<signal->data->length;i++){      
    
    //printf("%10.10e \n",hcross->data->data[i]);
    //}
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
        UINT4 upper=thisData->fHigh/injF->deltaF;
	for(SNR=0.0,j=thisData->fLow/injF->deltaF;j<upper;j++){
	  SNR+=pow(injF->data->data[j].re,2.0)/thisData->oneSidedNoisePowerSpectrum->data->data[j];
	  SNR+=pow(injF->data->data[j].im,2.0)/thisData->oneSidedNoisePowerSpectrum->data->data[j];
	}
        SNR*=4.0*injF->deltaF;
    }
    thisData->SNR=sqrt(SNR);
    NetworkSNR+=SNR;
    
    //if (thisData->SNR > previous_snr) {best_ifo_snr=highest_snr_index;    previous_snr=thisData->SNR;}
    //highest_snr_index++;

    if (!(BurstSNRpath==NULL)){ /* If the user provided a path with --snrpath store a file with injected SNRs */
    PrintBurstSNRsToFile(IFOdata , injEvent);
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
    fprintf(stdout,"Network SNR of event %d = %.4f\n",event,NetworkSNR);
    //if (NetworkSNR<8.0) {fprintf(stderr,"NetSNR below 8. Exiting...\n");exit(1);} 
    //if (NetworkSNR>50.0) {fprintf(stderr,"NetSNR above 50. Exiting...\n");exit(1);}

    thisData=IFOdata;
    //LALInferenceIFOData *IFOdataRed=NULL;
    //UINT4 Nifo=3;
    //IFOdataRed=calloc(sizeof(LALInferenceIFOData),Nifo-1);
    
    /*highest_snr_index=0;
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
    }*/
    
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
    LALInferenceAddVariable(vars, "hrss", &hrss, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
    LALInferenceAddVariable(vars, "polar_angle", &pol_angle, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
    LALInferenceAddVariable(vars, "eccentricity", &eccentricity, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
    LALInferenceAddVariable(vars, "polarisation", &(psi), LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
    LALInferenceAddVariable(vars, "declination", &dec, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
    LALInferenceAddVariable(vars, "rightascension", &ra, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
    LALInferenceAddVariable(vars, "loghrss", &loghrss, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);

}

static void PrintBurstSNRsToFile(LALInferenceIFOData *IFOdata ,SimBurst *inj_table){
    char SnrName[300];
    char ListOfIFOs[10]="";
    REAL8 NetSNR=0.0;
    
    LALInferenceIFOData *thisData=IFOdata;
    int nIFO=0;

    while(thisData){
         sprintf(ListOfIFOs,"%s%s",ListOfIFOs,thisData->name);
         thisData=thisData->next;
	nIFO++;
        }
    
    sprintf(SnrName,"%s/snr_%s_%10.1f.dat",BurstSNRpath,ListOfIFOs,(REAL8) inj_table->time_geocent_gps.gpsSeconds);
    FILE * snrout = fopen(SnrName,"w");
    if(!snrout){
	fprintf(stderr,"Unable to open the path %s for writing SNR files\n",BurstSNRpath);
	exit(1);
    }
    
    thisData=IFOdata; // restart from the first IFO
    while(thisData){
        fprintf(snrout,"%s:\t %4.2f\n",thisData->name,thisData->SNR);
        NetSNR+=(thisData->SNR*thisData->SNR);
        thisData=thisData->next;
    }		
    if (nIFO>1){  fprintf(snrout,"Network:\t");
    fprintf(snrout,"%4.2f\n",sqrt(NetSNR));}
    fclose(snrout);
}

void InjectSineGaussianFD(LALInferenceIFOData *IFOdata, SimBurst *inj_table, ProcessParamsTable *commandLine)
///*-------------- Inject in Frequency domain -----------------*/
{
    
     fprintf(stdout,"Injecting SineGaussian in the frequency domain\n");
          fprintf(stdout,"REMEMBER!!!!!! I HARD CODED h_plus=0 in LALSimBurst.c. Remember to restore  it .\n");
    /* Inject a gravitational wave into the data in the frequency domain */ 
    LALStatus status;
    memset(&status,0,sizeof(LALStatus));
    (void) commandLine;
    LALInferenceVariables *modelParams=NULL;
    LALInferenceIFOData * tmpdata=IFOdata;
    REAL8 Q =0.0;
    REAL8 hrss,loghrss = 0.0;
    REAL8 centre_frequency= 0.0;
    REAL8 polar_angle=0.0;
    REAL8 eccentricity=0.0;
    REAL8 latitude=0.0;
    REAL8 polarization=0.0;
    REAL8 injtime=0.0;
    REAL8 longitude;
	//LALInferenceIFOData *thisData=NULL;
    tmpdata->modelParams=XLALCalloc(1,sizeof(LALInferenceVariables));
	modelParams=tmpdata->modelParams;
    memset(modelParams,0,sizeof(LALInferenceVariables));
        
        Q=inj_table->q;
      centre_frequency=inj_table->frequency;
      hrss=inj_table->hrss;
     polarization=inj_table->psi;
      polar_angle=inj_table->pol_ellipse_angle;
      eccentricity=inj_table->pol_ellipse_e; // salvo
    loghrss=log(hrss);
injtime=inj_table->time_geocent_gps.gpsSeconds + 1e-9*inj_table->time_geocent_gps.gpsNanoSeconds;
latitude=inj_table->dec;
longitude=inj_table->ra;
//printf("----time before template call %10.10e\n",injtime);
    LALInferenceAddVariable(tmpdata->modelParams, "loghrss",&loghrss,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_LINEAR);
    LALInferenceAddVariable(tmpdata->modelParams, "Q",&Q,LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);  
    LALInferenceAddVariable(tmpdata->modelParams, "rightascension",&longitude,LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_CIRCULAR);  
    LALInferenceAddVariable(tmpdata->modelParams, "declination",&latitude,LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_CIRCULAR);  
    LALInferenceAddVariable(tmpdata->modelParams, "polarisation",&polarization,LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_CIRCULAR);  
    LALInferenceAddVariable(tmpdata->modelParams, "time",&injtime,LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
    LALInferenceAddVariable(tmpdata->modelParams, "polar_angle",&polar_angle,LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_CIRCULAR);
    LALInferenceAddVariable(tmpdata->modelParams, "eccentricity",&eccentricity,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_LINEAR);
    LALInferenceAddVariable(tmpdata->modelParams, "frequency",&centre_frequency,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_LINEAR);
    
      
    COMPLEX16FrequencySeries *freqModelhCross=NULL;
   freqModelhCross=XLALCreateCOMPLEX16FrequencySeries("freqDatahC",&(tmpdata->timeData->epoch),0.0,tmpdata->freqData->deltaF,&lalDimensionlessUnit,tmpdata->freqData->data->length);
    COMPLEX16FrequencySeries *freqModelhPlus=NULL;
    freqModelhPlus=XLALCreateCOMPLEX16FrequencySeries("freqDatahP",&(tmpdata->timeData->epoch),0.0,tmpdata->freqData->deltaF,&lalDimensionlessUnit,tmpdata->freqData->data->length);
    COMPLEX16FrequencySeries *freqTemplate=NULL;
    freqTemplate=XLALCreateCOMPLEX16FrequencySeries("freqTemplate",&(tmpdata->timeData->epoch),0.0,tmpdata->freqData->deltaF,&lalDimensionlessUnit,tmpdata->freqData->data->length);
    tmpdata->freqModelhPlus=freqModelhPlus;
    tmpdata->freqModelhCross=freqModelhCross;
    
      tmpdata->modelDomain = LALINFERENCE_DOMAIN_FREQUENCY;
    LALInferenceTemplateSineGaussianF(tmpdata);
    
     
    LALInferenceVariables *currentParams=IFOdata->modelParams;
       
  double Fplus, Fcross;
  double FplusScaled, FcrossScaled;
  REAL8 plainTemplateReal, plainTemplateImag;
  REAL8 templateReal, templateImag;
  int i, lower, upper;
  LALInferenceIFOData *dataPtr;
  double ra, dec, psi, gmst;
  LIGOTimeGPS GPSlal;
  double chisquared;
  double timedelay;  /* time delay b/w iterferometer & geocenter w.r.t. sky location */
  double timeshift;  /* time shift (not necessarily same as above)                   */
  double deltaT, deltaF, twopit, f, re, im;
 UINT4 j=0;
  REAL8 temp=0.0;
    REAL8 NetSNR=0.0;
  LALInferenceVariables intrinsicParams;

  /* determine source's sky location & orientation parameters: */
  ra        = *(REAL8*) LALInferenceGetVariable(currentParams, "rightascension"); /* radian      */
  dec       = *(REAL8*) LALInferenceGetVariable(currentParams, "declination");    /* radian      */
  psi       = *(REAL8*) LALInferenceGetVariable(currentParams, "polarisation");   /* radian      */
  /* figure out GMST: */
  //XLALGPSSetREAL8(&GPSlal, GPSdouble); //This is what used in the likelihood. It seems off by two seconds (should not make a big difference as the antenna patterns would not change much in such a short interval)
  XLALGPSSetREAL8(&GPSlal, injtime);
  //UandA.units    = MST_RAD;
  //UandA.accuracy = LALLEAPSEC_LOOSE;
  //LALGPStoGMST1(&status, &gmst, &GPSlal, &UandA);
  gmst=XLALGreenwichMeanSiderealTime(&GPSlal);
  intrinsicParams.head      = NULL;
  intrinsicParams.dimension = 0;
  LALInferenceCopyVariables(currentParams, &intrinsicParams);
  LALInferenceRemoveVariable(&intrinsicParams, "rightascension");
  LALInferenceRemoveVariable(&intrinsicParams, "declination");
  LALInferenceRemoveVariable(&intrinsicParams, "polarisation");
  LALInferenceRemoveVariable(&intrinsicParams, "time");
	
  /* loop over data (different interferometers): */
  dataPtr = IFOdata;
  for (j=0; j<freqTemplate->data->length; ++j){
        freqTemplate->data->data[j].re=freqTemplate->data->data[j].im=0.0;
    }
    
  while (dataPtr != NULL) {
     
      if (IFOdata->modelDomain == LALINFERENCE_DOMAIN_TIME) {
	  printf("There is a problem. You seem to be using a time domain model into the frequency domain injection function!. Exiting....\n"); 
      exit(1);
    }
      
    /*-- WF to inject is now in dataPtr->freqModelhPlus and dataPtr->freqModelhCross. --*/
    /* determine beam pattern response (F_plus and F_cross) for given Ifo: */
    XLALComputeDetAMResponse(&Fplus, &Fcross,
                             dataPtr->detector->response,
			     ra, dec, psi, gmst);
    /* signal arrival time (relative to geocenter); */
    timedelay = XLALTimeDelayFromEarthCenter(dataPtr->detector->location,
                                             ra, dec, &GPSlal);
    //printf("----time after template call %10.10e\n",(*(REAL8*) LALInferenceGetVariable(IFOdata->modelParams, "time")));
    dataPtr->injtime=injtime;
    /* (negative timedelay means signal arrives earlier at Ifo than at geocenter, etc.) */
    /* amount by which to time-shift template (not necessarily same as above "timedelay"): */
    timeshift =  (injtime - (*(REAL8*) LALInferenceGetVariable(IFOdata->modelParams, "time"))) + timedelay;
    twopit    = LAL_TWOPI * (timeshift);
    /* include distance (overall amplitude) effect in Fplus/Fcross: */
    FplusScaled  = Fplus ;
    FcrossScaled = Fcross;
    //printf("diff in inj %lf \n", (injtime - (*(REAL8*) LALInferenceGetVariable(IFOdata->modelParams, "time"))));
    dataPtr->fPlus = FplusScaled;
    dataPtr->fCross = FcrossScaled;
    dataPtr->timeshift = timeshift;
    
    char InjFileName[50];
    sprintf(InjFileName,"FD_injection_%s.dat",dataPtr->name);
    FILE *outInj=fopen(InjFileName,"w");
    REAL8 starttime=IFOdata->epoch.gpsSeconds+1.0e-9 *IFOdata->epoch.gpsNanoSeconds;
    //printf("time + shift= %lf \n", time+(injtime - (*(REAL8*) LALInferenceGetVariable(IFOdata->modelParams, "time"))));
    //printf("IFOdata time %lf, Time data time %lf\n",time,IFOdata->timeData->epoch.gpsSeconds+1.0e-9 *IFOdata->timeData->epoch.gpsNanoSeconds);
     char TInjFileName[50];
    sprintf(TInjFileName,"TD_injection_%s.dat",dataPtr->name);
    FILE *outTInj=fopen(TInjFileName,"w");
    
     /* determine frequency range & loop over frequency bins: */
    deltaT = dataPtr->timeData->deltaT;
    deltaF = 1.0 / (((double)dataPtr->timeData->data->length) * deltaT);
    REAL8 time_env_2sigma=Q / (LAL_TWOPI * centre_frequency);
    if (2.0*time_env_2sigma>1./deltaF)
        fprintf(stdout,"WARNING: 95 of the Gaussian envelop (%lf) is larger than seglen (%lf)!!\n",2.0*time_env_2sigma,1./deltaF);
    lower = (UINT4)ceil(dataPtr->fLow / deltaF);
    upper = (UINT4)floor(dataPtr->fHigh / deltaF);
     chisquared = 0.0;
    for (i=lower; i<=upper; ++i){
      /* derive template (involving location/orientation parameters) from given plus/cross waveforms: */
      
      plainTemplateReal =FplusScaled * IFOdata->freqModelhPlus->data->data[i].re  
                          + FcrossScaled * IFOdata->freqModelhCross->data->data[i].re; //SALVO
      plainTemplateImag = FplusScaled * IFOdata->freqModelhPlus->data->data[i].im  
                          + FcrossScaled * IFOdata->freqModelhCross->data->data[i].im;
                          
   
      f = ((double) i) * deltaF;
      /* real & imag parts of  exp(-2*pi*i*f*deltaT): */
      re = cos(twopit * f);
      im = - sin(twopit * f);
      templateReal = (plainTemplateReal*re - plainTemplateImag*im);
      templateImag = (plainTemplateReal*im + plainTemplateImag*re);
      freqTemplate->data->data[i].re=templateReal;
      freqTemplate->data->data[i].im=templateImag;
      dataPtr->freqData->data->data[i].re+=templateReal;
      dataPtr->freqData->data->data[i].im+=templateImag;
      temp = ((2.0/( deltaT*(double) dataPtr->timeData->data->length) * (templateReal*templateReal+templateImag*templateImag)) / dataPtr->oneSidedNoisePowerSpectrum->data->data[i]);
      chisquared  += temp;
      fprintf(outInj,"%lf %10.10e %10.10e\n",f,templateReal,templateImag);
    }
    printf("injected SNR %.1f in IFO %s\n",sqrt(2.0*chisquared),dataPtr->name);
    NetSNR+=2.0*chisquared;
    dataPtr->SNR=sqrt(2.0*chisquared);
     
    dataPtr = dataPtr->next;
    
     fclose(outInj);
  

    /* Calculate IFFT and print it to file */
    REAL8TimeSeries* timeData=NULL;
    //printf("doing IFFT. epoch %10.10f , dT %lf lenght %i \n",IFOdata->timeData->epoch.gpsSeconds+1.0e-9*IFOdata->timeData->epoch.gpsNanoSeconds,(REAL8)IFOdata->timeData->deltaT,IFOdata->timeData->data->length);
    timeData=(REAL8TimeSeries *) XLALCreateREAL8TimeSeries("name",&IFOdata->timeData->epoch,0.0,(REAL8)IFOdata->timeData->deltaT,&lalDimensionlessUnit,(size_t)IFOdata->timeData->data->length);
     
     XLALREAL8FreqTimeFFT(timeData,freqTemplate,IFOdata->freqToTimeFFTPlan);
     for (j=0;j<timeData->data->length;j++){
         
         fprintf(outTInj,"%10.10e %10.10e \n",starttime+j*deltaT,timeData->data->data[j]);
         
         
         }
     fclose(outTInj);
    //if(!timeData) XLAL_ERROR_NULL(XLAL_EFUNC);
    
  
  
  
  }
  

    LALInferenceDestroyVariables(&intrinsicParams);
    printf("injected Network SNR %.1f \n",sqrt(NetSNR)); 
  NetSNR=sqrt(NetSNR); 

    if (!(BurstSNRpath==NULL)){ /* If the user provided a path with --snrpath store a file with injected SNRs */
        PrintBurstSNRsToFile(IFOdata , inj_table);
    }
        
	XLALDestroyCOMPLEX16FrequencySeries(freqModelhCross);
    XLALDestroyCOMPLEX16FrequencySeries(freqModelhPlus);

}
