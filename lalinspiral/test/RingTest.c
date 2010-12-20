#include <stdlib.h>
#include <getopt.h>

#include <lal/LALStdlib.h>
#include <lal/LALStdio.h>
#include <lal/FrameCache.h>
#include <lal/FrameStream.h>
#include <lal/Units.h>
#include "LALInspiralMCMC.h"
#include "LALInspiralMCMCUser.h"
#include <lal/LIGOLwXMLInspiralRead.h>
#include <lal/Random.h>
#include <lal/TimeFreqFFT.h>
#include <lal/LALDetectors.h>
#include <lal/GeneratePPNInspiral.h>
#include <lal/SimulateCoherentGW.h>
#include <lal/LALStatusMacros.h>
#include <lal/LALNoiseModels.h>
#include <lal/Date.h>
#include <lal/LALInspiral.h>
#include <lal/GenerateInspiral.h>
#include <lal/FrequencySeries.h>
#include <lal/ResampleTimeSeries.h>
#include <lal/TimeSeries.h>
#include <lal/VectorOps.h>
#include <LALAppsVCSInfo.h>
#include <lalapps.h>
#include <lal/GeneratePPNAmpCorConsistency.h>
#include "nest_calc.h"

REAL8 MCMCLikelihoodMultiCoherentAmpCorTestGrid(REAL8 mc, REAL8 eta);
REAL8 MCMCLikelihoodMultiCoherentAmpCorGrid(REAL8 mc, REAL8 eta);
INT4 ampOrder=0;
INT4 SampleRate=1024;
INT4 nSegs=250;
REAL8 duration=2000.0;
UINT4 seglen=0;
UINT4 dataseed=10;

int main(void){
    static LALStatus status;
    LALMCMCInput inputMCMC;
    typedef void (NoiseFunc)(LALStatus *status,REAL8 *psd,REAL8 f);
    NoiseFunc *PSD=NULL;
    REAL8 mchirparray[100];
    REAL8 etaarray[100];
    REAL8 logL[100,100];
    LIGOTimeGPS datastart;
    REAL8 GPS=939936910.0;
    RandomParams *datarandparam=NULL;
    for (i=0; i<100; i++) {
        mchirparray[i]=9.5+i*(10.5-9.5)/100.0;
        etaarray[i]=0.15+i*(0.25-0.15)/100.0
    }
    
    datarandparam=XLALCreateRandomParams(dataseed);
    REAL8 segDur = duration/(REAL8)nSegs;
	seglen=(UINT4)(segDur*SampleRate);
    segDur = seglen/SampleRate;
	nSegs =(INT4)floor(duration/segDur);
    
    inputMCMC.deltaT=(REAL8 )(1.0/SampleRate);
    inputMCMC.deltaF = (REAL8)SampleRate/seglen;
    
    XLALGPSSetREAL8(&datastart,GPS);
    PSD = &LALAdvLIGOPsd; 
    scalefactor = 1E-49;
    inputMCMC.invspec=(REAL8FrequencySeries *)XLALCreateREAL8FrequencySeries("inverse spectrum",&datastart,0.0,(REAL8)(SampleRate)/seglen,&lalDimensionlessUnit,seglen/2 +1);
    for(j=0;j<inputMCMC.invspec->data->length;j++){ PSD(&status,&(inputMCMC.invspec->data->data[j]),j*inputMCMC.deltaF);}
    inputMCMC.stilde = (COMPLEX16FrequencySeries *)XLALCreateCOMPLEX16FrequencySeries("stilde",&datastart,0.0,inputMCMC.deltaF,&lalDimensionlessUnit,seglen/2 +1);
    memcpy(&(inputMCMC.stilde->epoch),&segmentStart,sizeof(LIGOTimeGPS));
    /* Create the fake data */
    for(j=0;j<inputMCMC.invspec->data->length;j++){
        inputMCMC.invspec->data->data[j]=1.0/(scalefactor*inputMCMC.invspec->data->data[j]);
        inputMCMC.stilde->data->data[j].re=XLALNormalDeviate(datarandparam)/(2.0*sqrt(inputMCMC.invspec->data->data[j]*inputMCMC.deltaF));
        inputMCMC.stilde->data->data[j].im=XLALNormalDeviate(datarandparam)/(2.0*sqrt(inputMCMC.invspec->data->data[j]*inputMCMC.deltaF));
    }
    inputMCMC.ampOrder=ampOrder;
    inputMCMC.phaseOrder=LAL_PNORDER_TWO;
    
    FILE *loglout;
    loglout=fopen("logl.out","w");
    for (i=0;i<100; i++) {
        for (j=0; j<100; j++) {
            logL[i,j]=MCMCLikelihoodMultiCoherentAmpCorTestGrid(mchirparray[i],etaarray[j])
            fprintf(outlogL,"%e\t %e\t %e\n",logL[i,j],mchirparray[i],etaarray[j]);
        }
    }
    fclose(loglout);
    
}    
REAL8 MCMCLikelihoodMultiCoherentAmpCorTestGrid(REAL8 mc,REAL8 eta){
    /* Calculate the likelihood for an amplitude-corrected waveform */
    /* This template is generated in the time domain */
    
    /* Test params inserted */
    REAL8 logL=0.0,chisq=0.0;
    REAL8 end_time,resp_r,resp_i,real,imag;
    UINT4 det_i=0,idx=0;
    UINT4 i=0;
    DetectorResponse det;
    static LALStatus status;
    CoherentGW coherent_gw;
    PPNConsistencyParamStruc PPNparams;
    LALDetAMResponse det_resp;
    REAL4TimeSeries *h_p_t=NULL,*h_c_t=NULL;
    COMPLEX8FrequencySeries *H_p_t=NULL, *H_c_t=NULL;
    size_t NFD = 0;
    memset(&PPNparams,0,sizeof(PPNparams));
    memset(&coherent_gw,0,sizeof(CoherentGW));
    memset(&status,0,sizeof(LALStatus));
    memset(&det,0,sizeof(DetectorResponse));
    /* Populate the structures */
    PPNparams.position.longitude=1.0;
    PPNparams.position.latitude=2.0;
    PPNparams.position.system=COORDINATESYSTEM_EQUATORIAL;
    PPNparams.psi=0.0;
    memcpy(&(PPNparams.epoch),&(inputMCMC->epoch),sizeof(LIGOTimeGPS));
    PPNparams.mTot=mc2mt(mc,eta);
    PPNparams.eta=eta;
    PPNparams.mTot_real8=mc2mt(mc,eta);
    PPNparams.eta_real8=eta;
    PPNparams.d=100.0*MpcInMeters;
    PPNparams.inc=1.0;
    PPNparams.phi=0.0;
    PPNparams.fStartIn=inputMCMC->fLow;
    PPNparams.fStopIn=0.5/inputMCMC->deltaT;
    PPNparams.deltaT=inputMCMC->deltaT;
    PPNparams.ampOrder = inputMCMC->ampOrder;
    
    LALPopulatePhasePNparams(&PPNparams,PhaseTestParam);
    
    /* GET TEST PHASE PARAMETER FROM MCMCSTRUCTURE */
    if (PhaseTestParam!=-1) {PPNparams.phasePNparams[PhaseTestParam] = XLALMCMCGetParameter(parameter,"phiTest");}
    
    /* Call LALGeneratePPNAmpCorConsistency */
    LALGeneratePPNAmpCorConsistency(&status,&coherent_gw,&PPNparams);
    
    if(status.statusCode)
    {
        REPORTSTATUS(&status);
        chisq=DBL_MAX;
        goto noWaveform;
    }
    
    
    /* Set the epoch so that the t_c is correct */
    end_time = XLALMCMCGetParameter(parameter,"time");
    
    REAL8 adj_epoch = end_time-PPNparams.tc;
    if(coherent_gw.h) XLALGPSSetREAL8(&(coherent_gw.h->epoch),adj_epoch);
    if(coherent_gw.a) XLALGPSSetREAL8(&(coherent_gw.a->epoch),adj_epoch);
    
    /* Inject h+ and hx into time domain signal of correct length */
    UINT4 NtimeDomain=inputMCMC->segment?inputMCMC->segment->data->length:2*(inputMCMC->stilde->data->length-1);
    h_p_t = XLALCreateREAL4TimeSeries("hplus",&inputMCMC->epoch,
                                      inputMCMC->fLow,inputMCMC->deltaT,
                                      &lalADCCountUnit,NtimeDomain);
    h_c_t = XLALCreateREAL4TimeSeries("hcross",&inputMCMC->epoch,
                                      inputMCMC->fLow,inputMCMC->deltaT,
                                      &lalADCCountUnit,NtimeDomain);
    if(!(h_p_t && h_c_t)){
        fprintf(stderr,"Unable to allocate signal buffer\n");
        exit(1);
    }
    
    /* Separate the + and x parts */
    for(i=0;i< (NtimeDomain<coherent_gw.h->data->length?NtimeDomain:coherent_gw.h->data->length) ;i++){
        h_p_t->data->data[i]=coherent_gw.h->data->data[2*i];
        h_c_t->data->data[i]=coherent_gw.h->data->data[2*i + 1];
    }
    for(;i<NtimeDomain;i++){
        h_p_t->data->data[i]=0.0;
        h_c_t->data->data[i]=0.0;
    }
    
    /* Get H+ and Hx in the Freq Domain */
    NFD=inputMCMC->stilde[det_i]->data->length;
    H_p_t = XLALCreateCOMPLEX8FrequencySeries("Hplus",&inputMCMC->epoch,0,(REAL8)inputMCMC->deltaF,&lalDimensionlessUnit,(size_t)NFD);
    H_c_t = XLALCreateCOMPLEX8FrequencySeries("Hcross",&inputMCMC->epoch,0,(REAL8)inputMCMC->deltaF,&lalDimensionlessUnit,(size_t)NFD);
    
    if(!(H_p_t && H_c_t)){
        fprintf(stderr,"Unable to allocate F-domain signal buffer\n");
        exit(1);
    }
    
    if(inputMCMC->likelihoodPlan==NULL) {LALCreateForwardREAL4FFTPlan(&status,&inputMCMC->likelihoodPlan,
                                                                      NtimeDomain,
                                                                      FFTW_PATIENT); fprintf(stderr,"Created FFTW plan\n");}
    LALTimeFreqRealFFT(&status,H_p_t,h_p_t,inputMCMC->likelihoodPlan);
    LALTimeFreqRealFFT(&status,H_c_t,h_c_t,inputMCMC->likelihoodPlan);
    
    XLALDestroyREAL4TimeSeries(h_p_t);
    XLALDestroyREAL4TimeSeries(h_c_t);
    
    /* The epoch of observation and the accuracy required ( we don't care about a few leap seconds) */
    LALSource source; /* The position and polarisation of the binary */
    source.equatorialCoords.longitude =1.0;
    source.equatorialCoords.latitude = 2.0;
    source.equatorialCoords.system = COORDINATESYSTEM_EQUATORIAL;
    source.orientation = 1.0;
    
    /* This also holds the source and the detector, LAL has two different structs for this! */
    LALDetAndSource det_source;
    det_source.pSource=&source;
    
    REAL8 TimeShiftToGC=XLALMCMCGetParameter(parameter,"time");
    
    TimeShiftToGC-=inputMCMC->epoch.gpsSeconds + 1.e-9*inputMCMC->epoch.gpsNanoSeconds;
    TimeShiftToGC-=PPNparams.tc;
    
    REAL8 TimeFromGC;
    /* Set up the detector */
    det.site=inputMCMC->detector;
    
    TimeFromGC = XLALTimeDelayFromEarthCenter(inputMCMC->detector->location, source.equatorialCoords.longitude, source.equatorialCoords.latitude, &(inputMCMC->epoch)); /* Compute time delay */
    /* Compute detector amplitude response */
    det_source.pDetector = (inputMCMC->detector); /* select detector */
    LALComputeDetAMResponse(&status,&det_resp,&det_source,&(inputMCMC->epoch)); /* Compute det_resp */
    /* No need to multiply by cos(iota) as GenerateAmpCorPPNInspiral() takes this into account */
    chisq=0.0;
    /* Calculate the logL */
    REAL8 deltaF = inputMCMC->stilde->deltaF;
    UINT4 lowBin = (UINT4)(inputMCMC->fLow / inputMCMC->stilde->deltaF);
    UINT4 highBin;
    /* Only compute sum up to maximum frquency of the waveform. PPNparams.fStop is the maximum of the 2*f_orb harmonic */
    
    REAL8 fMultiplier = (inputMCMC->ampOrder + 2.0)/2.0; /* The frequency of the highest harmonic as determined by ampOrder */
    highBin = (UINT4)(PPNparams.fStop * fMultiplier / inputMCMC->stilde->deltaF);
    
    if(highBin==0 || highBin>inputMCMC->stilde->data->length-1)
        highBin=inputMCMC->stilde->data->length-1;  /* AmpCor waveforms don't set the highest frequency of the highest harmonic */
    for(idx=lowBin;idx<=highBin;idx++){
        /* The phase shift angle, determined by the time FROM the geocentre to the Detector, plus the time */
        /* FROM the start of the segment to the t_c at geocentre */
        /* Phase shift is exp(-i*ang), but - signs on sin below take the negative into account */
        /* exp(-i*ang) = cos(ang) - sin(ang) */
        REAL8 ang = 2.0*LAL_PI*(TimeFromGC+TimeShiftToGC)*inputMCMC->stilde->deltaF*idx;
        /* Calculate rotated parts of the plus and cross */
        
        /* Negative signs on sins: see comment above for definition of ang */
        REAL4 plus_re,plus_im,cross_re,cross_im;
        plus_re = H_p_t->data->data[idx].re*cos(ang) + H_p_t->data->data[idx].im*sin(ang);
        plus_im = H_p_t->data->data[idx].im*cos(ang) - H_p_t->data->data[idx].re*sin(ang);
        cross_re = H_c_t->data->data[idx].re*cos(ang) + H_c_t->data->data[idx].im*sin(ang);
        cross_im = H_c_t->data->data[idx].im*cos(ang) - H_c_t->data->data[idx].re*sin(ang);
        
        /* Compute total real and imaginary responses */
        resp_r = (REAL8)( plus_re*det_resp.plus + cross_re*det_resp.cross );
        resp_i = (REAL8)( plus_im*det_resp.plus + cross_im*det_resp.cross );
        real=inputMCMC->stilde->data->data[idx].re - resp_r;
        imag=inputMCMC->stilde->data->data[idx].im - resp_i;
        
        /* Gaussian version */
        chisq+=(real*real + imag*imag)*inputMCMC->invspec->data->data[idx];
    }
    /* Add on the remaining sum, consulting the lookup table */
    if(highBin<inputMCMC->stilde->data->length-2 && highBin>lowBin) chisq+=topdown_sum->data[highBin+1];
    else if(highBin<=lowBin) chisq+=topdown_sum->data[highBin+1];
    chisq*=2.0*deltaF; /* for 2 sigma^2 on denominator, also in student-t version */
    
    logL-=chisq;
    
}
/* Destroy the response series */
if(coherent_gw.f) XLALDestroyREAL4TimeSeries(coherent_gw.f);
if(coherent_gw.phi) XLALDestroyREAL8TimeSeries(coherent_gw.phi);
if(coherent_gw.shift) XLALDestroyREAL4TimeSeries(coherent_gw.shift);
if(coherent_gw.h) {XLALDestroyREAL4VectorSequence(coherent_gw.h->data); LALFree(coherent_gw.h);}
if(coherent_gw.a) {XLALDestroyREAL4VectorSequence(coherent_gw.a->data); LALFree(coherent_gw.a);}
XLALDestroyCOMPLEX8FrequencySeries(H_p_t);
XLALDestroyCOMPLEX8FrequencySeries(H_c_t);

noWaveform:
/* return logL */
return(logL);

}

REAL8 MCMCLikelihoodMultiCoherentAmpCorGrid(REAL8 mc,REAL8 eta){
    /* Calculate the likelihood for an amplitude-corrected waveform */
    /* This template is generated in the time domain */
    REAL8 logL=0.0,chisq=0.0;
    REAL8 mc,eta,end_time,resp_r,resp_i,real,imag;
    UINT4 det_i=0,idx=0;
    UINT4 i=0;
    DetectorResponse det;
    static LALStatus status;
    CoherentGW coherent_gw;
    PPNParamStruc PPNparams;
    LALDetAMResponse det_resp;
    REAL4TimeSeries *h_p_t=NULL,*h_c_t=NULL;
    COMPLEX8FrequencySeries *H_p_t=NULL, *H_c_t=NULL;
    size_t NFD = 0;
    memset(&PPNparams,0,sizeof(PPNparams));
    memset(&coherent_gw,0,sizeof(CoherentGW));
    memset(&status,0,sizeof(LALStatus));
    memset(&det,0,sizeof(DetectorResponse));
    /* Populate the structures */
    if(XLALMCMCCheckParameter(parameter,"logM")) mc=exp(XLALMCMCGetParameter(parameter,"logM"));
    else mc=XLALMCMCGetParameter(parameter,"mchirp");
    eta=XLALMCMCGetParameter(parameter,"eta");
    PPNparams.position.longitude=XLALMCMCGetParameter(parameter,"long");
    PPNparams.position.latitude=XLALMCMCGetParameter(parameter,"lat");
    PPNparams.position.system=COORDINATESYSTEM_EQUATORIAL;
    PPNparams.psi=XLALMCMCGetParameter(parameter,"psi");
    memcpy(&(PPNparams.epoch),&(inputMCMC->epoch),sizeof(LIGOTimeGPS));
    PPNparams.mTot=mc2mt(mc,eta);
    PPNparams.eta=eta;
    if (XLALMCMCCheckParameter(parameter,"logdist")) PPNparams.d=exp(XLALMCMCGetParameter(parameter,"logdist"))*MpcInMeters;
    else PPNparams.d=XLALMCMCGetParameter(parameter,"distMpc")*MpcInMeters;
    PPNparams.inc=XLALMCMCGetParameter(parameter,"iota");
    PPNparams.phi=XLALMCMCGetParameter(parameter,"phi");
    PPNparams.fStartIn=inputMCMC->fLow;
    PPNparams.fStopIn=0.5/inputMCMC->deltaT;
    PPNparams.deltaT=inputMCMC->deltaT;
    PPNparams.ampOrder = inputMCMC->ampOrder;
    
    /* Call LALGeneratePPNAmpCorInspiral */
    LALGeneratePPNAmpCorInspiral(&status,&coherent_gw,&PPNparams);
    
    if(status.statusCode)
    {
        REPORTSTATUS(&status);
        chisq=DBL_MAX;
        goto noWaveform;
    }
    
    /* Set the epoch so that the t_c is correct */
    end_time = XLALMCMCGetParameter(parameter,"time");
    
    REAL8 adj_epoch = end_time-PPNparams.tc;
    if(coherent_gw.h) XLALGPSSetREAL8(&(coherent_gw.h->epoch),adj_epoch);
    if(coherent_gw.a) XLALGPSSetREAL8(&(coherent_gw.a->epoch),adj_epoch);
    
    /* Inject h+ and hx into time domain signal of correct length */
    UINT4 NtimeDomain=inputMCMC->segment[det_i]?inputMCMC->segment[det_i]->data->length:2*(inputMCMC->stilde[det_i]->data->length-1);
    h_p_t = XLALCreateREAL4TimeSeries("hplus",&inputMCMC->epoch,
                                      inputMCMC->fLow,inputMCMC->deltaT,
                                      &lalADCCountUnit,NtimeDomain);
    h_c_t = XLALCreateREAL4TimeSeries("hcross",&inputMCMC->epoch,
                                      inputMCMC->fLow,inputMCMC->deltaT,
                                      &lalADCCountUnit,NtimeDomain);
    if(!(h_p_t && h_c_t)){
        fprintf(stderr,"Unable to allocate signal buffer\n");
        exit(1);
    }
    
    /* Separate the + and x parts */
    for(i=0;i< (NtimeDomain<coherent_gw.h->data->length?NtimeDomain:coherent_gw.h->data->length) ;i++){
        h_p_t->data->data[i]=coherent_gw.h->data->data[2*i];
        h_c_t->data->data[i]=coherent_gw.h->data->data[2*i + 1];
    }
    for(;i<NtimeDomain;i++){
        h_p_t->data->data[i]=0.0;
        h_c_t->data->data[i]=0.0;
    }
    
    /* Get H+ and Hx in the Freq Domain */
    NFD=inputMCMC->stilde[det_i]->data->length;
    H_p_t = XLALCreateCOMPLEX8FrequencySeries("Hplus",&inputMCMC->epoch,0,(REAL8)inputMCMC->deltaF,&lalDimensionlessUnit,(size_t)NFD);
    H_c_t = XLALCreateCOMPLEX8FrequencySeries("Hcross",&inputMCMC->epoch,0,(REAL8)inputMCMC->deltaF,&lalDimensionlessUnit,(size_t)NFD);
    
    if(!(H_p_t && H_c_t)){
        fprintf(stderr,"Unable to allocate F-domain signal buffer\n");
        exit(1);
    }
    
    if(inputMCMC->likelihoodPlan==NULL) {LALCreateForwardREAL4FFTPlan(&status,&inputMCMC->likelihoodPlan,
                                                                      NtimeDomain,
                                                                      FFTW_PATIENT); fprintf(stderr,"Created FFTW plan\n");}
    LALTimeFreqRealFFT(&status,H_p_t,h_p_t,inputMCMC->likelihoodPlan);
    LALTimeFreqRealFFT(&status,H_c_t,h_c_t,inputMCMC->likelihoodPlan);
    
#if DEBUGMODEL !=0
    char tmodelname[100];
    sprintf(tmodelname,"tmodel_plus_%i.dat",det_i);
    modelout = fopen(tmodelname,"w");
    for(i=0;i<h_p_t->data->length;i++){
        fprintf(modelout,"%g %g\n",h_p_t->data->data[i],h_c_t->data->data[i]);
    }
    fclose(modelout);
#endif
    XLALDestroyREAL4TimeSeries(h_p_t);
    XLALDestroyREAL4TimeSeries(h_c_t);
    
    /* The epoch of observation and the accuracy required ( we don't care about a few leap seconds) */
    LALSource source; /* The position and polarisation of the binary */
    source.equatorialCoords.longitude = XLALMCMCGetParameter(parameter,"long");
    source.equatorialCoords.latitude = XLALMCMCGetParameter(parameter,"lat");
    source.equatorialCoords.system = COORDINATESYSTEM_EQUATORIAL;
    source.orientation = XLALMCMCGetParameter(parameter,"psi");
    
    /* This also holds the source and the detector, LAL has two different structs for this! */
    LALDetAndSource det_source;
    det_source.pSource=&source;
    
    REAL8 TimeShiftToGC=XLALMCMCGetParameter(parameter,"time");
    
    TimeShiftToGC-=inputMCMC->epoch.gpsSeconds + 1.e-9*inputMCMC->epoch.gpsNanoSeconds;
    TimeShiftToGC-=PPNparams.tc;
    
    
    /* For each IFO */
    for(det_i=0;det_i<inputMCMC->numberDataStreams;det_i++){
        REAL8 TimeFromGC;
        /* Set up the detector */
        det.site=inputMCMC->detector[det_i];
        /* Simulate the response */
#if DEBUGMODEL !=0
        char modelname[100];
        sprintf(modelname,"model_%i.dat",det_i);
        modelout = fopen(modelname,"w");
#endif
        
        TimeFromGC = XLALTimeDelayFromEarthCenter(inputMCMC->detector[det_i]->location, source.equatorialCoords.longitude, source.equatorialCoords.latitude, &(inputMCMC->epoch)); /* Compute time delay */
        /* Compute detector amplitude response */
        det_source.pDetector = (inputMCMC->detector[det_i]); /* select detector */
        LALComputeDetAMResponse(&status,&det_resp,&det_source,&(inputMCMC->epoch)); /* Compute det_resp */
        /* No need to multiply by cos(iota) as GenerateAmpCorPPNInspiral() takes this into account */
        chisq=0.0;
        /* Calculate the logL */
        REAL8 deltaF = inputMCMC->stilde[det_i]->deltaF;
        UINT4 lowBin = (UINT4)(inputMCMC->fLow / inputMCMC->stilde[det_i]->deltaF);
        UINT4 highBin;
        /* Only compute sum up to maximum frquency of the waveform. PPNparams.fStop is the maximum of the 2*f_orb harmonic */
        
        REAL8 fMultiplier = (inputMCMC->ampOrder + 2.0)/2.0; /* The frequency of the highest harmonic as determined by ampOrder */
        highBin = (UINT4)(PPNparams.fStop * fMultiplier / inputMCMC->stilde[det_i]->deltaF);
        
        if(highBin==0 || highBin>inputMCMC->stilde[det_i]->data->length-1)
            highBin=inputMCMC->stilde[det_i]->data->length-1;  /* AmpCor waveforms don't set the highest frequency of the highest harmonic */
        for(idx=lowBin;idx<=highBin;idx++){
            /* The phase shift angle, determined by the time FROM the geocentre to the Detector, plus the time */
            /* FROM the start of the segment to the t_c at geocentre */
            /* Phase shift is exp(-i*ang), but - signs on sin below take the negative into account */
            /* exp(-i*ang) = cos(ang) - sin(ang) */
            REAL8 ang = 2.0*LAL_PI*(TimeFromGC+TimeShiftToGC)*inputMCMC->stilde[det_i]->deltaF*idx;
            /* Calculate rotated parts of the plus and cross */
            
            /* Negative signs on sins: see comment above for definition of ang */
            REAL4 plus_re,plus_im,cross_re,cross_im;
            plus_re = H_p_t->data->data[idx].re*cos(ang) + H_p_t->data->data[idx].im*sin(ang);
            plus_im = H_p_t->data->data[idx].im*cos(ang) - H_p_t->data->data[idx].re*sin(ang);
            cross_re = H_c_t->data->data[idx].re*cos(ang) + H_c_t->data->data[idx].im*sin(ang);
            cross_im = H_c_t->data->data[idx].im*cos(ang) - H_c_t->data->data[idx].re*sin(ang);
            
            /* Compute total real and imaginary responses */
            resp_r = (REAL8)( plus_re*det_resp.plus + cross_re*det_resp.cross );
            resp_i = (REAL8)( plus_im*det_resp.plus + cross_im*det_resp.cross );
            real=inputMCMC->stilde[det_i]->data->data[idx].re - resp_r;
            imag=inputMCMC->stilde[det_i]->data->data[idx].im - resp_i;
            
            /* Gaussian version */
            chisq+=(real*real + imag*imag)*inputMCMC->invspec[det_i]->data->data[idx];
            
#if DEBUGMODEL !=0
            fprintf(modelout,"%lf %10.10e %10.10e %10.10e %10.10e %10.10e %10.10e\n",idx*deltaF,resp_r,resp_i,H_p_t->data->data[idx].re,H_p_t->data->data[idx].im,H_c_t->data->data[idx].re,H_c_t->data->data[idx].im);
#endif
        }
#if DEBUGMODEL !=0
        fclose(modelout);
#endif
        /* Add on the remaining sum, consulting the lookup table */
        if(highBin<inputMCMC->stilde[det_i]->data->length-2 && highBin>lowBin) chisq+=topdown_sum[det_i]->data[highBin+1];
        else if(highBin<=lowBin) chisq+=topdown_sum[det_i]->data[highBin+1];
        chisq*=2.0*deltaF; /* for 2 sigma^2 on denominator, also in student-t version */
        
        logL-=chisq;
        
    }
    /* Destroy the response series */
    if(coherent_gw.f) XLALDestroyREAL4TimeSeries(coherent_gw.f);
    if(coherent_gw.phi) XLALDestroyREAL8TimeSeries(coherent_gw.phi);
    if(coherent_gw.shift) XLALDestroyREAL4TimeSeries(coherent_gw.shift);
    if(coherent_gw.h) {XLALDestroyREAL4VectorSequence(coherent_gw.h->data); LALFree(coherent_gw.h);}
    if(coherent_gw.a) {XLALDestroyREAL4VectorSequence(coherent_gw.a->data); LALFree(coherent_gw.a);}
    XLALDestroyCOMPLEX8FrequencySeries(H_p_t);
    XLALDestroyCOMPLEX8FrequencySeries(H_c_t);
    
noWaveform:
    /* return logL */
    parameter->logLikelihood=logL;
    return(logL);
    
}
