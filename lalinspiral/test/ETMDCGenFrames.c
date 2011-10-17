/* *
* *  ETmdc.c Created by Tania Regimbau on 8/8/09.
* *  ET Noise added by Carl Rodriguez and B.S. Sathyaprakash 02/12/09.
* *
* */

/*
* simulate times series from the extragalactic
* population of BNS + ET Gaussian noise
*/


#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include <unistd.h>
#include <getopt.h>

//#include <lal/FrameL.h>

#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALConfig.h>
#include <lal/Units.h>
#include <lal/LALConstants.h>
#include <lal/LALDatatypes.h>
#include <lal/LALRCSID.h>
#include <lal/LALStatusMacros.h>
#include <lal/AVFactories.h>
#include <lal/TimeSeries.h>
#include <lal/PrintFTSeries.h>
#include <lal/RealFFT.h>
#include <lal/ReadFTSeries.h>
#include <lal/StreamInput.h>
#include <lal/PrintVector.h>
#include <lal/VectorOps.h>
#include <lal/FileIO.h>
#include <lal/FrameStream.h>
#include <lal/FrequencySeries.h>
#include <lal/TimeFreqFFT.h>
#include <lal/Random.h>
#include <lal/Interpolate.h>
#include <lal/Random.h>
#include <lal/LALInspiral.h>
#include <lal/DetectorSite.h>
#include <lal/TimeDelay.h>
#include <lal/DetResponse.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_fft_complex.h>
#include <gsl/gsl_fft_halfcomplex.h>
#include <gsl/gsl_statistics.h>
/*
#include <lal/RealFFT.h>
*/
#define Pi LAL_PI
#define Pi2 LAL_TWOPI
#define year LAL_YRSID_SI
#define Msolar LAL_MSUN_SI*1000.
#define Ts LAL_MTSUN_SI
#define c_l LAL_C_SI*10.
#define day LAL_DAYSID_SI
/* to enable the printing of FFTs, inv PSDs and waveforms set DEBUG to 1 */
#define DEBUG 0

NRCSID (ET_MdcC, "$Id: ET_Mdc.c,v 1.2 2010/06/18 11:07:52 tania Exp $");
RCSID ("$Id: ET_Mdcs.c,v 1.1 2009/06/18 11:07:52 tania Exp $");

/* cvs info */
#define CVS_ID "$Id: ET_Mdc.c,v 1.2 2010/06/18 11:07:52 tania Exp $"
#define CVS_REVISION "$Revision: 1.2 $"
#define CVS_DATE "$Date: 2010/06/18 11:07:52 $"
#define PROGRAM_NAME "ET_Mdc"

/* variables for getopt options parsing */
extern char *optarg;
extern int optind;

/* flag for getopt_long */
static int verbose_flag = 0;
static int ascii_flag = 0;
static int catalog_flag = 0;
static int stat_flag = 0;
static int null_flag = 0;
static int noise_flag = 0;
static int write_psd_flag = 0;

/* user parameters */
UINT4 seed1=12;
UINT4 seed2=12;
UINT4 noiseSeed = 10;

UINT8 startTime = 800000000;
REAL8 Tobs = 86400.; /* 1 day */
UINT4 job = 0;
UINT4 jobStart = 0;
UINT4 Njob = 1;
int Ndet = 1; /* 1-3 */

/* source */
REAL8 zMin=0.;
REAL8 zMax=6.;

/* detector */
REAL8 rate = 8192.; /*sampling */
REAL8 fMin = 10.;

Approximant approx=TaylorT4;
UINT4 order = 7;
UINT4 ampOrder = 0;

REAL8 localRate = 1.0; /* per Mpc^3 per Myr */

REAL8 mMin = 1.2;
REAL8 mMax = 3.;
REAL8 mMean = 1.4;
REAL8 mSigma = 0.53;

/* output directory */
CHAR outputDir[500] = "output";
CHAR waveout[100];
INT4 lalDebugLevel=33;

void parseOptions(INT4 argc, CHAR *argv[]);
void displayUsage(INT4 exitcode);
int ET_noise(gsl_rng *r, double TIME, long int N, double F_START, double SAMPLE_FREQ, double *noise_data);
int LIGOI_noise(gsl_rng *r, double TIME, long int N, double F_START, double SAMPLE_FREQ, double *noise_data);
void InvPSD(REAL8TimeSeries Data, LIGOTimeGPS *GPSStart, REAL8 duration, UINT4 SampleRate, UINT4 nSegs);
void PrintHT(REAL4Vector *s, FILE **f, REAL8 t0);

int main (INT4 argc, CHAR *argv[])
{
  static LALStatus status;
  UINT4 i, k;
  UINT4 kstart = 0;
  UINT4 indexMax = 0;
  UINT8 N, Nseg, Ntot, startTimeSeg, stopTimeSeg;
  REAL4 Tseg, deltat;
  REAL8 reject, p;
  REAL8 delta, logdelta, logzMax;
  REAL8 ra, decl, cosDecl, psi, cosI;
  REAL8 PzMax, z;
  REAL8 tau0;
  REAL8 duration, durationMax;
  UINT4 length;
  LIGOTimeGPS gpsStartTime, gpsT0;
  REAL8 tm, dt, t0;
  REAL8Vector *tc = NULL;
  REAL4Vector *hp, *hc, *h;
  InspiralTemplate params;
  REAL4 apFac, acFac;
  LALDetector dE1, dE2, dE3;
  REAL4TimeSeries *Fp = NULL, *Fc = NULL;
  REAL8TimeSeries serie1, serie2, serie3, serieNS;
  REAL8 meangw,vargw,skewgw,kurtgw;
  CHAR fileName[LALNameLength];	
  CHAR frameName1[LALNameLength];
  CHAR frameName2[LALNameLength];
  CHAR frameName3[LALNameLength];
  CHAR frameNameNS[LALNameLength];	
	
  /*pointers to input/output files*/
  FILE *pfout1 = NULL, *pfout2 = NULL, *pfout3 = NULL;	
	
  /* random generator */
  RandomParams *randParams1 = NULL, *randParams2 = NULL;
	
  /* parse command line options */
  parseOptions(argc, argv);
	
  job = jobStart+job;	
	
  /* output file name */

  snprintf(frameName1, LALNameLength,"%s/E1",outputDir);
  FrOutPar frSerie1 = {frameName1,"E1:STRAIN", ADCDataChannel, 1, 0, 0 };
  snprintf(frameName2, LALNameLength,"%s/E2",outputDir);
  FrOutPar frSerie2 = {frameName2,"E2:STRAIN", ADCDataChannel, 1, 0, 0 };
  snprintf(frameName3, LALNameLength,"%s/E3",outputDir);
  FrOutPar frSerie3 = {frameName3,"E3:STRAIN", ADCDataChannel, 1, 0, 0 };
  snprintf(frameNameNS, LALNameLength,"%s/NS",outputDir);
  FrOutPar frSerieNS = {frameNameNS,"NS:STRAIN", ADCDataChannel, 1, 0, 0 };
	
  if(catalog_flag)
  {
    
    snprintf(fileName, LALNameLength,"%s/catalog_%d.dat",outputDir,job);
    pfout1 = fopen(fileName,"w");
  }
	
  if(stat_flag)
  {
   snprintf(fileName, LALNameLength,"%s/stat_%d.dat",outputDir,job);
   pfout3 = fopen(fileName,"w");
	}
  
  if(ascii_flag)
  {
    snprintf(fileName, LALNameLength,"%s/serie_%d.dat",outputDir,job);
    pfout2 = fopen(fileName,"w");
  }
  
  /* initialize status pointer */
  status.statusPtr = NULL;
  
  /* response of the detectors */
  dE1.type = dE2.type = dE3.type = LALDETECTORTYPE_IFODIFF;
  dE1.location[0] = dE2.location[0] = dE3.location[0] = 4.546374099002599e6;
  dE1.location[1] = dE2.location[1] = dE3.location[1] = 8.42989697626334e5;
  dE1.location[2] = dE2.location[2] = dE3.location[2] = 4.378576962409281e6;
	
  dE1.response[0][0] = 0.166589852497480;
  dE1.response[1][1] = -0.248382035405337;
  dE1.response[2][2] = 0.081792182907857;
  dE1.response[0][1] = dE1.response[1][0] = -0.218849471235102 ;
  dE1.response[0][2] = dE1.response[2][0] = -0.129963871963915;
  dE1.response[1][2] = dE1.response[2][1] = 0.273214957676611;
  
  dE2.response[0][0] = -0.199221201378560;
  dE2.response[1][1] = 0.423356724499319;
  dE2.response[2][2]=-0.224135523120759;
  dE2.response[0][1] = dE2.response[1][0] = -0.070223802479191;
  dE2.response[0][2] = dE2.response[2][0] = 0.218900453442919;
  dE2.response[1][2] = dE2.response[2][1] = -0.008534697228688;

  dE3.response[0][0] = 0.032631348881079;
  dE3.response[1][1] = -0.174974689093981 ;
  dE3.response[2][2] = 0.142343340212902;
  dE3.response[0][1] = dE3.response[1][0] = 0.289073273714293;
  dE3.response[0][2] = dE3.response[2][0] = -0.088936581479004;
  dE3.response[1][2] = dE3.response[2][1] = -0.264680260447922;
  
	
  Tseg = Tobs/Njob;
  Nseg = floor(Tseg*rate);
  deltat = 1./rate;
  noiseSeed = noiseSeed +job;
	
  /* set start time structure */
  startTimeSeg = startTime + job*(UINT8)Tseg;
  stopTimeSeg = startTimeSeg + (UINT8)Tseg;
  gpsStartTime.gpsSeconds = startTimeSeg ;
  gpsStartTime.gpsNanoSeconds = 0.;
  
  
  /* initialize output time series */
  strncpy(serie1.name, "E1:STRAIN", LALNameLength);
  serie1.sampleUnits = lalStrainUnit;
  serie1.epoch = gpsStartTime;
  serie1.deltaT = 1./(REAL8)rate;
  serie1.f0 = fMin;
  serie1.data = XLALCreateREAL8Vector(Nseg);
  memset( serie1.data->data, 0, Nseg * sizeof( REAL8 ) );
  
  if (Ndet>1)
  {
    strncpy(serie2.name, "E2:STRAIN", LALNameLength);
    serie2.sampleUnits = lalStrainUnit;
    serie2.epoch = gpsStartTime;
    serie2.deltaT = 1./(REAL8)rate;
    serie2.f0 = fMin;
    serie2.data = XLALCreateREAL8Vector(Nseg);
    memset( serie2.data->data, 0, Nseg * sizeof( REAL8 ) );
  }
  
  if (Ndet>2)
  {
    strncpy(serie3.name, "E3:STRAIN", LALNameLength);
    serie3.sampleUnits = lalStrainUnit;
    serie3.epoch = gpsStartTime;
    serie3.deltaT = 1./(REAL8)rate;
    serie3.f0 = fMin;
    serie3.data = XLALCreateREAL8Vector(Nseg);
    memset( serie3.data->data, 0, Nseg * sizeof( REAL8 ) );
  }
  
  if(null_flag)
  {
	strncpy(serieNS.name, "NS:STRAIN", LALNameLength);
	serieNS.sampleUnits = lalStrainUnit;
	serieNS.epoch = gpsStartTime;
	serieNS.deltaT = 1./(REAL8)rate;
	serieNS.f0 = fMin;
	serieNS.data = XLALCreateREAL8Vector(Nseg);
	memset( serieNS.data->data, 0, Nseg * sizeof( REAL8 ) );
  }
  /* calculate some constants */
  tau0 = 646972.*pow(fMin,-8./3.);/*5*c^5/(256*Pi^(8/3)*G^(5/3))/Msolar^(5/3)*/
  
  /* set general waveform parameters*/
  memset( &params, 0, sizeof(params) );
  params.approximant=approx;
  //params.order = LAL_PNORDER_THREE_POINT_FIVE;
  //params.ampOrder = LAL_PNORDER_NEWTONIAN;
  params.order = order;
  params.ampOrder = ampOrder;	
  params.fLower = fMin;
  params.fCutoff = rate/2.-1.;
  params.tSampling=rate;
  params.massChoice = m1Andm2;
  params.startTime = 0.0;
  params.signalAmplitude = 1.0;
  
  /* maximum of the redshift probability distribution */
  if (zMax<1.6)
  {
    PzMax = -0.000429072589677 + (zMax * (-0.036349728568888 + (zMax * (0.860892111762314
      + ( zMax * (-0.740935488674010 + zMax * (0.265848831356864
      + zMax * (-0.050041573542298 + zMax * (0.005184554232421
      + zMax * (-0.000281450045300 + zMax *0.000006400690921))))))))));
  }
  else
  {
    PzMax = 0.42;
  }

  /*coalescence times (Poisson statistic)*/
  if (verbose_flag)
    printf("generate coalescence times\n");
  
  randParams1 = XLALCreateRandomParams(seed1);
  durationMax = 1.26*tau0; /* for a 1.-1. binary at z=0 */
  
  /* calculate average time interval between merger times
  (integrating merger rate up to zMax */
  logzMax=log10(zMax);
  logdelta = -0.039563*pow(logzMax,6.)-0.15282*pow(logzMax,5.)-0.017596*pow(logzMax,4.)
           + 0.67193*pow(logzMax,3.)+1.1347*pow(logzMax,2.)-2.3543*logzMax+ 2.0228;
  delta = pow(10.,logdelta) /localRate;
  
  /* draw coalescence times */
  i = -1;
  tm = (REAL8)startTime; /* start at first coalescence; */
  while (tm<=stopTimeSeg+durationMax)
  {
    dt = -delta*log(XLALUniformDeviate(randParams1));
    tm += dt;
    if (tm>=(REAL8)startTimeSeg)
    {
      i++;
      if (i==0)
      {
        tc = XLALCreateREAL8Vector(1);
      }
      else
      {
        tc = XLALResizeREAL8Vector(tc,i+1);
      }
      tc->data[i]=tm;
    }
    else
    {
      seed2++;
    }
  }
  N=i+1;
  if (verbose_flag)
  printf("generate %d sources\n", (UINT4)N);
  Ntot=0;

  /* generate sources */
  for (i=0;i<N;i++)
  {
    if((i%100)==0)
    {
      if (verbose_flag)
      printf("source %d... select parameters...\n",i);
    }
    
    randParams2 = XLALCreateRandomParams(seed2+i);
    
    /* draw redshift */
    do
    {
      z = (zMax-zMin)*(REAL8)XLALUniformDeviate(randParams2)+zMin;
      reject = PzMax*(REAL8)XLALUniformDeviate(randParams2);
      p = -0.000429072589677 + (z * (-0.036349728568888 + (z * (0.860892111762314
        + (z * (-0.740935488674010 + z * (0.265848831356864 + z * (-0.050041573542298
        + z * (0.005184554232421 + z * (-0.000281450045300 + z * 0.000006400690921))))))))));
    }
    while (reject>p);
	if(z<zMin){break;}
	  
    /* draw masses between 1.2-3, gaussianly distributed with mu=1.4 and sigma=0.53 */
    do
    {
      params.mass1 = (mMax-mMin)*(REAL8)XLALUniformDeviate(randParams2)+mMin;
      reject = (REAL8)XLALUniformDeviate(randParams2);
      p = exp(-1./(2.*mSigma*mSigma)*(params.mass1-mMean)*(params.mass1-mMean));
    }
    while (reject>p);
    params.mass1*=(1.+z);
    
    do
    {
      params.mass2 = (mMax-mMin)*(REAL8)XLALUniformDeviate(randParams2)+mMin;
      reject = (REAL8)XLALUniformDeviate(randParams2);
	  p = exp(-1./(2.*mSigma*mSigma)*(params.mass2-mMean)*(params.mass2-mMean));
    }
    while (reject>p);
    params.mass2*=(1.+z);
    
    /* calculate duration of the waveform and check if it is in our observation window */
    LALInspiralParameterCalc(&status, &params);
    duration = params.tC;
    //printf("%e | %e | %e | %d \n", tc->data[i], (REAL8)stopTimeSeg,params.tC, tc->data[i]<=stopTimeSeg+params.tC );
    if(tc->data[i]<=stopTimeSeg+params.tC)
    {
      Ntot++;
      /* select sky position in equatorial coordinates */
      ra = Pi2*(REAL8)XLALUniformDeviate(randParams2)-Pi;
      cosDecl=(REAL8)XLALUniformDeviate(randParams2);
      if (XLALUniformDeviate(randParams2)<0.5)
      {
        decl=-acos(cosDecl);
      }
      else
      {
        decl=acos(cosDecl);
      }
      
      
      /* select inclination */
      cosI=(REAL8)XLALUniformDeviate(randParams2)*2.-1.;
      
      /* select polarization */
      psi= Pi2*(REAL8)XLALUniformDeviate(randParams2);
      
      /* select phase at coalescence */
      params.startPhase = Pi2*(REAL8)XLALUniformDeviate(randParams2);
      
      /* write parameters to file */
      if (catalog_flag)
      fprintf(pfout1,"%d % d %f %f %f %f %f %f %f %f %f %f\n",
        i,(int)seed2+i,tc->data[i],duration, z,ra,decl,cosI,psi,params.startPhase,params.mass1,params.mass2);
      
      /* fit of luminosity distance in Mpc for h0=0.7, omega_m=0.3, omega_v=0.7 */
      params.distance = -2.89287707063171+(z*(4324.33492012756+(z*(3249.74193862773
                      +(z*(-1246.66339928289+z*(335.354613407693+z*(-56.1194965448065
                     +z*(5.20261234121263+z*(-0.203151569744028))))))))));
      
      /* calculate contribution of the source to the time serie */
      if((i%100)==0)
      {
        if (verbose_flag)
        printf("calculate waveform...\n");
      }
      
      /* Estimate the length of signal vector */
      LALInspiralWaveLength( &status, &length, params);
      
      length = (UINT4) duration*rate;
      
      hp = XLALCreateREAL4Vector(length);
      hc = XLALCreateREAL4Vector(length);

      memset( hp->data, 0, length * sizeof(REAL4) );
      memset( hc->data, 0, length * sizeof(REAL4) );
      
      /* Amplitude is not needed since distance is specified */
      /* But we do need to convert distance into SI units */
      params.distance = params.distance * 1.0e6 * LAL_PC_SI;
      /*LALInspiralRestrictedAmplitude(&status, &params);*/
      LALInspiralWaveTemplates(&status,hp, hc, &params);

      /* Now tC will have been replaced with the actual length */
      
      t0 = tc->data[i] - params.tC;
      gpsT0.gpsSeconds = floor(t0);
      gpsT0.gpsNanoSeconds = (UINT8)((t0-gpsT0.gpsSeconds)*1000000000);
      
      apFac = -(1.+cosI*cosI)/2.;
      acFac = - cosI;
      for (k=0;k<length;k++)
      {
        hp->data[k] *=apFac;
        hc->data[k] *=acFac;
      }
      /* first detector */
      //printf("%ld %f %d\n",length,t0,gpsT0.gpsSeconds);
      //printf("%f %f %f\n",ra,decl,psi);
      /* problem here: segmentation fault happens in some cases*/
      if ( XLALComputeDetAMResponseSeries(&Fp,&Fc,dE1.response,
             ra,decl,psi,&gpsT0,deltat,length) == XLAL_FAILURE )
      {
	    fprintf( stderr, "Error in response function\n" );
        exit( 1 );
      }
      
      h = XLALCreateREAL4Vector(length);
      if ( !h )
      {
        fprintf( stderr, "Memory allocation error!\n" );
        exit( 1 );
      }
      for (k=0;k<length;k++)
      {
        h->data[k] = hp->data[k]*(Fp)->data->data[k]+hc->data[k]*(Fc)->data->data[k];
        //printf("%e \n", h->data[k]);
      }
      
      if(t0>=startTimeSeg)
      {
        if (t0 < stopTimeSeg)
        {
          kstart = floor((t0-startTimeSeg)*rate);
          indexMax = length;
          if ( kstart + length > serie1.data->length )
          {
            indexMax = serie1.data->length - kstart;
          }
          for (k=0;k<indexMax;k++)
          {
            serie1.data->data[k+kstart] += h->data[k];
          }
        }
      }
      else
      {
        kstart = ((startTimeSeg-t0)*rate);
        if ( kstart < length )
        {
          if ( length - kstart  > serie1.data->length )
          {
            indexMax = serie1.data->length;
          }
          else
          {
            indexMax = h->length - kstart;
          }
          for (k=0;k<indexMax;k++)
          {
            serie1.data->data[k] += h->data[k+kstart];
          }
        }
      }
      if (DEBUG) {
        sprintf(waveout,"%s%d%s","gw_",i,".txt");
        FILE *hout=fopen(waveout,"w");
        PrintHT(h, &hout, t0);
        fclose(hout);
      }
      if (Ndet>1)
      {
        /* second detector */
        XLALComputeDetAMResponseSeries(&Fp,&Fc,dE2.response,
                     ra,decl,psi,&gpsT0,deltat,length);

        for (k=0;k<length;k++)
        {
          h->data[k] = hp->data[k]*(Fp)->data->data[k]+hc->data[k]*(Fc)->data->data[k];
        }
        
        if(t0>=startTimeSeg )
        {
          if ( t0 < stopTimeSeg )
          {
            for (k=0;k<indexMax;k++)
            {
              serie2.data->data[k+kstart] += h->data[k];
            }
          }
        }
        else if ( kstart < length )
        {
          for (k=0;k<indexMax;k++)
          {
            serie2.data->data[k] += h->data[k+kstart];
          }
        }
        XLALDestroyREAL4TimeSeries( Fp );
        XLALDestroyREAL4TimeSeries( Fc );
      }
      
      if (Ndet>2)
      {
        /* third detector */
        XLALComputeDetAMResponseSeries(&Fp,&Fc,dE3.response,
                        ra,decl,psi,&gpsT0,deltat,length);

        for (k=0;k<length;k++)
        {
          h->data[k] = hp->data[k]*(Fp)->data->data[k]+hc->data[k]*(Fc)->data->data[k];
        }
        if(t0>=startTimeSeg )
        {
          if ( t0 < stopTimeSeg )
          {
            for (k=0;k<indexMax;k++)
              serie3.data->data[k+kstart] += h->data[k];
          }
        }
        else if (kstart < length )
        {
          for (k=0;k<indexMax;k++)
            serie3.data->data[k] += h->data[k+kstart];
        }
        XLALDestroyREAL4TimeSeries( Fp );
        XLALDestroyREAL4TimeSeries( Fc );
      }
      XLALDestroyREAL4Vector(hp);
      XLALDestroyREAL4Vector(hc);
      XLALDestroyREAL4Vector(h);
      
    }
    XLALDestroyRandomParams( randParams2 );
	
	
  }

  	REAL8Vector * in1 = NULL;
	in1 = XLALCreateREAL8Vector(serie1.data->length);
	
	FILE *ht_out = fopen("ht.dat", "w");
    for(k=0;k<serie1.data->length;k++){
 		in1->data[k] = serie1.data->data[k];
 		fprintf(ht_out, "%e\t%e\n", k/((REAL8)rate), serie1.data->data[k]);
 		}
    if (DEBUG){
        REAL8FFTPlan * fftforwardplan = XLALCreateForwardREAL8FFTPlan(serie1.data->length, 0);
        COMPLEX16Vector * fftsinglewaveform = XLALCreateCOMPLEX16Vector(serie1.data->length/2+1);
        XLALREAL8ForwardFFT(fftsinglewaveform, in1, fftforwardplan);

        FILE *hf_out = fopen("hf.dat", "w");    
        for(k=0;k<fftsinglewaveform->length;k++){
            fprintf(hf_out, "%e\t%e\t%e\n", k*((REAL8)rate)/((REAL8)serie1.data->length), fftsinglewaveform->data[k].re, fftsinglewaveform->data[k].im);
            }
        XLALDestroyREAL8FFTPlan(fftforwardplan);
    }

  if (noise_flag)
  {
    gsl_rng *r;
    const gsl_rng_type * T;
    REAL8Vector *noise=NULL;
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);
    gsl_rng_set( r, noiseSeed );
    
    
    noise = XLALCreateREAL8Vector(Nseg);
    
    ET_noise(r, Tseg, Nseg, fMin, rate, noise->data);
    for(k=0;k<Nseg;k++)
    {
      serie1.data->data[k] += (REAL4) noise->data[k];
    }    
    
    if (Ndet>1)
    {
    ET_noise(r, Tseg, Nseg, fMin, rate, noise->data);
      for(k=0;k<Nseg;k++)
        serie2.data->data[k] += (REAL4) noise->data[k];
    }
    
    if (Ndet>2)
    {
    ET_noise(r, Tseg, Nseg, fMin, rate, noise->data);
      for(k=0;k<Nseg;k++)
        serie3.data->data[k] += (REAL4) noise->data[k];
    }
    
    XLALDestroyREAL8Vector(noise);
  }
  
  if (null_flag)
  {
	for(k=0;k<Nseg;k++)
	 serieNS.data->data[k] = serie1.data->data[k]+serie2.data->data[k]+serie3.data->data[k];
  }	
    if (DEBUG){
        REAL8Vector * in2 = NULL;
        in2 = XLALCreateREAL8Vector(serie1.data->length);
	
        FILE *serieT_out = fopen("serieT.dat", "w");
        for(k=0;k<serie1.data->length;k++){
            in2->data[k] = serieNS.data->data[k];
            fprintf(serieT_out, "%e\t%e\n", k/((REAL8)rate), serieNS.data->data[k]);
            }
        fclose(serieT_out);
        REAL8FFTPlan * fftforwardplan = XLALCreateForwardREAL8FFTPlan(serie1.data->length, 0);
        COMPLEX16Vector * fftsinglewaveform = XLALCreateCOMPLEX16Vector(serie1.data->length/2+1);
        XLALREAL8ForwardFFT(fftsinglewaveform, in2, fftforwardplan);

        FILE *serieF_out = fopen("serieF.dat", "w");
        for(k=0;k<fftsinglewaveform->length;k++){
            fprintf(serieF_out, "%e\t%e\t%e\n", k*((REAL8)rate)/((REAL8)serie1.data->length), fftsinglewaveform->data[k].re, fftsinglewaveform->data[k].im);
            }
        XLALDestroyREAL8FFTPlan(fftforwardplan);
        fprintf(stderr,"Computing Inverse PSD\n");
        InvPSD(serieNS, &serie1.epoch, Tobs, rate, 1);
        fprintf(stderr,"done\n");
        fclose(serieF_out);     
    }

  if (verbose_flag)
  printf("%d sources... write time serie to file\n",(UINT4)Ntot);
  /* write to file */
  if(ascii_flag)
  {
    for (k=0;k<Nseg;k++)
    {
      if (Ndet==1)
        fprintf(pfout2,"%f %e\n",
             k/rate,serie1.data->data[k]);
      else if (Ndet==2)
        fprintf(pfout2,"%f %e %e\n",
             k/rate,serie1.data->data[k],serie2.data->data[k]);
      else if (Ndet==3)
	  {
		if (null_flag)
		 fprintf(pfout2,"%f %e %e %e %e\n", 
		    k/rate,serie1.data->data[k],serie2.data->data[k],serie3.data->data[k],serieNS.data->data[k]);
		else
		 fprintf(pfout2,"%f %e %e %e\n", k/rate,serie1.data->data[k],serie2.data->data[k],serie3.data->data[k]);
	  }
	}}
  else
  {
      LALFrWriteREAL8TimeSeries(&status,&serie1,&frSerie1);
    if (Ndet>1)
      LALFrWriteREAL8TimeSeries(&status,&serie2,&frSerie2);
    if (Ndet>2)
      LALFrWriteREAL8TimeSeries(&status,&serie3,&frSerie3);
	if (null_flag)
	  LALFrWriteREAL8TimeSeries(&status,&serieNS,&frSerieNS);   
  }
	
  if(stat_flag)
  {	
   meangw= gsl_stats_mean(serie1.data->data, 1, serie1.data->length);
   vargw = gsl_stats_variance(serie1.data->data, 1, serie1.data->length);
   skewgw= gsl_stats_skew(serie1.data->data, 1, serie1.data->length);
   kurtgw = gsl_stats_kurtosis(serie1.data->data, 1, serie1.data->length);
   fprintf(pfout3,"%e %e %e %e\n",meangw,vargw,skewgw,kurtgw);
   if(verbose_flag)
   {
    printf("mean=%e variance=%e skewness=%e kurtosis=%e\n",meangw,vargw,skewgw,kurtgw);	  
   }
	   
  }
	
  if (verbose_flag)
    printf("clean up and exit\n");
  /* clean up */
  LALDDestroyVector(&status, &tc);
  XLALDestroyREAL8Vector(serie1.data);
  if (Ndet>1)
    XLALDestroyREAL8Vector(serie2.data);
  if (Ndet>2)
    XLALDestroyREAL8Vector(serie3.data);

  XLALDestroyRandomParams(randParams1);
  if (catalog_flag) 
  {
    fclose(pfout1); 
  }
  if (ascii_flag) 
  {
    fclose(pfout2); 
  }
  if (stat_flag) 
  {
   fclose(pfout3); 
  }
	
  
  LALCheckMemoryLeaks();
  return(0);
}


/* parse command line options */
void parseOptions(INT4 argc, CHAR *argv[])
{
  int c = -1;
  
  while(1)
  {
    static struct option long_options[] =
    {
      /* options that set a flag */
      {"verbose", no_argument, &verbose_flag, 1},
      {"ascii", no_argument, &ascii_flag, 1},
      {"catalog", no_argument, &catalog_flag, 1},
	  {"stat", no_argument, &stat_flag, 1},	
      {"noise", no_argument, &noise_flag, 1},
	  {"null", no_argument, &null_flag, 1},	
	  {"psd", no_argument, &write_psd_flag, 1},	
      /* options that don't set a flag */
      {"help", no_argument, 0, 'h'},
      {"seed1", required_argument, 0, 's'},
      {"seed2", required_argument, 0, 'S'},
      {"noise-seed", required_argument, 0, 'N'},
      {"job-number", required_argument, 0, 'j'},
	  {"start-job-number", required_argument, 0, 'J'},
      {"nodes", required_argument, 0, 'n'},
      {"number-detectors", required_argument, 0, 'D'},
      {"start-time", required_argument, 0, 't'},
      {"duration", required_argument, 0, 'd'},
      {"sampling-rate", required_argument, 0, 'r'},
      {"approximant", required_argument, 0, 'a'},
      {"phase-order", required_argument, 0, 'o'},
      {"amplitude-order", required_argument, 0, 'O'},
      {"frequency-min", required_argument, 0, 'f'},
      {"z-min", required_argument, 0, 'z'},
	  {"z-max", required_argument, 0, 'Z'},	
      {"local-rate", required_argument, 0, 'l'},
	  {"output-dir", required_argument, 0, 'R'},
	  {"mass-min", required_argument, 0, 'm'},
	  {"mass-max", required_argument, 0, 'M'},
	  {"mass-mean", required_argument, 0, 'p'},
	  {"mass-sigma", required_argument, 0, 'P'},
      {"version", no_argument, 0, 'v'},
      {0, 0, 0, 0}
    };
    
    /* getopt_long stores the option here */
    int option_index = 0;
    
    c = getopt_long(argc, argv,
    "hs:S:N:j:J:n:D:t:d:r:a:o:O:f:z:Z:l:R:m:M:p:P:v:",
    long_options, &option_index);
    
    if (c == -1)
    {
      /* end of options, break loop */
      break;
    }
    
    switch(c)
    {
      case 0:
      /* If this option set a flag, do nothing else now. */
      if (long_options[option_index].flag != 0)
      break;
      printf ("option %s", long_options[option_index].name);
      if (optarg)
      printf (" with arg %s", optarg);
      printf ("\n");
      break;
      
      case 'h':
      /* HELP!!! */
      displayUsage(0);
      break;
      
      case 's':
      /* seed1 */
      seed1 = atoi(optarg);
      break;
      
      case 'S':
      /* seed2 */
      seed2 = atoi(optarg);
      break;
      
      case 'N':
      /* noise seed */
      noiseSeed = atoi(optarg);
      break;
      
      case 'j':
      /* job number */
      job = atoi(optarg);
      break;
			
	  case 'J':
	  /* start job number */
	  jobStart = atoi(optarg);
	  break;
      
      case 'n':
      /* number of nodes */
      Njob = atoi(optarg);
      break;
      
      case 'D':
      /* number of detectors */
      Ndet = atoi(optarg);
      break;
      
      case 't':
      /* start time */
      startTime = atoi(optarg);
      break;
      
      case 'd':
      /* duration */
      Tobs = atof(optarg);
      break;
      
      case 'r':
      /* sampling rate */
      rate = atof(optarg);
      break;
      
      //case 'a':
      /* approximant */
      //	  approximant = atoi(optarg);
      //	  break;
      
      case 'o':
      /* phase order */
      order = atoi(optarg);
      break;
      
      case 'O':
      /* amplitude order */
      ampOrder = atoi(optarg);
      break;
      
      case 'f':
      /* minimal frequency */
      fMin = atof(optarg);
      break;
      
      case 'z':
      /* minimal redshift */
      zMin = atof(optarg);
      break;
			
	  case 'Z':
	  /* maximal redshift */
	  zMax = atof(optarg);
	  break;
      
      case 'l':
      /* local rate */
      localRate = atof(optarg);
      break;
      
      case 'm':
      /* minimal mass */
      mMin = atof(optarg);
      break;
			
      case 'M':
      /* maximal mass */
      mMax = atof(optarg);
      break;

      case 'p':
      /* mass mean */
      mMean = atof(optarg);
      break;

      case 'P':
      /* mass std */
      mSigma = atof(optarg);
      break;
	
	  case 'R':
	  /* directory for output files */
	  strncpy(outputDir, optarg, LALNameLength);
	  break;		
			
      case 'v':
      /* display version info and exit */
      fprintf(stdout, "simulation of extragalactic BNS\n" CVS_ID "\n");
      exit(0);
      break;
      
      default:
      displayUsage(1);
    }
  }
  
  if (optind < argc)
  {
    displayUsage(1);
  }
  
  if (noiseSeed == 0 )
  {
    fprintf( stderr, "--noise-seed must be set.\n" );
    displayUsage(1);
  }
  
  return;
}

/* display program usage */
void displayUsage(INT4 exitcode)
{
  fprintf(stderr, "Usage: pipeline [options]\n");
  fprintf(stderr, "Options:\n");
  fprintf(stderr, " -h                    print this message\n");
  
  fprintf(stderr, " --verbose             verbose mode\n");
  fprintf(stderr, " --ascii               write to ascii files\n");
  fprintf(stderr, " --catalog             write source parameters to files\n");
  fprintf(stderr, " --stat                calculate moments of the GW signal\n");	
  fprintf(stderr, " --psd                 write psd \n");
  fprintf(stderr, " --noise               write noise realization\n");\
  fprintf(stderr, " --null                calculate null stream\n");
  fprintf(stderr, " -s                    seed for coalescence times\n");
  fprintf(stderr, " -S                    initial seed for source parameters\n");
  fprintf(stderr, " -N                    seed for noise generation\n");
  fprintf(stderr, " -j                    job number\n");
  fprintf(stderr, " -J                    job number\n");
  fprintf(stderr, " -n                    number of nodes\n");
  fprintf(stderr, " -D                    number of detectors (1-3)\n");
  fprintf(stderr, " -t                    start time of the serie\n");
  fprintf(stderr, " -d                    duration of the time serie\n");
  fprintf(stderr, " -r                    sampling rate of the time serie\n");
  fprintf(stderr, " -a                    approximant\n");
  fprintf(stderr, " -o                    phase order\n");
  fprintf(stderr, " -O                    amplitude order\n");
  fprintf(stderr, " -f                    minimal frequency\n");
  fprintf(stderr, " -z                    maximal redshift\n");
  fprintf(stderr, " -Z                    minimal redshift\n");
  fprintf(stderr, " -l                    local rate per Myr^3 per Mpc (default 1)\n");
  fprintf(stderr, " -m                    minimal mass\n");
  fprintf(stderr, " -M                    maximal mass\n");
  fprintf(stderr, " -p                    Mean of the mass distribution\n");
  fprintf(stderr, " -P                    Std of the mass distribution\n");
  fprintf(stderr, " -R                    directory for output files\n");	
  fprintf(stderr, " -v                    display version\n");	
  exit(exitcode);
}
/*
ET Mock Data Challenge.
Signals are injected into colored gaussian noise generated
from the proposed ET sensitivity curve

Version 1.0 (August 3rd, 2009)

Author: Carl Rodriguez
*/

int ET_noise(gsl_rng *r, double TIME, long int N, double F_START, double SAMPLE_FREQ, double *noise_data)
{
	int i, Nby2;
	double Sn, rtS, S, f, df;
	FILE *f1 = NULL;

	/* To avoid artefacts, we will taper the start and end of the PSD */
	/* At the start, it is tapered from zero to lowTaper */
	/* At the end, it is tapered from highTaper to Nyquist */
	int lowTaper;
	int highTaper;
	double lowTaperGradient;
	double highTaperGradient;
	double sHigh;   /* S at the frequency beyond which we taper */
	const double F_HIGH = 2048.0; /*Freq beyond which we taper the PSD */
	
	const REAL8 c1 = 2.39e-27;
	const REAL8 c2 = 0.349;
	const REAL8 c3 = 1.76;
	const REAL8 c4 = 0.409;	
	const REAL8 p1 = -15.64;
	const REAL8 p2 = -2.145;
	const REAL8 p3 = -0.12;
	const REAL8 p4 = 1.10;	
	
	REAL8 xt;
	
	if(write_psd_flag)	
		f1 = fopen("psd.dat", "w");
	Nby2  = N/2;
	df = 1./TIME;
	
	noise_data[0] = 0.0;
	
	/* Calculate indices for tapering */
	lowTaper  = F_START / df;
	highTaper = F_HIGH / df;
	
	/* Calculate the gradients */
	xt = F_START/100.;
	rtS= 1.0e-25 * (c1 * pow( xt, p1 )+ c2 * pow( xt, p2 )+ c3 * pow( xt, p3 )+ c4 * pow( xt, p4 ));	
	Sn = rtS*rtS;
	lowTaperGradient = rtS / (double)lowTaper;
	
	xt = F_HIGH/100.;
	sHigh= 1.0e-25 * (c1 * pow( xt, p1 )+ c2 * pow( xt, p2 )+ c3 * pow( xt, p3 )+ c4 * pow( xt, p4 ));	
	Sn = sHigh*sHigh;	
	highTaperGradient = - sHigh / (double)(Nby2 - highTaper);
	
	for ( i = 1; i <=  Nby2; i++ )
	{
		f = i*df;
		noise_data[i] = 0.;
		noise_data[N-i] = 0.;
		
		if ( f < F_START )
		{
			rtS = (double)i * lowTaperGradient;
		}
		else if (f >= F_START && i <= highTaper )
		{
			xt = f/100.;
			rtS= 1.0e-25 * (c1 * pow( xt, p1 )+ c2 * pow( xt, p2 )+ c3 * pow( xt, p3 )+ c4 * pow( xt, p4 ));		
			Sn = rtS*rtS;
			if(write_psd_flag)	
				fprintf(f1, "%e %e\n", f, rtS);
		}
		else if ( f < SAMPLE_FREQ / 2 )
		{
			rtS = sHigh + (double)( i - highTaper ) * highTaperGradient;
		}
		S = rtS * sqrt(N*SAMPLE_FREQ/4.);
		noise_data[i]   = gsl_ran_gaussian(r, 1.) * S;
		noise_data[N-i] = gsl_ran_gaussian(r, 1.) * S;
	}
	
	gsl_fft_halfcomplex_radix2_inverse(noise_data, 1, N);
	if(write_psd_flag)
		fclose(f1);
	
	return 0;
}

int LIGOI_noise(gsl_rng *r, double TIME, long int N, double F_START, double SAMPLE_FREQ, double *noise_data)
{
  int i;
  double S0, Sn, S, f, x;
  double rtS;
  double df;
  int Nby2;
  
  /* To avoid artefacts, we will taper the start and end of the PSD */
  /* At the start, it is tapered from zero to lowTaper */
  /* At the end, it is tapered from highTaper to Nyquist */
  
  int lowTaper;
  int highTaper;
  double lowTaperGradient;
  double highTaperGradient;
  double sHigh;   /* S at the frequency beyond which we taper */
  const double F_HIGH = 2048.0; /*Freq beyond which we taper the PSD */
  
  
  FILE *f1 = NULL;
  
  if(write_psd_flag)
   f1 = fopen("psd.dat", "w");
  Nby2  = N/2;
  df = 1./TIME;
  S0 = 9e-46;
  
  noise_data[0] = 0.0;
  
  /* Calculate indices for tapering */
  lowTaper  = F_START / df;
  highTaper = F_HIGH / df;
  
  /* Calculate the gradients */
  x = lowTaper * df / 150.0;
  Sn = S0 * (pow(4.49 * x, -56) + 0.16 * pow(x, -4.52) + 0.52 + 0.32 * pow(x, 2));
  rtS = sqrt(Sn);
  lowTaperGradient = rtS / (double)lowTaper;
  
  x     = highTaper * df / 150.0;
  Sn    = S0 * (pow(4.49 * x, -56) + 0.16 * pow(x, -4.52) + 0.52 + 0.32 * pow(x, 2));
  sHigh   = sqrt(Sn);
  highTaperGradient = - sHigh / (double)(Nby2 - highTaper);
  
  for ( i = 1; i <=  Nby2; i++ )
  {
    f = i*df;
    noise_data[i] = 0.;
    noise_data[N-i] = 0.;
    
    if ( f < F_START )
    {
      rtS = (double)i * lowTaperGradient;
    }
    else if (f >= F_START && i <= highTaper )
    {
      x = f / 150.;
      Sn = S0 * (pow(4.49 * x, -56) + 0.16 * pow(x, -4.52) + 0.52 + 0.32 * pow(x, 2));
      rtS = sqrt(Sn);
	  if(write_psd_flag)	
       fprintf(f1, "%e %e\n", f, rtS);
    }
    else if ( f < SAMPLE_FREQ / 2 )
    {
      rtS = sHigh + (double)( i - highTaper ) * highTaperGradient;
    }
    S = rtS * sqrt(N*SAMPLE_FREQ/4.);
    noise_data[i]   = gsl_ran_gaussian(r, 1.) * S;
    noise_data[N-i] = gsl_ran_gaussian(r, 1.) * S;
  }
  
  
  gsl_fft_halfcomplex_radix2_inverse(noise_data, 1, N);
  if(write_psd_flag)	
   fclose(f1);
  
  return 0;
}

void InvPSD(REAL8TimeSeries Data, LIGOTimeGPS *GPSStart, REAL8 duration, UINT4 SampleRate, UINT4 nSegs)
{

    int check = 0;
    UINT4 j=0;
    
    REAL8 segDur = duration/(REAL8)nSegs;
    UINT4 seglen = (UINT4)(segDur*SampleRate); // GET THE LENGTH OF A SEGMENT
    segDur = seglen/SampleRate; // UPDATE WITH THE TRUNCATED SEGLEN
  	nSegs = (INT4)floor(duration/segDur); // UPDATE NSEGS TO ACCOMODATE SAMPLERATE
    UINT4 stride = seglen; /* Overlap the padding */
    //REAL8 strideDur = stride / SampleRate;
    REAL8 deltaF=(REAL8)SampleRate/seglen;
    REAL8 end_freq=10.0; /* cutoff frequency */ 
    
    REAL8Window  *windowplan = XLALCreateTukeyREAL8Window(seglen,0.1*(REAL8)8.0*SampleRate/(REAL8)seglen);;
    REAL8FFTPlan *fwdplan = XLALCreateForwardREAL8FFTPlan( seglen, 1 );
    REAL8FFTPlan *revplan = XLALCreateReverseREAL8FFTPlan( seglen, 1 );
    REAL8FrequencySeries *inverse_spectrum = XLALCreateREAL8FrequencySeries("inverse spectrum",GPSStart,0.0,deltaF,&lalDimensionlessUnit,seglen/2+1);
    
    fprintf(stdout,"... Computing power spectrum, seglen %i stride %i\n",seglen,stride);
    
    // CALCULATE THE INVERSE SPECTRUM 
		check = XLALREAL8AverageSpectrumMedian(inverse_spectrum,&Data,(UINT4)seglen,(UINT4)stride,windowplan,fwdplan);
    if (check) {fprintf(stderr,"Failed computing the XLALREAL8AverageSpectrumMedian \n");exit(-1);}
    
    // TRUNCATE INVERSE SPECTRUM
    check = XLALREAL8SpectrumInvertTruncate(inverse_spectrum, end_freq, seglen, (seglen-stride)/4, fwdplan, revplan ); 
    if (check) {fprintf(stderr,"Failed computing XLALREAL8SpectrumInvertTruncate \n");exit(-1);}

    FILE *psdout;
    psdout=fopen("inv_spec.dat","w");
    
    // PRINT FILE
    for(j=0;j<inverse_spectrum->data->length;j++) {
       fprintf(psdout,"%10.10lf %10.10e\n",j*deltaF,inverse_spectrum->data->data[j]); 
    }
    fclose(psdout);

    return;
}

void PrintHT(REAL4Vector *s, FILE **f, REAL8 t0){
    UINT4 k=0;
	t0-=800000000.;
    for(k=0;k<s->length;k++){
 		fprintf(*f, "%e\t%e\n", t0+k/((REAL8)rate), s->data[k]);
    }
    return;
}