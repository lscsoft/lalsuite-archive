#include <stdio.h>
#include <complex.h>
#include <math.h>
#include <string.h>

#include <gsl/gsl_const.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_odeiv.h>

#include <lal/LALSimInspiral.h>
#include <lal/LALSimIMR.h>
#include <lal/LALConstants.h>
#include <lal/LALDict.h>
#include <lal/LALStdlib.h>
#include <lal/LALInference.h>
#include <lal/LALInferenceReadData.h>
#include <lal/LALInferenceLikelihood.h>
#include <lal/Sequence.h>
#include <lal/TimeSeries.h>
#include <lal/FrequencySeries.h>
#include <lal/TimeFreqFFT.h>
#include <lal/RealFFT.h>
#include <lal/Units.h>
#include <lal/SphericalHarmonics.h>
#include <lal/LIGOMetadataTables.h>

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

#define APPROX_NAME_SIZE 64
#define FILE_NAME_SIZE 256
#define PARAM_NAME_SIZE 32

typedef struct tagParameter 
{
  REAL8 m1;
  REAL8 m2;
  REAL8 d;
  REAL8 iota;
  REAL8 pol;
  REAL8 orb;
  REAL8 coal;
  REAL8 ra;
  REAL8 dec;
  REAL8 s1x;
  REAL8 s1y;
  REAL8 s1z;
  REAL8 s2x;
  REAL8 s2y;
  REAL8 s2z;
  REAL8 ecc;
  REAL8 f_ecc;
  INT4 ecc_order;
  char approx1[APPROX_NAME_SIZE];
  char approx2[APPROX_NAME_SIZE];
  INT4 ampOrder1;
  INT4 ampOrder2;
  INT4 phaseOrder1;
  INT4 phaseOrder2;
  char outfile[FILE_NAME_SIZE];
  INT4 dimension; // 1d or 2d
  char p1n[PARAM_NAME_SIZE];
  REAL8 p1s;
  REAL8 p1e;
  REAL8 p1d;
  char p2n[PARAM_NAME_SIZE];
  REAL8 p2s;
  REAL8 p2e;
  REAL8 p2d;
} Parameter;

int parseParameters(ProcessParamsTable *procParams, Parameter *inParams);
int printParameters(FILE *out, Parameter *inParams);
REAL8 LALInferenceComputeFrequencyDomainOverlapNormalised(LALInferenceIFOData * dataPtr, 
                                   COMPLEX16Vector * freqData1, 
                                   COMPLEX16Vector * freqData2);
void generateGNUfile(Parameter *inParams);

int main(int argc, char *argv[])
/*===========================================================================================
/ overlap calculation for two template wave function
/ Parameters
/  --m1 : central mass 1 in M_sun
/  --m2 : central mass 2 in M_sun
/  --d  : distance in Mpc
/  --iota : inclination in radian
/  --pol : polarization in radian
/  --orb : orbital phase in radian
/  --coal : coalescence time in seconds
/  --ra : right ascension in radian
/  --dec : declination in radian
/  --s1x : x component of spin 1
/  --s1y : y component of spin 1
/  --s1z : z component of spin 1
/  --s2x : x component of spin 2
/  --s2y : y component of spin 2
/  --s2z : z component of spin 2
/  --ecc : ecc value
/  --f_ecc : f_ecc value
/  --ecc_order : ecc_order value
/  --approx1 : first approximant, waveform template
/  --approx2 : second approximant, waveform template it could be same as approx1
/  --ampOrder1 : first amplitude PN order (1 -> 0.5PN,  2 -> 1PN)
/  --ampOrder2 : second amplitude PN order (1 -> 0.5PN,  2 -> 1PN)
/  --phaseOrder1 : first phase PN order (1 -> 0.5PN,  2 -> 1PN)
/  --phaseOrder2 : second phase PN order (1 -> 0.5PN,  2 -> 1PN)
/  --outfile : output file name string
/  --1d : 1 dimensional graph (one parameter vary)
/  --2d : 2 dimensional graph (two parameters vary) 1d and 2d should be exclusive
/  --p1n : first parameter name
/  --p1s : first parameter start value
/  --p1e : first parameter end value
/  --p1d : first parameter delta value
/  --p2n : second parameter name
/  --p2s : second parameter start value
/  --p2e : second parameter end value
/  --p2d : second parameter delta value
=========================================================================================*/
{
  LALInferenceRunState *irs = NULL;
  //LALInferenceIFOData *ifoPtr, *ifoStart;
  COMPLEX16FrequencySeries *hptilde1=NULL;
  COMPLEX16FrequencySeries *hptilde2=NULL;
  COMPLEX16FrequencySeries *hctilde1=NULL;
  COMPLEX16FrequencySeries *hctilde2=NULL;
  //COMPLEX16 *data, *data1;
  ProcessParamsTable *procParams = NULL;
  //ProcessParamsTable *ppt = NULL;
  Parameter inParams ;
  int ret;
  REAL8 deltaF, phiRef = 0.0, m1=0.0, m2=0.0, S1z = 0.0, S2z = 0.0;
  REAL8 deltaT, f_min=0.0, f_max = 0.0, r=0.0;
  INT4 phaseO, amplitudeO, lower, upper, ecc_order;
  REAL8 ecc, f_ecc;
  char file_name[256];
  REAL8 overlap, overlapNorm, p1, p1v;
  //LALInferenceVariables currentParams;
  FILE *outf=NULL;
  LALDict *LALparams = XLALCreateDict();
  memset(&inParams, 0x00, sizeof(Parameter));
  /* Read command line and parse */
  //printf("main DEBUG -5\n");
  procParams = LALInferenceParseCommandLine(argc, argv);
  //printf("main DEBUG -4\n");
  parseParameters(procParams, &inParams);
  //printf("main DEBUG -3\n");
  printParameters(stdout, &inParams);
  //printf("main DEBUG -2\n");

  irs = XLALCalloc(1, sizeof(LALInferenceRunState));
  fprintf(stdout, "LALInferenceReadData(): started.\n");
  irs->commandLine = procParams;
  irs->data = LALInferenceReadData(irs->commandLine);
  fprintf(stdout, "LALInferenceReadData(): finished.\n");
  deltaT = irs->data->timeData->deltaT;
  deltaF = 1.0 / (((double)irs->data->timeData->data->length) * deltaT);
  lower = ceil(irs->data->fLow / deltaF);
  upper = floor(irs->data->fHigh / deltaF);
  fprintf(stdout, "deltaT = %g, deltaF = %g, length = %d\n", irs->data->timeData->deltaT, deltaF, irs->data->timeData->data->length);
  fprintf(stdout, "fLow = %g, fHigh = %g\n", irs->data->fLow, irs->data->fHigh);
  fprintf(stdout, "lower = %d, upper = %d\n", lower, upper);
  //for(i=lower; i <= upper; i++)
  //{
  //  printf("i=%d, psd = %g\n", i, irs->data->oneSidedNoisePowerSpectrum->data->data[i]);
  //} 
  m1 = inParams.m1*LAL_MSUN_SI; // convert to kg
  m2 = inParams.m2*LAL_MSUN_SI;
  f_min = irs->data->fLow;
  r = inParams.d*3.08567758e22; // convert to m
  phaseO = inParams.phaseOrder1;
  amplitudeO = inParams.ampOrder1;
  //inclination = inParams.iota;
  //S1x = inParams.s1x;
  //S1y = inParams.s1y;
  S1z = inParams.s1z;
  //S2x = inParams.s2x;
  //S2y = inParams.s2y;
  S2z = inParams.s2z;
  ecc = inParams.ecc;
  ecc_order = inParams.ecc_order;
  f_ecc = inParams.f_ecc;
  XLALSimInspiralWaveformParamsInsertPNPhaseOrder(LALparams, phaseO);
  XLALSimInspiralWaveformParamsInsertPNAmplitudeOrder(LALparams, amplitudeO);
  XLALSimInspiralWaveformParamsInsertPNEccentricityOrder(LALparams, ecc_order);
  XLALSimInspiralWaveformParamsInsertEccentricityFreq(LALparams, f_ecc);
  XLALSimInspiralWaveformParamsInsertTidalLambda1(LALparams, 0.);
  XLALSimInspiralWaveformParamsInsertTidalLambda2(LALparams, 0.);

  if(strstr(inParams.approx1, "TaylorF2Ecc"))
  {
    ret = XLALSimInspiralTaylorF2Ecc(&hptilde1, phiRef, deltaF, 
                    m1, m2, S1z, S2z,
                    f_min, f_max, 0, r, 1, 1, 0, 0, ecc, ecc_order, f_ecc, 0, 0, phaseO, amplitudeO, NULL); // qm1=1, qm2=1 hard fixed
    //ret = XLALSimInspiralTaylorF2AmpPlus(&hptilde1, phiRef, deltaF, inclination,
    //                m1, m2, S1x, S1y, S1z, S2x, S2y, S2z,
    //                LNhatx, LNhaty, LNhatz, f_min, f_max, r, phaseO, amplitudeO);
  }
  else if(strstr(inParams.approx1, "TaylorF2"))
  {
    ret = XLALSimInspiralTaylorF2(&hptilde1, phiRef, deltaF, 
                    m1, m2, S1z, S2z,
                    f_min, f_max, 0, r, LALparams); // qm1=1, qm2=1 hard fixed
  }
  else 
  {
    ret = XLALSimInspiralTaylorF2Ecc(&hptilde1, phiRef, deltaF, 
                    m1, m2, S1z, S2z,
                    f_min, f_max, 0, r, 1, 1, 0, 0, ecc, ecc_order, f_ecc, 0, 0, phaseO, amplitudeO, NULL); // qm1=1, qm2=1 hard fixed
    //ret = XLALSimInspiralTaylorF2AmpPlus(&hptilde1, phiRef, deltaF, inclination,
    //                m1, m2, S1x, S1y, S1z, S2x, S2y, S2z,
    //                LNhatx, LNhaty, LNhatz, f_min, f_max, r, phaseO, amplitudeO);
  }
  if ((unsigned int)upper > hptilde1->data->length)
  {
    upper = hptilde1->data->length - 1;
    irs->data->fHigh = upper * deltaF + 0.5*deltaF; //reduce fHigh value
  }
  fprintf(stdout, "hptilde1->data->length = %d\n", hptilde1->data->length);
  //for(i=lower; i <= upper; i++)
  //{
  //  printf("i=%d, hptilde1 = %g, %g\n", i, creal(hptilde1->data->data[i]), cimag(hptilde1->data->data[i]));
  //} 
  sprintf(file_name, "%s.dat", inParams.outfile);
  outf = fopen(file_name, "w");
  if(strstr(inParams.p1n, "m1")) // mass1 changed
  {
    for(p1 = inParams.p1s; p1 < inParams.p1e + 0.5*inParams.p1d; p1 += inParams.p1d)
    {
      p1v = p1*LAL_MSUN_SI; // convert to kg
      if(strstr(inParams.approx2, "TaylorF2Ecc"))
      {
        ret = XLALSimInspiralTaylorF2Ecc(&hptilde2, phiRef, deltaF, 
                    p1v, m2, S1z, S2z,
                    f_min, f_max, 0, r, 1, 1, 0, 0, ecc, ecc_order, f_ecc, 0, 0, phaseO, amplitudeO, NULL); // qm1=qm2=1 hard code
        //ret = XLALSimInspiralTaylorF2AmpPlus(&hptilde2, phiRef, deltaF, inclination,
        //            p1v, m2,S1x, S1y, S1z, S2x, S2y, S2z,
        //            LNhatx, LNhaty, LNhatz, f_min, f_max, r, phaseO, amplitudeO);
      }
      else if(strstr(inParams.approx2, "TaylorF2"))
      {
        ret = XLALSimInspiralTaylorF2(&hptilde2, phiRef, deltaF, 
                    p1v, m2, S1z, S2z,
                    f_min, f_max, 0, r, LALparams); // qm1=qm2=1 hard code
      }
      else
      {
        ret = XLALSimInspiralTaylorF2Ecc(&hptilde2, phiRef, deltaF, 
                    p1v, m2, S1z, S2z,
                    f_min, f_max, 0, r, 1, 1, 0, 0, ecc, ecc_order, f_ecc, 0, 0, phaseO, amplitudeO, NULL); // qm1=qm2=1 hard code
        //ret = XLALSimInspiralTaylorF2AmpPlus(&hptilde2, phiRef, deltaF, inclination,
        //            p1v, m2,S1x, S1y, S1z, S2x, S2y, S2z,
        //            LNhatx, LNhaty, LNhatz, f_min, f_max, r, phaseO, amplitudeO);
      }
      overlap = LALInferenceComputeFrequencyDomainOverlap(irs->data, hptilde1->data, hptilde2->data);
      overlapNorm = LALInferenceComputeFrequencyDomainOverlapNormalised(irs->data, hptilde1->data, hptilde2->data);
      //loglikelihood = LALInferenceFreqDomainLogLikelihood(&currentParams, runstate->data, LALInferenceTemplateXLALSimInspiralChooseWaveform));
      fprintf(outf, "%20.8e %20.8e %20.8e\n", p1, overlap, overlapNorm);
      if ( hptilde2 ) XLALDestroyCOMPLEX16FrequencySeries(hptilde2);
      hptilde2 = NULL;
    }
  }
  else if(strstr(inParams.p1n, "s1z")) // s1z changed
  {
    for(p1 = inParams.p1s; p1 < inParams.p1e + 0.5*inParams.p1d; p1 += inParams.p1d)
    {
      p1v = p1; 
      if(strstr(inParams.approx2, "TaylorF2Ecc"))
      {
        ret = XLALSimInspiralTaylorF2Ecc(&hptilde2, phiRef, deltaF, 
                    m1, m2, p1v, S2z,
                    f_min, f_max, 0, r, 1, 1, 0, 0, ecc, ecc_order, f_ecc, 0, 0, phaseO, amplitudeO, NULL); // qm1=qm2 = 1 hard code
        //ret = XLALSimInspiralTaylorF2AmpPlus(&hptilde2, phiRef, deltaF, inclination,
        //            m1, m2,S1x, S1y, p1v, S2x, S2y, S2z,
        //            LNhatx, LNhaty, LNhatz, f_min, f_max, r, phaseO, amplitudeO);
      }
      else if(strstr(inParams.approx2, "TaylorF2"))
      {
        ret = XLALSimInspiralTaylorF2(&hptilde2, phiRef, deltaF, 
                    m1, m2, p1v, S2z,
                    f_min, f_max, 0, r, LALparams); // qm1=qm2 = 1 hard code
      }
      else
      {
        ret = XLALSimInspiralTaylorF2Ecc(&hptilde2, phiRef, deltaF, 
                    m1, m2, p1v, S2z,
                    f_min, f_max, 0, r, 1, 1, 0, 0, ecc, ecc_order, f_ecc, 0, 0, phaseO, amplitudeO, NULL); // qm1=qm2 = 1 hard code
        //ret = XLALSimInspiralTaylorF2AmpPlus(&hptilde2, phiRef, deltaF, inclination,
        //            m1, m2,S1x, S1y, p1v, S2x, S2y, S2z,
        //            LNhatx, LNhaty, LNhatz, f_min, f_max, r, phaseO, amplitudeO);
      }
      overlap = LALInferenceComputeFrequencyDomainOverlap(irs->data, hptilde1->data, hptilde2->data);
      overlapNorm = LALInferenceComputeFrequencyDomainOverlapNormalised(irs->data, hptilde1->data, hptilde2->data);
      //loglikelihood = LALInferenceFreqDomainLogLikelihood(&currentParams, runstate->data, LALInferenceTemplateXLALSimInspiralChooseWaveform));
      fprintf(outf, "%20.8e %20.8e %20.8e\n", p1, overlap, overlapNorm);
      if ( hptilde2 ) XLALDestroyCOMPLEX16FrequencySeries(hptilde2);
      hptilde2 = NULL;
    }
  }
  else if(strstr(inParams.p1n, "m2")) // mass2 changed
  {
    for(p1 = inParams.p1s; p1 < inParams.p1e + 0.5*inParams.p1d; p1 += inParams.p1d)
    {
      p1v = p1*LAL_MSUN_SI; // convert to kg
      ret = XLALSimInspiralTaylorF2Ecc(&hptilde2, phiRef, deltaF, 
                    m1, p1v, S1z, S2z,
                    f_min, f_max, 0, r, 1, 1, 0, 0, ecc, ecc_order, f_ecc, 0, 0, phaseO, amplitudeO, NULL); // qm1=qm2 = 1 hard code
      //ret = XLALSimInspiralTaylorF2AmpPlus(&hptilde2, phiRef, deltaF, inclination,
      //              m1, p1v, S1x, S1y, S1z, S2x, S2y, S2z,
      //              LNhatx, LNhaty, LNhatz, f_min, f_max, r, phaseO, amplitudeO);
      overlap = LALInferenceComputeFrequencyDomainOverlap(irs->data, hptilde1->data, hptilde2->data);
      overlapNorm = LALInferenceComputeFrequencyDomainOverlapNormalised(irs->data, hptilde1->data, hptilde2->data);
      //loglikelihood = LALInferenceFreqDomainLogLikelihood(&currentParams, runstate->data, LALInferenceTemplateXLALSimInspiralChooseWaveform));
      fprintf(outf, "%20.8e %20.8e %20.8e\n", p1, overlap, overlapNorm);
      if ( hptilde2 ) XLALDestroyCOMPLEX16FrequencySeries(hptilde2);
      hptilde2 = NULL;
    }
  }
  else if(strstr(inParams.p1n, "ecc")) // ecc changed
  {
    for(p1 = inParams.p1s; p1 < inParams.p1e + 0.5*inParams.p1d; p1 += inParams.p1d)
    {
      p1v = p1; 
      ret = XLALSimInspiralTaylorF2Ecc(&hptilde2, phiRef, deltaF, 
                    m1, m2, S1z, S2z,
                    f_min, f_max, 0, r, 1, 1, 0, 0, p1v, ecc_order, f_ecc, 0, 0, phaseO, amplitudeO, NULL); // qm1=qm2 = 1 hard code
      //ret = XLALSimInspiralTaylorF2AmpPlus(&hptilde2, phiRef, deltaF, inclination,
      //              m1, p1v, S1x, S1y, S1z, S2x, S2y, S2z,
      //              LNhatx, LNhaty, LNhatz, f_min, f_max, r, phaseO, amplitudeO);
      overlap = LALInferenceComputeFrequencyDomainOverlap(irs->data, hptilde1->data, hptilde2->data);
      overlapNorm = LALInferenceComputeFrequencyDomainOverlapNormalised(irs->data, hptilde1->data, hptilde2->data);
      //loglikelihood = LALInferenceFreqDomainLogLikelihood(&currentParams, runstate->data, LALInferenceTemplateXLALSimInspiralChooseWaveform));
      fprintf(outf, "%20.8e %20.8e %20.8e\n", p1, overlap, overlapNorm);
      if ( hptilde2 ) XLALDestroyCOMPLEX16FrequencySeries(hptilde2);
      hptilde2 = NULL;
    }
  }
  fclose(outf);
  generateGNUfile(&inParams);
/*
  ret = XLALSimInspiralTaylorF2(&hptildeF2, phiRef, deltaF, m1, m2, S1z, S2z, f_min, f_ref, f_max, r, 1, 1, 0, 0, 0, 0,10, 0, 0, phaseO, amplitudeO, NULL);// qm1=qm2 = 1 hard code
  hctildeF2 = XLALCreateCOMPLEX16FrequencySeries("FD hcross",
                    &((hptildeF2)->epoch), (hptildeF2)->f0, (hptildeF2)->deltaF,
                    &((hptildeF2)->sampleUnits), (hptildeF2)->data->length);
  cfac = cos(inclination);
  pfac = 0.5 * (1. + cfac*cfac);
  //printf("End of calling\n");

  for(j = 0; j < (hptildeF2)->data->length; j++) {
    (hctildeF2)->data->data[j] = I*cfac*(hptildeF2)->data->data[j];
    (hptildeF2)->data->data[j] *= pfac;
  }

*/
  if( irs != NULL)
    XLALFree(irs);
  irs = NULL;
  if ( hptilde1 ) XLALDestroyCOMPLEX16FrequencySeries(hptilde1);
  if ( hptilde2 ) XLALDestroyCOMPLEX16FrequencySeries(hptilde2);
  if ( hctilde1 ) XLALDestroyCOMPLEX16FrequencySeries(hctilde1);
  if ( hctilde2 ) XLALDestroyCOMPLEX16FrequencySeries(hctilde2);
  if ( LALparams ) XLALDestroyDict(LALparams);
  printf("End of Ovelap program\n");
  return ret;
}

int parseParameters(ProcessParamsTable *procParams, Parameter *inParams)
{
  ProcessParamsTable *ppt = NULL;
  ppt = LALInferenceGetProcParamVal(procParams, "--m1");
  if (ppt) {
    inParams->m1 = atof(ppt->value);
  }
  ppt = LALInferenceGetProcParamVal(procParams, "--m2");
  if (ppt) {
    inParams->m2 = atof(ppt->value);
  }
  ppt = LALInferenceGetProcParamVal(procParams, "--d");
  if (ppt) {
    inParams->d = atof(ppt->value);
  }
  ppt = LALInferenceGetProcParamVal(procParams, "--iota");
  if (ppt) {
    inParams->iota = atof(ppt->value);
  }
  ppt = LALInferenceGetProcParamVal(procParams, "--pol");
  if (ppt) {
    inParams->pol = atof(ppt->value);
  }
  ppt = LALInferenceGetProcParamVal(procParams, "--orb");
  if (ppt) {
    inParams->orb = atof(ppt->value);
  }
  ppt = LALInferenceGetProcParamVal(procParams, "--coal");
  if (ppt) {
    inParams->coal = atof(ppt->value);
  }
  ppt = LALInferenceGetProcParamVal(procParams, "--ra");
  if (ppt) {
    inParams->ra = atof(ppt->value);
  }
  //printf("DEBUG 1\n");
  ppt = LALInferenceGetProcParamVal(procParams, "--s1x");
  inParams->s1x = 0.0;
  if (ppt) {
    inParams->s1x = atof(ppt->value);
  }
  ppt = LALInferenceGetProcParamVal(procParams, "--s1y");
  inParams->s1y = 0.0;
  if (ppt) {
    inParams->s1y = atof(ppt->value);
  }
  ppt = LALInferenceGetProcParamVal(procParams, "--s1z");
  inParams->s1z = 0.0;
  if (ppt) {
    inParams->s1z = atof(ppt->value);
  }
  ppt = LALInferenceGetProcParamVal(procParams, "--s2x");
  inParams->s2x = 0.0;
  if (ppt) {
    inParams->s2x = atof(ppt->value);
  }
  ppt = LALInferenceGetProcParamVal(procParams, "--s2y");
  inParams->s2y = 0.0;
  if (ppt) {
    inParams->s2y = atof(ppt->value);
  }
  ppt = LALInferenceGetProcParamVal(procParams, "--s2z");
  inParams->s2z = 0.0;
  if (ppt) {
    inParams->s2z = atof(ppt->value);
  }
  ppt = LALInferenceGetProcParamVal(procParams, "--ecc");
  inParams->ecc = 0.0;
  if (ppt) {
    inParams->ecc = atof(ppt->value);
  }
  ppt = LALInferenceGetProcParamVal(procParams, "--f_ecc");
  inParams->f_ecc = 0.0;
  if (ppt) {
    inParams->f_ecc = atof(ppt->value);
  }
  ppt = LALInferenceGetProcParamVal(procParams, "--ecc_order");
  inParams->ecc_order = -1;
  if (ppt) {
    inParams->ecc_order = atoi(ppt->value);
  }
  ppt = LALInferenceGetProcParamVal(procParams, "--dec");
  if (ppt) {
    inParams->dec = atof(ppt->value);
  }
  //printf("DEBUG 2\n");
  ppt = LALInferenceGetProcParamVal(procParams, "--approx1");
  if (ppt) {
    strcpy(inParams->approx1, ppt->value);
  }
  //printf("DEBUG 3\n");
  ppt = LALInferenceGetProcParamVal(procParams, "--approx2");
  if (ppt) {
    strcpy(inParams->approx2, ppt->value);
  }
  //printf("DEBUG 4\n");
  ppt = LALInferenceGetProcParamVal(procParams, "--ampOrder1");
  if (ppt) {
    inParams->ampOrder1 = atoi(ppt->value);
  }
  //printf("DEBUG 5\n");
  ppt = LALInferenceGetProcParamVal(procParams, "--ampOrder2");
  if (ppt) {
    inParams->ampOrder2 = atoi(ppt->value);
  }
  //printf("DEBUG 6\n");
  ppt = LALInferenceGetProcParamVal(procParams, "--phaseOrder1");
  if (ppt) {
    inParams->phaseOrder1 = atoi(ppt->value);
  }
  //printf("DEBUG 7\n");
  ppt = LALInferenceGetProcParamVal(procParams, "--phaseOrder2");
  if (ppt) {
    inParams->phaseOrder2 = atoi(ppt->value);
  }
  //printf("DEBUG 8\n");
  ppt = LALInferenceGetProcParamVal(procParams, "--outfile");
  if (ppt) {
    //printf("ppt->value = %s\n", ppt->value);
    strcpy(inParams->outfile, ppt->value);
  }
  //printf("outfile = %s\n", inParams->outfile);
  ppt = LALInferenceGetProcParamVal(procParams, "--2d");
  if (ppt) {
    inParams->dimension = 2;
  }
  //printf("DEBUG 9\n");
  ppt = LALInferenceGetProcParamVal(procParams, "--1d");
  if (ppt) {
    inParams->dimension = 1;
  }
  ppt = LALInferenceGetProcParamVal(procParams, "--p1n");
  if (ppt) {
    strcpy(inParams->p1n, ppt->value);
  }
  ppt = LALInferenceGetProcParamVal(procParams, "--p1s");
  if (ppt) {
    inParams->p1s = atof(ppt->value);
  }
  ppt = LALInferenceGetProcParamVal(procParams, "--p1e");
  if (ppt) {
    inParams->p1e = atof(ppt->value);
  }
  ppt = LALInferenceGetProcParamVal(procParams, "--p1d");
  if (ppt) {
    inParams->p1d = atof(ppt->value);
  }
  ppt = LALInferenceGetProcParamVal(procParams, "--p2n");
  if (ppt) {
    strcpy(inParams->p2n, ppt->value);
  }
  ppt = LALInferenceGetProcParamVal(procParams, "--p2s");
  if (ppt) {
    inParams->p2s = atof(ppt->value);
  }
  ppt = LALInferenceGetProcParamVal(procParams, "--p2e");
  if (ppt) {
    inParams->p2e = atof(ppt->value);
  }
  ppt = LALInferenceGetProcParamVal(procParams, "--p2d");
  if (ppt) {
    inParams->p2d = atof(ppt->value);
  }
  //printf("DEBUG 100\n");
  return 0;
}
int printParameters(FILE *out, Parameter *inParams)
{
  fprintf(out, "m1 = %.3fM_sun\n", inParams->m1);
  fprintf(out, "m2 = %.3fM_sun\n", inParams->m2);
  return 0;
}
REAL8 LALInferenceComputeFrequencyDomainOverlapNormalised(LALInferenceIFOData * dataPtr, 
                                   COMPLEX16Vector * freqData1, 
                                   COMPLEX16Vector * freqData2){
  REAL8 rho1, rho2, rho0, rhoN;
  rho1 = LALInferenceComputeFrequencyDomainOverlap(dataPtr, freqData1, freqData1);
  rho2 = LALInferenceComputeFrequencyDomainOverlap(dataPtr, freqData2, freqData2);
  rho0 = LALInferenceComputeFrequencyDomainOverlap(dataPtr, freqData1, freqData2);
  rhoN = rho0/(sqrt(rho1)*sqrt(rho2));
  return rhoN;
}
void generateGNUfile(Parameter *inParams)
{
  FILE *outf;
  char gnu_file[FILE_NAME_SIZE];
  sprintf(gnu_file, "%s.gnu", inParams->outfile);
  outf = fopen(gnu_file, "w");
  fprintf(outf, "reset\n");
  fprintf(outf, "set title \"Overlaps %s-%s w.r.t. Parameters\"\n", inParams->approx1, inParams->approx2);
  if(inParams->dimension == 1) // only 1d graph is allowed currently
  {
    fprintf(outf, "set grid x y2\n");
    fprintf(outf, "set xlabel \"%s\"\n", inParams->p1n);
    fprintf(outf, "set xrange [%e:%e]\n", inParams->p1s*0.9, inParams->p1e*1.1);
    fprintf(outf, "set x2range [%e:%e]\n", inParams->p1s*0.9, inParams->p1e*1.1);
    fprintf(outf, "set ylabel \"Overlap\"\n");
    fprintf(outf, "set y2label \"Normalized Overlap\"\n");
    fprintf(outf, "set ytics nomirror\n");
    fprintf(outf, "set y2tics\n");
    fprintf(outf, "set tics out\n");
    fprintf(outf, "set autoscale y\n");
    fprintf(outf, "set autoscale y2\n");
    fprintf(outf, "set terminal postscript landscape colour\n");
    fprintf(outf, "set output \"%s_abs.eps\"\n", inParams->outfile);
    fprintf(outf, "plot \"%s.dat\" using ($1):(abs($2)) with lines axes x1y1 title \"Unnormalised\", \"%s.dat\" using ($1):(abs($3)) with lines axes x2y2 title \"Normalised\"\n", 
		inParams->outfile, inParams->outfile);
    fprintf(outf, "set output \"%s.eps\"\n", inParams->outfile);
    fprintf(outf, "plot \"%s.dat\" using ($1):($2) with lines axes x1y1 title \"Unnormalised\", \"%s.dat\" using ($1):($3) with lines axes x2y2 title \"Normalised\"\n", 
		inParams->outfile, inParams->outfile);
    fprintf(outf, "set terminal qt\n");
    fprintf(outf, "plot \"%s.dat\" using ($1):(abs($2)) with lines axes x1y1 title \"Unnormalised\", \"%s.dat\" using ($1):(abs($3)) with lines axes x2y2 title \"Normalised\"\n", 
		inParams->outfile, inParams->outfile);
  }
  fclose(outf);
}
