
/*
 *  LALInferencePriorVolumes.c
 *
 *  Copyright (C) 2016 Salvatore Vitale
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
#include <math.h>
#include <gsl/gsl_integration.h>
#include <lal/LALInference.h>
#include <lal/LALInferencePriorVolumes.h>
#include <lal/LALInferencePrior.h>
#include <lal/LALCosmologyCalculator.h>

typedef struct {
  double mc;
  LALInferenceVariables *prior;
} innerParams;

typedef struct {
  LALInferenceVariables *prior;
} outerParams;

static double copt=-0.1, zopt=0.5;
static double chirp_to_comp_jacobian(double mc,double mzc);
static double mass_indicator(double mc, double mzc,LALInferenceVariables *priorParams);
static double integrand(double mzc,void *params);
static double inner_integral(double mc, void *params);
static double mass_outer_integral(LALInferenceVariables *priorArgs);
static double loudness_volume(LALInferenceRunState *state);

static double chirp_to_comp_jacobian(double mc,double mzc){

  double kopt=5.0/(3.0*copt -2.0*(copt*zopt));
  double delta_opt=sqrt(1.0 - 4.0 * pow(mc, kopt*copt*(1.0+zopt)) * pow(mzc, -kopt));
  double m1 = 0.5* pow(mc, -kopt * (copt*zopt)) * pow(mzc, 3.0*kopt/5.0) * (1.0 + delta_opt);
  double m2 = 0.5* pow(mc, -kopt * (copt*zopt)) * pow(mzc, 3.0*kopt/5.0) * (1.0 - delta_opt);
  return 5.0*m1*m2*(m1+m2)/((2*(zopt * copt)-3*copt)*mc*mzc*(m1-m2));
}


static double mass_indicator(double mc, double mzc,LALInferenceVariables *priorParams){
  /*
   * 
   * This function gets mchirp and mzc
   * Calculates component and total mass and check if all parameters
   * are within their priors.
   * 
   * Return 1 or 0 
   * 
   * t work with mzc
   * 
   * */
  double kopt=5.0/(3.0*copt -2.0*(copt*zopt));
  double delta_opt=sqrt(1.0 - 4.0 * pow(mc, kopt*copt*(1.0+zopt)) * pow(mzc, -kopt));
  double m1 = 0.5* pow(mc, -kopt * (copt*zopt)) * pow(mzc, 3.0*kopt/5.0) * (1.0 + delta_opt);
  double m2 = 0.5* pow(mc, -kopt * (copt*zopt)) * pow(mzc, 3.0*kopt/5.0) * (1.0 - delta_opt);
printf("\n inside mass_indicator\n");
printf("\n m1=%f\n",m1);
printf("\n m2=%f\n",m2);
  /* Check for individual mass priors */
  if(LALInferenceCheckVariable(priorParams,"mass1_min"))
    if(LALInferenceGetREAL8Variable(priorParams,"mass1_min") > m1)
      {printf("\n1 is working\n");return 0.0;}
  if(LALInferenceCheckVariable(priorParams,"mass1_max"))
    if(LALInferenceGetREAL8Variable(priorParams,"mass1_max") < m1)
      {printf("\n2 is working\n");return 0.0;}
  if(LALInferenceCheckVariable(priorParams,"mass2_min"))
    if(LALInferenceGetREAL8Variable(priorParams,"mass2_min") > m2)
      {printf("\n3 is working\n");return 0.0;}
  if(LALInferenceCheckVariable(priorParams,"mass2_max"))
    if(LALInferenceGetREAL8Variable(priorParams,"mass2_max") < m2)
      {printf("\n4 is working\n");return 0.0;}
  if(LALInferenceCheckVariable(priorParams,"MTotMax"))
    if(*(REAL8 *)LALInferenceGetVariable(priorParams,"MTotMax") < m1+m2)
      {printf("\n5 is working\n");return 0.0;}
  if(LALInferenceCheckVariable(priorParams,"MTotMin"))
    if(*(REAL8 *)LALInferenceGetVariable(priorParams,"MTotMin") > m1+m2)
      {printf("\n6 is working\n");return 0.0;}
  if (LALInferenceCheckVariable(priorParams,"mzc_min")){
    double mzc_min,mzc_max;
    LALInferenceGetMinMaxPrior(priorParams, "mzc", &mzc_min, &mzc_max);
    if (mzc<mzc_min || mzc>mzc_max)       {printf("\n7 is working\n");return 0.0;}
  }
  if (LALInferenceCheckVariable(priorParams,"chirpmass_min")){
    double mc_min,mc_max;
    LALInferenceGetMinMaxPrior(priorParams, "chirpmass", &mc_min, &mc_max);
    if (mc<mc_min || mc>mc_max) {printf("\n8 is working\n");return 0.0;}
  }
  return 1.0;

}



static double integrand(double mzc,void *params){
  
  /* Integrand for the dobule integral over Mc and mzc
   * 
   * This is the jacobian within the support of the prior, zero otherwise
   * 
   * */
  innerParams iData = *(innerParams *)params;
  double mc= iData.mc;
  //printf("\nintegrand val  =\t%f\n", chirp_to_comp_jacobian(mc,mzc));
  printf("\nintegrand val  =\t%f\n", mass_indicator(mc,mzc, iData.prior));
  return chirp_to_comp_jacobian(mc,mzc)*mass_indicator(mc,mzc, iData.prior);
  
}


static double inner_integral(double mc, void *params){
  
  // mzc comes through params
  printf("\ninner integral fn invoked . .\n");
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (10000);
  double result, error;
  const double epsabs = 1e-4;
  const double epsrel = 1e-4;
  const size_t wsSize = 10000;
  gsl_function F;
  F.function = &integrand;
  innerParams iParams;
  F.params = &iParams;
  outerParams oParams=*(outerParams *) params;
  iParams.mc=mc;
  iParams.prior= (LALInferenceVariables *) oParams.prior;
  
  double mzc_min,mzc_max;
  if (LALInferenceCheckVariable(iParams.prior,"mzc_min")){
    
    LALInferenceGetMinMaxPrior(iParams.prior, "mzc", &mzc_min, &mzc_max);
    
  }
  else{
    fprintf(stderr,"ERROR: mzc doesn't seem to be a valid param. Exiting\n");
    exit(1);
    }
  // TODO: make it work with  and q
  //printf("\nline before gsl call\n");
  //printf("error=%f\n",error);
  printf("\nb4 mzc\n");
  int status = gsl_integration_qags (&F, mzc_min, mzc_max,epsabs, epsrel, wsSize,
                        w, &result, &error);
  printf("\naftr mzc\n");
  if (status)
        XLAL_ERROR_REAL8(XLAL_EFUNC | XLAL_EDATA, "Bad data; GSL integration failed.");
  gsl_integration_workspace_free(w);
  printf("\ninner integral fn at return . .\n");
  return result;  
  }

  
static double mass_outer_integral(LALInferenceVariables *priorArgs){
  printf("\nouter integral fn invoked . .\n");
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (10000);
  double result, error;
  const double epsabs = 1e-4;
  const double epsrel = 1e-4;
  const size_t wsSize = 10000;

  gsl_function F;
  F.function = &inner_integral;
  outerParams oParams;
  F.params = &oParams;
  oParams.prior=priorArgs;
  
  double mc_min,mc_max;
  if (LALInferenceCheckVariable(priorArgs,"chirpmass_min")){
    
    LALInferenceGetMinMaxPrior(priorArgs, "chirpmass", &mc_min, &mc_max);
  }
  else{
    fprintf(stderr,"ERROR: chirpmass doesn't seem to be a valid param. Exiting\n");
    exit(1);
    }
  /* this integrates on chirpmass */
  printf("\nb4 mc\n");
  int status = gsl_integration_qags (&F, mc_min, mc_max, epsabs, epsrel, wsSize, w, &result, &error); 
  printf("\naftr mc\n");
  
  if (status)
        XLAL_ERROR_REAL8(XLAL_EFUNC | XLAL_EDATA, "Bad data; GSL integration failed.");

  gsl_integration_workspace_free(w);
  printf("\nouter integral fn at return . .\n");
  return result;
  }

typedef REAL8 (*LALInferenceLoudnessPriorFunction) (double x,LALCosmologicalParameters *omega);
  
typedef struct {
  LALInferenceVariables *priorargs;
  LALInferenceLoudnessPriorFunction priorfunc;
  LALCosmologicalParameters *omega;
} loudnessParams;

static double distance_prior(double d,LALCosmologicalParameters *omega){
  (void) omega;
  return d*d;
}
static double logdistance_prior(double ld,LALCosmologicalParameters *omega){
  (void) omega;
  return ld*ld*ld;
}
static double redshift_prior(double z,LALCosmologicalParameters *omega){
  return XLALUniformComovingVolumeDensity(z,omega);
}
  
static double loudness_integrand(double x,void *params){
  
  loudnessParams lParams=*(loudnessParams *) params;
  LALInferenceLoudnessPriorFunction priorf= lParams.priorfunc;
  LALCosmologicalParameters *omega=lParams.omega;
  return priorf(x,omega);
}
  
static double loudness_volume(LALInferenceRunState *state){


  LALInferenceVariables *priorArgs=state->priorArgs;
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (10000);
  double result, error;
  const double epsabs = 1e-4;
  const double epsrel = 1e-4;
  const size_t wsSize = 10000;

  gsl_function F;
  F.function = &loudness_integrand;
  loudnessParams lParams;
  F.params = &lParams;
  double intmin,intmax;
  if (LALInferenceCheckVariable(priorArgs,"redshift_min")){
    lParams.priorfunc=&redshift_prior;
    // WILL need to enable this after we push the cosmology patch, for the moment set to null
    //lParams.omega=state->omega;
    lParams.omega=NULL;

    LALInferenceGetMinMaxPrior(priorArgs, "redshift", &intmin, &intmax);
  }
  else if (LALInferenceCheckVariable(priorArgs,"distance_min")){
    lParams.priorfunc=&distance_prior;
    lParams.omega=NULL;
    LALInferenceGetMinMaxPrior(priorArgs, "distance", &intmin, &intmax);
  }
  else if (LALInferenceCheckVariable(priorArgs,"logdistance_min")){
    lParams.priorfunc=&logdistance_prior;
    lParams.omega=NULL;
    LALInferenceGetMinMaxPrior(priorArgs, "logdistance", &intmin, &intmax);
  }
  else XLAL_ERROR_REAL8(XLAL_EINVAL,"No known distance parameter found, unable to proceed\n");

  lParams.priorargs=priorArgs;
  printf("b4 3rd");
  int status = gsl_integration_qags (&F, intmin, intmax, epsabs, epsrel, wsSize, w, &result, &error); 
  printf("aftr 3rd");
  
  if (status)
        XLAL_ERROR_REAL8(XLAL_EFUNC | XLAL_EDATA, "Bad data; GSL integration failed.");

  gsl_integration_workspace_free(w);
  return result;
  
}
  
double LALInferenceMassPriorVolume(LALInferenceRunState *state){

  LALInferenceVariables *priorArgs=state->priorArgs;
  return mass_outer_integral(priorArgs);
}

double LALInferenceMassDistancePriorVolume(LALInferenceRunState *state){
  
  return LALInferenceMassPriorVolume(state)*loudness_volume(state);
}
