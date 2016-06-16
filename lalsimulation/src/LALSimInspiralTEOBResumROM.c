/*
 *  Copyright (C) 2014 Michael Puerrer, John Veitch
 *  Reduced Order Model for SEOBNR
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

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <time.h>
#include <stdbool.h>
#include <alloca.h>
#include <string.h>
#include <libgen.h>


#include <gsl/gsl_errno.h>
#include <gsl/gsl_chebyshev.h>

// #include <gsl/gsl_blas.h>
// #include <gsl/gsl_min.h>
#include <lal/Units.h>
#include <lal/SeqFactories.h>
#include <lal/LALConstants.h>
#include <lal/XLALError.h>
#include <lal/FrequencySeries.h>
#include <lal/Date.h>
#include <lal/StringInput.h>
#include <lal/Sequence.h>
#include <lal/LALStdio.h>
#include <lal/FileIO.h>

#include <lal/LALSimInspiral.h>
#include <lal/LALSimIMR.h>

//gsl_bspline needs to be imported here, because this was not done in ROMUtilities!
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_spline.h>
//needed for read_vector
#include "LALSimIMRSEOBNRROMUtilities.c"

#include <lal/LALConfig.h>
#ifdef LAL_PTHREAD_LOCK
#include <pthread.h>
#endif


/******************************************************************
 * The double-spin SEOBNRv2 ROM consists of a number of submodels *
 * for subdomains of the whole parameter space:                   *
 * These submodels may change in the future.                      *
 *                                                                *
 * "Core" submodel 1                                              *
 * B-spline points: 54x24x24                                      *
 * Frequency points: {133, 139}                                   *
 *                                                                *
 * "Near q=1" submodel 2                                          *
 * B-spline points: 11x60x60                                      *
 * Frequency points: {200, 187}                                   *
 *                                                                *
 * "High q, high chi1" submodel 3                                 *
 * B-spline points: 24x30x24                                      *
 * frequency points: {133, 139}                                   *
 *                                                                *
 *****************************************************************/


/********* Input data for spline basis points **************/

//Prepending G to recognise it as a global variable
#define Gntimes 69072 //(will be reduced in the future)
#define Gnamp 21
#define Gnphase 7
#define Gnq 16
#define Gnlambda1 16
#define Gnlambda2 16

#ifdef LAL_PTHREAD_LOCK
static pthread_once_t TEOBResumROM_is_initialized = PTHREAD_ONCE_INIT;
#endif

static const double Gparams_min[] = {0.5,50.,50.}; //qmin,lambda1min,lambda2min
static const double Gparams_max[] = {1.0,5000.,5000.}; //qmax,lambda1max,lambda2max

/*************** type definitions ******************/

typedef struct tagTEOBResumROMdataDS_coeff
{
  gsl_vector* c_amp;
  gsl_vector* c_phi;
} TEOBResumROMdataDS_coeff;

struct tagTEOBResumROMdataDS_submodel
{
  gsl_vector* cvec_amp;      // amplitude projection coefficients
  gsl_vector* cvec_phi;      // phase projection coefficients
  gsl_matrix *Bamp;          // Reduced SVD basis for amplitude
  gsl_matrix *Bphi;          // Reduced SVD basis for phase
  gsl_vector* times;  // AMplitude prefactor coefficient
  const double *params_min;
  const double *params_max;
  int n_amp;                 // Number frequency points for amplitude
  int n_phi;                 // Number of frequency points for phase
  int nq, nl1, nl2, ntimes;         // Number of points in eta, chi1, chi2
};
typedef struct tagTEOBResumROMdataDS_submodel TEOBResumROMdataDS_submodel;

struct tagTEOBResumROMdataDS
{
  UINT4 setup;
  TEOBResumROMdataDS_submodel* sub1;
  TEOBResumROMdataDS_submodel* sub2;
  TEOBResumROMdataDS_submodel* sub3;
};
typedef struct tagTEOBResumROMdataDS TEOBResumROMdataDS;

static TEOBResumROMdataDS __lalsim_TEOBResumROMDS_data;

typedef int (*load_dataPtr)(const char*, gsl_vector *, gsl_vector *, gsl_matrix *, gsl_matrix *, gsl_vector *);

// typedef struct tagChebyshevData
// {
//   gsl_cheb_series *csx ;
//   gsl_cheb_series *csy ;
//   gsl_cheb_series *csz ;
// } ChebyshevData;

/**************** Internal functions **********************/

static void TEOBResumROM_Init_LALDATA(void);
static int TEOBResumROM_Init(const char dir[]);
static bool TEOBResumROM_IsSetup(void);

static int TEOBResumROMdataDS_Init(TEOBResumROMdataDS *romdata, const char dir[]);
static void TEOBResumROMdataDS_Cleanup(TEOBResumROMdataDS *romdata);

static int TEOBResumROMdataDS_Init_submodel(
  TEOBResumROMdataDS_submodel **submodel,
  const int n_amp,
  const int n_phi,
  const int nq,
  const int nl1,
  const int nl2,
  const int ntimes,
  const double *params_min,
  const double *params_max,
  const char dir[],
  load_dataPtr load_data
);


// static void ChebyshevData_Init(
//   ChebyshevData **chebydata);
//
// static void ChebyshevData_Destroy(
//   ChebyshevData *chebydata);

static double gsl_cheb_evaluate_polynomial(int n, double x);
static double gsl_cheb_eval_3d(gsl_vector *c_ijk, int nx, int ny, int nz, double x, double y, double z);
static int chebyshev_interpolation3d(
  double q,
  double lambda1,
  double lambda2,
  int nx, int ny, int nz,
  gsl_vector *cvec_amp,
  gsl_vector *cvec_phi,
  //gsl_vector *times,
  int nk_amp,
  int nk_phi,
  gsl_vector *interp_amp,
  gsl_vector *interp_phi);

static void TEOBResumROMdataDS_Cleanup_submodel(TEOBResumROMdataDS_submodel *submodel);

static int TEOBResumROMCore_test(
  REAL8TimeSeries **hPlus,
  REAL8TimeSeries **hCross,
  double phiRef,
  double deltaT,
  double fRef,
  double distance,
  double inclination,
  double Mtot_sec,
  double eta,
  double lambda1,
  double lambda2
);

static int TEOBResumROMCore(
  REAL8TimeSeries **hPlus,
  REAL8TimeSeries **hCross,
  double phiRef,
  double deltaT,
  double fRef,
  double distance,
  double inclination,
  double Mtot_sec,
  double eta,
  double lambda1,
  double lambda2
);

// static void SEOBNRROMdataDS_coeff_Init(SEOBNRROMdataDS_coeff **romdatacoeff, int nk_amp, int nk_phi);
// static void SEOBNRROMdataDS_coeff_Cleanup(SEOBNRROMdataDS_coeff *romdatacoeff);

// static size_t NextPow2(const size_t n);

static int load_data_romeos(const char dir[], gsl_vector *cvec_amp, gsl_vector *cvec_phi, gsl_matrix *Bamp, gsl_matrix *Bphi, gsl_vector *times);

/********************* Definitions begin here ********************/

/** Setup SEOBNRv2ROMDoubleSpin model using data files installed in dir
 */
static int TEOBResumROM_Init(const char dir[]) {
  if(__lalsim_TEOBResumROMDS_data.setup) {
    XLALPrintError("Error: DSTEOBResumROMdata was already set up!");
    XLAL_ERROR(XLAL_EFAILED);
  }

  TEOBResumROMdataDS_Init(&__lalsim_TEOBResumROMDS_data, dir);

  if(__lalsim_TEOBResumROMDS_data.setup) {
    return(XLAL_SUCCESS);
  }
  else {
    return(XLAL_EFAILED);
  }
}

/** Helper function to check if the SEOBNRv2ROMDoubleSpin model has been initialised */
static bool TEOBResumROM_IsSetup(void) {
  if(__lalsim_TEOBResumROMDS_data.setup)
    return true;
  else
    return false;
}

// Read binary ROM data for basis functions and coefficients for submodel 1
static int load_data_romeos(const char dir[], gsl_vector *cvec_amp, gsl_vector *cvec_phi, gsl_matrix *Bamp, gsl_matrix *Bphi, gsl_vector *times) {
  // Load binary data for amplitude and phase spline coefficients and reduced bases as computed in Mathematica
  // "Core" submodel 1
  // B-spline points: 54x24x24
  // Frequency points: {133, 139}
  int ret = XLAL_SUCCESS;
  ret |= read_vector(dir, "TEOBResumROM_Amp_ciall.dat", cvec_amp);
  ret |= read_vector(dir, "TEOBResumROM_Phase_ciall.dat", cvec_phi);
  ret |= read_matrix(dir, "TEOBResumROM_Bamp_matrix.dat", Bamp);
  ret |= read_matrix(dir, "TEOBResumROM_Bphase_matrix.dat", Bphi);
  ret |= read_vector(dir, "TEOBResumROM_times.dat", times);
  return(ret);
}

/* Set up a new ROM submodel, using data contained in dir */
static int TEOBResumROMdataDS_Init_submodel(
  TEOBResumROMdataDS_submodel **submodel,
  const int n_amp,
  const int n_phi,
  const int nq,
  const int nl1,
  const int nl2,
  const int ntimes,
  const double *params_min,
  const double *params_max,
  const char dir[],
  load_dataPtr load_data
) {
  int ret = XLAL_FAILURE;

  if(!submodel) exit(1);
  /* Create storage for submodel structures */
  if (!*submodel)
    *submodel = XLALCalloc(1,sizeof(TEOBResumROMdataDS_submodel));
  else
    TEOBResumROMdataDS_Cleanup_submodel(*submodel);

  int N = nq*nl1*nl2; // Total number of points over parameter space = size of the data matrix for one SVD-mode

  // Initalize actual ROM data
  (*submodel)->cvec_amp = gsl_vector_alloc(N*n_amp);
  (*submodel)->cvec_phi = gsl_vector_alloc(N*n_phi);
  (*submodel)->Bamp = gsl_matrix_alloc(n_amp, ntimes);
  (*submodel)->Bphi = gsl_matrix_alloc(n_phi, ntimes);
  (*submodel)->times = gsl_vector_alloc(ntimes);

  // Load ROM data for this submodel
  ret=load_data(dir, (*submodel)->cvec_amp, (*submodel)->cvec_phi, (*submodel)->Bamp, (*submodel)->Bphi, (*submodel)->times);

  // Initialize other members
  (*submodel)->n_amp = n_amp;
  (*submodel)->n_phi = n_phi;
  (*submodel)->nq = nq;
  (*submodel)->nl1 = nl1;
  (*submodel)->nl2 = nl2;
  (*submodel)->ntimes = ntimes;

  (*submodel)->params_min = params_min;
  (*submodel)->params_max = params_max;

  return ret;
}

/* Deallocate contents of the given TEOBResumROMdataDS_submodel structure */
static void TEOBResumROMdataDS_Cleanup_submodel(TEOBResumROMdataDS_submodel *submodel) {
  if(submodel->cvec_amp) gsl_vector_free(submodel->cvec_amp);
  if(submodel->cvec_phi) gsl_vector_free(submodel->cvec_phi);
  if(submodel->Bamp) gsl_matrix_free(submodel->Bamp);
  if(submodel->Bphi) gsl_matrix_free(submodel->Bphi);
  if(submodel->times) gsl_vector_free(submodel->times);
}

/* Set up a new ROM model, using data contained in dir */
int TEOBResumROMdataDS_Init(TEOBResumROMdataDS *romdata, const char dir[]) {
  int ret = XLAL_FAILURE;

  /* Create storage for structures */
  if(romdata->setup) {
    XLALPrintError("WARNING: You tried to setup the TEOBResumROM model that was already initialised. Ignoring\n");
    return (XLAL_FAILURE);
  }

  gsl_set_error_handler(&err_handler);

  load_dataPtr load_data = &load_data_romeos;
  ret = TEOBResumROMdataDS_Init_submodel(&(romdata)->sub1, Gnamp, Gnphase, Gnq, Gnlambda1, Gnlambda2, Gntimes, Gparams_min, Gparams_max, dir, load_data);
  if (ret==XLAL_SUCCESS) XLALPrintInfo("%s : submodel 1 loaded sucessfully.\n", __func__);

  if(XLAL_SUCCESS==ret)
    romdata->setup=1;
  else
    TEOBResumROMdataDS_Cleanup(romdata);

  return (ret);
}

/* Deallocate contents of the given TEOBResumROMdataDS structure */
static void TEOBResumROMdataDS_Cleanup(TEOBResumROMdataDS *romdata) {
  TEOBResumROMdataDS_Cleanup_submodel((romdata)->sub1);
  XLALFree((romdata)->sub1);
  (romdata)->sub1 = NULL;
  romdata->setup=0;
}

/* Structure for internal use */
// static void TEOBResumROMdataDS_coeff_Init(TEOBResumROMdataDS_coeff **romdatacoeff, int nk_amp, int nk_phi) {
//
//   if(!romdatacoeff) exit(1);
//   /* Create storage for structures */
//   if(!*romdatacoeff)
//     *romdatacoeff=XLALCalloc(1,sizeof(TEOBResumROMdataDS_coeff));
//   else
//     TEOBResumROMdataDS_coeff_Cleanup(*romdatacoeff);
//
//   (*romdatacoeff)->c_amp = gsl_vector_alloc(nk_amp);
//   (*romdatacoeff)->c_phi = gsl_vector_alloc(nk_phi);
// }
//
// /* Deallocate contents of the given TEOBResumROMdataDS_coeff structure */
// static void TEOBResumROMdataDS_coeff_Cleanup(TEOBResumROMdataDS_coeff *romdatacoeff) {
//   if(romdatacoeff->c_amp) gsl_vector_free(romdatacoeff->c_amp);
//   if(romdatacoeff->c_phi) gsl_vector_free(romdatacoeff->c_phi);
//   XLALFree(romdatacoeff);
// }

/* Return the closest higher power of 2  */
// Note: NextPow(2^k) = 2^k for integer values k.
// static size_t NextPow2(const size_t n) {
//   return 1 << (size_t) ceil(log2(n));
// }

// static void ChebyshevData_Init(
//   ChebyshevData **chebydata
// )
// {
//   if(!chebydata) exit(1);
//   if(*chebydata) ChebyshevData_Destroy(*chebydata);
//
//   gsl_cheb_series *csx = gsl_cheb_alloc(3);
//   gsl_cheb_series *csy = gsl_cheb_alloc(3);
//   gsl_cheb_series *csz = gsl_cheb_alloc(3);
//
//   (*chebydata)=XLALCalloc(1,sizeof(ChebyshevData));
// }
//
// static void ChebyshevData_Destroy(ChebyshevData *chebydata)
// {
//   if(!chebydata) return;
//   if(chebydata->csx) gsl_bspline_free(chebydata->csx);
//   if(chebydata->csy) gsl_bspline_free(chebydata->csy);
//   if(chebydata->csz) gsl_bspline_free(chebydata->csz);
//   XLALFree(chebydata);
// }



/* To evaluate T_n(x) this function will call gsl_cheb_eval(csx,x)
 * gsl_cheb_eval(csx,x) <=> f(x) = (c_0 / 2) + \sum_{n=1} c_n T_n(x)
 *   Define f^N(x) := gsl_cheb_eval(csxN,x), where csxN is a list of ones of length n
 * We can then calculate any T_n(x) as:
 *   T_n(x) = f^n(x) - f^{n-1}(x)
 */
static double gsl_cheb_evaluate_polynomial(int n, double x){

  double Tnx = 0.0 ;
  //T_0(x)
  if (n==0){
    Tnx = 1.0 ;
  }
  //T_1(x)
  else if (n==1){
    Tnx = x ;
  }
  //T_2(x)
  else if (n==2){
    Tnx = 2.0*x*x - 1.0 ;
  }
  //T_n(x)
  else {
    gsl_cheb_series *csx1 = gsl_cheb_alloc(n);
    csx1->a = -1.0 ;
    csx1->b = 1.0 ;
    int i=0;
    for (i = 0; i < n; i++){
      csx1->c[i]=0.0;
    }
    csx1->c[n]=1.0;

    //f^n
    Tnx = gsl_cheb_eval(csx1, x);
    gsl_cheb_free(csx1);
  }

  return Tnx ;

}


/* Wrapper function to evaluate 3d chebyshev polynomial.
 * p(x,y,z) = \sum_{i,j,k} c_{ijk} T_i(x) T_j(y) T_k(z)
 */
static double gsl_cheb_eval_3d(gsl_vector *c_ijk, int nx, int ny, int nz, double x, double y, double z){

  double sum=0.0;
  int i,j,k;
  int index=0;
  double Tix=0.,Tjy=0.,Tkz=0.;

  for (i=0;i<nx;i++){
    Tix=gsl_cheb_evaluate_polynomial(i,x);
    for (j=0;j<ny;j++){
      Tjy=gsl_cheb_evaluate_polynomial(j,y);
      for (k=0;k<nz;k++){
        Tkz=gsl_cheb_evaluate_polynomial(k,z);
        sum+=gsl_vector_get(c_ijk,index)*Tix*Tjy*Tkz;
        index+=1;
      }
    }
  }

  return sum ;

}

/*
 * This function calculates the phase and amplitude interpolants at a time Tj.
 *  pj(q,l1,l2) = sum_l sum_m sum_n b_lmn Tl(q) Tm(l1) Tn(l2)
 * GSL can only evaluate a 1D chebyshev polynomial with gsl_cheb_eval_n, which calculates
 *  f(x) = (c_0 / 2) + \sum_{n=1}^(N_coeffs) c_n T_n(x)
 * This means that p(x,y,z) will have to be evaluated as
 *  p(x,y,z) = f'(x,csx,nx) f'(y,csy,ny) f'(z,csz,nz)
 * where csx are the chebyshev coefficients for parameter x , nx its size and f' is defined as
 *  f'(x;csx,nx) = f(x,csx,nx) + csx[0]/2.0
 */
/* NOTE: when we calculate p(x,y,z) = \sum_{i,j,k} c_{ijk} T_i(q) T_j(l1) T_k(l2)
 *       (done here with gsl_cheb_eval_3d for each k in nk_amp and nk_phi),
 *       we get an array for amp and phi with length nk_amp and nk_phi resp.
 *       All this is just for a single Tj. What do we do with these?
 */
static int chebyshev_interpolation3d(
  double q,
  double lambda1,
  double lambda2,
  int nx, int ny, int nz,
  gsl_vector *cvec_amp,
  gsl_vector *cvec_phi,
  int nk_amp,
  int nk_phi,
  gsl_vector *interp_amp,   //return: A(T_j;q,lambda1,lambda2)    <=> p_j(q,lambda1,lambda2)
  gsl_vector *interp_phi)   //return: \Phi(T_j;q,lambda1,lambda2) <=> p_j(q,lambda1,lambda2)
{

  double sum = 0.0;
  int k = 0;
  int N=nx*ny*nz ;

  //TODO: ranges should be given as function argument
  //Rescale so that x,y and z are always between -1 and 1
  double xrescale = (q-0.5*(Gparams_max[0]+Gparams_min[0])) / (0.5*(Gparams_max[0]-Gparams_min[0]));
  double yrescale = (lambda1-0.5*(Gparams_max[1]+Gparams_min[1])) / (0.5*(Gparams_max[1]-Gparams_min[1]));
  double zrescale = (lambda2-0.5*(Gparams_max[2]+Gparams_min[2])) / (0.5*(Gparams_max[2]-Gparams_min[2]));

  /*-- Calculate interp_amp --*/
  for (k=0; k<nk_amp; k++) { // For each empirical node
    gsl_vector v = gsl_vector_subvector(cvec_amp, k*N, N).vector; // Pick out the coefficient matrix corresponding to the k-th node.
    sum = gsl_cheb_eval_3d(&v, nx, ny, nz, xrescale, yrescale, zrescale) ;
    gsl_vector_set(interp_amp, k, sum); //write p_k(x,y,z)
  }

  /*-- Calculate interp_phi --*/
  for (k=0; k<nk_phi; k++) { // For each empirical node
    gsl_vector v = gsl_vector_subvector(cvec_phi, k*N, N).vector; // Pick out the coefficient matrix corresponding to the k-th node.
    sum = gsl_cheb_eval_3d(&v, nx, ny, nz, xrescale, yrescale, zrescale) ;
    gsl_vector_set(interp_phi, k, sum); //write p_k(x,y,z)
  }

  return 0;

}

static int TEOBResumROMCore_test(
  REAL8TimeSeries **hPlus,
  REAL8TimeSeries **hCross,
  double phiRef, // orbital reference phase
  double deltaT,
  double fRef,
  double distance,
  double inclination,
  double Mtot, // in Msol
  double eta,
  double lambda1,
  double lambda2
)
{

  /* Check output arrays */
  if(!hPlus || !hCross)
    XLAL_ERROR(XLAL_EFAULT);
  TEOBResumROMdataDS *romdata=&__lalsim_TEOBResumROMDS_data;
  if(*hPlus || *hCross)
  {
    XLALPrintError("(*hPlus) and (*hCross) are supposed to be NULL, but got %p and %p",(*hPlus),(*hCross));
    XLAL_ERROR(XLAL_EFAULT);
  }
  int retcode=0;

  REAL8TimeSeries *hp;
  REAL8TimeSeries *hc;

  /* Select ROM submodel */
  TEOBResumROMdataDS_submodel *submodel;
  submodel = romdata->sub1;

  fprintf(stdout,"--- EOSROMEOSCore input parameters ---\n");
  fprintf(stdout,"  phiRef       = %.2f\n",phiRef);
  fprintf(stdout,"  fRef         = %.2f\n",fRef);
  fprintf(stdout,"  distance     = %.2f Mpc\n",distance/(LAL_PC_SI*1000000.));
  fprintf(stdout,"  inclination  = %.2f\n",inclination);
  fprintf(stdout,"  Mtot         = %.2f\n",Mtot);
  fprintf(stdout,"  eta          = %.2f\n",eta);
  fprintf(stdout,"  lambda1      = %.2e\n",lambda1);
  fprintf(stdout,"  lambda2      = %.2e\n",lambda2);
  fprintf(stdout,"  deltaT       = %.2f\n\n",deltaT);

  fprintf(stdout,"--- submodel check ---\n");
  fprintf(stdout,"  cvec_amp size = %d\n",(int) (submodel->cvec_amp->size));
  fprintf(stdout,"  cvec_phi size = %d\n",(int) (submodel->cvec_phi->size));
  fprintf(stdout,"  Bamp size     = %dx%d\n",(int) (submodel->Bamp->size1),(int) (submodel->Bamp->size2));
  fprintf(stdout,"  Bphi size     = %dx%d\n",(int) (submodel->Bphi->size1),(int) (submodel->Bphi->size2));
  fprintf(stdout,"  times size    = %d\n",(int) (submodel->times->size));
  fprintf(stdout,"  q-min         = %.2f\n",submodel->params_min[0]);
  fprintf(stdout,"  q-max         = %.2f\n",submodel->params_max[0]);
  fprintf(stdout,"  lambda1-min   = %.2f\n",submodel->params_min[1]);
  fprintf(stdout,"  lambda1-max   = %.2f\n",submodel->params_max[1]);
  fprintf(stdout,"  lambda2-min   = %.2f\n",submodel->params_min[2]);
  fprintf(stdout,"  lambda2-max   = %.2f\n",submodel->params_max[2]);

  double x = sqrt(1.0-4.0*eta) ;
  double q = (1-x)/(1+x);
  fprintf(stdout,"  q             = %.2f\n\n",q);

  //Allocate space for the nodes
  gsl_vector *amp_at_nodes = gsl_vector_alloc(submodel->times->size);
  gsl_vector *phi_at_nodes = gsl_vector_alloc(submodel->times->size);

  double *amp_interp = calloc(Gntimes,sizeof(double));
  double *phi_interp = calloc(Gntimes,sizeof(double));
  double *freqs = calloc(Gntimes,sizeof(double));
  double *physical_times = calloc(Gntimes,sizeof(double));

  fprintf(stdout,"--- calculating chebyshev interpolants (A(T_j),Phi(T_j)) ---\n");
  retcode = chebyshev_interpolation3d(q,lambda1,lambda2,
                                      Gnq, Gnlambda1, Gnlambda2,
                                      submodel->cvec_amp,submodel->cvec_phi,
                                      Gnamp,Gnphase,amp_at_nodes,phi_at_nodes);

  fprintf(stdout,"\n--- calculating A(t) and Phi(t) at nodes ---\n");
  double BjAmp_tn=0.0;
  double BjPhi_tn=0.0;
  int n,j;
  double c3 = LAL_C_SI*LAL_C_SI*LAL_C_SI ;
  double time_to_phys = LAL_G_SI*Mtot*LAL_MSUN_SI/c3 ;
  FILE *out1 = fopen("./ampinterp_phiinterp_times.txt","w");
  for (n=0;n<Gntimes;n++){
    BjAmp_tn=0.0 ;
    BjPhi_tn=0.0 ;
    for (j=0;j<Gnamp;j++){
      BjAmp_tn+=gsl_vector_get(amp_at_nodes,j)*gsl_matrix_get(submodel->Bamp,j,n);
    }
    for (j=0;j<Gnphase;j++){
      BjPhi_tn+=gsl_vector_get(phi_at_nodes,j)*gsl_matrix_get(submodel->Bphi,j,n);
    }
    //convert times in to seconds
    physical_times[n]=gsl_vector_get(submodel->times,n)*time_to_phys;
    amp_interp[n]=BjAmp_tn;
    phi_interp[n]=BjPhi_tn;
    fprintf(out1,"%.9e %.9e %.9e\n",amp_interp[n],phi_interp[n],physical_times[n]);
  }
  fclose(out1);



  fprintf(stdout,"--- Resampling A(t) and Phi(t) to arbitrary deltaT ---\n");
  gsl_interp_accel *acc = gsl_interp_accel_alloc();
  gsl_spline *ampoft_spline = gsl_spline_alloc (gsl_interp_cspline, Gntimes);
  gsl_spline *phioft_spline = gsl_spline_alloc (gsl_interp_cspline, Gntimes);

  fprintf(stdout,"    -.- initializing splines\n");
  gsl_spline_init(ampoft_spline, physical_times, amp_interp, Gntimes);
  gsl_spline_init(phioft_spline, physical_times, phi_interp, Gntimes);

  double der ;
  //int i_end_mono ;
  fprintf(stdout,"    -.- calculating frequencies at nodes (derivative of phi(t)/2pi)\n");
  int i_end_mono = Gntimes;
  FILE *out2 = fopen("./times_freqs.txt","w");
  for (n=0;n<Gntimes;n++) {
    der = gsl_spline_eval_deriv (phioft_spline, physical_times[n], acc);
    freqs[n] = 0.5*der/LAL_PI;//omegaoft(time_phys)/(2*np.pi)
    fprintf(out2,"%.9e %.9e\n",physical_times[n],freqs[n]);
    //determine up to where f is monotonically increasing
    if (n > 0) {
      if (freqs[n] < freqs[n-1]) i_end_mono = n ;
    }
  }
  fclose(out2);
  fprintf(stdout,"         * i_end_mono %d\n",i_end_mono);

  fprintf(stdout,"    -.- creating t(f) spline\n");
  gsl_spline *toffreq_spline = gsl_spline_alloc (gsl_interp_cspline, i_end_mono);
  //construct t(f)
  gsl_spline_init(toffreq_spline, freqs, physical_times, i_end_mono);

  fprintf(stdout,"    -.- calculate parameters to resample with even spacing\n");
  double tstart = gsl_spline_eval(toffreq_spline, fRef, acc);
  fprintf(stdout,"         * tstart     = %.2f\n",tstart);
  int Ntimes_res = (int) ceil((physical_times[Gntimes-1]-tstart)/deltaT);
  fprintf(stdout,"         * Ntimes_res = %d\n",Ntimes_res);
  double *times_res = calloc(Ntimes_res,sizeof(double));
  double *amp_res = calloc(Ntimes_res,sizeof(double));
  double *phi_res = calloc(Ntimes_res,sizeof(double));
  double t=tstart;

  //for scaling the amplitude
  double h22_to_h = 4.0*eta*sqrt(5.0/LAL_PI)/8.0;
  double amp_units = LAL_G_SI*Mtot*LAL_MSUN_SI/(LAL_C_SI*LAL_C_SI*distance) ;

  //Adjust for inclination angle [0,pi]
  double cosi = cos(inclination);
  double inc_plus = (1.0+cosi*cosi)/2.0;
  double inc_cross = cosi;


  fprintf(stdout,"--- Generate h+(t) and hx(t) ---\n");

  //XLALGPSAdd(&tC, -1 / deltaF);  /* coalesce at t=0 */
  LIGOTimeGPS tC = LIGOTIMEGPSZERO;
  //XLALGPSAdd(&tC, -1.0*j*deltaT);
  /* Allocate hplus and hcross */
  hp = XLALCreateREAL8TimeSeries("hplus: TD waveform", &tC, 0.0, deltaT, &lalStrainUnit, Ntimes_res);
  if (!hp) XLAL_ERROR(XLAL_EFUNC);
  memset(hp->data->data, 0, Ntimes_res * sizeof(REAL8));

  hc = XLALCreateREAL8TimeSeries("hcross: TD waveform", &tC, 0.0, deltaT, &lalStrainUnit, Ntimes_res);
  if (!hc) XLAL_ERROR(XLAL_EFUNC);
  memset(hc->data->data, 0, Ntimes_res * sizeof(REAL8));

  FILE *out3 = fopen("./test_out.txt","w");
  times_res[0] = t ;
  amp_res[0] = gsl_spline_eval(ampoft_spline, t, acc)*amp_units*h22_to_h;
  double phi0 = gsl_spline_eval(phioft_spline, t, acc);
  phi_res[0] = 0.0;
  hp->data->data[0] = inc_plus*amp_res[0]*cos(phi_res[0]);
  hc->data->data[0] = inc_cross*amp_res[0]*sin(phi_res[0]);
  fprintf(out3,"%.9e %.9e %.9e %.9e %.9e\n",t,amp_res[0],phi_res[0],hp->data->data[0],hc->data->data[0]);
  t+=deltaT;
  for (n=1;n<Ntimes_res;n++) {
    times_res[n] = t;
    amp_res[n] = gsl_spline_eval(ampoft_spline, t, acc)*amp_units*h22_to_h;
    //Zero the phase at the beginning (-phi0)
    phi_res[n] = gsl_spline_eval(phioft_spline, t, acc)-phi0;

    hp->data->data[n] = inc_plus*amp_res[n]*cos(phi_res[n]);
    hc->data->data[n] = inc_cross*amp_res[n]*sin(phi_res[n]);

    fprintf(out3,"%.9e %.9e %.9e %.9e %.9e\n",t,amp_res[n],phi_res[n],hp->data->data[n],hc->data->data[n]);

    t+=deltaT;
  }
  fclose(out3);

  *hPlus = hp;
  *hCross = hc;

  gsl_spline_free (ampoft_spline);
  gsl_spline_free (phioft_spline);
  gsl_spline_free (toffreq_spline);
  gsl_interp_accel_free (acc);

  gsl_vector_free(amp_at_nodes);
  gsl_vector_free(phi_at_nodes);

  free(amp_interp);
  free(phi_interp);
  free(freqs);
  free(physical_times);
  free(times_res);
  free(amp_res);
  free(phi_res);

  if (retcode==0){
    return(XLAL_SUCCESS);
  } else {
    return(retcode);
  }

}

static int TEOBResumROMCore(
  REAL8TimeSeries **hPlus,
  REAL8TimeSeries **hCross,
  double phiRef, // orbital reference phase NOTE: unused
  double deltaT,
  double fRef,
  double distance,
  double inclination,
  double Mtot, // in Msol
  double eta,
  double lambda1, //in range 50 - 5000
  double lambda2  //in range 50 - 5000
)
{

  //NOTE: silly
  inclination = inclination + phiRef - phiRef ;

  /* Check output arrays */
  if(!hPlus || !hCross)
    XLAL_ERROR(XLAL_EFAULT);
  TEOBResumROMdataDS *romdata=&__lalsim_TEOBResumROMDS_data;
  if(*hPlus || *hCross)
  {
    XLALPrintError("(*hPlus) and (*hCross) are supposed to be NULL, but got %p and %p",(*hPlus),(*hCross));
    XLAL_ERROR(XLAL_EFAULT);
  }
  int retcode=0;

  REAL8TimeSeries *hp;
  REAL8TimeSeries *hc;

  /* Select ROM submodel */
  TEOBResumROMdataDS_submodel *submodel;
  submodel = romdata->sub1;

  double x = sqrt(1.0-4.0*eta) ;
  double q = (1-x)/(1+x);

  //Allocate space for the nodes
  gsl_vector *amp_at_nodes = gsl_vector_alloc(submodel->times->size);
  gsl_vector *phi_at_nodes = gsl_vector_alloc(submodel->times->size);

  double *amp_interp = calloc(Gntimes,sizeof(double));
  double *phi_interp = calloc(Gntimes,sizeof(double));
  double *freqs = calloc(Gntimes,sizeof(double));
  double *physical_times = calloc(Gntimes,sizeof(double));

  //calculate chebyshev interpolants (A(T_j),Phi(T_j))
  retcode = chebyshev_interpolation3d(q,lambda1,lambda2,
                                      Gnq, Gnlambda1, Gnlambda2,
                                      submodel->cvec_amp,submodel->cvec_phi,
                                      Gnamp,Gnphase,amp_at_nodes,phi_at_nodes);

  //calculate A(T_j) and Phi(T_j) at nodes
  double BjAmp_tn=0.0;
  double BjPhi_tn=0.0;
  int n,j;
  double c3 = LAL_C_SI*LAL_C_SI*LAL_C_SI ;
  double time_to_phys = LAL_G_SI*Mtot*LAL_MSUN_SI/c3 ;
  for (n=0;n<Gntimes;n++){
    BjAmp_tn=0.0 ;
    BjPhi_tn=0.0 ;
    for (j=0;j<Gnamp;j++){
      BjAmp_tn+=gsl_vector_get(amp_at_nodes,j)*gsl_matrix_get(submodel->Bamp,j,n);
    }
    for (j=0;j<Gnphase;j++){
      BjPhi_tn+=gsl_vector_get(phi_at_nodes,j)*gsl_matrix_get(submodel->Bphi,j,n);
    }
    //convert times in to seconds
    physical_times[n]=gsl_vector_get(submodel->times,n)*time_to_phys;
    amp_interp[n]=BjAmp_tn;
    phi_interp[n]=BjPhi_tn;
  }



  //Resampling A(t) and Phi(t) to arbitrary deltaT ---\n");
  gsl_interp_accel *acc = gsl_interp_accel_alloc();
  gsl_spline *ampoft_spline = gsl_spline_alloc (gsl_interp_cspline, Gntimes);
  gsl_spline *phioft_spline = gsl_spline_alloc (gsl_interp_cspline, Gntimes);

  //initializing splines
  gsl_spline_init(ampoft_spline, physical_times, amp_interp, Gntimes);
  gsl_spline_init(phioft_spline, physical_times, phi_interp, Gntimes);

  double der ;
  //calculate frequencies at nodes (derivative of phi(t)/2pi)
  int i_end_mono = Gntimes;
  for (n=0;n<Gntimes;n++) {
    der = gsl_spline_eval_deriv (phioft_spline, physical_times[n], acc);
    freqs[n] = 0.5*der/LAL_PI;//omegaoft(time_phys)/(2*np.pi)
    //determine up to where f is monotonically increasing
    if (n > 0) {
      if (freqs[n] < freqs[n-1]) i_end_mono = n ;
    }
  }

  //create t(f) spline
  gsl_spline *toffreq_spline = gsl_spline_alloc (gsl_interp_cspline, i_end_mono);
  gsl_spline_init(toffreq_spline, freqs, physical_times, i_end_mono);

  //calculate parameters to resample with even spacing
  double tstart = gsl_spline_eval(toffreq_spline, fRef, acc);
  int Ntimes_res = (int) ceil((physical_times[Gntimes-1]-tstart)/deltaT);
  double *times_res = calloc(Ntimes_res,sizeof(double));
  double *amp_res = calloc(Ntimes_res,sizeof(double));
  double *phi_res = calloc(Ntimes_res,sizeof(double));
  double t=tstart;

  //for scaling the amplitude
  double h22_to_h = 4.0*eta*sqrt(5.0/LAL_PI)/8.0;
  double amp_units = LAL_G_SI*Mtot*LAL_MSUN_SI/(LAL_C_SI*LAL_C_SI*distance) ;

  //Adjust for inclination angle [0,pi]
  double cosi = cos(inclination);
  double inc_plus = (1.0+cosi*cosi)/2.0;
  double inc_cross = cosi;


  //Generate h+(t) and hx(t)

  //XLALGPSAdd(&tC, -1 / deltaF);  /* coalesce at t=0 */
  LIGOTimeGPS tC = LIGOTIMEGPSZERO;
  //XLALGPSAdd(&tC, tstart);
  //XLALGPSAdd(&tC, -1.0*j*deltaT);
  /* Allocate hplus and hcross */
  hp = XLALCreateREAL8TimeSeries("hplus: TD waveform", &tC, 0.0, deltaT, &lalStrainUnit, Ntimes_res);
  if (!hp) XLAL_ERROR(XLAL_EFUNC);
  memset(hp->data->data, 0, Ntimes_res * sizeof(REAL8));

  hc = XLALCreateREAL8TimeSeries("hcross: TD waveform", &tC, 0.0, deltaT, &lalStrainUnit, Ntimes_res);
  if (!hc) XLAL_ERROR(XLAL_EFUNC);
  memset(hc->data->data, 0, Ntimes_res * sizeof(REAL8));

  times_res[0] = t ;
  amp_res[0] = gsl_spline_eval(ampoft_spline, t, acc)*amp_units*h22_to_h;
  double phi0 = gsl_spline_eval(phioft_spline, t, acc);
  phi_res[0] = 0.0;
  hp->data->data[0] = inc_plus*amp_res[0]*cos(phi_res[0]);
  hc->data->data[0] = inc_cross*amp_res[0]*sin(phi_res[0]);
  t+=deltaT;
  for (n=1;n<Ntimes_res;n++) {
    times_res[n] = t;
    amp_res[n] = gsl_spline_eval(ampoft_spline, t, acc)*amp_units*h22_to_h;
    //Zero the phase at the beginning (-phi0)
    phi_res[n] = gsl_spline_eval(phioft_spline, t, acc)-phi0;

    hp->data->data[n] = inc_plus*amp_res[n]*cos(phi_res[n]);
    hc->data->data[n] = inc_cross*amp_res[n]*sin(phi_res[n]);

    t+=deltaT;
  }

  *hPlus = hp;
  *hCross = hc;

  gsl_spline_free (ampoft_spline);
  gsl_spline_free (phioft_spline);
  gsl_spline_free (toffreq_spline);
  gsl_interp_accel_free (acc);

  gsl_vector_free(amp_at_nodes);
  gsl_vector_free(phi_at_nodes);

  free(amp_interp);
  free(phi_interp);
  free(freqs);
  free(physical_times);
  free(times_res);
  free(amp_res);
  free(phi_res);

  if (retcode==0){
    return(XLAL_SUCCESS);
  } else {
    return(retcode);
  }

}

/**
 * @addtogroup LALSimIMRSEOBNRROM_c
 *
 * @{
 *
 * @name SEOBNRv2 Reduced Order Model (Double Spin)
 *
 * @author Michael Puerrer, John Veitch
 *
 * @brief C code for SEOBNRv2 reduced order model (double spin version).
 * See CQG 31 195010, 2014, arXiv:1402.4146 for the basic approach.
 * Further details in PRD 93, 064041, 2016, arXiv:1512.02248.
 *
 * This is a frequency domain model that approximates the time domain SEOBNRv2 model.
 *
 * The binary data files are available at https://dcc.ligo.org/T1400701-v1.
 * Put the untared data into a location in your LAL_DATA_PATH.
 *
 * @note Note that due to its construction the iFFT of the ROM has a small (~ 20 M) offset
 * in the peak time that scales with total mass as compared to the time-domain SEOBNRv2 model.
 *
 * @note Parameter ranges:
 *   * 0.01 <= eta <= 0.25
 *   * -1 <= chi_i <= 0.99
 *   * Mtot >= 12Msun
 *
 *  Aligned component spins chi1, chi2.
 *  Symmetric mass-ratio eta = m1*m2/(m1+m2)^2.
 *  Total mass Mtot.
 *
 * @{
 */


//Function to test data reading
int XLALSimInspiralTEOBResumROM(
  REAL8TimeSeries **hPlus, /**< Output: Frequency-domain waveform h+ */
  REAL8TimeSeries **hCross, /**< Output: Frequency-domain waveform hx */
  REAL8 phiRef,                                 /**< Orbital phase at reference frequency*/
  REAL8 deltaT,                                 /**< Sampling frequency (Hz) */
  REAL8 fLow,                                   /**< Starting GW frequency (Hz) */
  REAL8 fRef,                                   /**< Reference frequency (Hz); 0 defaults to fLow */
  REAL8 distance,                               /**< Distance of source (m) */
  REAL8 inclination,                            /**< Inclination of source (rad) */
  REAL8 m1SI,                                   /**< Mass of companion 1 (kg) */
  REAL8 m2SI,                                   /**< Mass of companion 2 (kg) */
  REAL8 lambda1,                                /**< dimensionless tidal deformability of body 1 */
  REAL8 lambda2)                                /**< dimensionless tidal deformability of body 1 */
{
  /* Internally we need m1 > m2, so change around if this is not the case */
  if (m1SI < m2SI) {
    // Swap m1 and m2
    double m1temp = m1SI;
    double l1temp = lambda1;
    m1SI = m2SI;
    lambda1 = lambda2;
    m2SI = m1temp;
    lambda2 = l1temp;
  }

  /* Get masses in terms of solar mass */
  double mass1 = m1SI / LAL_MSUN_SI;
  double mass2 = m2SI / LAL_MSUN_SI;
  double Mtot = mass1+mass2;
  double eta = mass1 * mass2 / (Mtot*Mtot);    /* Symmetric mass-ratio */

//   REAL8 m1sec = (m1SI/LAL_MSUN_SI)*LAL_MTSUN_SI ;
//   REAL8 m2sec = (m2SI/LAL_MSUN_SI)*LAL_MTSUN_SI ;
//   REAL8 m1sec5=m1sec*m1sec*m1sec*m1sec*m1sec ;
//   REAL8 m2sec5=m2sec*m2sec*m2sec*m2sec*m2sec ;
//   lambda1*=m1sec5 ;
//   lambda2*=m2sec5 ;

  if (fRef==0.0) fRef=fLow;

  // Load ROM data if not loaded already
  fprintf(stdout,"initializing with TEOBResumROM_Init_LALDATA()\n");
  #ifdef LAL_PTHREAD_LOCK
  (void) pthread_once(&TEOBResumROM_is_initialized, TEOBResumROM_Init_LALDATA);
  #else
  TEOBResumROM_Init_LALDATA();
  #endif

  if(!TEOBResumROM_IsSetup()) XLAL_ERROR(XLAL_EFAILED,"Error setting up TEOBResumROM data - check your $LAL_DATA_PATH\n");

  int retcode = TEOBResumROMCore(hPlus,hCross, phiRef, deltaT, fRef, distance, inclination, Mtot, eta, lambda1, lambda2);

  return(retcode);
}



//Function to test data reading
int XLALSimInspiralTEOBResumROM_test(
  REAL8TimeSeries **hPlus, /**< Output: Frequency-domain waveform h+ */
  REAL8TimeSeries **hCross, /**< Output: Frequency-domain waveform hx */
  REAL8 phiRef,                                 /**< Orbital phase at reference frequency*/
  REAL8 deltaT,                                 /**< Sampling frequency (Hz) */
  REAL8 fLow,                                   /**< Starting GW frequency (Hz) */
  REAL8 fRef,                                   /**< Reference frequency (Hz); 0 defaults to fLow */
  REAL8 distance,                               /**< Distance of source (m) */
  REAL8 inclination,                            /**< Inclination of source (rad) */
  REAL8 m1SI,                                   /**< Mass of companion 1 (kg) */
  REAL8 m2SI,                                   /**< Mass of companion 2 (kg) */
  REAL8 lambda1,                                /**< Dimensionless aligned component spin 1 */
  REAL8 lambda2)                                /**< Dimensionless aligned component spin 2 */
{
  /* Internally we need m1 > m2, so change around if this is not the case */
  if (m1SI < m2SI) {
    // Swap m1 and m2
    double m1temp = m1SI;
    double l1temp = lambda1;
    m1SI = m2SI;
    lambda1 = lambda2;
    m2SI = m1temp;
    lambda2 = l1temp;
  }

  /* Get masses in terms of solar mass */
  double mass1 = m1SI / LAL_MSUN_SI;
  double mass2 = m2SI / LAL_MSUN_SI;
  double Mtot = mass1+mass2;
  double eta = mass1 * mass2 / (Mtot*Mtot);    /* Symmetric mass-ratio */

  REAL8 m1sec = (m1SI/LAL_MSUN_SI)*LAL_MTSUN_SI ;
  REAL8 m2sec = (m2SI/LAL_MSUN_SI)*LAL_MTSUN_SI ;
  REAL8 m1sec5=m1sec*m1sec*m1sec*m1sec*m1sec ;
  REAL8 m2sec5=m2sec*m2sec*m2sec*m2sec*m2sec ;
  lambda1*=m1sec5 ;
  lambda2*=m2sec5 ;

  if (fRef==0.0) fRef=fLow;

  // Load ROM data if not loaded already
  fprintf(stdout,"initializing with TEOBResumROM_Init_LALDATA()\n");
  #ifdef LAL_PTHREAD_LOCK
  (void) pthread_once(&TEOBResumROM_is_initialized, TEOBResumROM_Init_LALDATA);
  #else
  TEOBResumROM_Init_LALDATA();
  #endif

  if(!TEOBResumROM_IsSetup()) XLAL_ERROR(XLAL_EFAILED,"Error setting up TEOBResumROM data - check your $LAL_DATA_PATH\n");

  // Use fLow, fHigh, deltaF to compute freqs sequence
  // Instead of building a full sequency we only transfer the boundaries and let
  // the internal core function do the rest (and properly take care of corner cases).

  int retcode = TEOBResumROMCore_test(hPlus,hCross, phiRef, deltaT, fRef, distance, inclination, Mtot, eta, lambda1, lambda2);

  TEOBResumROMdataDS_Cleanup(&__lalsim_TEOBResumROMDS_data);

  return(retcode);
}



/** Setup SEOBNRv2ROMDoubleSpin model using data files installed in $LAL_DATA_PATH
 */
static void TEOBResumROM_Init_LALDATA(void)
{
  if (TEOBResumROM_IsSetup()) return;

  // If we find one ROM datafile in a directory listed in LAL_DATA_PATH,
  // then we expect the remaining datafiles to also be there.
  char datafile[] = "TEOBResumROM_Phase_ciall.dat";

  char *path = XLALFileResolvePathLong(datafile, PKG_DATA_DIR);
  if (path==NULL)
    XLAL_ERROR_VOID(XLAL_EIO, "Unable to resolve data file %s in $LAL_DATA_PATH\n", datafile);
  char *dir = dirname(path);
  int ret = TEOBResumROM_Init(dir);
  XLALFree(path);

  if(ret!=XLAL_SUCCESS)
    XLAL_ERROR_VOID(XLAL_FAILURE, "Unable to find TEOBResumROM data files in $LAL_DATA_PATH\n");
}
