#include <lal/LALSimInspiral.h>
#include <lal/LALSimIMR.h>
#include <gsl/gsl_spline.h>
#include <math.h>

#include "LALSimIMREOBNRv2.h"

#ifndef _LALSIMIMRSPINEOB_H
#define _LALSIMIMRSPINEOB_H

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif



/**
 * Set the total number of multipoles
 * */
#define MAX_NUM_MODES 7



struct
SpinEOBModes
{
  INT4 lmModes[MAX_NUM_MODES][2];
  REAL8TimeSeries *hlposm[MAX_NUM_MODES];
  REAL8TimeSeries *hlnegm[MAX_NUM_MODES];

  struct SpinEOBModes *next;
};

/**
 * Parameters for the spinning EOB model, used in calculating the Hamiltonian.
 * The Hamiltonian is given in Barausse and Buonanno (http://arxiv.org/pdf/0912.3517)
 * The parameters correspond to the following as found in the paper:
 * KK - K found in the equations for \f$\Delta_i\f$ (Eqn 5.77-5.81)
 * k0 - \f$\Delta_0\f$ given in Eqn 5.77
 * k1 - \f$\Delta_1\f$ given in Eqn 5.78
 * k2 - \f$\Delta_2\f$ given in Eqn 5.79
 * k3 - \f$\Delta_3\f$ given in Eqn 5.80
 * k4 - \f$\Delta_4\f$ given in Eqn 5.81
 * b3 - \f$\omega^{fd}_2\f$ given in Eqn 5.40
 * bb3 - \f$\omega^{fd}_1\f$ given in Eqn 5.40
 * d1 - SO calibration parameter of SEOBNRv1
 * d1v2 - SO calibration parameter of SEOBNRv2
 * dheffSS - SS calibration parameter of SEOBNRv1
 * dheffSSv2 - SS calibration parameter of SEOBNRv2
 */

typedef struct
tagSpinEOBHCoeffs
{
  double KK;
  double k0;
  double k1;
  double k2;
  double k3;
  double k4;
  double k5;
  double k5l;
  double b3;
  double bb3;
  double d1;
  double d1v2;
  double dheffSS;
  double dheffSSv2;
  UINT4    SpinAlignedEOBversion;
  int      updateHCoeffs;
}
SpinEOBHCoeffs;

/**
 * Parameters for the spinning EOB model.
 * 1) eobParams contains parameters common to nonspin and spin EOBNR models,
 * including mass ratio, masses, pre-computed coefficients for potential, flux and waveforms,
 * NQC coefficients and Newtonian multiple prefixes.
 * 2) seobCoeffs contans parameters for calculating the spin EOB Hamiltonian.
 * 3) s1Vec and s2Vec are individual spins, in unit of total mass.
 * 4) sigmaStar and sigmaKerr are effective spins of the test-particle and background.
 * 5) a is the spin value being used for test-particle limit spin terms.
 * 6) alignedSpins and tortoise are controling flags.
 */

typedef struct
tagSpinEOBParams
{
  EOBParams               *eobParams;
  SpinEOBHCoeffs          *seobCoeffs;
  EOBNonQCCoeffs          *nqcCoeffs;
  REAL8Vector             *s1Vec;
  REAL8Vector             *s2Vec;
  REAL8Vector             *sigmaStar;
  REAL8Vector             *sigmaKerr;
  REAL8                   a;
  REAL8                   chi1;
  REAL8                   chi2;
  int                     alignedSpins;
  int                     tortoise;
}
SpinEOBParams;

/* We need to encapsulate the data for the GSL derivative function */
typedef
struct tagHcapDerivParams
{
   const REAL8   *values;
   SpinEOBParams *params;
   UINT4         varyParam;
}
HcapDerivParams;

/* We need to encapsulate the data for calculating spherical 2nd derivatives */
typedef
struct tagHcapSphDeriv2Params
{
  const REAL8     *sphValues;
  SpinEOBParams   *params;
  UINT4           varyParam1;
  UINT4           varyParam2;
}
HcapSphDeriv2Params;

/* We need to encapsulate the data for the GSL derivative function */
typedef
struct tagPrecEulerAnglesIntegration
{
   gsl_spline *alpha_spline;
   gsl_spline *beta_spline;
   gsl_interp_accel *alpha_acc;
   gsl_interp_accel *beta_acc;
}
PrecEulerAnglesIntegration;


//
REAL8 XLALCalculateEOBA( 
        const REAL8 r,                     /**<< Orbital separation (in units of total mass M) */
        EOBACoefficients * restrict coeffs /**<< Pre-computed coefficients for the A function */
 );




REAL8  XLALEffectiveHamiltonian( const REAL8 eta,          /**<< Symmetric mass ratio */
                                const REAL8 r,            /**<< Orbital separation */
                                const REAL8 pr,           /**<< Tortoise co-ordinate */
                                const REAL8 pp,           /**<< Momentum pphi */
                                EOBACoefficients *aCoeffs /**<< Pre-computed coefficients in A function */
                              );



typedef struct
{
   /* coefficients in the Pade expression of new energy function */
   REAL8 ePaN, ePa1, ePa2, ePa3;
   /* coefficients in the Taylor expansion of usual energy function */
   REAL8 ETaN, ETa1, ETa2, ETa3, ETa5, ETa6;
   /* coefficients in the Taylor expansion of the derivative of the
    usual energy function */
   REAL8 dETaN, dETa1, dETa2, dETa3, dETa5, dETa6;

   /* Taylor expansion coefficients of energy flux*/
   REAL8 FTaN, FTa1, FTa2, FTa3, FTa4, FTa5, FTa6, FTa7, FTa8, FTl6, FTl8, FTa10, FTa12;
   /* Coefficients of the corresponding P-approximant */
   REAL8 fPaN, fPa1, fPa2, fPa3, fPa4, fPa5, fPa6, fPa7, fPa8;

   /* symmetric mass ratio, total mass, component masses */
   REAL8 eta, totalmass;

   /* initial and final values of frequency, time, velocity; lso
    values of velocity and frequency; final phase. */
   REAL8 vlso;

   /* last stable orbit and pole defined by various Taylor and P-approximants */
   REAL8 vlsoT0, vlsoT2, vlsoT4, vlsoT6;
   REAL8 vlsoP0, vlsoP2, vlsoP4, vlsoP6;
   REAL8 vlsoPP;
   REAL8 vpoleP4, vpoleP6;
   REAL8 vpolePP;
}  expnCoeffsdEnergyFlux;


typedef REAL8 (*EnergyFunction)(
   REAL8 v,
   expnCoeffsdEnergyFlux *ak);


typedef REAL8 (*FluxFunction)(
   REAL8 v,
   expnCoeffsdEnergyFlux *ak);


typedef struct
{
   EnergyFunction dEnergy;
   FluxFunction flux;
   expnCoeffsdEnergyFlux *coeffs;
} TofVIntegrandIn;


typedef struct
{
   REAL8 t;
   REAL8 v0;
   REAL8 t0;
   REAL8 vlso;
   REAL8 totalmass;
   EnergyFunction dEnergy;
   FluxFunction flux;
   expnCoeffsdEnergyFlux *coeffs;
} TofVIn;

/**
 * Structure consisting SEOBNR parameters that can be used by gsl root finders
 */
typedef
struct tagSEOBRootParams
{
  REAL8          values[12]; /**<< Dynamical variables, x, y, z, px, py, pz, S1x, S1y, S1z, S2x, S2y and S2z */
  SpinEOBParams *params;     /**<< Spin EOB parameters -- physical, pre-computed, etc. */
  REAL8          omega;      /**<< Orbital frequency */
}
SEOBRootParams;




 REAL8  XLALSimInspiralep2(REAL8 v, expnCoeffsdEnergyFlux *ak);

 REAL8  XLALSimInspiralep4(REAL8 v, expnCoeffsdEnergyFlux *ak);

 REAL8  XLALSimInspiralep6(REAL8 v, expnCoeffsdEnergyFlux *ak);

 REAL8  XLALSimInspiraldEp2(REAL8 v, expnCoeffsdEnergyFlux *ak);

 REAL8  XLALSimInspiraldEp4(REAL8 v, expnCoeffsdEnergyFlux *ak);

 REAL8  XLALSimInspiraldEp6(REAL8 v, expnCoeffsdEnergyFlux *ak);

 REAL8  XLALSimInspiralFp3(REAL8 v, expnCoeffsdEnergyFlux *ak);

 REAL8  XLALSimInspiralFp4(REAL8 v, expnCoeffsdEnergyFlux *ak);

 REAL8  XLALSimInspiralFp5(REAL8 v, expnCoeffsdEnergyFlux *ak);

 REAL8  XLALSimInspiralFp6(REAL8 v, expnCoeffsdEnergyFlux *ak);

 REAL8  XLALSimInspiralFp7(REAL8 v, expnCoeffsdEnergyFlux *ak);

 REAL8  XLALSimInspiralFp8PP(REAL8 v, expnCoeffsdEnergyFlux *ak);




int XLALSimIMRSpinEOBInitialConditions(
                      REAL8Vector   *initConds, /**<< OUTPUT, Initial dynamical variables */
                      const REAL8    mass1,     /**<< mass 1 */
                      const REAL8    mass2,     /**<< mass 2 */
                      const REAL8    fMin,      /**<< Initial frequency (given) */
                      const REAL8    inc,       /**<< Inclination */
                      const REAL8    spin1[],   /**<< Initial spin vector 1 */
                      const REAL8    spin2[],   /**<< Initial spin vector 2 */
                      SpinEOBParams *params     /**<< Spin EOB parameters */
                      );




 int  XLALSimIMREOBComputeNewtonMultipolePrefixes(
                NewtonMultipolePrefixes *prefix, /**<< OUTPUT Structure containing the coeffs */
                const REAL8             m1,      /**<< Mass of first component */
                const REAL8             m2       /**<< Nass of second component */
                );


  int 
XLALSimIMREOBCalculateNewtonianMultipole(
                 COMPLEX16 *multipole, /**<< OUTPUT, Newtonian multipole */
                 REAL8 x,              /**<< Dimensionless parameter \f$\equiv v^2\f$ */
                  REAL8 r,       /**<< Orbital separation (units of total mass M) */
                 REAL8 phi,            /**<< Orbital phase (in radians) */
                 UINT4  l,             /**<< Mode l */
                 INT4  m,              /**<< Mode m */
                 EOBParams *params     /**<< Pre-computed coefficients, parameters, etc. */
                 );


  int
XLALSimIMRSpinEOBCalculateNewtonianMultipole(
                 COMPLEX16 *multipole, /**<< OUTPUT, Newtonian multipole */
                 REAL8 x,              /**<< Dimensionless parameter \f$\equiv v^2\f$ */
                  REAL8 r,       /**<< Orbital separation (units of total mass M */
                 REAL8 phi,            /**<< Orbital phase (in radians) */
                 UINT4  l,             /**<< Mode l */
                 INT4  m,              /**<< Mode m */
                 EOBParams *params     /**<< Pre-computed coefficients, parameters, etc. */
                 );


  int
XLALSimIMRSpinEOBFluxCalculateNewtonianMultipole(
                 COMPLEX16 *multipole, /**<< OUTPUT, Newtonian multipole */
                 REAL8 x,              /**<< Dimensionless parameter \f$\equiv v^2\f$ */
                  REAL8 r,       /**<< Orbital separation (units of total mass M */
                  REAL8 phi,     /**<< Orbital phase (in radians) */
                 UINT4  l,             /**<< Mode l */
                 INT4  m,              /**<< Mode m */
                 EOBParams *params     /**<< Pre-computed coefficients, parameters, etc. */
                 );





 REAL8   XLALSimIMREOBGetNRPeakDeltaT( 
                         INT4 l,    /**<< Mode l */ 
                         INT4 m,    /**<< Mode m */
                         REAL8 eta  /**<< Symmetric mass ratio */
                         );

  int XLALSimIMREOBGetCalibratedNQCCoeffs( 
                                EOBNonQCCoeffs *coeffs, /**<< OUTPUT, Structure for NQC coeffs */
                                INT4            l,      /**<< Mode l */
                                INT4            m,      /**<< Mode m */
                                REAL8           eta     /**<< Symmetric mass ratio */
                                );


  int  XLALSimIMREOBNonQCCorrection(
                      COMPLEX16      * restrict nqc,    /**<< OUTPUT, The NQC correction */
                      REAL8Vector    * restrict values, /**<< Dynamics r, phi, pr, pphi */
                      const REAL8               omega,  /**<< Angular frequency */
                      EOBNonQCCoeffs * restrict coeffs  /**<< NQC coefficients */
                     );


  int  XLALSimIMRSpinEOBNonQCCorrection(
                      COMPLEX16      * restrict nqc,    /**<< OUTPUT, The NQC correction */
                      REAL8Vector    * restrict values, /**<< Dynamics r, phi, pr, pphi */
                      const REAL8               omega,  /**<< Angular frequency */
                      EOBNonQCCoeffs * restrict coeffs  /**<< NQC coefficients */
                     );


  int XLALSimIMREOBCalculateNQCCoefficients(
                 EOBNonQCCoeffs * restrict coeffs,    /**<< OUTPUT, NQC coefficients */
                 REAL8Vector    * restrict amplitude, /**<< Waveform amplitude, func of time */
                 REAL8Vector    * restrict phase,     /**<< Waveform phase(rad), func of time */
                 REAL8Vector    * restrict q1,        /**<< Function of dynamics (see DCC doc) */
                 REAL8Vector    * restrict q2,        /**<< Function of dynamics (see DCC doc) */
                 REAL8Vector    * restrict q3,        /**<< Function of dynamics (see DCC doc) */
                 REAL8Vector    * restrict p1,        /**<< Function of dynamics (see DCC doc) */
                 REAL8Vector    * restrict p2,        /**<< Function of dynamics (see DCC doc) */
                 INT4                      l,         /**<< Mode l */
                 INT4                      m,         /**<< Mode m */
                 REAL8                     timePeak,  /**<< Time of peak orbital frequency */
                 REAL8                     deltaT,    /**<< Sampling interval */
                 REAL8                     eta        /**<< Symmetric mass ratio */
                 );



  int XLALSimIMRGetEOBCalibratedSpinNQC( EOBNonQCCoeffs * restrict coeffs,                                   INT4  l,                                   INT4  m,                                   REAL8 eta,                                   REAL8 a );


  int XLALSimIMRGetEOBCalibratedSpinNQCv2( EOBNonQCCoeffs * restrict coeffs,                                   INT4  l,                                   INT4  m,                                   REAL8 eta,                                   REAL8 a );


  int XLALSimIMRGetEOBCalibratedSpinNQCv2chiAmax( EOBNonQCCoeffs * restrict coeffs,                                                             INT4  l,                                                             INT4  m,                                                             REAL8 eta,                                                             REAL8 a );


  int XLALSimIMRGetEOBCalibratedSpinNQCv2chiAmed( EOBNonQCCoeffs * restrict coeffs,                                                             INT4  l,                                                             INT4  m,                                                             REAL8 eta,                                                             REAL8 a );



  int XLALSimIMRGetEOBCalibratedSpinNQCv2chiAmin( EOBNonQCCoeffs * restrict coeffs,                                                             INT4  l,                                                             INT4  m,                                                             REAL8 eta,                                                             REAL8 a );


  int XLALSimIMRGetEOBCalibratedSpinNQC3D(
                        EOBNonQCCoeffs * restrict coeffs,
                        INT4  l,
                        INT4  m,
                        REAL8 m1,
                        REAL8 m2,
                        REAL8 a,
                        REAL8 chiAin );


  int XLALSimIMRSpinEOBCalculateNQCCoefficients(
                 REAL8Vector    * restrict amplitude,   /**<< Waveform amplitude, func of time */
                 REAL8Vector    * restrict phase,       /**<< Waveform phase(rad), func of time */
                 REAL8Vector    * restrict rVec,        /**<< Position-vector, function of time */
                 REAL8Vector    * restrict prVec,       /**<< Momentum vector, function of time */
                 REAL8Vector    * restrict orbOmegaVec, /**<< Orbital frequency, func of time */
                 INT4                      l,           /**<< Mode index l */
                 INT4                      m,           /**<< Mode index m */
                 REAL8                     timePeak,    /**<< Time of peak orbital frequency */
                 REAL8                     deltaT,      /**<< Sampling interval */
                 REAL8                     m1,          /**<< Component mass 1 */
                 REAL8                     m2,          /**<< Component mass 2 */
                 REAL8                     a,           /**<< Normalized spin of deformed-Kerr */
                 REAL8                     chiA,        /**<< Assymmetric dimensionless spin combination */
                 REAL8                     chiS,        /**<< Symmetric dimensionless spin combination */
                 EOBNonQCCoeffs * restrict coeffs,      /**<< OUTPUT, NQC coefficients */
                 UINT4                     SpinAlignedEOBversion  /**<< 1 for SEOBNRv1, 2 for SEOBNRv2 */
);



/* */
int XLALCalculateEOBACoefficients(
          EOBACoefficients * const coeffs, /**<< A coefficients (populated in function) */
          const REAL8              eta     /**<< Symmetric mass ratio */
          );


/* */
REAL8 XLALCalculateEOBdAdr( 
        const REAL8 r,                     /**<< Orbital separation (in units of total mass M) */
        EOBACoefficients * restrict coeffs /**<< Pre-computed coefficients for the A function */
     );




  int XLALSimIMREOBCalcFacWaveformCoefficients(
          FacWaveformCoeffs * const coeffs, /**<< Structure containing coefficients (populated in function) */
          const REAL8               eta     /**<< Symmetric mass ratio */
          );

  int XLALSimIMREOBModifyFacWaveformCoefficients( 
                FacWaveformCoeffs * const coeffs, /**<< Structure containing coefficients */
                const REAL8 eta                   /**<< Symmetric mass ratio */
        );


 REAL8  
nonKeplerianCoefficient(
        REAL8Vector * restrict values, /**<< Dynamics r, phi, pr, pphi */
        const REAL8       eta,         /**<< Symmetric mass ratio */
        EOBACoefficients *coeffs       /**<< Pre-computed A coefficients */
        );


  int  XLALSimIMREOBGetFactorizedWaveform( 
            COMPLEX16   * restrict hlm,    /**<< The value of hlm (populated by the function) */
            REAL8Vector * restrict values, /**<< Vector containing dynamics r, phi, pr, pphi for a given point */
            const REAL8 v,                 /**<< Velocity (in geometric units) */
            const INT4  l,                 /**<< Mode l */
            const INT4  m,                 /**<< Mode m */
            EOBParams   * restrict params  /**<< Structure containing pre-computed coefficients, etc. */
            );

 REAL8   XLALSimIMREOBFactorizedFlux(
        REAL8Vector  *values, /**<< Dynamics r, phi, pr, pphi */
        const REAL8  omega,   /**<< Angular frequency omega */
        EOBParams    *ak,     /**<< Structure containing pre-computed parameters */
        const INT4   lMax     /**<< Maximum l to include when calculating flux (between 2 and 8) */
        );

 INT4   XLALSimIMREOBHybridRingdownWave(
  REAL8Vector          *rdwave1,   /**<< OUTPUT, Real part of ringdown waveform */
  REAL8Vector          *rdwave2,   /**<< OUTPUT, Imag part of ringdown waveform */
  const REAL8           dt,        /**<< Sampling interval */
  const REAL8           mass1,     /**<< First component mass (in Solar masses) */
  const REAL8           mass2,     /**<< Second component mass (in Solar masses) */
  REAL8VectorSequence  *inspwave1, /**<< Values and derivs of real part inspiral waveform */
  REAL8VectorSequence  *inspwave2, /**<< Values and derivs of imag part inspiral waveform */
  COMPLEX16Vector      *modefreqs, /**<< Complex freqs of ringdown (scaled by total mass) */
  REAL8Vector          *matchrange /**<< Times which determine the comb of ringdown attachment */
  );


 INT4   XLALGenerateHybridWaveDerivatives (
	REAL8Vector	*rwave,      /**<< OUTPUT, values of the waveform at comb points */
	REAL8Vector	*dwave,      /**<< OUTPUT, 1st deriv of the waveform at comb points */
	REAL8Vector	*ddwave,     /**<< OUTPUT, 2nd deriv of the waveform at comb points */
    REAL8Vector	*timeVec,    /**<< Vector containing the time */
	REAL8Vector	*wave,       /**<< Last part of inspiral waveform */
	REAL8Vector	*matchrange, /**<< Times which determine the size of the comb */
    REAL8           dt,          /**<< Sample time step */
    REAL8           mass1,       /**<< First component mass (in Solar masses) */
    REAL8           mass2        /**<< Second component mass (in Solar masses) */
	);

 INT4   XLALSimIMREOBHybridAttachRingdown(
  REAL8Vector *signal1,    /**<< OUTPUT, Real of inspiral waveform to which we attach ringdown */
  REAL8Vector *signal2,    /**<< OUTPUT, Imag of inspiral waveform to which we attach ringdown */
  const INT4   l,          /**<< Current mode l */
  const INT4   m,          /**<< Current mode m */
  const REAL8  dt,         /**<< Sample time step (in seconds) */
  const REAL8  mass1,      /**<< First component mass (in Solar masses) */
  const REAL8  mass2,      /**<< Second component mass (in Solar masses) */
  const REAL8  spin1x,     /**<<The spin of the first object; only needed for spin waveforms */
  const REAL8  spin1y,     /**<<The spin of the first object; only needed for spin waveforms */
  const REAL8  spin1z,     /**<<The spin of the first object; only needed for spin waveforms */
  const REAL8  spin2x,     /**<<The spin of the second object; only needed for spin waveforms */
  const REAL8  spin2y,     /**<<The spin of the second object; only needed for spin waveforms */
  const REAL8  spin2z,     /**<<The spin of the second object; only needed for spin waveforms */
  REAL8Vector *timeVec,    /**<< Vector containing the time values */
  REAL8Vector *matchrange, /**<< Time values chosen as points for performing comb matching */
  Approximant  approximant /**<<The waveform approximant being used */
  );

REAL8 XLALKronecker( const INT4 i, const INT4 j );


/**
 * Calculates the dot product of two vectors
 */
static inline
REAL8
CalculateDotProduct( const REAL8 a[], const REAL8 b[] )
{
  return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}


/**
 * calculates the ith component of the cross product of two vectors,
 * where i is in the range 0-2.
 */
static inline
REAL8
CalculateCrossProduct( const INT4 i, const REAL8 a[], const REAL8 b[] )
{
  return a[(i+1)%3]*b[(i+2)%3] - a[(i+2)%3]*b[(i+1)%3];
}


/**
 * Normalizes the given vector
 */
static inline
int
NormalizeVector( REAL8 a[] )
{
  REAL8 norm = sqrt( a[0]*a[0] + a[1]*a[1] + a[2]*a[2] );

  a[0] /= norm;
  a[1] /= norm;
  a[2] /= norm;

  return XLAL_SUCCESS;
}





                 



#endif /* _LALSIMIMRSPINEOB_H */
