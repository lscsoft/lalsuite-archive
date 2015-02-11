#include <lal/LALSimInspiral.h>
#include <lal/LALSimIMR.h>

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

#include <math.h>

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

UNUSED  static
REAL8 XLALCalculateEOBA( 
        const REAL8 r,                     /**<< Orbital separation (in units of total mass M) */
        EOBACoefficients * restrict coeffs /**<< Pre-computed coefficients for the A function */
 );


static INT4 UNUSED XLALSimIMRSpinEOBGetSpinFactorizedWaveform(
                                COMPLEX16             * restrict hlm,
                                REAL8Vector           * restrict values,
                                const REAL8           v,
                                const REAL8           Hreal,
                                const INT4            l,
                                const INT4            m,
                                SpinEOBParams         * restrict params
                                );

static INT4 UNUSED XLALSimIMRSpinEOBGetPrecSpinFactorizedWaveform( 
                 COMPLEX16         * restrict hlm,    /**< OUTPUT, hlm waveforms */
                 REAL8Vector       * restrict values, /**< dyanmical variables: (r,\phi,p_r,p_\phi) */
                 REAL8Vector       * restrict cartvalues, /**< dyanmical variables */
                 const REAL8         v,               /**< velocity */
                 const REAL8         Hreal,           /**< real Hamiltonian */
                 const INT4          l,               /**< l mode index */
                 const INT4          m,               /**< m mode index */
                 SpinEOBParams     * restrict params  /**< Spin EOB parameters */
                 );
static INT4 UNUSED  XLALSimIMRSpinEOBFluxGetSpinFactorizedWaveform(
                                COMPLEX16             * restrict hlm,
                                REAL8Vector           * restrict values,
                                const REAL8           v,
                                const REAL8           Hreal,
                                const INT4            l,
                                const INT4            m,
                                SpinEOBParams         * restrict params
                                );

static INT4 UNUSED XLALSimIMRSpinEOBFluxGetPrecSpinFactorizedWaveform( 
                 COMPLEX16         * restrict hlm,    /**< OUTPUT, hlm waveforms */
                 REAL8Vector       * restrict values, /**< dyanmical variables: (r,\phi,p_r,p_\phi) */
                 REAL8Vector       * restrict cartvalues, /**< dyanmical variables */
                 const REAL8         v,               /**< velocity */
                 const REAL8         Hreal,           /**< real Hamiltonian */
                 const INT4          l,               /**< l mode index */
                 const INT4          m,               /**< m mode index */
                 SpinEOBParams     * restrict params  /**< Spin EOB parameters */
                 );


static
REAL8 UNUSED XLALEffectiveHamiltonian( const REAL8 eta,          /**<< Symmetric mass ratio */
                                const REAL8 r,            /**<< Orbital separation */
                                const REAL8 pr,           /**<< Tortoise co-ordinate */
                                const REAL8 pp,           /**<< Momentum pphi */
                                EOBACoefficients *aCoeffs /**<< Pre-computed coefficients in A function */
                              );



REAL8 UNUSED  GSLSpinHamiltonianWrapper( double x, void *params );

int UNUSED  XLALSpinHcapNumericalDerivative(
                          double                t,
                          const REAL8           values[],
                          REAL8                 dvalues[],
                          void                  *funcParams
                               ) UNUSED;

int UNUSED XLALSpinHcapNumericalDerivativeNoFlux(
                          double                t,
                          const REAL8           values[],
                          REAL8                 dvalues[],
                          void                  *funcParams
                               ) UNUSED;

REAL8 UNUSED XLALSpinHcapNumDerivWRTParam(
                       const INT4 paramIdx,
                       const REAL8 values[],
                       SpinEOBParams *params
                       ) UNUSED;


static inline REAL8 UNUSED XLALCalculateA5( REAL8 eta );

static inline REAL8 UNUSED  XLALCalculateA6( REAL8 eta );

static
REAL8 UNUSED  XLALCalculateEOBD( REAL8    r,
                         REAL8	eta) UNUSED;


/**
 * Calculates the a5 parameter in the A potential function in EOBNRv2
 */
static inline
REAL8 UNUSED  XLALCalculateA5( const REAL8 eta /**<< Symmetric mass ratio */
                     )
{
  return - 5.82827 - 143.486 * eta + 447.045 * eta * eta;
}

/**
 * Calculates the a6 parameter in the A potential function in EOBNRv2
 */
static inline
REAL8 UNUSED XLALCalculateA6( const REAL8 UNUSED eta /**<< Symmetric mass ratio */
                     )
{
  return 184.0;
}


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

static REAL8 UNUSED XLALSimInspiralEt0(REAL8 v, expnCoeffsdEnergyFlux *ak);

static REAL8 UNUSED XLALSimInspiralEt2(REAL8 v, expnCoeffsdEnergyFlux *ak);

static REAL8 UNUSED XLALSimInspiralEt4(REAL8 v, expnCoeffsdEnergyFlux *ak);

static REAL8 UNUSED XLALSimInspiralEt6(REAL8 v, expnCoeffsdEnergyFlux *ak);

static REAL8 UNUSED XLALSimInspiraldEt0(REAL8 v, expnCoeffsdEnergyFlux *ak);

static REAL8 UNUSED XLALSimInspiraldEt2(REAL8 v, expnCoeffsdEnergyFlux *ak);

static REAL8 UNUSED XLALSimInspiraldEt4(REAL8 v, expnCoeffsdEnergyFlux *ak);

static REAL8 UNUSED XLALSimInspiraldEt6(REAL8 v, expnCoeffsdEnergyFlux *ak);

static REAL8 UNUSED XLALSimInspiralFt0(REAL8 v, expnCoeffsdEnergyFlux *ak);

static REAL8 UNUSED XLALSimInspiralFt2(REAL8 v, expnCoeffsdEnergyFlux *ak);

static REAL8 UNUSED XLALSimInspiralFt3(REAL8 v, expnCoeffsdEnergyFlux *ak);

static REAL8 UNUSED XLALSimInspiralFt4(REAL8 v, expnCoeffsdEnergyFlux *ak);

static REAL8 UNUSED XLALSimInspiralFt5(REAL8 v, expnCoeffsdEnergyFlux *ak);

static REAL8 UNUSED XLALSimInspiralFt6(REAL8 v, expnCoeffsdEnergyFlux *ak);

static REAL8 UNUSED XLALSimInspiralFt7(REAL8 v, expnCoeffsdEnergyFlux *ak);


static REAL8 UNUSED XLALSimInspiralep2(REAL8 v, expnCoeffsdEnergyFlux *ak);

static REAL8 UNUSED XLALSimInspiralep4(REAL8 v, expnCoeffsdEnergyFlux *ak);

static REAL8 UNUSED XLALSimInspiralep6(REAL8 v, expnCoeffsdEnergyFlux *ak);

static REAL8 UNUSED XLALSimInspiraldEp2(REAL8 v, expnCoeffsdEnergyFlux *ak);

static REAL8 UNUSED XLALSimInspiraldEp4(REAL8 v, expnCoeffsdEnergyFlux *ak);

static REAL8 UNUSED XLALSimInspiraldEp6(REAL8 v, expnCoeffsdEnergyFlux *ak);

static REAL8 UNUSED XLALSimInspiralFp3(REAL8 v, expnCoeffsdEnergyFlux *ak);

static REAL8 UNUSED XLALSimInspiralFp4(REAL8 v, expnCoeffsdEnergyFlux *ak);

static REAL8 UNUSED XLALSimInspiralFp5(REAL8 v, expnCoeffsdEnergyFlux *ak);

static REAL8 UNUSED XLALSimInspiralFp6(REAL8 v, expnCoeffsdEnergyFlux *ak);

static REAL8 UNUSED XLALSimInspiralFp7(REAL8 v, expnCoeffsdEnergyFlux *ak);

static REAL8 UNUSED XLALSimInspiralFp8PP(REAL8 v, expnCoeffsdEnergyFlux *ak);

static REAL8 UNUSED 
XLALAssociatedLegendreXIsZero( const int l,
                               const int m
                             );

static int UNUSED 
XLALScalarSphHarmThetaPiBy2(COMPLEX16 *y,
                         INT4 l,
                         INT4  m,
                         REAL8 phi);

static int UNUSED 
XLALAbsScalarSphHarmThetaPiBy2(COMPLEX16 *y,
                         INT4 l,
                         INT4  m);

static int UNUSED 
CalculateThisMultipolePrefix(
               COMPLEX16 *prefix,
               const REAL8 m1,
               const REAL8 m2,
               const INT4 l,
               const INT4 m );

static int UNUSED XLALSimIMREOBComputeNewtonMultipolePrefixes(
                NewtonMultipolePrefixes *prefix, /**<< OUTPUT Structure containing the coeffs */
                const REAL8             m1,      /**<< Mass of first component */
                const REAL8             m2       /**<< Nass of second component */
                );


UNUSED static int 
XLALSimIMREOBCalculateNewtonianMultipole(
                 COMPLEX16 *multipole, /**<< OUTPUT, Newtonian multipole */
                 REAL8 x,              /**<< Dimensionless parameter \f$\equiv v^2\f$ */
                 UNUSED REAL8 r,       /**<< Orbital separation (units of total mass M) */
                 REAL8 phi,            /**<< Orbital phase (in radians) */
                 UINT4  l,             /**<< Mode l */
                 INT4  m,              /**<< Mode m */
                 EOBParams *params     /**<< Pre-computed coefficients, parameters, etc. */
                 );


UNUSED static int
XLALSimIMRSpinEOBCalculateNewtonianMultipole(
                 COMPLEX16 *multipole, /**<< OUTPUT, Newtonian multipole */
                 REAL8 x,              /**<< Dimensionless parameter \f$\equiv v^2\f$ */
                 UNUSED REAL8 r,       /**<< Orbital separation (units of total mass M */
                 REAL8 phi,            /**<< Orbital phase (in radians) */
                 UINT4  l,             /**<< Mode l */
                 INT4  m,              /**<< Mode m */
                 EOBParams *params     /**<< Pre-computed coefficients, parameters, etc. */
                 );


UNUSED static int
XLALSimIMRSpinEOBFluxCalculateNewtonianMultipole(
                 COMPLEX16 *multipole, /**<< OUTPUT, Newtonian multipole */
                 REAL8 x,              /**<< Dimensionless parameter \f$\equiv v^2\f$ */
                 UNUSED REAL8 r,       /**<< Orbital separation (units of total mass M */
                 UNUSED REAL8 phi,     /**<< Orbital phase (in radians) */
                 UINT4  l,             /**<< Mode l */
                 INT4  m,              /**<< Mode m */
                 EOBParams *params     /**<< Pre-computed coefficients, parameters, etc. */
                 );

static int UNUSED 
XLALScalarSphHarmThetaPiBy2(
                 COMPLEX16 *y, /**<< OUTPUT, Ylm(0,phi) */
                 INT4 l,       /**<< Mode l */
                 INT4  m,      /**<< Mode m */
                 REAL8 phi     /**<< Orbital phase (in radians) */
                 );


static int UNUSED 
XLALAbsScalarSphHarmThetaPiBy2(
                 COMPLEX16 *y, /**<< OUTPUT, Ylm(0,phi) */
                 INT4 l,       /**<< Mode l */
                 INT4  m      /**<< Mode m */
                 );


static REAL8 UNUSED 
XLALAssociatedLegendreXIsZero( const int l,
                               const int m );


static int UNUSED 
CalculateThisMultipolePrefix(
                 COMPLEX16 *prefix, /**<< OUTPUT, Prefix value */
                 const REAL8 m1,    /**<< mass 1 */
                 const REAL8 m2,    /**<< mass 2 */
                 const INT4 l,      /**<< Mode l */
                 const INT4 m       /**<< Mode m */
                 );


static REAL8  UNUSED XLALSimIMREOBGetNRPeakDeltaT( 
                         INT4 l,    /**<< Mode l */ 
                         INT4 m,    /**<< Mode m */
                         REAL8 eta  /**<< Symmetric mass ratio */
                         );

UNUSED static int XLALSimIMREOBGetCalibratedNQCCoeffs( 
                                EOBNonQCCoeffs *coeffs, /**<< OUTPUT, Structure for NQC coeffs */
                                INT4            l,      /**<< Mode l */
                                INT4            m,      /**<< Mode m */
                                REAL8           eta     /**<< Symmetric mass ratio */
                                );


UNUSED static int  XLALSimIMREOBNonQCCorrection(
                      COMPLEX16      * restrict nqc,    /**<< OUTPUT, The NQC correction */
                      REAL8Vector    * restrict values, /**<< Dynamics r, phi, pr, pphi */
                      const REAL8               omega,  /**<< Angular frequency */
                      EOBNonQCCoeffs * restrict coeffs  /**<< NQC coefficients */
                     );


UNUSED static int  XLALSimIMRSpinEOBNonQCCorrection(
                      COMPLEX16      * restrict nqc,    /**<< OUTPUT, The NQC correction */
                      REAL8Vector    * restrict values, /**<< Dynamics r, phi, pr, pphi */
                      const REAL8               omega,  /**<< Angular frequency */
                      EOBNonQCCoeffs * restrict coeffs  /**<< NQC coefficients */
                     );


UNUSED static int XLALSimIMREOBCalculateNQCCoefficients(
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



UNUSED static int XLALSimIMRGetEOBCalibratedSpinNQC( EOBNonQCCoeffs * restrict coeffs,                                   INT4 UNUSED l,                                   INT4 UNUSED m,                                   REAL8 eta,                                   REAL8 a );


UNUSED static int XLALSimIMRGetEOBCalibratedSpinNQCv2( EOBNonQCCoeffs * restrict coeffs,                                   INT4 UNUSED l,                                   INT4 UNUSED m,                                   REAL8 eta,                                   REAL8 a );


UNUSED static int XLALSimIMRGetEOBCalibratedSpinNQCv2chiAmax( EOBNonQCCoeffs * restrict coeffs,                                                             INT4 UNUSED l,                                                             INT4 UNUSED m,                                                             REAL8 eta,                                                             REAL8 a );


UNUSED static int XLALSimIMRGetEOBCalibratedSpinNQCv2chiAmed( EOBNonQCCoeffs * restrict coeffs,                                                             INT4 UNUSED l,                                                             INT4 UNUSED m,                                                             REAL8 eta,                                                             REAL8 a );



UNUSED static int XLALSimIMRGetEOBCalibratedSpinNQCv2chiAmin( EOBNonQCCoeffs * restrict coeffs,                                                             INT4 UNUSED l,                                                             INT4 UNUSED m,                                                             REAL8 eta,                                                             REAL8 a );


UNUSED static int XLALSimIMRGetEOBCalibratedSpinNQC3D(
                        EOBNonQCCoeffs * restrict coeffs,
                        INT4 UNUSED l,
                        INT4 UNUSED m,
                        REAL8 m1,
                        REAL8 m2,
                        REAL8 a,
                        REAL8 chiAin );


UNUSED static int XLALSimIMRSpinEOBCalculateNQCCoefficients(
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

UNUSED static
int XLALCalculateEOBACoefficients(
          EOBACoefficients * const coeffs, /**<< A coefficients (populated in function) */
          const REAL8              eta     /**<< Symmetric mass ratio */
          );


UNUSED static
REAL8 XLALCalculateEOBdAdr( 
        const REAL8 r,                     /**<< Orbital separation (in units of total mass M) */
        EOBACoefficients * restrict coeffs /**<< Pre-computed coefficients for the A function */
     );



static REAL8 UNUSED  XLALCalculateEOBD( 
        REAL8   r, /**<< Orbital separation (in units of total mass M) */
        REAL8 eta  /**<< Symmetric mass ratio */       
        );


UNUSED static int XLALSimIMREOBCalcFacWaveformCoefficients(
          FacWaveformCoeffs * const coeffs, /**<< Structure containing coefficients (populated in function) */
          const REAL8               eta     /**<< Symmetric mass ratio */
          );

UNUSED static int XLALSimIMREOBModifyFacWaveformCoefficients( 
                FacWaveformCoeffs * const coeffs, /**<< Structure containing coefficients */
                const REAL8 eta                   /**<< Symmetric mass ratio */
        );


static REAL8 UNUSED 
nonKeplerianCoefficient(
        REAL8Vector * restrict values, /**<< Dynamics r, phi, pr, pphi */
        const REAL8       eta,         /**<< Symmetric mass ratio */
        EOBACoefficients *coeffs       /**<< Pre-computed A coefficients */
        );


UNUSED static int  XLALSimIMREOBGetFactorizedWaveform( 
            COMPLEX16   * restrict hlm,    /**<< The value of hlm (populated by the function) */
            REAL8Vector * restrict values, /**<< Vector containing dynamics r, phi, pr, pphi for a given point */
            const REAL8 v,                 /**<< Velocity (in geometric units) */
            const INT4  l,                 /**<< Mode l */
            const INT4  m,                 /**<< Mode m */
            EOBParams   * restrict params  /**<< Structure containing pre-computed coefficients, etc. */
            );

static REAL8 UNUSED  XLALSimIMREOBFactorizedFlux(
        REAL8Vector  *values, /**<< Dynamics r, phi, pr, pphi */
        const REAL8  omega,   /**<< Angular frequency omega */
        EOBParams    *ak,     /**<< Structure containing pre-computed parameters */
        const INT4   lMax     /**<< Maximum l to include when calculating flux (between 2 and 8) */
        );

static INT4 UNUSED  XLALSimIMREOBHybridRingdownWave(
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


static INT4 UNUSED  XLALGenerateHybridWaveDerivatives (
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

static INT4 UNUSED  XLALSimIMREOBHybridAttachRingdown(
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










                 





































/*
int XLALSimIMREOBCalcSpinFacWaveformCoefficients(
          FacWaveformCoeffs * const coeffs,
          const REAL8               m1,
          const REAL8               m2,
          const REAL8               eta,
          const REAL8               a,
          const REAL8               chiS,
          const REAL8               chiA,
          const UINT4               SpinAlignedEOBversion
          );*/
#endif /* _LALSIMIMRSPINEOB_H */
