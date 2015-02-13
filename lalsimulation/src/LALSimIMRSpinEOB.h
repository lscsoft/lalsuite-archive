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


 INT4  XLALSimIMRSpinEOBGetSpinFactorizedWaveform(
                                COMPLEX16             * restrict hlm,
                                REAL8Vector           * restrict values,
                                const REAL8           v,
                                const REAL8           Hreal,
                                const INT4            l,
                                const INT4            m,
                                SpinEOBParams         * restrict params
                                );

 INT4  XLALSimIMRSpinEOBGetPrecSpinFactorizedWaveform( 
                 COMPLEX16         * restrict hlm,    /**< OUTPUT, hlm waveforms */
                 REAL8Vector       * restrict values, /**< dyanmical variables: (r,\phi,p_r,p_\phi) */
                 REAL8Vector       * restrict cartvalues, /**< dyanmical variables */
                 const REAL8         v,               /**< velocity */
                 const REAL8         Hreal,           /**< real Hamiltonian */
                 const INT4          l,               /**< l mode index */
                 const INT4          m,               /**< m mode index */
                 SpinEOBParams     * restrict params  /**< Spin EOB parameters */
                 );
 INT4   XLALSimIMRSpinEOBFluxGetSpinFactorizedWaveform(
                                COMPLEX16             * restrict hlm,
                                REAL8Vector           * restrict values,
                                const REAL8           v,
                                const REAL8           Hreal,
                                const INT4            l,
                                const INT4            m,
                                SpinEOBParams         * restrict params
                                );

 INT4  XLALSimIMRSpinEOBFluxGetPrecSpinFactorizedWaveform( 
                 COMPLEX16         * restrict hlm,    /**< OUTPUT, hlm waveforms */
                 REAL8Vector       * restrict values, /**< dyanmical variables: (r,\phi,p_r,p_\phi) */
                 REAL8Vector       * restrict cartvalues, /**< dyanmical variables */
                 const REAL8         v,               /**< velocity */
                 const REAL8         Hreal,           /**< real Hamiltonian */
                 const INT4          l,               /**< l mode index */
                 const INT4          m,               /**< m mode index */
                 SpinEOBParams     * restrict params  /**< Spin EOB parameters */
                 );



REAL8  XLALEffectiveHamiltonian( const REAL8 eta,          /**<< Symmetric mass ratio */
                                const REAL8 r,            /**<< Orbital separation */
                                const REAL8 pr,           /**<< Tortoise co-ordinate */
                                const REAL8 pp,           /**<< Momentum pphi */
                                EOBACoefficients *aCoeffs /**<< Pre-computed coefficients in A function */
                              );



REAL8   GSLSpinHamiltonianWrapper( double x, void *params );

int   XLALSpinHcapNumericalDerivative(
                          double                t,
                          const REAL8           values[],
                          REAL8                 dvalues[],
                          void                  *funcParams
                               ) ;

int  XLALSpinHcapNumericalDerivativeNoFlux(
                          double                t,
                          const REAL8           values[],
                          REAL8                 dvalues[],
                          void                  *funcParams
                               ) ;

REAL8  XLALSpinHcapNumDerivWRTParam(
                       const INT4 paramIdx,
                       const REAL8 values[],
                       SpinEOBParams *params
                       ) ;


 inline REAL8  XLALCalculateA5( REAL8 eta );

 inline REAL8   XLALCalculateA6( REAL8 eta );


/**
 * Calculates the a5 parameter in the A potential function in EOBNRv2
 */
 inline
REAL8   XLALCalculateA5( const REAL8 eta /**<< Symmetric mass ratio */
                     )
{
  return - 5.82827 - 143.486 * eta + 447.045 * eta * eta;
}

/**
 * Calculates the a6 parameter in the A potential function in EOBNRv2
 */
 inline
REAL8  XLALCalculateA6( UNUSED const REAL8 eta /**<< Symmetric mass ratio */)
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

 REAL8  
XLALAssociatedLegendreXIsZero( const int l,
                               const int m
                             );

 int  
XLALScalarSphHarmThetaPiBy2(COMPLEX16 *y,
                         INT4 l,
                         INT4  m,
                         REAL8 phi);

 int  
XLALAbsScalarSphHarmThetaPiBy2(COMPLEX16 *y,
                         INT4 l,
                         INT4  m);

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

int XLALSpinAlignedHcapDerivative(
            double                t,
            const REAL8           values[],
            REAL8                 dvalues[],
            void                  *funcParams
 );

 int  
CalculateThisMultipolePrefix(
               COMPLEX16 *prefix,
               const REAL8 m1,
               const REAL8 m2,
               const INT4 l,
               const INT4 m );

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

 int  
XLALScalarSphHarmThetaPiBy2(
                 COMPLEX16 *y, /**<< OUTPUT, Ylm(0,phi) */
                 INT4 l,       /**<< Mode l */
                 INT4  m,      /**<< Mode m */
                 REAL8 phi     /**<< Orbital phase (in radians) */
                 );


 int  
XLALAbsScalarSphHarmThetaPiBy2(
                 COMPLEX16 *y, /**<< OUTPUT, Ylm(0,phi) */
                 INT4 l,       /**<< Mode l */
                 INT4  m      /**<< Mode m */
                 );


 REAL8  
XLALAssociatedLegendreXIsZero( const int l,
                               const int m );


 int  
CalculateThisMultipolePrefix(
                 COMPLEX16 *prefix, /**<< OUTPUT, Prefix value */
                 const REAL8 m1,    /**<< mass 1 */
                 const REAL8 m2,    /**<< mass 2 */
                 const INT4 l,      /**<< Mode l */
                 const INT4 m       /**<< Mode m */
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




/**
 * This function calculates the DeltaR potential function in the spin EOB Hamiltonian
 */
REAL8 XLALSimIMRSpinEOBHamiltonianDeltaR(
        SpinEOBHCoeffs *coeffs, /**<< Pre-computed coefficients which appear in the function */
        const REAL8    r,       /**<< Current orbital radius (in units of total mass) */
        const REAL8    eta,     /**<< Symmetric mass ratio */
        const REAL8    a        /**<< Normalized deformed Kerr spin */
        );

REAL8 XLALSimIMRSpinEOBHamiltonian(
               const REAL8    eta,
               REAL8Vector    * restrict x,
               REAL8Vector    * restrict p,
               REAL8Vector    * restrict s1Vec,
               REAL8Vector    * restrict s2Vec,
               REAL8Vector    * restrict sigmaKerr,
               REAL8Vector    * restrict sigmaStar,
               int                       tortoise,
               SpinEOBHCoeffs *coeffs);

int XLALSimIMRCalculateSpinEOBHCoeffs(
        SpinEOBHCoeffs *coeffs,
        const REAL8    eta,
        const REAL8    a,
        const UINT4    SpinAlignedEOBversion
        );

REAL8 XLALSimIMRSpinEOBHamiltonianDeltaT( 
        SpinEOBHCoeffs *coeffs,
        const REAL8    r,
        const REAL8    eta,
        const REAL8    a
        );

REAL8 XLALSimIMRSpinAlignedEOBCalcOmega(
                      const REAL8          values[],
                      SpinEOBParams        *funcParams
                      );

REAL8 XLALSimIMRSpinAlignedEOBNonKeplerCoeff(
                      const REAL8           values[],
                      SpinEOBParams         *funcParams
                      );

double GSLSpinAlignedHamiltonianWrapper( double x, void *params );

REAL8  inner_product( const REAL8 values1[], 
                             const REAL8 values2[]
                             );

/* static */
REAL8* cross_product( const REAL8 values1[],
                              const REAL8 values2[],
                              REAL8 result[] 
                              );

REAL8 XLALSimIMRSpinEOBCalcOmega(
                      const REAL8          values[],
                      SpinEOBParams        *funcParams
                      );

REAL8 XLALSimIMRSpinEOBNonKeplerCoeff(
                      const REAL8           values[],
                      SpinEOBParams         *funcParams
                      );

int XLALSpinHcapRvecDerivative(
                 double      t,         /**<<  */
                 const  REAL8      values[],  /**<< Dynamical variables */
                 REAL8             dvalues[], /**<< Time derivatives of variables (returned) */
                 void             *funcParams /**<< EOB parameters */
                               ) ;
                              
double GSLSpinHamiltonianWrapperFordHdpphi( double x, void *params );

double GSLSpinHamiltonianWrapperForRvecDerivs( double x, void *params );





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



/* */
REAL8 XLALCalculateEOBD( 
        REAL8   r, /**<< Orbital separation (in units of total mass M) */
        REAL8 eta  /**<< Symmetric mass ratio */       
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

int XLALSimIMRSpinEOBCalculateSigmaKerr(
                                   REAL8Vector *sigmaKerr,
                                   REAL8        mass1,
                                   REAL8        mass2,
                                   REAL8Vector *s1,
                                   REAL8Vector *s2 );

int XLALSimIMRSpinEOBCalculateSigmaStar(
                                   REAL8Vector *sigmaStar,
                                   REAL8        mass1,
                                   REAL8        mass2,
                                   REAL8Vector *s1,
                                   REAL8Vector *s2 );


int XLALSimIMREOBCalcSpinFacWaveformCoefficients(
          FacWaveformCoeffs * const coeffs, /**< OUTPUT, pre-computed waveform coefficients */
          const REAL8               m1,     /**< mass 1 */
          const REAL8               m2,     /**< mass 2 */
          const REAL8               eta,    /**< symmetric mass ratio */
          const REAL8               a,      /**< Kerr spin parameter for test-particle terms */
          const REAL8               chiS,   /**< (chi1+chi2)/2 */
          const REAL8               chiA,   /**< (chi1-chi2)/2 */
          const UINT4               SpinAlignedEOBversion  /**< 1 for SEOBNRv1; 2 for SEOBNRv2 */
          );




/***** Inline functions  **********/

/**
 * The time difference between the orbital peak and the peak amplitude
 * of the mode in question (currently only 2,2 implemented ).
 * Eq. 33 of Taracchini et al. PRD 86, 024011 (2012).
 */
 static 
inline REAL8 XLALSimIMREOBGetNRSpinPeakDeltaT( 
                 INT4 l,           /**<< Mode l */
                 INT4 m,           /**<< Mode m */
                 REAL8  eta, /**<< Symmetric mass ratio */
                 REAL8 a           /**<< Dimensionless spin */
                 )
{

  switch ( l )
  {
    case 2:
      switch ( m )
      {
        case 2:
          /* DeltaT22 defined here is a minus sign different from Eq. (33) of Taracchini et al. */
          if ( a <= 0.0 )
          {
            return 2.5;
          }
          else
          {
            return (2.5 + 1.77*a*a*a*a/(0.43655*0.43655*0.43655*0.43655)/(1.0-2.0*eta)/(1.0-2.0*eta)/(1.0-2.0*eta)/(1.0-2.0*eta));
          }
          break;
        default:
          XLAL_ERROR_REAL8( XLAL_EINVAL );
      }
      break;
    default:
      XLAL_ERROR_REAL8( XLAL_EINVAL );
  }

  /* We should never get here, but I expect a compiler whinge without it... */
  XLALPrintError( "XLAL Error %s - We should never get here!!\n", __func__ );
  XLAL_ERROR_REAL8( XLAL_EINVAL );
}

/**
 * Peak amplitude predicted by fitting NR results (currently only 2,2 available).
 * Tables IV and V and Eq. 42 of Taracchini et al. PRD 86, 024011 (2012).
 */
 static inline REAL8 GetNRSpinPeakAmplitude( INT4 UNUSED l, INT4 UNUSED m, REAL8  eta, REAL8 UNUSED a )
{
  /* Fit for HOMs missing */
  return 1.3547468629743946*eta + 0.9187885481024214*eta*eta;
}

/**
 * Peak amplitude curvature predicted by fitting NR results (currently only 2,2 available).
 * Tables IV and V and Eq. 42 of Taracchini et al. PRD 86, 024011 (2012).
 */
 static inline REAL8 GetNRSpinPeakADDot( INT4 UNUSED  l, INT4  UNUSED m, REAL8  eta, REAL8  a )
{
  /* Fit for HOMs missing */
  return eta*(-0.0024971911410897156 + (-0.006128515435641139 + 0.01732656*a/(2.0-4.0*eta))*eta);
}

/**
 * Peak frequency predicted by fitting NR results (currently only 2,2 available).
 * Tables IV and V and Eq. 42 of Taracchini et al. PRD 86, 024011 (2012).
 */
 static inline REAL8 GetNRSpinPeakOmega( INT4 UNUSED  l, INT4 UNUSED  m, REAL8  eta, REAL8 a )
{
  /* Fit for HOMs missing */
  return 0.27581190323955274 + 0.19347381066059993*eta
       - 0.08898338208573725*log(1.0 - a/(1.0-2.0*eta))
       + eta*eta*(1.78832*(0.2690779744133912 + a/(2.0-4.0*eta))*(1.2056469070395925
       + a/(2.0-4.0*eta)) + 1.423734113371796*log(1.0 - a/(1.0-2.0*eta)));
}

/**
 * Peak frequency slope predicted by fitting NR results (currently only 2,2 available).
 * Tables IV and V and Eq. 42 of Taracchini et al. PRD 86, 024011 (2012).
 */
 static inline REAL8 GetNRSpinPeakOmegaDot( INT4 UNUSED  l, INT4 UNUSED  m, REAL8  eta, REAL8  a )
{
  /* Fit for HOMs missing */
  return 0.006075014646800278 + 0.012040017219351778*eta
       + (0.0007353536801336875 + 0.0015592659912461832*a/(1.0-2.0*eta))*log(1.0-a/(1.0-2.0*eta))
       + eta*eta*(0.03575969677378844 + (-0.011765658882139 - 0.02494825585993893*a/(1.0-2.0*eta))
       * log(1.0 - a/(1.0-2.0*eta)));
}


/**
 * The time difference between the orbital peak and the peak amplitude
 * of the mode in question (currently only 2,2 implemented ).
 */
 static inline REAL8 XLALSimIMREOBGetNRSpinPeakDeltaTv2(
                 UNUSED INT4  l,    /**<< Mode l */
                 UNUSED INT4  m,    /**<< Mode m */
                 REAL8  m1,  /**<< mass 1 */
                 REAL8  m2,  /**<< mass 22 */
                 REAL8 chi1,       /**<< Dimensionless spin1 */
                 REAL8 chi2        /**<< Dimensionless spin2 */
                 )
{
  REAL8 chi, chichi;
  REAL8 eta = m1*m2 / ((m1+m2)*(m1+m2));
  chi    = (chi1+chi2)/2. + (chi1-chi2)/2.*((m1-m2)/(m1+m2))/(1-2.*eta);
    
  chichi = chi*chi;
  if ( chi > 0.8 )
  {
    return (0.75*eta*chi + sqrt(1.-4.*eta)) * (57.1755-48.0564*chi);
  }
  else if ( chi > 0.0 )
  {
    return (0.75*eta*chi + sqrt(1.-4.*eta)) * (2.5+10.*chichi+24.*chichi*chichi);
  }
  else
  {
    return 2.5 + (1.+2.5*chi) * (-2.5+2.5*sqrt(1.-4.*eta));
  }
}


/**
 * Peak frequency predicted by fitting NR results (currently only 2,2 available).
 * Take from unpublished SEOBNRv2 results.
 */
 static inline REAL8 GetNRSpinPeakOmegav2( INT4 UNUSED  l, INT4 UNUSED  m, REAL8  eta, REAL8 a )
{
  REAL8 chi = a/(1.0 - 2.0*eta);
  REAL8 eta2= eta*eta;
  if ( eta > 50./51./51. )
  {
    return 0.43747541927878864 + (-0.10933208665273314 - 0.007325831113333813*chi)*log(
            4.500844771420863 - 9.681916048928946*eta +
            chi*(-4.254886879579986 + 11.513558950322647*eta));
  }
  else
  {
    return 1.5609526077704716 - 122.25721149839733 * eta +
           3586.2192688666914 * eta2 - 13869.506144441548 * eta*eta2 +
           (eta - 0.25) * (1651.5823693445805 * (-0.01922337588094282 + eta) *
           (-0.01922337536857659 + eta) + 66.87492814925524 * chi *
           (0.0003695381704106058 - 0.03844675124951941 * eta + eta2)) *
           log(5600.67382718678 - 5555.824895398546 * chi) +
           (-1412.8186461833657 + 67.66455403259023 * chi) * (eta - 0.001) *
           (0.0003695381704106056 - 0.038446751249519406 * eta + eta2) *
           log(0.5680439481719505 - 0.36813967358200156 * chi) +
           0.012328326527732041 * log(4.500844771420863 - 9.681916048928946 * eta + 
           chi * (-4.254886879579986 + 11.513558950322647 * eta)) +
           0.0008260634258180991 * chi * log(4.500844771420863 -9.681916048928946 * eta +
           chi * (-4.254886879579986 + 11.513558950322647 * eta)) -
           12.6575493872956 * eta * log(4.500844771420863 -
           9.681916048928946 * eta + chi * (-4.254886879579986 + 11.513558950322647 * eta)) -
           0.8481231078533651 * chi * eta * log(4.500844771420863 - 
           9.681916048928946 * eta + chi * (-4.254886879579986 + 11.513558950322647 * eta)) +
           329.2228595635586 * eta2 * log(4.500844771420863 - 9.681916048928946 * eta +
           chi * (-4.254886879579986 + 11.513558950322647 * eta)) +
           22.05968203526603 * chi * eta2 * log(4.500844771420863 -9.681916048928946 * eta +
           chi * (-4.254886879579986 + 11.513558950322647 * eta));
  }
}

/**
 * Peak frequency slope predicted by fitting NR results (currently only 2,2 available).
 * Take from unpublished SEOBNRv2 results.
 */
 static inline REAL8 GetNRSpinPeakOmegaDotv2( INT4 UNUSED  l, INT4 UNUSED  m, REAL8  eta, REAL8  a )
{
  REAL8 chi = a/(1.0-2.0*eta);
  REAL8 eta2= eta*eta;
  /* Fit for HOMs missing */
  if (chi < 0.8 )
  {
    return -0.07086074186161867 * chi * (-0.26367236731979804 + eta) *
           (-0.0010019969893089581 + eta) + 0.2893863668183948 *
           (-0.16845695144529893 + eta) * (0.23032241797163952 + eta) +
           (0.004086861548547749 - 0.06538978477676398 * eta2 +
           chi * (0.0006334026884930817 - 0.010134443015889307 * eta2)) *
           log(68.47466578101876 - 58.30148755701496 * chi);
  }
  else
  {
    return -0.10069512275335238 * (-0.46107388514323044 + eta) *
           (0.2832795481380979 + eta) + 0.2614619716504706 * chi *
           (-0.24838163750494138 + eta) * (0.320112993649413 + eta) +
           chi * chi * (0.010000160002560042 - 0.16000256004096067 * eta2);
  }
}


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
