/*
* Copyright (C) 2016 Lionel London
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


/* .................... */
/* HEADER SECTION       */
/* .................... */
#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <complex.h>
#include <stdlib.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_spline.h>
#include <lal/Units.h>
#include <lal/SeqFactories.h>
#include <lal/Date.h>
#include <lal/StringInput.h>
#include <lal/Sequence.h>
#include <lal/FileIO.h>

#include <lal/XLALError.h>

#include <lal/LALConstants.h>
#include <lal/LALStdio.h>
#include <lal/LALSimSphHarmSeries.h>
#include <lal/LALStdlib.h>
#include <lal/LALSimInspiral.h>
#include <lal/Window.h>
#include <lal/TimeSeries.h>
#include <lal/FrequencySeries.h>

#include "LALSimRingdownMMRDNS.h"
#include <lal/LALSimInspiralTestGRParams.h>


/*
* Domain mapping for dimnesionless BH spin
*/
REAL8 XLALKAPPA( double jf, int l, int m ){
  /* */
  /* if ( jf > 1.0 ) XLAL_ERROR(XLAL_EDOM, "Spin (dimensionless Kerr parameter) must not be greater than 1.0\n"); */
  /**/
  double alpha = log( 2.0 - jf ) * 0.91023922662; // 1/log(3)=1.09861228867
  double beta  = 1.0 / ( 2.0 + l-abs(m) );
  return pow( alpha , beta );
}

/*
* Dimensionless QNM Frequencies: Note that name encodes date of writing
*/
COMPLEX16 XLALcomplexOmega( double kappa,  /* Domain mapping for  remnant BH's spin (Dimensionless) */
                          int l,        /* Polar eigenvalue */
                          int input_m,  /* Azimuthal eigenvalue*/
                          int n ) {     /* Overtone Number*/

  /* Predefine powers to increase efficiency*/
  double kappa2 = kappa  * kappa;
  double kappa3 = kappa2 * kappa;
  double kappa4 = kappa3 * kappa;

  /* NOTE that |m| will be used to determine the fit to use, and if input_m < 0, then a conjugate will be taken*/
  int m = abs(input_m);

  /**/
  complex double j = _Complex_I;

  /* Initialize the answer*/
  double complex ans = 0.0;

  /* Use If-Else ladder to determine which mode function to evaluate*/
  if ( 2==l && 2==m && 0==n  ){

    /* Fit for (l,m,n) == (2,2,0). This is a zero-damped mode in the extremal Kerr limit.*/
    ans = 1.0 + kappa * (  1.557847  *cexp(2.903124*j) +
                           1.95097051*cexp(5.920970*j)*kappa +
                           2.09971716*cexp(2.760585*j)*kappa2 +
                           1.41094660*cexp(5.914340*j)*kappa3 +
                           0.41063923*cexp(2.795235*j)*kappa4  );

    } else if ( 2==l && 2==m && 1==n ) {

      /* Fit for (l,m,n) == (2,2,1). This is a zero-damped mode in the extremal Kerr limit.*/
      ans = 1.0 + kappa * (  1.870939*cexp(2.511247*j) +
                             2.71924916*cexp(5.424999*j)*kappa +
                             3.05648030*cexp(2.285698*j)*kappa2 +
                             2.05309677*cexp(5.486202*j)*kappa3 +
                             0.59549897*cexp(2.422525*j)*kappa4  );

    } else if ( 3==l && 2==m && 0==n ) {

      /* Define extra powers as needed*/
      double kappa5 = kappa4 * kappa;
      double kappa6 = kappa5 * kappa;

      /* Fit for (l,m,n) == (3,2,0). This is NOT a zero-damped mode in the extremal Kerr limit.*/
      ans = 1.022464*cexp(0.004870*j) +
            0.24731213*cexp(0.665292*j)*kappa +
            1.70468239*cexp(3.138283*j)*kappa2 +
            0.94604882*cexp(0.163247*j)*kappa3 +
            1.53189884*cexp(5.703573*j)*kappa4 +
            2.28052668*cexp(2.685231*j)*kappa5 +
            0.92150314*cexp(5.841704*j)*kappa6;

    } else if ( 4==l && 4==m && 0==n ) {

      /* Fit for (l,m,n) == (4,4,0). This is a zero-damped mode in the extremal Kerr limit.*/
      ans = 2.0 + kappa * (  2.658908*cexp(3.002787*j) +
                             2.97825567*cexp(6.050955*j)*kappa +
                             3.21842350*cexp(2.877514*j)*kappa2 +
                             2.12764967*cexp(5.989669*j)*kappa3 +
                             0.60338186*cexp(2.830031*j)*kappa4  );

    } else if ( 2==l && 1==m && 0==n ) {

      /* Define extra powers as needed*/
      double kappa5 = kappa4 * kappa;
      double kappa6 = kappa5 * kappa;

      /* Fit for (l,m,n) == (2,1,0). This is NOT a zero-damped mode in the extremal Kerr limit.*/
      ans = 0.589113*cexp(0.043525*j) +
            0.18896353*cexp(2.289868*j)*kappa +
            1.15012965*cexp(5.810057*j)*kappa2 +
            6.04585476*cexp(2.741967*j)*kappa3 +
            11.12627777*cexp(5.844130*j)*kappa4 +
            9.34711461*cexp(2.669372*j)*kappa5 +
            3.03838318*cexp(5.791518*j)*kappa6;

    } else if ( 3==l && 3==m && 0==n ) {

      /* Fit for (l,m,n) == (3,3,0). This is a zero-damped mode in the extremal Kerr limit.*/
      ans = 1.5 + kappa*(  2.095657*cexp(2.964973*j) +
                           2.46964352*cexp(5.996734*j)*kappa +
                           2.66552551*cexp(2.817591*j)*kappa2 +
                           1.75836443*cexp(5.932693*j)*kappa3 +
                           0.49905688*cexp(2.781658*j)*kappa4  );

    } else if ( 3==l && 3==m && 1==n ) {

      /* Fit for (l,m,n) == (3,3,1). This is a zero-damped mode in the extremal Kerr limit.*/
      ans = 1.5 + kappa*(  2.339070*cexp(2.649692*j) +
                           3.13988786*cexp(5.552467*j)*kappa +
                           3.59156756*cexp(2.347192*j)*kappa2 +
                           2.44895997*cexp(5.443504*j)*kappa3 +
                           0.70040804*cexp(2.283046*j)*kappa4  );

    } else if ( 4==l && 3==m && 0==n ) {

      /* Fit for (l,m,n) == (4,3,0). This is a zero-damped mode in the extremal Kerr limit.*/
      ans = 1.5 + kappa*(  0.205046*cexp(0.595328*j) +
                           3.10333396*cexp(3.016200*j)*kappa +
                           4.23612166*cexp(6.038842*j)*kappa2 +
                           3.02890198*cexp(2.826239*j)*kappa3 +
                           0.90843949*cexp(5.915164*j)*kappa4  );

    } else if ( 5==l && 5==m && 0==n ) {

      /* Fit for (l,m,n) == (5,5,0). This is a zero-damped mode in the extremal Kerr limit. */
      ans = 2.5 + kappa*(  3.240455*cexp(3.027869*j) +
                           3.49056455*cexp(6.088814*j)*kappa +
                           3.74704093*cexp(2.921153*j)*kappa2 +
                           2.47252790*cexp(6.036510*j)*kappa3 +
                           0.69936568*cexp(2.876564*j)*kappa4  );

    } else {

      /**/
      ans  = 0.0;

  } /* END of IF-ELSE Train for QNM cases */

  /* If m<0, then take the *Negative* conjugate */
  if ( input_m < 0 ) {
    /**/
    ans = -conj( ans );
  }

  return ans;

} /* END of complexOmega */


/*
* -------------------------------------------------------------------------------- *
* Low level models: QNM Frequencies, Separation Constants and Spheroidal Harmonics
* -------------------------------------------------------------------------------- *
*/

/*
* QNM Separation Constants: Note that name encodes date of writing
*/
COMPLEX16 XLALseparationConstant( double kappa,  /* Domain mapping for remnant BH's spin (Dimensionless) */
                          int l,        /* Polar eigenvalue */
                          int input_m,  /* Azimuthal eigenvalue */
                          int n ) {     /* Overtone Number */

  /* Predefine powers to increase efficiency */
  double kappa2 = kappa  * kappa;
  double kappa3 = kappa2 * kappa;
  double kappa4 = kappa3 * kappa;
  double kappa5 = kappa4 * kappa;
  double kappa6 = kappa5 * kappa;
  double kappa7 = kappa6 * kappa;

  /* NOTE that |m| will be used to determine the fit to use, and if input_m < 0, then a conjugate will be taken */
  int m = abs(input_m);

  /**/
  complex double j = _Complex_I;

  /* Initialize the answer */
  double complex ans = 0.0;

  /* Use If-Else ladder to determine which mode function to evaluate */
  if ( 2==l && 2==m && 0==n  ){

    /* Fit for (l,m,n) == (2,2,0). This is a zero-damped mode in the extremal Kerr limit. */
    ans = 0.55262405 + 6.54272463*cexp(0.24443847*j)*kappa + 5.94664565*cexp(3.88409012*j)*kappa2 + 5.39298183*cexp(1.01651284*j)*kappa3 + 3.58701474*cexp(4.53395559*j)*kappa4 + 1.36858235*cexp(1.57079633*j)*kappa5 + 0.18520700*cexp(4.71238898*j)*kappa6 ;

    } else if ( 2==l && 2==m && 1==n ) {

      /* Fit for (l,m,n) == (2,2,1). This is a zero-damped mode in the extremal Kerr limit. */
      ans = 0.55229247 + 7.94074969*cexp(0.64081239*j)*kappa + 12.55567057*cexp(4.41980669*j)*kappa2 + 13.68518711*cexp(1.48039237*j)*kappa3 + 10.43884041*cexp(4.72599435*j)*kappa4 + 4.20731453*cexp(1.57079633*j)*kappa5 + 0.76232588*cexp(4.71238898*j)*kappa6;

    } else if ( 3==l && 2==m && 0==n ) {

      /* Fit for (l,m,n) == (3,2,0). This is NOT a zero-damped mode in the extremal Kerr limit.  */
      ans = 8.18542769*cexp(6.27603422*j) + 1.55192720*cexp(1.79088081*j)*kappa + 8.94654695*cexp(5.18681710*j)*kappa2 + 28.66050158*cexp(1.63658858*j)*kappa3 + 60.77789497*cexp(4.72114050*j)*kappa4 + 72.13239907*cexp(1.57079633*j)*kappa5 + 45.38115278*cexp(4.71238898*j)*kappa6 + 11.84706755*cexp(1.57079633*j)*kappa7;

    } else if ( 4==l && 4==m && 0==n ) {

      /* Fit for (l,m,n) == (4,4,0). This is a zero-damped mode in the extremal Kerr limit. */
      ans = 13.05294185 + 9.23462388*cexp(0.14179514*j)*kappa + 7.09045393*cexp(3.69184561*j)*(kappa2) + 6.46711175*cexp(0.89254551*j)*(kappa3) + 4.96905278*cexp(4.43853588*j)*(kappa4) + 2.62299932*cexp(1.57079633*j)*(kappa5) + 0.58168681*cexp(4.71238898*j)*(kappa6);

    } else if ( 2==l && 1==m && 0==n ) {

      /* Fit for (l,m,n) == (2,1,0). This is NOT a zero-damped mode in the extremal Kerr limit. */
      ans = 3.10089518*cexp(6.25822093*j) + 2.69208437*cexp(1.95853947*j)*kappa + 16.58575360*cexp(4.98423605*j)*kappa2 + 57.84090876*cexp(1.63720921*j)*kappa3 + 118.21761290*cexp(4.72674943*j)*kappa4 + 135.93985738*cexp(1.57079633*j)*kappa5 + 82.81742189*cexp(4.71238898*j)*kappa6 + 20.85173245*cexp(1.57079633*j)*kappa7;

    } else if ( 3==l && 3==m && 0==n ) {

      /* Fit for (l,m,n) == (3,3,0). This is a zero-damped mode in the extremal Kerr limit.*/
      ans = 5.70465254 + 7.94433155*cexp(0.18039136*j)*kappa + 6.55099749*cexp(3.77926384*j)*kappa2 + 6.31422768*cexp(0.93863733*j)*kappa3 + 4.81214531*cexp(4.46906976*j)*kappa4 + 2.38927043*cexp(1.57079633*j)*kappa5 + 0.48077965*cexp(4.71238898*j)*kappa6;

    } else if ( 3==l && 3==m && 1==n ) {

      /* Fit for (l,m,n) == (3,3,1). This is a zero-damped mode in the extremal Kerr limit. */
      ans = 5.70318420 + 8.94926548*cexp(0.49834140*j)*kappa + 12.70528736*cexp(4.31772419*j)*kappa2 + 15.63533560*cexp(1.39390017*j)*kappa3 + 14.19057659*cexp(4.66913674*j)*kappa4 + 7.33238119*cexp(1.57079633*j)*kappa5 + 1.53701758*cexp(4.71238898*j)*kappa6;

    } else if ( 4==l && 3==m && 0==n ) {

      /* Fit for (l,m,n) == (4,3,0). This is a zero-damped mode in the extremal Kerr limit. */
      ans = 15.28866348 + 0.75297352*cexp(0.22048290*j)*kappa + 3.64936150*cexp(0.61644055*j)*kappa2 + 8.02530641*cexp(4.82756576*j)*kappa3 + 12.47205664*cexp(1.67334685*j)*kappa4 + 10.30282199*cexp(4.71238898*j)*kappa5 + 3.52885679*cexp(1.57079633*j)*kappa6;

    } else if ( 5==l && 5==m && 0==n ) {

      /* Fit for (l,m,n) == (5,5,0). This is a zero-damped mode in the extremal Kerr limit. */
      ans = 22.52292196 + 10.44137664*cexp(0.11607502*j)*kappa + 7.79707643*cexp(3.61247422*j)*kappa2 + 6.59989026*cexp(0.83792606*j)*kappa3 + 4.90367451*cexp(4.40545635*j)*kappa4 + 2.59913853*cexp(1.57079633*j)*kappa5 + 0.58985077*cexp(4.71238898*j)*kappa6;

    } else {

      /**/
      ans  = 0.0;

  } /* END of IF-ELSE Train for QNM cases */

  /* If m<0, then take the conjugate */
  if ( input_m < 0 ) {
    /**/
    ans = conj( ans );
  }

  return ans;

} /* END of separationConstant */


/*
* Spheroidal Harmonic Normalization Constants: Note that name encodes date of writing
*/
static double spheroidalHarmonicNormalization( double kappa, int l, int input_m, int n );
static double spheroidalHarmonicNormalization( double kappa,  /* Domain mapping for remnant BH's spin (Dimensionless)*/
                          int l,        /* Polar eigenvalue*/
                          int input_m,  /* Azimuthal eigenvalue*/
                          int n ) {     /* Overtone Number*/

  /* Predefine powers to increase efficiency*/
  double kappa2 = kappa  * kappa;
  double kappa3 = kappa2 * kappa;
  double kappa4 = kappa3 * kappa;
  double kappa5 = kappa4 * kappa;
  double kappa6 = kappa5 * kappa;

  /* NOTE that |m| will be used to determine the fit to use, and if input_m < 0, then a conjugate will be taken*/
  int m = abs(input_m);

  /* Initialize the answer*/
  double complex ans = 0.0;

  /* Use If-Else ladder to determine which mode function to evaluate */
  if ( 2==l && 2==m && 0==n  ){

    /* Fit for (l,m,n) == (2,2,0). This is a zero-damped mode in the extremal Kerr limit.*/
    ans = 7.86366171 - 3.61447483*kappa + 3.48996689*kappa2 - 2.29347705*kappa3 + 0.74425069*kappa4 ;

    } else if ( 2==l && 2==m && 1==n ) {

      /* Fit for (l,m,n) == (2,2,1). This is a zero-damped mode in the extremal Kerr limit. */
      ans = 7.86298703 - 3.59872285*kappa + 2.88459437*kappa2 - 0.92740734*kappa3 - 0.04445478*kappa4;

    } else if ( 3==l && 2==m && 0==n ) {

      /* Fit for (l,m,n) == (3,2,0). This is NOT a zero-damped mode in the extremal Kerr limit.*/
      ans = 0.74845717 - 0.08157463*kappa + 1.03748092*kappa2 - 3.27926931*kappa3 + 7.24584503*kappa4 - 7.41316799*kappa5 + 3.06056035*kappa6;

    } else if ( 4==l && 4==m && 0==n ) {

      /* Fit for (l,m,n) == (4,4,0). This is a zero-damped mode in the extremal Kerr limit.*/
      ans = 1.75389888 + 1.00111258*kappa + 1.55498487*kappa2 - 1.22344804*kappa3 + 1.64621074*kappa4;

    } else if ( 2==l && 1==m && 0==n ) {

      /* Fit for (l,m,n) == (2,1,0). This is NOT a zero-damped mode in the extremal Kerr limit.*/
      ans = 3.04393302 - 0.06877527*kappa + 0.87671129*kappa2 - 3.92206769*kappa3 + 8.59631959*kappa4 - 8.52199526*kappa5 + 3.31150324*kappa6;

    } else if ( 3==l && 3==m && 0==n ) {

      /* Fit for (l,m,n) == (3,3,0). This is a zero-damped mode in the extremal Kerr limit.*/
      ans = 3.51631915 + 0.16499714*kappa + 1.30114387*kappa2 - 0.83622153*kappa3 + 0.82020713*kappa4;

    } else if ( 3==l && 3==m && 1==n ) {

      /* Fit for (l,m,n) == (3,3,1). This is a zero-damped mode in the extremal Kerr limit.*/
      ans = 3.51530809 + 0.19285707*kappa + 0.96814190*kappa2 - 0.00547882*kappa3 + 0.24982172*kappa4;

    } else if ( 4==l && 3==m && 0==n ) {

      /* Fit for (l,m,n) == (4,3,0). This is a zero-damped mode in the extremal Kerr limit.*/
      ans = 0.39542385 - 0.09918352*kappa + 1.52850262*kappa2 - 5.09932727*kappa3 + 10.95647104*kappa4 - 10.99914124*kappa5 + 4.52212985*kappa6;

    } else if ( 5==l && 5==m && 0==n ) {

      /* Fit for (l,m,n) == (5,5,0). This is a zero-damped mode in the extremal Kerr limit. */
      ans = 0.91349889 + 0.89568178*kappa + 2.54404526*kappa2 - 2.82437113*kappa3 + 3.28143852*kappa4;

    } else {

      /**/
      ans  = 0.0;

  } /* END of IF-ELSE Train for QNM cases */

  return ans;

} /* END of spheroidalHarmonicNormalization */


/*
* Final Spin (Dimensionless)
*/
UNUSED static double finalSpinFit( double eta );
UNUSED static double finalSpinFit( double eta ) {
  /* Implement Final Spin Fit from arXiv:1404.3197 */
  return eta*( 3.4339 - eta*( 3.7988 + eta*(5.7733 - 6.3780*eta) ) );
}

/*
* Final Mass (Dimensionless: Relative to M=1 initially)
*/
UNUSED static double finalMassFit( double eta );
UNUSED static double finalMassFit( double eta ) {
  /* Implement Final Mass Fit from arXiv:1404.3197 */
  return 1.0 - eta*(0.046297 - eta*(0.71006 + eta*(1.5028 - eta*(4.0124 - eta*0.28448))));
}


/* ------------------------------------------------
          Angular parameter functions
 ------------------------------------------------ */
static double K1( int m, int s );
static double K1( int m, int s ){
  return 0.5*abs(m-s);
}
static double K2( int m, int s );
static double K2( int m, int s ){
  return 0.5*abs(m+s);
}
static complex double ALPHA_RD( int m, int s, int p );
static complex double ALPHA_RD( int m, int s, int p ){
  /**/
  double k1 = 0.5*abs(m-s);
  return -2.0*(p+1.0)*(p+2.0*k1+1.0);
}
static complex double BETA_RD( int m, int s, int p, complex double aw, complex double A_lm );
static complex double BETA_RD( int m, int s, int p, complex double aw, complex double A_lm ){
  /**/
  double k1 = K1(m,s);
  double k2 = K2(m,s);
  return  p*(p-1.0)+2.0*p*(k1+k2+1.0-2.0*aw) - ( 2.0*aw*(2.0*k1+s+1.0)-(k1+k2) * (k1+k2+1.0) ) - ( aw*aw + s*(s+1.0) + A_lm);
}
static complex double GAMMA_RD( int m, int s, int p, complex double aw );
static complex double GAMMA_RD( int m, int s, int p, complex double aw ){
  /**/
  double k1 = K1(m,s);
  double k2 = K2(m,s);
  return 2.0*aw*(p+k1+k2+s);
}

/*
* QNM Ampliutde models for MMRDNS
* NOTE that the terms here differ from 1404.3197v3 for accurate relative phases
*/
COMPLEX16 XLALMMRDNSAmplitudeOverOmegaSquared( double eta, int l, int input_m, int n );
COMPLEX16 XLALMMRDNSAmplitudeOverOmegaSquared( double eta, int l, int input_m, int n ){

  /* Initialize the answer */
  double complex ans = 0.0;

  /* NOTE that |m| will be used to determine the fit to use, and if input_m < 0, then a conjugate will be taken */
  int m = abs(input_m);

  /**/
  double eta2 = eta*eta;
  double eta3 = eta2*eta;
  double eta4 = eta3*eta;
  double eta5 = eta4*eta;

  /*  Evaluate QNM amplitude models for input l m n */
  if ( l==2 && m==2 && n==0 ) {

    /*A2201*/
    ans = 0.95846504*cexp(2.99318408*I)*eta + 0.47588079*cexp(0.82658128*I)*(eta2) + 1.23853419*cexp(2.30528861*I)*(eta3);

  } else if ( l==2 && m==2 && n==1 ) {

    /*A2211*/
    ans = 0.12750415*cexp(0.05809736*I)*eta + 1.18823931*cexp(1.51798243*I)*(eta2) + 8.27086561*cexp(4.42014780*I)*(eta3) + 26.23294960*cexp(1.16782950*I)*(eta4);

  } else if ( l==2 && m==1 && n==0 ) {

    /*A2101*/
    ans = sqrt(1-4*eta)  * (  0.47952344*cexp(5.96556090*I)*eta + 1.17357614*cexp(3.97472217*I)*(eta2) + 1.23033028*cexp(2.17322465*I)*(eta3)  );

  } else if ( l==3 && m==3 && n==0 ) {

    /*A3301*/
    ans = sqrt(1-4*eta)  * (  0.42472339*cexp(4.54734400*I)*eta + 1.47423728*cexp(2.70187807*I)*(eta2) + 4.31385024*cexp(5.12815819*I)*(eta3) + 15.72642073*cexp(2.25473854*I)*(eta4)  );

  } else if ( l==3 && m==3 && n==1 ) {

    /*A3311*/
    ans = sqrt(1-4*eta)  * (  0.14797161*cexp(2.03957081*I)*eta + 1.48738894*cexp(5.89538621*I)*(eta2) + 10.16366839*cexp(3.28354928*I)*(eta3) + 29.47859659*cexp(0.81061521*I)*(eta4)  );

  } else if ( l==3 && m==2 && n==0 ) {

    /*A3201*/
    ans = 0.19573228*cexp(0.54325509*I)*eta + 1.58299638*cexp(4.24509590*I)*(eta2) + 5.03380859*cexp(1.71003281*I)*(eta3) + 3.73662711*cexp(5.14735754*I)*(eta4);

  } else if ( l==4 && m==4 && n==0 ) {

    /*A4401*/
    ans = 0.25309908*cexp(5.16320109*I)*eta + 2.40404787*cexp(2.46899414*I)*(eta2) + 14.72733952*cexp(5.56235208*I)*(eta3) + 67.36237809*cexp(2.19824119*I)*(eta4) + 126.58579931*cexp(5.41735031*I)*(eta5);

  } else if ( l==4 && m==3 && n==0 ) {

    /*A4301*/
    ans = sqrt(1-4*eta)  * (  0.09383417*cexp(2.30765661*I)*eta + 0.82734483*cexp(6.10053234*I)*(eta2) + 3.33846327*cexp(3.87329126*I)*(eta3) + 4.66386840*cexp(1.75165690*I)*(eta4)  );

  } else if ( l==5 && m==5 && n==0 ) {

    /*A5501*/
    ans = sqrt(1-4*eta)  * (  0.15477314*cexp(1.06752431*I)*eta + 1.50914172*cexp(4.54983062*I)*(eta2) + 8.93331690*cexp(1.28981042*I)*(eta3) + 42.34309620*cexp(4.10035598*I)*(eta4) + 89.19466498*cexp(1.02508947*I)*(eta5)  );

  }

  /* If m<0, then take the conjugate */
  if ( input_m < 0 ) {
    /**/
    ans = conj( ans );
  }

  return ans;

}

/*
* Spheroical Harmonic Functions (Leaver's Formulation circa 1986/85)
*/
COMPLEX16 XLALSpinWeightedSpheroidalHarmonic( double jf,           /* Spin of remnant */
                   int l, int m, int n, /* QNM indices */
                   double theta,        /* polar angle */
                   double phi          /* azimuthal angle */
                 ) {

  /* Set spin weight */
  const int s = -2;

  /* Use tabulated cw and sc values from the core package*/
  double complex cw, sc, aw;
  double kappa = XLALKAPPA(jf,l,m);

  cw = XLALcomplexOmega( kappa, l, m, n );
  sc = XLALseparationConstant( kappa, l, m, n );

  /* Define dimensionless deformation parameter */
  aw = jf*cw;

  /* ------------------------------------------------
      Calculate the angular eigenfunction
   ------------------------------------------------ */

  /* Variable map for theta */
  double u = cos(theta);

  /* the non-sum part of eq 18 */
  complex double X = 1.0;
  X = X * cexp(aw*u) * pow(1.0+u, K1(m,s) );
  X = X * pow(1.0-u,K2(m,s));

  /* NOTE that here we apply the normalization constant */
  X = X / spheroidalHarmonicNormalization(kappa,l,m,n);

  /* initial series values */
  double a0 = 1.0;

  double a1 = -BETA_RD( m, s, 0, aw, sc )/ALPHA_RD( m, s, 0 );

  /* the sum part */
  complex double Y = a0;
  complex double S;
  complex double dY;
  double xx;
  complex double a2;
  int done = 0;
  int k, j;
  Y = Y + a1*(1.0+u);
  k = 1;
  int kmax = 2e3;
  double et = 1e-8;
  while ( ! done ) {
    k += 1;
    j = k-1;
    a2 = -1.0*( BETA_RD( m, s, j, aw, sc )*a1 + GAMMA_RD(m,s,j,aw)*a0 ) / ALPHA_RD(m,s,j);
    dY = pow(a2*(1.0+u),k);
    Y += dY;
    xx = cabs( dY );

    done = (k>=l) && ( (xx<et && k>30) || k>kmax );
    done = done || xx<et;
    a0 = a1;
    a1 = a2;
  }

  /* together now */
  S = X*Y*cexp( _Complex_I * m * phi );

  /* Use the same sign convention as spherical harmonics
  e.g http://en.wikipedia.org/wiki/Spin-weighted_spherical_harmonics#Calculating */
  double C = 1.0;
  C = C * pow(-1,fmax(-m,-s)) * pow(-1,l);
  S = C * S;

  /**/
  return S;

} /* END of Spheroical Harmonic function */

/* Mode included within the model */
UNUSED static const int MMRDNS_MODES[9][3] =  { {2,2,0},
                                         {2,2,1},
                                         {3,3,0},
                                         {3,3,1},
                                         {4,4,0},
                                         {5,5,0},
                                         {2,1,0},
                                         {3,2,0},
                                         {4,3,0}
                                       };

/**/
int XLALSimRingdownMMRDNSTD(
        UNUSED REAL8TimeSeries **hplus,        /**< OUTPUT h_+ vector */
        UNUSED REAL8TimeSeries **hcross,       /**< OUTPUT h_x vector */
        UNUSED REAL8 phiRef,                   /**< orbital phase at reference pt. */
        UNUSED REAL8 inclination,              /**< inclination angle */
        UNUSED REAL8 deltaT,                   /**< sampling interval (s) */
        UNUSED REAL8 m1,                       /**< mass of companion 1 (kg) */
        UNUSED REAL8 m2,                       /**< mass of companion 2 (kg) */
        UNUSED REAL8 r,                        /**< distance of source (m) */
        UNUSED const LALSimInspiralTestGRParam *nonGRparams /**< linked list containing the extra testing GR parameters */
        ){

        /* Declarations */
        const LIGOTimeGPS T0=LIGOTIMEGPSZERO;
        double kappa = XLALKAPPA( 0.68, 2, 2 );
        UNUSED double w22 = creal( XLALcomplexOmega(kappa,2,2,0) );
        UNUSED double tau22 = cimag( XLALcomplexOmega(kappa,2,2,0) );
        int waveform_length = 1000;


        /* Create the plus and cross polarization vectors */
        *hplus = XLALCreateREAL8TimeSeries("hplus", &T0, 0.0, deltaT, &lalStrainUnit, waveform_length);
        *hcross = XLALCreateREAL8TimeSeries("hcross", &T0, 0.0, deltaT, &lalStrainUnit, waveform_length);


  return 0;
}


/* XLALSimRingdownSingleModeTD: Time domain waveformgenerator for single QNM without angular dependence (i.e. this function generates multipole moments only ). In  */
int XLALSimRingdownGenerateSingleModeTD(
  UNUSED COMPLEX16TimeSeries **hk,            /**< OUTPUT vector for QNM time series */
  UNUSED const LIGOTimeGPS T0,                /**< Start time  */
  UNUSED REAL8 deltaT,                        /**< sampling interval (s) */
  UNUSED REAL8 Nsamples,                      /**< Number of time domain samples */
  UNUSED complex double Ak,                   /**< COMPLEX QNM Strain Amplitude */
  UNUSED complex double CWk                   /**< COMPLEX QNM Frequency */
){

  /*  */
  *hk = XLALCreateCOMPLEX16TimeSeries( "SingleRDTDMode", &T0, 0.0, deltaT, &lalStrainUnit, Nsamples);
  memset( (*hk)->data->data, 0, Nsamples*sizeof(COMPLEX16) ) ;

  /**/
  COMPLEX16 cwdt = 0;

  /**/
  cwdt = 0; /* CWk*deltaT; */

  /**/
  int k;
  for ( k=0; k < Nsamples; ++k ){

    /**/
    COMPLEX16 h;
    h = Ak * cexp(-I * cwdt * k);
    (*hk)->data->data[k] = h;

  }

  return 0;

}


/* XLALSimRingdownMMRDNSFD: Frequency domain waveformgenerator for all QNM with angular dependence */
int XLALSimRingdownMMRDNSFD(
        COMPLEX16FrequencySeries **hptilde,          /**< OUTPUT FD h_+ polarization */
        COMPLEX16FrequencySeries **hctilde,          /**< OUTPUT FD h_x polarization */
        const REAL8 deltaF,                          /**< Frequency resolution (Hz) */
        const REAL8 fStart,                          /**< Start GW frequency (Hz) */
        const REAL8 fEnd,                            /**< Highest GW frequency (Hz) */
        REAL8 Mf,                                    /**< Final BH Mass (kg) */
        REAL8 jf,                                    /**< Final BH dimensionaless spin */
        REAL8 eta,                                   /**< Symmetric mass ratio of two companions */
        REAL8 iota,                                  /**< inclination angle (in rad) */
        REAL8 phi_offset,                            /**< intrinsic phase offset */
        REAL8 r,                                     /**< distance of source (m) */
        LALSimInspiralTestGRParam *nonGRparams       /**< testing GR parameters */
    ){
        /* Perform some initial checks */
        if (Mf <= 0) XLAL_ERROR(XLAL_EDOM);
        if (jf >= 1) XLAL_ERROR(XLAL_EDOM);
        if (eta > 0.25 || eta <= 0) XLAL_ERROR(XLAL_EDOM);
        if (fStart <= 0) XLAL_ERROR(XLAL_EDOM);
        if (r <= 0) XLAL_ERROR(XLAL_EDOM);

        /* Declarations */
        UINT4 jStart, jMax;
        REAL8 f_max;
        LIGOTimeGPS tC = {0,0};
        //XLALGPSAdd(&tC, -1 / deltaF);

        REAL8 dfreq220=0., dfreq221=0., dfreq330=0., dfreq331=0., dfreq440=0., dfreq550=0., dfreq210=0., dfreq320=0., dfreq430=0.;
        REAL8 dtau220=0.,  dtau221=0.,  dtau330=0.,  dtau331=0.,  dtau440=0.,  dtau550=0.,  dtau210=0.,  dtau320=0.,  dtau430=0.;

        COMPLEX16FrequencySeries *hp220, *hp221, *hp330, *hp331, *hp440, *hp550, *hp210, *hp320, *hp430;
        COMPLEX16FrequencySeries *hc220, *hc221, *hc330, *hc331, *hc440, *hc550, *hc210, *hc320, *hc430;

        /* Get nonGRparams */
        char *nonGRParamName = malloc(512*sizeof(char));

        sprintf(nonGRParamName,"dfreq220") ;
        if (XLALSimInspiralTestGRParamExists(nonGRparams,nonGRParamName) )
          dfreq220 = XLALSimInspiralGetTestGRParam(nonGRparams,nonGRParamName) ;
        sprintf(nonGRParamName,"dtau220") ;
        if (XLALSimInspiralTestGRParamExists(nonGRparams,nonGRParamName) )
          dtau220 = XLALSimInspiralGetTestGRParam(nonGRparams,nonGRParamName) ;
        sprintf(nonGRParamName,"dfreq221") ;
        if (XLALSimInspiralTestGRParamExists(nonGRparams,nonGRParamName) )
          dfreq221 = XLALSimInspiralGetTestGRParam(nonGRparams,nonGRParamName) ;
        sprintf(nonGRParamName,"dtau221") ;
        if (XLALSimInspiralTestGRParamExists(nonGRparams,nonGRParamName) )
          dtau221 = XLALSimInspiralGetTestGRParam(nonGRparams,nonGRParamName) ;
        sprintf(nonGRParamName,"dfreq330") ;
        if (XLALSimInspiralTestGRParamExists(nonGRparams,nonGRParamName) )
          dfreq330 = XLALSimInspiralGetTestGRParam(nonGRparams,nonGRParamName) ;
        sprintf(nonGRParamName,"dtau330") ;
        if (XLALSimInspiralTestGRParamExists(nonGRparams,nonGRParamName) )
          dtau330 = XLALSimInspiralGetTestGRParam(nonGRparams,nonGRParamName) ;
        sprintf(nonGRParamName,"dfreq331") ;
        if (XLALSimInspiralTestGRParamExists(nonGRparams,nonGRParamName) )
          dfreq331 = XLALSimInspiralGetTestGRParam(nonGRparams,nonGRParamName) ;
        sprintf(nonGRParamName,"dtau331") ;
        if (XLALSimInspiralTestGRParamExists(nonGRparams,nonGRParamName) )
          dtau331 = XLALSimInspiralGetTestGRParam(nonGRparams,nonGRParamName) ;
        sprintf(nonGRParamName,"dfreq440") ;
        if (XLALSimInspiralTestGRParamExists(nonGRparams,nonGRParamName) )
          dfreq440 = XLALSimInspiralGetTestGRParam(nonGRparams,nonGRParamName) ;
        sprintf(nonGRParamName,"dtau440") ;
        if (XLALSimInspiralTestGRParamExists(nonGRparams,nonGRParamName) )
          dtau440 = XLALSimInspiralGetTestGRParam(nonGRparams,nonGRParamName) ;
        sprintf(nonGRParamName,"dfreq550") ;
        if (XLALSimInspiralTestGRParamExists(nonGRparams,nonGRParamName) )
          dfreq550 = XLALSimInspiralGetTestGRParam(nonGRparams,nonGRParamName) ;
        sprintf(nonGRParamName,"dtau550") ;
        if (XLALSimInspiralTestGRParamExists(nonGRparams,nonGRParamName) )
          dtau550 = XLALSimInspiralGetTestGRParam(nonGRparams,nonGRParamName) ;
        sprintf(nonGRParamName,"dfreq210") ;
        if (XLALSimInspiralTestGRParamExists(nonGRparams,nonGRParamName) )
          dfreq210 = XLALSimInspiralGetTestGRParam(nonGRparams,nonGRParamName) ;
        sprintf(nonGRParamName,"dtau210") ;
        if (XLALSimInspiralTestGRParamExists(nonGRparams,nonGRParamName) )
          dtau210 = XLALSimInspiralGetTestGRParam(nonGRparams,nonGRParamName) ;
        sprintf(nonGRParamName,"dfreq320") ;
        if (XLALSimInspiralTestGRParamExists(nonGRparams,nonGRParamName) )
          dfreq320 = XLALSimInspiralGetTestGRParam(nonGRparams,nonGRParamName) ;
        sprintf(nonGRParamName,"dtau320") ;
        if (XLALSimInspiralTestGRParamExists(nonGRparams,nonGRParamName) )
          dtau320 = XLALSimInspiralGetTestGRParam(nonGRparams,nonGRParamName) ;
        sprintf(nonGRParamName,"dfreq430") ;
        if (XLALSimInspiralTestGRParamExists(nonGRparams,nonGRParamName) )
          dfreq430 = XLALSimInspiralGetTestGRParam(nonGRparams,nonGRParamName) ;
        sprintf(nonGRParamName,"dtau430") ;
        if (XLALSimInspiralTestGRParamExists(nonGRparams,nonGRParamName) )
          dtau430 = XLALSimInspiralGetTestGRParam(nonGRparams,nonGRParamName) ;

        /* allocate hptilde and hctilde */
        if ( fEnd == 0. ){
          f_max = 2048.0;
        } else {
          f_max = fEnd;
        }
        jMax = (UINT4)(f_max/deltaF + 1);

        XLALSimRingdownGenerateSingleModeFD(&hp220, &hc220, deltaF, fStart, fEnd, Mf, jf, eta, iota, phi_offset, 2, 2, 0, r, dfreq220, dtau220);
        XLALSimRingdownGenerateSingleModeFD(&hp221, &hc221, deltaF, fStart, fEnd, Mf, jf, eta, iota, phi_offset, 2, 2, 1, r, dfreq221, dtau221);
        XLALSimRingdownGenerateSingleModeFD(&hp330, &hc330, deltaF, fStart, fEnd, Mf, jf, eta, iota, phi_offset, 3, 3, 0, r, dfreq330, dtau330);
        XLALSimRingdownGenerateSingleModeFD(&hp331, &hc331, deltaF, fStart, fEnd, Mf, jf, eta, iota, phi_offset, 3, 3, 1, r, dfreq331, dtau331);
        XLALSimRingdownGenerateSingleModeFD(&hp440, &hc440, deltaF, fStart, fEnd, Mf, jf, eta, iota, phi_offset, 4, 4, 0, r, dfreq440, dtau440);
        XLALSimRingdownGenerateSingleModeFD(&hp550, &hc550, deltaF, fStart, fEnd, Mf, jf, eta, iota, phi_offset, 5, 5, 0, r, dfreq550, dtau550);
        XLALSimRingdownGenerateSingleModeFD(&hp210, &hc210, deltaF, fStart, fEnd, Mf, jf, eta, iota, phi_offset, 2, 1, 0, r, dfreq210, dtau210);
        XLALSimRingdownGenerateSingleModeFD(&hp320, &hc320, deltaF, fStart, fEnd, Mf, jf, eta, iota, phi_offset, 3, 2, 0, r, dfreq320, dtau320);
        XLALSimRingdownGenerateSingleModeFD(&hp430, &hc430, deltaF, fStart, fEnd, Mf, jf, eta, iota, phi_offset, 4, 3, 0, r, dfreq430, dtau430);

        *hptilde = XLALCreateCOMPLEX16FrequencySeries("hptilde: FD waveform", &tC, 0.0, deltaF, &lalStrainUnit, jMax);
        if (!(*hptilde)) XLAL_ERROR(XLAL_EFUNC);
        memset((*hptilde)->data->data, 0, jMax * sizeof(COMPLEX16));
        XLALUnitMultiply(&(*hptilde)->sampleUnits, &(*hptilde)->sampleUnits, &lalSecondUnit);
        *hctilde = XLALCreateCOMPLEX16FrequencySeries("hctilde: FD waveform", &tC, 0.0, deltaF, &lalStrainUnit, jMax);
        if (!(*hctilde)) XLAL_ERROR(XLAL_EFUNC);
        memset((*hctilde)->data->data, 0, jMax * sizeof(COMPLEX16));
        XLALUnitMultiply(&(*hctilde)->sampleUnits, &(*hctilde)->sampleUnits, &lalSecondUnit);

        jStart   = (UINT4) ceil(fStart / deltaF);

        for ( UINT4 j=jStart ; j<jMax ; j++ ) {
          (*hptilde)->data->data[j] = hp220->data->data[j];
          (*hptilde)->data->data[j] += hp221->data->data[j];
          (*hptilde)->data->data[j] += hp330->data->data[j];
          (*hptilde)->data->data[j] += hp331->data->data[j];
          (*hptilde)->data->data[j] += hp440->data->data[j];
          (*hptilde)->data->data[j] += hp550->data->data[j];
          (*hptilde)->data->data[j] += hp210->data->data[j];
          (*hptilde)->data->data[j] += hp320->data->data[j];
          (*hptilde)->data->data[j] += hp430->data->data[j];

          (*hctilde)->data->data[j] = hc220->data->data[j];
          (*hctilde)->data->data[j] += hc221->data->data[j];
          (*hctilde)->data->data[j] += hc330->data->data[j];
          (*hctilde)->data->data[j] += hc331->data->data[j];
          (*hctilde)->data->data[j] += hc440->data->data[j];
          (*hctilde)->data->data[j] += hc550->data->data[j];
          (*hctilde)->data->data[j] += hc210->data->data[j];
          (*hctilde)->data->data[j] += hc320->data->data[j];
          (*hctilde)->data->data[j] += hc430->data->data[j];
        }

        /* Destroy the COMPLEX16FrequencySeries object */
        if (hp220) XLALDestroyCOMPLEX16FrequencySeries(hp220);
        if (hp221) XLALDestroyCOMPLEX16FrequencySeries(hp221);
        if (hp330) XLALDestroyCOMPLEX16FrequencySeries(hp330);
        if (hp331) XLALDestroyCOMPLEX16FrequencySeries(hp331);
        if (hp440) XLALDestroyCOMPLEX16FrequencySeries(hp440);
        if (hp550) XLALDestroyCOMPLEX16FrequencySeries(hp550);
        if (hp210) XLALDestroyCOMPLEX16FrequencySeries(hp210);
        if (hp320) XLALDestroyCOMPLEX16FrequencySeries(hp320);
        if (hp430) XLALDestroyCOMPLEX16FrequencySeries(hp430);

        if (hc220) XLALDestroyCOMPLEX16FrequencySeries(hc220);
        if (hc221) XLALDestroyCOMPLEX16FrequencySeries(hc221);
        if (hc330) XLALDestroyCOMPLEX16FrequencySeries(hc330);
        if (hc331) XLALDestroyCOMPLEX16FrequencySeries(hc331);
        if (hc440) XLALDestroyCOMPLEX16FrequencySeries(hc440);
        if (hc550) XLALDestroyCOMPLEX16FrequencySeries(hc550);
        if (hc210) XLALDestroyCOMPLEX16FrequencySeries(hc210);
        if (hc320) XLALDestroyCOMPLEX16FrequencySeries(hc320);
        if (hc430) XLALDestroyCOMPLEX16FrequencySeries(hc430);

  return XLAL_SUCCESS;
}

/* XLALSimRingdownGenerateSingleModeFD: Frequency domain waveformgenerator for single QNM with angular dependence */
int XLALSimRingdownGenerateSingleModeFD(
        COMPLEX16FrequencySeries **hptilde_lmn,      /**< OUTPUT FD h_+ polarization */
        COMPLEX16FrequencySeries **hctilde_lmn,      /**< OUTPUT FD h_x polarization */
        const REAL8 deltaF,                          /**< Frequency resolution (Hz) */
        const REAL8 fStart,                          /**< Start GW frequency (Hz) */
        const REAL8 fEnd,                            /**< Highest GW frequency (Hz) */
        REAL8 Mf,                                    /**< Final BH Mass (kg) */
        REAL8 jf,                                    /**< Final BH dimensionaless spin */
        REAL8 eta,                                   /**< Symmetric mass ratio of two companions */
        REAL8 iota,                                  /**< inclination angle (in rad) */
        REAL8 phi_offset,                            /**< intrinsic phase offset (in rad) */
        UINT4 l,                                     /**< Polar eigenvalue */
        UINT4 m,                                     /**< Azimuthal eigenvalue */
        UINT4 n,                                     /**< Overtone Number */
        REAL8 r,                                     /**< distance of source (m) */
        REAL8 dfreq,                                 /**< relative shift in the real frequency parameter */
        REAL8 dtau                                   /**< relative shift in the damping time parameter */
        ){

        /* Perform some initial checks */
        if (Mf <= 0) XLAL_ERROR(XLAL_EDOM);
        if (jf >= 1) XLAL_ERROR(XLAL_EDOM);
        if (eta > 0.25 || eta <= 0) XLAL_ERROR(XLAL_EDOM);
        if (fStart <= 0) XLAL_ERROR(XLAL_EDOM);
        if (r <= 0) XLAL_ERROR(XLAL_EDOM);

        /* Declarations */
        UINT4 jStart, jMax;
        COMPLEX16 h_f, hconj_mf;
        REAL8 f, f_max;
        LIGOTimeGPS tC = {0,0};
        //XLALGPSAdd(&tC, -1 / deltaF);

        REAL8 kappa   = XLALKAPPA( jf, l, m );
        REAL8 Mf_sec  = Mf*LAL_MTSUN_SI/LAL_MSUN_SI;
        REAL8 r_sec   = r/LAL_C_SI;

        COMPLEX16 A_lmn, S_lmn, Omega_lmn, Prefactor;

        /* Mode Component Calculation*/
        Omega_lmn = XLALcomplexOmega(kappa, l, m, n)/Mf_sec;
        A_lmn = XLALMMRDNSAmplitudeOverOmegaSquared(eta, l, m, n);
        Omega_lmn = creal(Omega_lmn)*(1.+dfreq) + I * cimag(Omega_lmn)/(1.+dtau);
        S_lmn = XLALSpinWeightedSpheroidalHarmonic(jf, l, m, n, iota, 0.0);
        Prefactor = cexp(I*phi_offset)*(Mf_sec/r_sec)*(A_lmn*S_lmn)*(-I);

        /* allocate htilde_p and htilde_c */
        /* The COMPLEX16FrequencySeries has to be destroyed by whom created it. */
        if ( fEnd == 0. ){
          f_max = 2048.0;
        } else {
          f_max = fEnd;
        }
 
        jMax = (UINT4)(f_max/deltaF + 1);

        *hptilde_lmn = XLALCreateCOMPLEX16FrequencySeries("hptilde_lmn: FD waveform", &tC, 0.0, deltaF, &lalStrainUnit, jMax);
        if (!(*hptilde_lmn)) XLAL_ERROR(XLAL_EFUNC);
        memset((*hptilde_lmn)->data->data, 0, jMax * sizeof(COMPLEX16));
        XLALUnitMultiply(&(*hptilde_lmn)->sampleUnits, &(*hptilde_lmn)->sampleUnits, &lalSecondUnit);
        *hctilde_lmn = XLALCreateCOMPLEX16FrequencySeries("hctilde_lmn: FD waveform", &tC, 0.0, deltaF, &lalStrainUnit, jMax);
        if (!(*hctilde_lmn)) XLAL_ERROR(XLAL_EFUNC);
        memset((*hctilde_lmn)->data->data, 0, jMax * sizeof(COMPLEX16));
        XLALUnitMultiply(&(*hctilde_lmn)->sampleUnits, &(*hctilde_lmn)->sampleUnits, &lalSecondUnit);

        jStart   = (UINT4) ceil(fStart / deltaF);
        f        = jStart*deltaF;
        h_f      = 0.0;
        hconj_mf = 0.0;

        for ( UINT4 j=jStart ; j<jMax ; j++ ) {
        h_f      =      Prefactor/(Omega_lmn-LAL_TWOPI*f)*cexp(I*Omega_lmn*10.0*Mf_sec-I*LAL_TWOPI*f*10.0*Mf_sec);
        hconj_mf = conj(Prefactor/(Omega_lmn+LAL_TWOPI*f)*cexp(I*Omega_lmn*10.0*Mf_sec+I*LAL_TWOPI*f*10.0*Mf_sec));

        (*hptilde_lmn)->data->data[j] = 0.5 * (h_f + hconj_mf);
        (*hctilde_lmn)->data->data[j] = 0.5 * I * (h_f - hconj_mf);

        h_f      = 0.0;
        hconj_mf = 0.0;
        f+=deltaF;

        }

  return XLAL_SUCCESS;
}


int XLALSimRingdownMMRDNS_time(
        REAL8TimeSeries **hplus,                     /**< plus-polarization waveform [returned] */
        REAL8TimeSeries **hcross,                    /**< cross-polarization waveform [returned] */
        const LIGOTimeGPS *T0,                       /**< start time of ringdown => NEEDS TO BE CHECKED! */
        REAL8 deltaT,                                /**< sampling interval (s) */
        REAL8 Mf,                                    /**< Final BH Mass (kg) */
        REAL8 jf,                                    /**< Final BH dimensionaless spin */
        REAL8 eta,                                   /**< Symmetric mass ratio of two companions */
        REAL8 iota,                                  /**< inclination angle (in rad) */
        REAL8 phi_offset,                            /**< intrinsic phase offset */
        REAL8 r,                                     /**< distance of source (m) */
        LALSimInspiralTestGRParam *nonGRparams       /**< testing GR parameters */
    ){
        /* Perform some initial checks */
        if (Mf <= 0) XLAL_ERROR(XLAL_EDOM);
        if (jf >= 1) XLAL_ERROR(XLAL_EDOM);
        if (eta > 0.25 || eta <= 0) XLAL_ERROR(XLAL_EDOM);
        if (r <= 0) XLAL_ERROR(XLAL_EDOM);

        /* Declarations */
        REAL8 dfreq220=0., dfreq221=0., dfreq330=0., dfreq331=0., dfreq440=0., dfreq550=0., dfreq210=0., dfreq320=0., dfreq430=0.;
        REAL8 dtau220=0.,  dtau221=0.,  dtau330=0.,  dtau331=0.,  dtau440=0.,  dtau550=0.,  dtau210=0.,  dtau320=0.,  dtau430=0.;
        COMPLEX16TimeSeries *h220, *h221, *h330, *h331, *h440, *h550, *h210, *h320, *h430;

        /* Get nonGRparams */
        char *nonGRParamName = malloc(512*sizeof(char));

        sprintf(nonGRParamName,"dfreq220") ;
        if (XLALSimInspiralTestGRParamExists(nonGRparams,nonGRParamName) )
          dfreq220 = XLALSimInspiralGetTestGRParam(nonGRparams,nonGRParamName) ;
        sprintf(nonGRParamName,"dtau220") ;
        if (XLALSimInspiralTestGRParamExists(nonGRparams,nonGRParamName) )
          dtau220 = XLALSimInspiralGetTestGRParam(nonGRparams,nonGRParamName) ;
        sprintf(nonGRParamName,"dfreq221") ;
        if (XLALSimInspiralTestGRParamExists(nonGRparams,nonGRParamName) )
          dfreq221 = XLALSimInspiralGetTestGRParam(nonGRparams,nonGRParamName) ;
        sprintf(nonGRParamName,"dtau221") ;
        if (XLALSimInspiralTestGRParamExists(nonGRparams,nonGRParamName) )
          dtau221 = XLALSimInspiralGetTestGRParam(nonGRparams,nonGRParamName) ;
        sprintf(nonGRParamName,"dfreq330") ;
        if (XLALSimInspiralTestGRParamExists(nonGRparams,nonGRParamName) )
          dfreq330 = XLALSimInspiralGetTestGRParam(nonGRparams,nonGRParamName) ;
        sprintf(nonGRParamName,"dtau330") ;
        if (XLALSimInspiralTestGRParamExists(nonGRparams,nonGRParamName) )
          dtau330 = XLALSimInspiralGetTestGRParam(nonGRparams,nonGRParamName) ;
        sprintf(nonGRParamName,"dfreq331") ;
        if (XLALSimInspiralTestGRParamExists(nonGRparams,nonGRParamName) )
          dfreq331 = XLALSimInspiralGetTestGRParam(nonGRparams,nonGRParamName) ;
        sprintf(nonGRParamName,"dtau331") ;
        if (XLALSimInspiralTestGRParamExists(nonGRparams,nonGRParamName) )
          dtau331 = XLALSimInspiralGetTestGRParam(nonGRparams,nonGRParamName) ;
        sprintf(nonGRParamName,"dfreq440") ;
        if (XLALSimInspiralTestGRParamExists(nonGRparams,nonGRParamName) )
          dfreq440 = XLALSimInspiralGetTestGRParam(nonGRparams,nonGRParamName) ;
        sprintf(nonGRParamName,"dtau440") ;
        if (XLALSimInspiralTestGRParamExists(nonGRparams,nonGRParamName) )
          dtau440 = XLALSimInspiralGetTestGRParam(nonGRparams,nonGRParamName) ;
        sprintf(nonGRParamName,"dfreq550") ;
        if (XLALSimInspiralTestGRParamExists(nonGRparams,nonGRParamName) )
          dfreq550 = XLALSimInspiralGetTestGRParam(nonGRparams,nonGRParamName) ;
        sprintf(nonGRParamName,"dtau550") ;
        if (XLALSimInspiralTestGRParamExists(nonGRparams,nonGRParamName) )
          dtau550 = XLALSimInspiralGetTestGRParam(nonGRparams,nonGRParamName) ;
        sprintf(nonGRParamName,"dfreq210") ;
        if (XLALSimInspiralTestGRParamExists(nonGRparams,nonGRParamName) )
          dfreq210 = XLALSimInspiralGetTestGRParam(nonGRparams,nonGRParamName) ;
        sprintf(nonGRParamName,"dtau210") ;
        if (XLALSimInspiralTestGRParamExists(nonGRparams,nonGRParamName) )
          dtau210 = XLALSimInspiralGetTestGRParam(nonGRparams,nonGRParamName) ;
        sprintf(nonGRParamName,"dfreq320") ;
        if (XLALSimInspiralTestGRParamExists(nonGRparams,nonGRParamName) )
          dfreq320 = XLALSimInspiralGetTestGRParam(nonGRparams,nonGRParamName) ;
        sprintf(nonGRParamName,"dtau320") ;
        if (XLALSimInspiralTestGRParamExists(nonGRparams,nonGRParamName) )
          dtau320 = XLALSimInspiralGetTestGRParam(nonGRparams,nonGRParamName) ;
        sprintf(nonGRParamName,"dfreq430") ;
        if (XLALSimInspiralTestGRParamExists(nonGRparams,nonGRParamName) )
          dfreq430 = XLALSimInspiralGetTestGRParam(nonGRparams,nonGRParamName) ;
        sprintf(nonGRParamName,"dtau430") ;
        if (XLALSimInspiralTestGRParamExists(nonGRparams,nonGRParamName) )
          dtau430 = XLALSimInspiralGetTestGRParam(nonGRparams,nonGRParamName) ;

        /* time runs from Mstart to Mend after merger */
        REAL8 Tstart = 10.0*Mf*LAL_MTSUN_SI/LAL_MSUN_SI;
        REAL8 Tend = 60.0*Mf*LAL_MTSUN_SI/LAL_MSUN_SI;
        UINT4 Nsamples = ceil((Tend-Tstart)/deltaT);

        /* Compute the modes seperately */
        XLALSimRingdownGenerateSingleModeMMRDNS_time(&h220, T0, deltaT, Mf, jf, eta, iota, phi_offset, 2, 2, 0, r, dfreq220, dtau220,Nsamples,Tstart);
        XLALSimRingdownGenerateSingleModeMMRDNS_time(&h221, T0, deltaT, Mf, jf, eta, iota, phi_offset, 2, 2, 1, r, dfreq221, dtau221,Nsamples,Tstart);
        XLALSimRingdownGenerateSingleModeMMRDNS_time(&h330, T0, deltaT, Mf, jf, eta, iota, phi_offset, 3, 3, 0, r, dfreq330, dtau330,Nsamples,Tstart);
        XLALSimRingdownGenerateSingleModeMMRDNS_time(&h331, T0, deltaT, Mf, jf, eta, iota, phi_offset, 3, 3, 1, r, dfreq331, dtau331,Nsamples,Tstart);
        XLALSimRingdownGenerateSingleModeMMRDNS_time(&h440, T0, deltaT, Mf, jf, eta, iota, phi_offset, 4, 4, 0, r, dfreq440, dtau440,Nsamples,Tstart);
        XLALSimRingdownGenerateSingleModeMMRDNS_time(&h550, T0, deltaT, Mf, jf, eta, iota, phi_offset, 5, 5, 0, r, dfreq550, dtau550,Nsamples,Tstart);
        XLALSimRingdownGenerateSingleModeMMRDNS_time(&h210, T0, deltaT, Mf, jf, eta, iota, phi_offset, 2, 1, 0, r, dfreq210, dtau210,Nsamples,Tstart);
        XLALSimRingdownGenerateSingleModeMMRDNS_time(&h320, T0, deltaT, Mf, jf, eta, iota, phi_offset, 3, 2, 0, r, dfreq320, dtau320,Nsamples,Tstart);
        XLALSimRingdownGenerateSingleModeMMRDNS_time(&h430, T0, deltaT, Mf, jf, eta, iota, phi_offset, 4, 3, 0, r, dfreq430, dtau430,Nsamples,Tstart);

        /* Add the modes to get the final waveform and  get cross and plus polarization */
        *hplus = XLALCreateREAL8TimeSeries( "hplus: TD waveform", T0, 0.0, deltaT, &lalStrainUnit, Nsamples);
        if (!(*hplus)) XLAL_ERROR(XLAL_EFUNC);
        memset((*hplus)->data->data, 0, Nsamples * sizeof(REAL8));
        *hcross = XLALCreateREAL8TimeSeries( "hcross: TD waveform", T0, 0.0, deltaT, &lalStrainUnit, Nsamples); 
        if (!(*hcross)) XLAL_ERROR(XLAL_EFUNC);
        memset((*hcross)->data->data, 0, Nsamples * sizeof(REAL8));


        /* prepare window */
        REAL8Window *window_rd;
        REAL8 start = 0; //TODO: adjust window
        window_rd = XLALCreatePlanckREAL8Window(Nsamples,start,Nsamples*deltaT-2.0,1.0/deltaT);


        COMPLEX16 h_val = 0.0;
        for ( UINT4 i=0 ; i<Nsamples ; i++ ) {
          h_val = h220->data->data[i];
          h_val += h221->data->data[i];
          h_val += h330->data->data[i];
          h_val += h331->data->data[i];
          h_val += h440->data->data[i];
          h_val += h550->data->data[i];
          h_val += h210->data->data[i];
          h_val += h320->data->data[i];
          h_val += h430->data->data[i];

          (*hplus)->data->data[i] = creal(h_val)*(window_rd->data->data)[i];
          (*hcross)->data->data[i] = -cimag(h_val)*(window_rd->data->data)[i];
          h_val = 0.0;
        }


        /* Destroy the COMPLEX16TimeSeries object */
        if (h220) XLALDestroyCOMPLEX16TimeSeries(h220);
        if (h221) XLALDestroyCOMPLEX16TimeSeries(h221);
        if (h330) XLALDestroyCOMPLEX16TimeSeries(h330);
        if (h331) XLALDestroyCOMPLEX16TimeSeries(h331);
        if (h440) XLALDestroyCOMPLEX16TimeSeries(h440);
        if (h550) XLALDestroyCOMPLEX16TimeSeries(h550);
        if (h210) XLALDestroyCOMPLEX16TimeSeries(h210);
        if (h320) XLALDestroyCOMPLEX16TimeSeries(h320);
        if (h430) XLALDestroyCOMPLEX16TimeSeries(h430);

        return 0;
        
}


int XLALSimRingdownGenerateSingleModeMMRDNS_time(
        COMPLEX16TimeSeries **htilde_lmn,            /**< OUTPUT TD waveform mode lmn */
        const LIGOTimeGPS *T0,                       /**< start time of ringdown */
        REAL8 deltaT,                                /**< sampling interval (s) */
        REAL8 Mf,                                    /**< Final BH Mass (kg) */
        REAL8 jf,                                    /**< Final BH dimensionaless spin */
        REAL8 eta,                                   /**< Symmetric mass ratio of two companions */
        REAL8 iota,                                  /**< inclination angle (in rad) */
        REAL8 phi_offset,                            /**< intrinsic phase offset (in rad) */
        UINT4 l,                                     /**< Polar eigenvalue */
        UINT4 m,                                     /**< Azimuthal eigenvalue */
        UINT4 n,                                     /**< Overtone Number */
        REAL8 r,                                     /**< distance of source (m) */
        REAL8 dfreq,                                 /**< relative shift in the real frequency parameter */
        REAL8 dtau,                                  /**< relative shift in the damping time parameter */
        UINT4 Nsamples,                              /**< waveform length */
        REAL8 Tstart                                 /**< starting time of waveform (10M at zero) */
        ){

        /* Declarations */
        COMPLEX16 h_lmn = 0;
        COMPLEX16 prefactor = 0;
        REAL8 kappa   = XLALKAPPA( jf, l, m );
        REAL8 Mf_sec  = Mf*LAL_MTSUN_SI/LAL_MSUN_SI;
        REAL8 r_sec   = r/LAL_C_SI;
        REAL8 t = 0;
        COMPLEX16 A_lmn = 0;
        COMPLEX16 S_lmn = 0;
        COMPLEX16 Omega_lmn = 0;

        /* Mode Component Calculation*/
        Omega_lmn = XLALcomplexOmega(kappa, l, m, n)/Mf_sec;
        A_lmn = XLALMMRDNSAmplitudeOverOmegaSquared(eta, l, m, n);
        Omega_lmn = creal(Omega_lmn)*(1.+dfreq) + I * cimag(Omega_lmn)/(1.+dtau);
        S_lmn = XLALSpinWeightedSpheroidalHarmonic(jf, l, m, n, iota, 0.0);
        prefactor = -A_lmn*S_lmn*Mf_sec/(r_sec)*cexp(I*phi_offset);

        /* allocate htilde_lmn */
        *htilde_lmn = XLALCreateCOMPLEX16TimeSeries("htilde_lmn: TD waveform", T0, 0.0, deltaT, &lalStrainUnit, Nsamples);
        if (!(*htilde_lmn)) XLAL_ERROR(XLAL_EFUNC);
        memset((*htilde_lmn)->data->data, 0, Nsamples * sizeof(COMPLEX16));
        
        /* 10M is zero in model -> shift 10M to 0 on time axis*/
        Tstart -= 10.0*Mf_sec;

        /* fill waveform */
        for ( UINT4 i=0 ; i<Nsamples ; i++ ) {
            t = Tstart+i*deltaT;
            h_lmn = prefactor*cexp(I*Omega_lmn*t);
            (*htilde_lmn)->data->data[i] = h_lmn;
            h_lmn = 0.0;
            t = 0.0;
        }

  return XLAL_SUCCESS;

}
