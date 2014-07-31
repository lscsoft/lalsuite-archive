/*
*  Copyright (C) 2014 Chinmay Kalaghatgi, P. Ajith
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

#include <lal/AVFactories.h>
#include <stdlib.h>
#include <math.h>
#include <lal/LALConstants.h>
#include <lal/LALDatatypes.h>
#include <lal/LALSimInspiral.h>
#include <lal/Units.h>
#include <lal/XLALError.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_permutation.h>

/* Computed in Mathematica with N[\[Pi]^2, 35] */
#define Pi_p2 9.8696044010893586188344909998761511
#define Pi_p3 31.006276680299820175476315067101395
#define Pi_p4 97.409091034002437236440332688705111
#define Pi_p5 306.01968478528145326274131004343561
#define Pi_p4by3 4.6011511144704899609836928301589280
#define Pi_p5by3 6.7388085956981412852178979102191822
#define Pi_p10by3 45.411541289455155012188835024560824
#define Pi_p7by3 14.454942539276981063196177096150014
#define Pi_p1by3 1.4645918875615232630201425272637904
#define Pi_p2by3 2.1450293971110256000774441009412356
#define Pi_11by3 66.509374974201175524807943232081902
#define tenbyPi_p2by3 2.1638812222639821850478371484217674
#define twoPi_p2by3 3.4050219214767547368607032470366207

/**
 * Function to compute the metric elements using waveform derivatives
 */
static REAL8 MetricCoeffs(REAL8Vector *A, REAL8Vector *dPsii, REAL8Vector *dPsij,
        REAL8Vector *dAi, REAL8Vector*dAj, REAL8Vector *Sh, REAL8 hSqr, REAL8 df) {
    size_t k = A->length;
    REAL8 gij = 0.;
    for (;k--;) {
        gij += (A->data[k]*A->data[k]*dPsii->data[k]*dPsij->data[k]
                + dAi->data[k]*dAj->data[k])/Sh->data[k];
    }
    return 4.*df*gij/(2.*hSqr);
}

/**
 * Frequency domain amplitude of the TaylorF2 Reduced Spin waveforms
 */
static REAL8 XLALSimIMRPhenomBAofF(
    const REAL8 mc_msun,   /**< chirp mass (solar masses) */
    const REAL8 eta,   /**< symmetric mass ratio  */
    const REAL8 chi,   /**< reduced-spin parameter */
    const REAL8 f   /**< Fourier frequency (Hz) */
) {
    const REAL8 mc = mc_msun * LAL_MTSUN_SI;
    const REAL8 etap2 = eta*eta;
    const REAL8 etap3 = etap2*eta;
    const REAL8 chip2 = chi*chi;
    const REAL8 mcp2  = mc*mc;

    REAL8 eta_three_fifths = pow(eta, 0.6);
    REAL8 f_mc_eta = (f*mc)/eta_three_fifths;
    REAL8 f_mc_eta_one_third = cbrt(f_mc_eta);

    return (sqrt(0.8333333333333334)*pow(mc/eta_three_fifths,0.8333333333333334)*
        sqrt(eta)*(1 + (Pi_p4by3*
         f_mc_eta_one_third * f_mc_eta_one_third *
         (743 + 924*eta))/672. -
        (Pi_p10by3 *
         f_mc_eta * f_mc_eta * f_mc_eta_one_third*
         (5111593 + 8088752*eta + 151088*etap2))/2.709504e6 +
        (f*mc*LAL_PI*(-2*LAL_PI + (113*chi)/24.))/eta_three_fifths +
        (Pi_p5by3 *
         f_mc_eta * f_mc_eta / f_mc_eta_one_third *
         (12*LAL_PI*(-4757 + 4788*eta) -
           (113*(502429 - 591368*eta + 1680*etap2)*chi)/
            (-113 + 76*eta)))/16128. +
        (Pi_p5 * Pi_p1by3 *
         f_mc_eta * f_mc_eta_one_third *
         (7266251 + 9532152*eta + 9730224*etap2 +
           (3243530304*(-81 + 4*eta)*chip2)/
            pow(113 - 76*eta,2)))/8.128512e6 +
        (f*f*mcp2*Pi_p2*
         (-58.601030974347324 + (10*Pi_p2)/3. +
           (3526813753*eta)/2.7869184e7 -
           (451*Pi_p2*eta)/96. -
           (1041557*etap2)/258048. +
           (67999*etap3)/82944. +
           (856*(3*LAL_GAMMA + log((64*f*mc*LAL_PI)/eta_three_fifths)))/315.))/(eta_three_fifths * eta_three_fifths)))/
        (2.*pow(f,1.1666666666666667)*Pi_p4by3);
}

/**
 * Derivative of the amplitude with respect to \f$\chi\f$ (reduced-spin parameter)
 */
static REAL8 XLALSimIMRPhenomBDerivAChi(
        const REAL8 mc_msun,   /**< chirp mass (solar masses) */
        const REAL8 eta,  /**< symmetric mass ratio  */
        const REAL8 chi,  /**< reduced-spin parameter */
        const REAL8 f     /**< Fourier frequency (Hz) */
) {
    const REAL8 mc = mc_msun * LAL_MTSUN_SI;
    const REAL8 etap2 = eta*eta;
    REAL8 etap35 = pow(eta, 0.6);

    return (113*sqrt(0.8333333333333334)*Pi_p1by3*
        pow(mc/etap35,0.8333333333333334)*sqrt(eta)*
        ((672*f*mc)/etap35 -
        (Pi_p1by3*Pi_p1by3*
         pow((f*mc)/etap35,1.6666666666666667)*
         (502429 - 591368*eta + 1680*etap2))/
        (-113 + 76*eta) + (113904*Pi_p1by3*
         pow((f*mc)/etap35,1.3333333333333333)*
         (-81 + 4*eta)*chi)/pow(113 - 76*eta,2)))/
    (32256.*pow(f,1.1666666666666667));
}

/**
 * Derivative of the amplitude with respect to \f$\eta\f$ (symm. mass ratio)
 */
static REAL8 XLALSimIMRPhenomBDerivAEta(
    const REAL8 mc_msun,   /**< chirp mass (M_sun) */
    const REAL8 eta,  /**< symmetric mass ratio  */
    const REAL8 chi,  /**< reduced-spin parameter */
    const REAL8 f     /**< Fourier frequency (Hz) */
) {
    const REAL8 mc = mc_msun * LAL_MTSUN_SI;
    const REAL8 etap2 = eta*eta;
    const REAL8 etap3 = etap2*eta;
    const REAL8 etap4 = etap3*eta;
    const REAL8 etap5 = etap4*eta;
    const REAL8 chip2 = chi*chi;
    const REAL8 mcp2  = mc*mc;
    const REAL8 eta_fac = (-113 + 76 * eta) * (-113 + 76 * eta) * (-113 + 76 * eta);
    REAL8 eta_three_fifths = pow(eta, 0.6);
    REAL8 f_mc_eta_one_third = cbrt((f*mc)/eta_three_fifths);

    return (pow(mc/eta_three_fifths,0.8333333333333334)*
        (2235340800*f_mc_eta_one_third * f_mc_eta_one_third*
        eta_three_fifths*eta_three_fifths*eta_fac*(-743 + 1386*eta) +
        f*f*mcp2*Pi_p4by3*
        eta_fac *
        (257959397806029 - 36738273116160*LAL_GAMMA -
         95047630643350*eta - 12126223216800*etap2 +
         5541700903200*etap3 +
         7823692800*Pi_p2*(-1920 + 451*eta) -
         1940400*Pi_p4by3*
          f_mc_eta_one_third*
          (-5111593 - 2311072*eta + 64752*etap2)) +
        369600*f*mc*Pi_p1by3*eta_three_fifths*
        (12192768*LAL_PI*eta_fac +
         35962920*Pi_p5by3*
          f_mc_eta_one_third * f_mc_eta_one_third *
          eta_fac -
         28703808*eta_fac*chi -
         71190*Pi_p2by3*
          f_mc_eta_one_third * f_mc_eta_one_third *
          (-6415515901 + 12944580756*eta -
            10861276272*etap2 + 3401313728*etap3)*
          chi + Pi_p1by3*
          f_mc_eta_one_third *
          (34636018030144*etap3 -
            27532505500416*etap4 +
            6407002215936*etap5 -
            1442897*(-7266251 + 20575296*chip2) -
            21696*etap2*(-4887203 + 102257316*chip2) +
            76614*eta*(-320998087 + 907387488*chip2))) -
        12246091038720*f*f*mcp2*Pi_p4by3*
        eta_fac*log((64*f*mc*LAL_PI)/eta_three_fifths)))/
        (1.5021490176e12*sqrt(30)*pow(f,1.1666666666666667)*pow(eta,1.7)*
        eta_fac);
}

/**
 * Derivative of the amplitude with respect to the chirp mass
 */
static REAL8 XLALSimIMRPhenomBDerivAMChirp(
    const REAL8 mc_msun,   /**< chirp mass (M_sun) */
    const REAL8 eta,  /**< symmetric mass ratio  */
    const REAL8 chi,  /**< reduced-spin parameter */
    const REAL8 f     /**< Fourier frequency (Hz) */
) {
    const REAL8 mc = mc_msun * LAL_MTSUN_SI;
    const REAL8 etap2 = eta*eta;
    const REAL8 etap3 = etap2*eta;
    const REAL8 etap4 = etap3*eta;
    const REAL8 chip2 = chi*chi;
    const REAL8 mcp2  = mc*mc;
    REAL8 eta_three_fifths = pow(eta, 0.6);

    return (sqrt(0.8333333333333334)*(5 -
      (17*f*f*mcp2*Pi_p2*Pi_p2*(-320 + 451*eta))/
       (96.*eta_three_fifths*eta_three_fifths) +
      (3*Pi_p2by3*
         pow((f*mc)/eta_three_fifths,0.6666666666666666)*
         (743 + 924*eta))/224. +
      (5*Pi_p3/Pi_p1by3*
         pow((f*mc)/eta_three_fifths,1.6666666666666667)*
         (-4757 + 4788*eta))/448. -
      (19*pow(LAL_PI,3.3333333333333335)*
         pow((f*mc)/eta_three_fifths,2.3333333333333335)*
         (5111593 + 8088752*eta + 151088*etap2))/2.709504e6 +
      (f*mc*Pi_p2*(-33047278387200*eta_three_fifths +
           f*mc*(-1471974996766431 + 208183547658240*LAL_GAMMA +
              3231619441873900*eta - 103072897342800*etap2 +
              20935314523200*etap3)))/
       (1.5021490176e12*eta_three_fifths*eta_three_fifths) +
      (1243*f*mc*LAL_PI*chi)/(24.*eta_three_fifths) -
      (565*Pi_p5by3*
         pow((f*mc)/eta_three_fifths,1.6666666666666667)*
         (502429 - 591368*eta + 1680*etap2)*chi)/
       (5376.*(-113 + 76*eta)) +
      (13*Pi_p4by3*
         pow((f*mc)/eta_three_fifths,1.3333333333333333)*
         (2490853280*etap2 - 112068617472*etap3 +
           56201773824*etap4 +
           1808*eta*(-1708561 + 7175952*chip2) -
           12769*(-7266251 + 20575296*chip2)))/
       (8.128512e6*pow(113 - 76*eta,2)) +
      (14552*f*f*mcp2*Pi_p2*
         log((64*f*mc*LAL_PI)/eta_three_fifths))/(315.*eta_three_fifths*eta_three_fifths)))/
    (12.*pow(f,1.1666666666666667)*Pi_p2by3*
    pow(mc/eta_three_fifths,0.16666666666666666)*pow(eta,0.1));
}

/**
 * Derivative of the phasae with respect to \f$\chi\f$ (reduced spin parameter)
 */
static REAL8 XLALSimIMRPhenomBDerivPsiChi(
    const REAL8 mc_msun,   /**< chirp mass (M_sun) */
    const REAL8 eta,  /**< symmetric mass ratio  */
    const REAL8 chi,  /**< reduced-spin parameter */
    const REAL8 f     /**< Fourier frequency (Hz) */
) {
    const REAL8 mc = mc_msun * LAL_MTSUN_SI;
    const REAL8 etap2 = eta*eta;
    REAL8 eta_three_fifths = pow(eta, 0.6);

    return (113*((756*f*mc)/eta_three_fifths +
      (320355*Pi_p1by3*
         pow((f*mc)/eta_three_fifths,1.3333333333333333)*
         (-81 + 4*eta)*chi)/pow(113 - 76*eta,2) -
      (5*Pi_p2by3*
         pow((f*mc)/eta_three_fifths,1.6666666666666667)*
         (-146597 + 135856*eta + 17136*etap2)*
         (1 + log(f / 40.)))/
       (-113 + 76*eta)))/
    (96768.*Pi_p2by3*
    pow((f*mc)/eta_three_fifths,1.6666666666666667)*eta);
}

/**
 * Derivative of the phasae with respect to \f$\eta\f$ (symmetric mass ratio)
 */
static REAL8 XLALSimIMRPhenomBDerivPsiEta(
    const REAL8 mc_msun,   /**< chirp mass (M_sun) */
    const REAL8 eta,  /**< symmetric mass ratio  */
    const REAL8 chi,  /**< reduced-spin parameter */
    const REAL8 f     /**< Fourier frequency (Hz) */
) {
    const REAL8 mc = mc_msun * LAL_MTSUN_SI;
    const REAL8 etap2 = eta*eta;
    const REAL8 etap3 = etap2*eta;
    const REAL8 etap4 = etap3*eta;
    const REAL8 etap5 = etap4*eta;
    const REAL8 chip2 = chi*chi;
    const REAL8 mcp2  = mc*mc;
    REAL8 eta_three_fifths = pow(eta, 0.6);
    const REAL8 eta_fac = -113 + 76 * eta;
    const REAL8 fmc_fac = cbrt(f * mc / eta_three_fifths);

    return (-77616000*f*mc*LAL_PI*eta_fac*
     (23187*LAL_PI*eta_fac*eta_fac -
       113*(16565461 - 22282744*eta + 12261424*etap2)*chi)*
     log(fmc_fac/
       (2.*cbrt(5)*
         cbrt(mc/eta_three_fifths))) +
    (31046400*fmc_fac*fmc_fac*
        eta_three_fifths * eta_three_fifths *eta_fac*eta_fac*eta_fac*(-743 + 1386*eta) -
       f*f*mcp2*LAL_PI * Pi_p1by3 *
        eta_fac * eta_fac * eta_fac *
        (33984313019673 - 4592284139520*LAL_GAMMA -
          12118079538950*eta - 413215941600*etap2 +
          2083465692000*etap3 +
          977961600*Pi_p2*(-3072 + 451*eta) +
          323400*LAL_PI*Pi_p1by3*
           fmc_fac*
           (15419335 + 3633744*eta + 2132496*etap2)) -
       18480*f*mc*Pi_p1by3*eta_three_fifths*
        (-6096384*LAL_PI*eta_fac * eta_fac * eta_fac +
          32461800*Pi_p2/Pi_p1by3*
           fmc_fac*fmc_fac*
           eta_fac * eta_fac * eta_fac +
          14351904*eta_fac * eta_fac * eta_fac*chi -
          158200*LAL_PI/Pi_p1by3*
           fmc_fac*fmc_fac*
           (-1871897093 + 3776925108*eta -
             3079029456*etap2 + 931868224*etap3)*
           chi - 5*Pi_p1by3*
           fmc_fac*
           (14990425815136*etap3 -
             12186233587584*etap4 +
             2866657264128*etap5 -
             1442897*(-3058673 + 5143824*chip2) +
             612912*eta*(-17749451 + 28355859*chip2) -
             2712*etap2*(-202619251 + 204514632*chip2)))
        + 1530761379840*f*f*mcp2*LAL_PI*Pi_p1by3*
        eta_fac * eta_fac * eta_fac*log((64*f*mc*LAL_PI)/eta_three_fifths))/
     (fmc_fac * fmc_fac *eta_three_fifths))/
    (5.007163392e11*f*mc*LAL_PI*etap2*eta_fac * eta_fac * eta_fac);
}

/**
 * Derivative of the phase with respect to the chirp mass
 */
static REAL8 XLALSimIMRPhenomBDerivPsiMChirp(
    const REAL8 mc_msun,   /**< chirp mass (M_sun) */
    const REAL8 eta,  /**< symmetric mass ratio  */
    const REAL8 chi,  /**< reduced-spin parameter */
    const REAL8 f     /**< Fourier frequency (Hz) */
) {
    const REAL8 mc = mc_msun * LAL_MTSUN_SI;
    const REAL8 etap2 = eta*eta;
    const REAL8 etap3 = etap2*eta;
    const REAL8 etap4 = etap3*eta;
    const REAL8 chip2 = chi*chi;
    const REAL8 mcp2  = mc*mc;
    const REAL8 mcp3  = mcp2*mc;
    REAL8 eta_three_fifths = pow(eta, 0.6);

    return (cbrt((f*mc)/eta_three_fifths)*
    (-93139200*pow(113 - 76*eta,2)*eta_three_fifths * eta_three_fifths *
       (252 + Pi_p2by3*
          pow((f*mc)/eta_three_fifths,0.6666666666666666)*
          (743 + 924*eta)) +
      f*f*mcp2*LAL_PI*LAL_PI*pow(113 - 76*eta,2)*
       (10052469856691 - 1530761379840*LAL_GAMMA -
         24236159077900*eta + 206607970800*etap2 -
         462992376000*etap3 +
         1955923200*Pi_p2*(-512 + 451*eta) -
         184800*Pi_p4by3*
          cbrt((f*mc)/eta_three_fifths)*
          (-15419335 - 12718104*eta + 4975824*etap2)) +
      9240*f*mc*LAL_PI*eta_three_fifths*
       (16257024*LAL_PI*pow(113 - 76*eta,2) -
         38271744*pow(113 - 76*eta,2)*chi -
         5*cbrt(LAL_PI)*
          cbrt((f*mc)/eta_three_fifths)*
          (-20737091296*etap2 - 43167841920*etap3 +
            25146116352*etap4 +
            904*eta*(19183315 + 3587976*chip2) -
            12769*(-3058673 + 5143824*chip2))) -
      510253793280*f*f*mcp2*LAL_PI*LAL_PI*
       pow(113 - 76*eta,2)*log((64*f*mc*LAL_PI)/eta_three_fifths)))/
    (6.0085960704e11*f*f*mcp3*pow(LAL_PI,1.6666666666666667)*
    pow(113 - 76*eta,2)*eta);
}


/**
 * Compute the template-space metric of "reduced-spin" PN templates in
 * Mchirp-eta-chi parameter space.
 */
int XLALSimIMRPhenomBMetricMChirpEtaChi(
    REAL8 *gamma00,  /**< template metric coeff. 00 in mChirp-eta-chi */
    REAL8 *gamma01,  /**< template metric coeff. 01/10 in mChirp-eta-chi */
    REAL8 *gamma02,  /**< template metric coeff. 02/20 in mChirp-eta-chi */
    REAL8 *gamma11,  /**< template metric coeff. 11 in mChirp-eta-chi */
    REAL8 *gamma12,  /**< template metric coeff. 12/21 in mChirp-eta-chi */
    REAL8 *gamma22,  /**< template metric coeff. 22 in mChirp-eta-chi */
    const REAL8 mc,     /**< chirp mass (in solar mass) */
    const REAL8 eta,    /**< symmetric mass ratio */
    const REAL8 chi,    /**< reduced-spin parameter */
    const REAL8 fLow,   /**< low-frequency cutoff (Hz) */
    const REAL8FrequencySeries *Sh  /**< PSD in strain per root Hertz */
) {
    REAL8Vector *A=NULL, *dAMc=NULL, *dAEta=NULL;
    REAL8Vector *dAChi=NULL, *dAT0=NULL, *dAPhi0=NULL, *dPsiMc=NULL;
    REAL8Vector *dPsiEta=NULL, *dPsiChi=NULL, *dPsiT0=NULL, *dPsiPhi0=NULL;

    const REAL8 df = Sh->deltaF;
    REAL8 hSqr = 0.;
    int s = 0;

    /* compute the Schwarzschild ISCO frequency */
    REAL8 fCut = 1. / (LAL_PI*mc*pow(eta, -0.6) * sqrt(6 * 6 * 6) * LAL_MTSUN_SI);

    /* create a view of the PSD between fLow and fCut */
    size_t nBins = (fCut - fLow) / df;
    size_t k = nBins;
    REAL8Vector Shdata = {nBins, Sh->data->data + (size_t) (fLow / df)}; /* copy the Vector, including its pointer to the actual data */
    /* drop low-frequency samples */
    Shdata.length = nBins;  /* drop high-frequency samples */

    /* allocate memory for various vectors */
    A = XLALCreateREAL8Vector(nBins);
    dAMc = XLALCreateREAL8Vector(nBins);
    dAEta = XLALCreateREAL8Vector(nBins);
    dAChi = XLALCreateREAL8Vector(nBins);
    dAT0 = XLALCreateREAL8Vector(nBins);
    dAPhi0 = XLALCreateREAL8Vector(nBins);
    dPsiMc = XLALCreateREAL8Vector(nBins);
    dPsiEta = XLALCreateREAL8Vector(nBins);
    dPsiChi = XLALCreateREAL8Vector(nBins);
    dPsiT0 = XLALCreateREAL8Vector(nBins);
    dPsiPhi0 = XLALCreateREAL8Vector(nBins);

    /* create a frequency vector from fLow to fCut with frequency resolution df */
    for (;k--;) {
        const REAL8 f = fLow + k * df;

        /* compute the amplitude of the frequency-domain waveform */
        A->data[k] = XLALSimIMRPhenomBAofF(mc, eta, chi, f);

        /* compute the square of the template norm */
        hSqr += A->data[k] * A->data[k] / Shdata.data[k];

        /* compute the waveform deratives with respect to the parameters */
        dAMc->data[k] = XLALSimIMRPhenomBDerivAMChirp(mc, eta, chi, f);
        dAEta->data[k] = XLALSimIMRPhenomBDerivAEta(mc, eta, chi, f);
        dAChi->data[k] = XLALSimIMRPhenomBDerivAChi(mc, eta, chi, f);
        dAT0->data[k] = 0.;
        dAPhi0->data[k] = 0.;
        dPsiMc->data[k] = XLALSimIMRPhenomBDerivPsiMChirp(mc, eta, chi, f);
        dPsiEta->data[k] = XLALSimIMRPhenomBDerivPsiEta(mc, eta, chi, f);
        dPsiChi->data[k] = XLALSimIMRPhenomBDerivPsiChi(mc, eta, chi, f);
        dPsiT0->data[k] = LAL_TWOPI * f;
        dPsiPhi0->data[k] = 1.;
    }
    hSqr *= 4 * df;

    /* allocate memory, and initialize the Fisher matrix */
    gsl_matrix * g = gsl_matrix_calloc (5, 5);

    /* compute the components of the Fisher matrix in coordinates mc, eta, chi, t0, phi0 */
    gsl_matrix_set (g, 0,0, MetricCoeffs(A, dPsiMc, dPsiMc, dAMc, dAMc, &Shdata, hSqr, df));
    gsl_matrix_set (g, 0,1, MetricCoeffs(A, dPsiMc, dPsiEta, dAMc, dAEta, &Shdata, hSqr, df));
    gsl_matrix_set (g, 0,2, MetricCoeffs(A, dPsiMc, dPsiChi, dAMc, dAChi, &Shdata, hSqr, df));
    gsl_matrix_set (g, 0,3, MetricCoeffs(A, dPsiMc, dPsiT0, dAMc, dAT0, &Shdata, hSqr, df));
    gsl_matrix_set (g, 0,4, MetricCoeffs(A, dPsiMc, dPsiPhi0, dAMc, dAPhi0, &Shdata, hSqr, df));

    gsl_matrix_set (g, 1,0, gsl_matrix_get(g, 0,1));
    gsl_matrix_set (g, 1,1, MetricCoeffs(A, dPsiEta, dPsiEta, dAEta, dAEta, &Shdata, hSqr, df));
    gsl_matrix_set (g, 1,2, MetricCoeffs(A, dPsiEta, dPsiChi, dAEta, dAChi, &Shdata, hSqr, df));
    gsl_matrix_set (g, 1,3, MetricCoeffs(A, dPsiEta, dPsiT0, dAEta, dAT0, &Shdata, hSqr, df));
    gsl_matrix_set (g, 1,4, MetricCoeffs(A, dPsiEta, dPsiPhi0, dAEta, dAPhi0, &Shdata, hSqr, df));

    gsl_matrix_set (g, 2,0, gsl_matrix_get(g, 0,2));
    gsl_matrix_set (g, 2,1, gsl_matrix_get(g, 1,2));
    gsl_matrix_set (g, 2,2, MetricCoeffs(A, dPsiChi, dPsiChi, dAChi, dAChi, &Shdata, hSqr, df));
    gsl_matrix_set (g, 2,3, MetricCoeffs(A, dPsiChi, dPsiT0, dAChi, dAT0, &Shdata, hSqr, df));
    gsl_matrix_set (g, 2,4, MetricCoeffs(A, dPsiChi, dPsiPhi0, dAChi, dAPhi0, &Shdata, hSqr, df));

    gsl_matrix_set (g, 3,0, gsl_matrix_get(g, 0,3));
    gsl_matrix_set (g, 3,1, gsl_matrix_get(g, 1,3));
    gsl_matrix_set (g, 3,2, gsl_matrix_get(g, 2,3));
    gsl_matrix_set (g, 3,3, MetricCoeffs(A, dPsiT0, dPsiT0, dAT0, dAT0, &Shdata, hSqr, df));
    gsl_matrix_set (g, 3,4, MetricCoeffs(A, dPsiT0, dPsiPhi0, dAT0, dAPhi0, &Shdata, hSqr, df));

    gsl_matrix_set (g, 4,0, gsl_matrix_get(g, 0,4));
    gsl_matrix_set (g, 4,1, gsl_matrix_get(g, 1,4));
    gsl_matrix_set (g, 4,2, gsl_matrix_get(g, 2,4));
    gsl_matrix_set (g, 4,3, gsl_matrix_get(g, 3,4));
    gsl_matrix_set (g, 4,4, MetricCoeffs(A, dPsiPhi0, dPsiPhi0, dAPhi0, dAPhi0, &Shdata, hSqr, df));

    /* free the memory */
    XLALDestroyREAL8Vector(A);
    XLALDestroyREAL8Vector(dAMc);
    XLALDestroyREAL8Vector(dAEta);
    XLALDestroyREAL8Vector(dAChi);
    XLALDestroyREAL8Vector(dAT0);
    XLALDestroyREAL8Vector(dAPhi0);
    XLALDestroyREAL8Vector(dPsiMc);
    XLALDestroyREAL8Vector(dPsiEta);
    XLALDestroyREAL8Vector(dPsiChi);
    XLALDestroyREAL8Vector(dPsiT0);
    XLALDestroyREAL8Vector(dPsiPhi0);

    {
    /* Form submatrices g1, g2, g3, g4, defined as:
     *              g = [ g1 g2
     *                    g4 g3 ]                           */
    gsl_matrix_view g1v = gsl_matrix_submatrix (g, 0, 0, 3, 3);
    gsl_matrix_view g2v = gsl_matrix_submatrix (g, 0, 3, 3, 2);
    gsl_matrix_view g3v = gsl_matrix_submatrix (g, 3, 3, 2, 2);
    gsl_matrix_view g4v = gsl_matrix_submatrix (g, 3, 0, 2, 3);

    gsl_matrix * g1 = gsl_matrix_calloc (3, 3);
    gsl_matrix * g2 = gsl_matrix_calloc (3, 2);
    gsl_matrix * g3 = gsl_matrix_calloc (2, 2);
    gsl_matrix * g4 = gsl_matrix_calloc (2, 3);
    gsl_matrix * g3invg4 = gsl_matrix_calloc (2, 3);
    gsl_matrix * g2g3invg4 = gsl_matrix_calloc (3, 3);

    /* Project out the t0 and phi0 dimensions: gamma =  g1 - g2 g3^{-1} g4 */
    gsl_matrix * g3inv = gsl_matrix_calloc (2, 2);
    gsl_permutation * p = gsl_permutation_calloc (2);

    gsl_matrix_memcpy (g1, &g1v.matrix);
    gsl_matrix_memcpy (g2, &g2v.matrix);
    gsl_matrix_memcpy (g3, &g3v.matrix);
    gsl_matrix_memcpy (g4, &g4v.matrix);
    gsl_matrix_free (g);

    gsl_linalg_LU_decomp (g3, p, &s);
    gsl_linalg_LU_invert (g3, p, g3inv);
    gsl_permutation_free (p);
    gsl_matrix_free (g3);

    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, g3inv, g4,  0.0, g3invg4);
    gsl_matrix_free (g4);
    gsl_matrix_free (g3inv);
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, g2, g3invg4,  0.0, g2g3invg4);
    gsl_matrix_free (g2);
    gsl_matrix_free (g3invg4);

    gsl_matrix_sub (g1, g2g3invg4);
    gsl_matrix_free (g2g3invg4);

    *gamma00 = gsl_matrix_get(g1, 0, 0);
    *gamma01 = gsl_matrix_get(g1, 0, 1);
    *gamma02 = gsl_matrix_get(g1, 0, 2);
    *gamma11 = gsl_matrix_get(g1, 1, 1);
    *gamma12 = gsl_matrix_get(g1, 1, 2);
    *gamma22 = gsl_matrix_get(g1, 2, 2);
    gsl_matrix_free (g1);
    }

    return XLAL_SUCCESS;
}

