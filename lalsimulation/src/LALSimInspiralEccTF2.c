/*
 *  Copyright (C) 2015 Maria Haney, A. Gopakumar
 *
 *  Copyright (C) 2007 Jolien Creighton, B.S. Sathyaprakash, Thomas Cokelaer
 *  Copyright (C) 2012 Leo Singer, Evan Ochsner, Les Wade, Alex Nitz
 *  Assembled from code found in:
 *    - LALInspiralStationaryPhaseApproximation2.c
 *    - LALInspiralChooseModel.c
 *    - LALInspiralSetup.c
 *    - LALSimInspiralTaylorF2ReducedSpin.c
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

#include <stdlib.h>
#include <math.h>
#include <lal/Date.h>
#include <lal/FrequencySeries.h>
#include <lal/LALConstants.h>
#include <lal/Sequence.h>
#include <lal/LALDatatypes.h>
#include <lal/LALSimInspiral.h>
#include <lal/Units.h>
#include <lal/XLALError.h>
#include "LALSimInspiralPNCoefficients.c"

#ifndef _OPENMP
#define omp ignore
#endif


#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif



/*
 * Leading-order eccentric corrections to TaylorF2 phasing up to 3PN phase order
 * 3PN extension of Eqs. A3 and A4 in LIGO-P1500148-v1
 */
static void UNUSED
XLALSimInspiralPNPhasing_eccF2(
        PNPhasingSeries *pfv19by3, /**< \todo UNDOCUMENTED */
        PNPhasingSeries *pfv25by3, /**< \todo UNDOCUMENTED */
        PNPhasingSeries *pfv28by3, /**< \todo UNDOCUMENTED */
        PNPhasingSeries *pfv31by3, /**< \todo UNDOCUMENTED */
        PNPhasingSeries *pfv34by3, /**< \todo UNDOCUMENTED */
        PNPhasingSeries *pfv37by3, /**< \todo UNDOCUMENTED */
        const REAL8 m1, /**< Mass of body 1, in Msol */
        const REAL8 m2, /**< Mass of body 2, in Msol */
        const REAL8 e_min, /**< Initial eccentricity at f0 */
        const REAL8 f_min /**< Initial frequency */
        )
{
    const REAL8 e0 = e_min;
    const REAL8 f0 = f_min;

    const REAL8 mtot = m1 + m2;
    const REAL8 eta = m1*m2/mtot/mtot;
    const REAL8 m_sec = mtot * LAL_MTSUN_SI;  /* total mass in seconds */
    const REAL8 piM = LAL_PI * m_sec;
    const REAL8 v0 = cbrt(piM*f0);

    const REAL8 eta2 = eta * eta;
    const REAL8 eta3 = eta2 * eta;
    const REAL8 v02 = v0 * v0;
    const REAL8 v03 = v02 * v0;
    const REAL8 v04 = v02 * v02;
    const REAL8 v05 = v04 * v0;
    const REAL8 v06 = v03 * v03;

    const REAL8 pfeN = 3.L/(128.L * eta) * e0 * e0 * pow(v0, 19.L/3.L);

    memset(pfv19by3, 0, sizeof(PNPhasingSeries));
    memset(pfv25by3, 0, sizeof(PNPhasingSeries));
    memset(pfv28by3, 0, sizeof(PNPhasingSeries));
    memset(pfv31by3, 0, sizeof(PNPhasingSeries));
    memset(pfv34by3, 0, sizeof(PNPhasingSeries));
    memset(pfv37by3, 0, sizeof(PNPhasingSeries));

    /* Leading-order eccentric corrections to non-spin phasing terms */

    /* v^(-19/3) coefficients at every phase order */

    pfv19by3->v[0] = -2355.L/1462.L;
    pfv19by3->v[2] = -2045665.L/348096.L - (128365.L/12432.L) * eta;
    pfv19by3->v[3] = (65561.L/4080.L) * LAL_PI;
    pfv19by3->v[4] = -111064865.L/14141952.L - 165068815.L/4124736.L * eta - 10688155.L/294624.L * eta2;
    pfv19by3->v[5] = (15803101.L/229824.L * eta + 3873451.L/100548.L) * LAL_PI;
    pfv19by3->vlogv[6] = -734341.L/16800.L;
    pfv19by3->v[6] = (-21508213.L/276480.L + 103115.L/6144.L * eta) * LAL_PI * LAL_PI 
			- 409265200567.L/585252864.L * eta - 4726688461.L/34836480.L * eta2 
			- 9663919.L/50400.L * log(2.) - 69237581.L/746496.L * eta3 
			- 734341.L/16800.L * LAL_GAMMA + 4602177.L/44800.L * log(3.) 
			+ 59648637301056877.L/112661176320000.L;

    /* v^(-25/3) coefficients at every phase order */

    pfv25by3->v[2] = ((154645.L/17544.L) * eta - 2223905.L/491232.L) * v02;
    pfv25by3->v[4] = (-5795368945.L/350880768.L + 25287905.L/447552.L * eta2 + 4917245.L/1566432.L * eta) * v02;
    pfv25by3->v[5] = (-12915517.L/146880.L * eta + 185734313.L/4112640.L) * v02 * LAL_PI;
    pfv25by3->v[6] = (-314646762545.L/14255087616.L + 11585856665.L/98993664.L * eta2 
			- 1733730575525.L/24946403328.L * eta + 2105566535.L/10606464.L * eta3) * v02;

    /* v^(-28/3) coefficients at every phase order */

    pfv28by3->v[3] = -(295945.L/35088.L) * LAL_PI * v03;
    pfv28by3->v[5] = (-771215705.L/25062912.L - 48393605.L/895104.L * eta) * v03 * LAL_PI;
    pfv28by3->v[6] = 24716497.L/293760.L * LAL_PI * LAL_PI * v03;

    /* v^(-31/3) coefficients at every phase order */

    pfv31by3->v[4] = (2133117395.L/1485485568.L - 14251675.L/631584.L * eta2 + 6623045.L/260064.L * eta) * v04;
    pfv31by3->v[6] = (-2330466575.L/16111872.L * eta3 + 5558781650755.L/1061063442432.L 
			+ 12168668755.L/150377472.L * eta2 + 3869704471075.L/37895122944.L * eta) * v04;

    /* v^(-34/3) coefficients at every phase order */

    pfv34by3->v[5] = (149064749.L/2210544.L * eta - 7063901.L/520128.L) * v05 * LAL_PI;

    /* v^(-37/3) coefficients at every phase order */

    pfv37by3->v[6] = (-96423905.L/5052672.L + 748105.L/561408.L * eta) * v06 * LAL_PI * LAL_PI 
			+ (2603845.L/61404.L * log(v0) + 2425890995.L/68211072.L * eta3 
			+ 2101401811655.L/53477480448.L * eta + 2603845.L/61404.L * LAL_GAMMA 
			- 50292983575.L/636636672.L * eta2 - 4101692791967447.L/16471063977984.L 
			+ 1898287.L/184212.L * log(2.) + 12246471.L/163744.L * log(3.)) * v06;

    /* At the very end, multiply everything in the series by pfeN */

     for(int ii = 0; ii <= PN_PHASING_SERIES_MAX_ORDER; ii++)
     {
         pfv19by3->v[ii] *= pfeN;
         pfv19by3->vlogv[ii] *= pfeN;
         pfv19by3->vlogvsq[ii] *= pfeN;
     }

     for(int ii = 0; ii <= PN_PHASING_SERIES_MAX_ORDER; ii++)
     {
         pfv25by3->v[ii] *= pfeN;
         pfv25by3->vlogv[ii] *= pfeN;
         pfv25by3->vlogvsq[ii] *= pfeN;
     }

     for(int ii = 0; ii <= PN_PHASING_SERIES_MAX_ORDER; ii++)
     {
         pfv28by3->v[ii] *= pfeN;
         pfv28by3->vlogv[ii] *= pfeN;
         pfv28by3->vlogvsq[ii] *= pfeN;
     }

     for(int ii = 0; ii <= PN_PHASING_SERIES_MAX_ORDER; ii++)
     {
         pfv31by3->v[ii] *= pfeN;
         pfv31by3->vlogv[ii] *= pfeN;
         pfv31by3->vlogvsq[ii] *= pfeN;
     }

     for(int ii = 0; ii <= PN_PHASING_SERIES_MAX_ORDER; ii++)
     {
         pfv34by3->v[ii] *= pfeN;
         pfv34by3->vlogv[ii] *= pfeN;
         pfv34by3->vlogvsq[ii] *= pfeN;
     }

     for(int ii = 0; ii <= PN_PHASING_SERIES_MAX_ORDER; ii++)
     {
         pfv37by3->v[ii] *= pfeN;
         pfv37by3->vlogv[ii] *= pfeN;
         pfv37by3->vlogvsq[ii] *= pfeN;
     }
}

/** \brief Returns structure containing phasing coefficients of `eccentric` TaylorF2 for given
 *  *  physical parameters.
 *   */
int XLALSimInspiralEccTF2AlignedPhasing(
        PNPhasingSeries **pn,   /**< circular phasing coefficients (output) */
	PNPhasingSeries **pfv19by3, /**< eccentric phasing coefficients (output) */
        PNPhasingSeries **pfv25by3, /**< eccentric phasing coefficients (output) */
        PNPhasingSeries **pfv28by3, /**< eccentric phasing coefficients (output) */
        PNPhasingSeries **pfv31by3, /**< eccentric phasing coefficients (output) */
        PNPhasingSeries **pfv34by3, /**< eccentric phasing coefficients (output) */
        PNPhasingSeries **pfv37by3, /**< eccentric phasing coefficients (output) */
        const REAL8 m1,         /**< mass of body 1 */
        const REAL8 m2,         /**< mass of body 2 */
        const REAL8 chi1,       /**< aligned spin parameter of body 1 */
        const REAL8 chi2,       /**< aligned spin parameter of body 2 */
        const REAL8 qm_def1,    /**< quadrupole-monopole parameter of body 1 (set 1 for BH) */
        const REAL8 qm_def2,    /**< quadrupole-monopole parameter of body 2 (set 1 for BH) */
	const REAL8 e_min, /**< Initial eccentricity at f0 */
        const REAL8 f_min, /**< Initial frequency */
        const LALSimInspiralSpinOrder spinO  /**< PN order for spin contributions */
        )
{
    PNPhasingSeries *pfcir;
    PNPhasingSeries *pfeccA;
    PNPhasingSeries *pfeccB;
    PNPhasingSeries *pfeccC;
    PNPhasingSeries *pfeccD;
    PNPhasingSeries *pfeccE;
    PNPhasingSeries *pfeccF;

    if (!pn) XLAL_ERROR(XLAL_EFAULT);
    if (*pn) XLAL_ERROR(XLAL_EFAULT);
    if (!pfv19by3) XLAL_ERROR(XLAL_EFAULT);
    if (*pfv19by3) XLAL_ERROR(XLAL_EFAULT);
    if (!pfv25by3) XLAL_ERROR(XLAL_EFAULT);
    if (*pfv25by3) XLAL_ERROR(XLAL_EFAULT);
    if (!pfv28by3) XLAL_ERROR(XLAL_EFAULT);
    if (*pfv28by3) XLAL_ERROR(XLAL_EFAULT);
    if (!pfv31by3) XLAL_ERROR(XLAL_EFAULT);
    if (*pfv31by3) XLAL_ERROR(XLAL_EFAULT);
    if (!pfv34by3) XLAL_ERROR(XLAL_EFAULT);
    if (*pfv34by3) XLAL_ERROR(XLAL_EFAULT);
    if (!pfv37by3) XLAL_ERROR(XLAL_EFAULT);
    if (*pfv37by3) XLAL_ERROR(XLAL_EFAULT);

    pfcir = (PNPhasingSeries *) LALMalloc(sizeof(PNPhasingSeries));
    pfeccA = (PNPhasingSeries *) LALMalloc(sizeof(PNPhasingSeries));
    pfeccB = (PNPhasingSeries *) LALMalloc(sizeof(PNPhasingSeries));
    pfeccC = (PNPhasingSeries *) LALMalloc(sizeof(PNPhasingSeries));
    pfeccD = (PNPhasingSeries *) LALMalloc(sizeof(PNPhasingSeries));
    pfeccE = (PNPhasingSeries *) LALMalloc(sizeof(PNPhasingSeries));
    pfeccF = (PNPhasingSeries *) LALMalloc(sizeof(PNPhasingSeries));

    XLALSimInspiralPNPhasing_F2(pfcir, m1, m2, chi1, chi2, chi1*chi1, chi2*chi2, chi1*chi2, qm_def1, qm_def2, spinO);
    XLALSimInspiralPNPhasing_eccF2(pfeccA, pfeccB, pfeccC, pfeccD, pfeccE, pfeccF, m1, m2, e_min, f_min);

    *pn = pfcir;
    *pfv19by3 = pfeccA;
    *pfv25by3 = pfeccB;
    *pfv28by3 = pfeccC;
    *pfv31by3 = pfeccD;
    *pfv34by3 = pfeccE;
    *pfv37by3 = pfeccF;
    

    return XLAL_SUCCESS;
}

int XLALSimInspiralEccTF2Core(
        COMPLEX16FrequencySeries **htilde_out, /**< FD waveform */
	const REAL8Sequence *freqs,            /**< frequency points at which to evaluate the waveform (Hz) */
        const REAL8 phi_ref,                   /**< reference orbital phase (rad) */
        const REAL8 m1_SI,                     /**< mass of companion 1 (kg) */
        const REAL8 m2_SI,                     /**< mass of companion 2 (kg) */
        const REAL8 S1z,                       /**<  z component of the spin of companion 1 */
        const REAL8 S2z,                       /**<  z component of the spin of companion 2  */
        const REAL8 f_ref,                     /**< Reference GW frequency (Hz) - if 0 reference point is coalescence */
	const REAL8 shft,		       /**< time shift to be applied to frequency-domain phase (sec)*/
        const REAL8 r,                         /**< distance of source (m) */
        const REAL8 quadparam1,                /**< quadrupole deformation parameter of body 1 (dimensionless, 1 for BH) */
        const REAL8 quadparam2,                /**< quadrupole deformation parameter of body 2 (dimensionless, 1 for BH) */
        const REAL8 lambda1,                   /**< (tidal deformation of body 1)/(mass of body 1)^5 */
        const REAL8 lambda2,                   /**< (tidal deformation of body 2)/(mass of body 2)^5 */
	const REAL8 f_min,		       /**< initial GW frequency in Hz */
	const REAL8 e_min,		       /**< initial eccentricity at f_min */
        const LALSimInspiralSpinOrder spinO,  /**< twice PN order of spin effects */
        const LALSimInspiralTidalOrder tideO,  /**< flag to control tidal effects */
        const INT4 phaseO,                     /**< twice PN phase order */
        const INT4 amplitudeO                  /**< twice PN amplitude order */
        )
{

    if (!htilde_out) XLAL_ERROR(XLAL_EFAULT);
    if (!freqs) XLAL_ERROR(XLAL_EFAULT);
    /* external: SI; internal: solar masses */
    const REAL8 m1 = m1_SI / LAL_MSUN_SI;
    const REAL8 m2 = m2_SI / LAL_MSUN_SI;
    const REAL8 m = m1 + m2;
    const REAL8 m_sec = m * LAL_MTSUN_SI;  /* total mass in seconds */
    const REAL8 eta = m1 * m2 / (m * m);
    const REAL8 piM = LAL_PI * m_sec;
    const REAL8 m1OverM = m1 / m;
    const REAL8 m2OverM = m2 / m;
    REAL8 amp0;
    size_t i;
    COMPLEX16 *data = NULL;
    LIGOTimeGPS tC = {0, 0};
    INT4 iStart = 0;

    COMPLEX16FrequencySeries *htilde = NULL;

    if (*htilde_out) { //case when htilde_out has been allocated in XLALSimInspiralTaylorF2
	    htilde = *htilde_out;
	    iStart = htilde->data->length - freqs->length; //index shift to fill pre-allocated data
	    if(iStart < 0) XLAL_ERROR(XLAL_EFAULT);
    }
    else { //otherwise allocate memory here
	    htilde = XLALCreateCOMPLEX16FrequencySeries("htilde: FD waveform", &tC, freqs->data[0], 0., &lalStrainUnit, freqs->length);
	    if (!htilde) XLAL_ERROR(XLAL_EFUNC);
	    XLALUnitMultiply(&htilde->sampleUnits, &htilde->sampleUnits, &lalSecondUnit);
    }

    /* circular phasing coefficients */
    PNPhasingSeries pfa;
    XLALSimInspiralPNPhasing_F2(&pfa, m1, m2, S1z, S2z, S1z*S1z, S2z*S2z, S1z*S2z, quadparam1, quadparam2, spinO);

    REAL8 pfaN = 0.;
    REAL8 pfa2 = 0.; REAL8 pfa3 = 0.; REAL8 pfa4 = 0.;
    REAL8 pfa5 = 0.; REAL8 pfl5 = 0.;
    REAL8 pfa6 = 0.; REAL8 pfl6 = 0.;
    REAL8 pfa7 = 0.;

    switch (phaseO)
    {
        case -1:
        case 7:
            pfa7 = pfa.v[7];
        case 6:
            pfa6 = pfa.v[6];
            pfl6 = pfa.vlogv[6];
        case 5:
            pfa5 = pfa.v[5];
            pfl5 = pfa.vlogv[5];
        case 4:
            pfa4 = pfa.v[4];
        case 3:
            pfa3 = pfa.v[3];
        case 2:
            pfa2 = pfa.v[2];
        case 0:
            pfaN = pfa.v[0];
            break;
	case 1:
	    XLALPrintWarning( "There is no 0.5PN phase coefficient, returning Newtonian-order phase.\n" );
            pfaN = pfa.v[0];
            break;
        default:
            XLAL_ERROR(XLAL_ETYPE, "Invalid phase PN order %d", phaseO);
    }
    
    /* phasing coefficients of leading-order ecentricity corrections */
    PNPhasingSeries pfv19by3, pfv25by3, pfv28by3, pfv31by3, pfv34by3, pfv37by3;
    XLALSimInspiralPNPhasing_eccF2(&pfv19by3, &pfv25by3, &pfv28by3, &pfv31by3, &pfv34by3, &pfv37by3, m1, m2, e_min, f_min);

    REAL8 pfeNv19by3 = 0.; 
    REAL8 pfe2v19by3 = 0.;	REAL8 pfe2v25by3 = 0.;
    REAL8 pfe3v19by3 = 0.;				REAL8 pfe3v28by3 = 0.;	
    REAL8 pfe4v19by3 = 0.;	REAL8 pfe4v25by3 = 0.;				REAL8 pfe4v31by3 = 0.;	
    REAL8 pfe5v19by3 = 0.;	REAL8 pfe5v25by3 = 0.;	REAL8 pfe5v28by3 = 0.;				REAL8 pfe5v34by3 = 0.;	
    REAL8 pfe6v19by3 = 0.;	REAL8 pfe6v25by3 = 0.;	REAL8 pfe6v28by3 = 0.;	REAL8 pfe6v31by3 = 0.;				REAL8 pfe6v37by3 = 0.;
    REAL8 pfl6v19by3 = 0.;
    REAL8 pfe7 = 0.;

    switch (phaseO)
    {
        case -1:
        case 7:
            pfe7 = 0.;
        case 6:
            pfe6v19by3 = pfv19by3.v[6];
            pfl6v19by3 = pfv19by3.vlogv[6];
	    pfe6v25by3 = pfv25by3.v[6];
	    pfe6v28by3 = pfv28by3.v[6];
	    pfe6v31by3 = pfv31by3.v[6];
	    pfe6v37by3 = pfv37by3.v[6];
        case 5:
            pfe5v19by3 = pfv19by3.v[5];
            pfe5v25by3 = pfv25by3.v[5];
	    pfe5v28by3 = pfv28by3.v[5];
	    pfe5v34by3 = pfv34by3.v[5];
        case 4:
            pfe4v19by3 = pfv19by3.v[4];
	    pfe4v25by3 = pfv25by3.v[4];
	    pfe4v31by3 = pfv31by3.v[4];
        case 3:
            pfe3v19by3 = pfv19by3.v[3];
	    pfe3v28by3 = pfv28by3.v[3];
        case 2:
            pfe2v19by3 = pfv19by3.v[2];
	    pfe2v25by3 = pfv25by3.v[2];
        case 0:
            pfeNv19by3 = pfv19by3.v[0];
            break;
	case 1:
	    XLALPrintWarning( "There is no 0.5PN phase coefficient, returning Newtonian-order phase.\n" );
            pfeNv19by3 = pfv19by3.v[0];
            break;
        default:
            XLAL_ERROR(XLAL_ETYPE, "Invalid phase PN order %d", phaseO);
    }
    

    /* Validate expansion order arguments.
     * This must be done here instead of in the OpenMP parallel loop
     * because when OpenMP parallelization is turned on, early exits
     * from loops (via return or break statements) are not permitted.
     */

    /* Validate amplitude PN order. */
    switch (amplitudeO)
    {
        case -1:
        case 7:
        case 6:
        case 5:
        case 4:
        case 3:
        case 2:
        case 0:
            break;
        default:
            XLAL_ERROR(XLAL_ETYPE, "Invalid amplitude PN order %d", amplitudeO);
    }

    /* Generate tidal terms separately.
     * Enums specifying tidal order are in LALSimInspiralWaveformFlags.h
     */
    REAL8 pft10 = 0.;
    REAL8 pft12 = 0.;
    switch( tideO )
    {
	    case LAL_SIM_INSPIRAL_TIDAL_ORDER_ALL:
        case LAL_SIM_INSPIRAL_TIDAL_ORDER_6PN:
	    pft12 = pfaN * (lambda1*XLALSimInspiralTaylorF2Phasing_12PNTidalCoeff(m1OverM) + lambda2*XLALSimInspiralTaylorF2Phasing_12PNTidalCoeff(m2OverM) );
        case LAL_SIM_INSPIRAL_TIDAL_ORDER_5PN:
            pft10 = pfaN * ( lambda1*XLALSimInspiralTaylorF2Phasing_10PNTidalCoeff(m1OverM) + lambda2*XLALSimInspiralTaylorF2Phasing_10PNTidalCoeff(m2OverM) );
        case LAL_SIM_INSPIRAL_TIDAL_ORDER_0PN:
            break;
        default:
            XLAL_ERROR(XLAL_EINVAL, "Invalid tidal PN order %d", tideO);
    }

    /* The flux and energy coefficients below are used to compute SPA amplitude corrections */

    /* flux coefficients */
    const REAL8 FTaN = XLALSimInspiralPNFlux_0PNCoeff(eta);
    const REAL8 FTa2 = XLALSimInspiralPNFlux_2PNCoeff(eta);
    const REAL8 FTa3 = XLALSimInspiralPNFlux_3PNCoeff(eta);
    const REAL8 FTa4 = XLALSimInspiralPNFlux_4PNCoeff(eta);
    const REAL8 FTa5 = XLALSimInspiralPNFlux_5PNCoeff(eta);
    const REAL8 FTl6 = XLALSimInspiralPNFlux_6PNLogCoeff(eta);
    const REAL8 FTa6 = XLALSimInspiralPNFlux_6PNCoeff(eta);
    const REAL8 FTa7 = XLALSimInspiralPNFlux_7PNCoeff(eta);

    /* energy coefficients */
    const REAL8 dETaN = 2. * XLALSimInspiralPNEnergy_0PNCoeff(eta);
    const REAL8 dETa1 = 2. * XLALSimInspiralPNEnergy_2PNCoeff(eta);
    const REAL8 dETa2 = 3. * XLALSimInspiralPNEnergy_4PNCoeff(eta);
    const REAL8 dETa3 = 4. * XLALSimInspiralPNEnergy_6PNCoeff(eta);


    /* Perform some initial checks */
    if (m1_SI <= 0) XLAL_ERROR(XLAL_EDOM);
    if (m2_SI <= 0) XLAL_ERROR(XLAL_EDOM);
    if (f_ref < 0) XLAL_ERROR(XLAL_EDOM);
    if (r <= 0) XLAL_ERROR(XLAL_EDOM);

    /* extrinsic parameters */
    amp0 = -4. * m1 * m2 / r * LAL_MRSUN_SI * LAL_MTSUN_SI * sqrt(LAL_PI/12.L);

    data = htilde->data->data;

    /* Compute the SPA phase at the reference point
     * N.B. f_ref == 0 means we define the reference time/phase at "coalescence"
     * when the frequency approaches infinity. In that case,
     * the integrals Eq. 3.15 of arXiv:0907.0700 vanish when evaluated at
     * f_ref == infinity. If f_ref is finite, we must compute the SPA phase
     * evaluated at f_ref, store it as ref_phasing and subtract it off.
     */
    REAL8 ref_phasing = 0.;
    if( f_ref != 0. ) {
        const REAL8 vref = cbrt(piM*f_ref);
        const REAL8 logvref = log(vref);
        const REAL8 v2ref = vref * vref;
        const REAL8 v3ref = vref * v2ref;
        const REAL8 v4ref = vref * v3ref;
        const REAL8 v5ref = vref * v4ref;
        const REAL8 v6ref = vref * v5ref;
        const REAL8 v7ref = vref * v6ref;
        const REAL8 v8ref = vref * v7ref;
        const REAL8 v9ref = vref * v8ref;
        const REAL8 v10ref = vref * v9ref;
        const REAL8 v12ref = v2ref * v10ref;
	/* orders of v in leading-oder eccentricity corrections */
	const REAL8 v19by3ref = pow(vref, 19./3.);
	const REAL8 v25by3ref = v19by3ref * v2ref;
	const REAL8 v28by3ref = v19by3ref * v3ref;
	const REAL8 v31by3ref = v19by3ref * v4ref;
	const REAL8 v34by3ref = v19by3ref * v5ref;
	const REAL8 v37by3ref = v19by3ref * v6ref;
        ref_phasing += pfa7 * v7ref;
	ref_phasing += pfe7 * v7ref; /* 3.5PN e0^2 corrections (0) */
        ref_phasing += (pfa6 + pfl6 * logvref) * v6ref;
	ref_phasing += (pfe6v19by3 / v19by3ref + pfe6v25by3 / v25by3ref + pfe6v28by3 / v28by3ref 
			+ pfe6v31by3 / v31by3ref + pfe6v37by3 / v37by3ref 
			+ pfl6v19by3 / v19by3ref * logvref) * v6ref; /* 3PN e0^2 corrections */
        ref_phasing += (pfa5 + pfl5 * logvref) * v5ref;
	ref_phasing += (pfe5v19by3 / v19by3ref + pfe5v25by3 / v25by3ref 
			+ pfe5v28by3 / v28by3ref + pfe5v34by3 / v34by3ref) * v5ref; /* 2.5PN e0^2 corrections */
        ref_phasing += pfa4 * v4ref;
	ref_phasing += (pfe4v19by3 / v19by3ref + pfe4v25by3 / v25by3ref 
			+ pfe4v31by3 / v31by3ref) * v4ref; /* 2PN e0^2 corrections */
        ref_phasing += pfa3 * v3ref;
	ref_phasing += (pfe3v19by3 / v19by3ref + pfe3v28by3 / v28by3ref) * v3ref; /* 1.5PN e0^2 corrections */
        ref_phasing += pfa2 * v2ref;
	ref_phasing += (pfe2v19by3 / v19by3ref + pfe2v25by3 / v25by3ref) * v2ref; /* 1PN e0^2 corrections */
        ref_phasing += pfaN;
	ref_phasing += pfeNv19by3 / v19by3ref; /* Newt. e0^2 corrections */

        /* Tidal terms in reference phasing */
        ref_phasing += pft12 * v12ref;
        ref_phasing += pft10 * v10ref;

        ref_phasing /= v5ref;
    } /* End of if(f_ref != 0) block */

    #pragma omp parallel for
    for (i = 0; i < freqs->length; i++) {
        const REAL8 f = freqs->data[i];
        const REAL8 v = cbrt(piM*f);
        const REAL8 logv = log(v);
        const REAL8 v2 = v * v;
        const REAL8 v3 = v * v2;
        const REAL8 v4 = v * v3;
        const REAL8 v5 = v * v4;
        const REAL8 v6 = v * v5;
        const REAL8 v7 = v * v6;
        const REAL8 v8 = v * v7;
        const REAL8 v9 = v * v8;
        const REAL8 v10 = v * v9;
        const REAL8 v12 = v2 * v10;
	/* orders of v in leading-order eccentricity corrections */
	const REAL8 v19by3 = pow(v, 19./3.);
	const REAL8 v25by3 = v19by3 * v2;
	const REAL8 v28by3 = v19by3 * v3;
	const REAL8 v31by3 = v19by3 * v4;
	const REAL8 v34by3 = v19by3 * v5;
	const REAL8 v37by3 = v19by3 * v6;
        REAL8 phasing = 0.;
        REAL8 dEnergy = 0.;
        REAL8 flux = 0.;
        REAL8 amp;

        phasing += pfa7 * v7;
	phasing += pfe7 * v7; /* 3.5PN e0^2 corrections (0) */
        phasing += (pfa6 + pfl6 * logv) * v6;
	phasing += (pfe6v19by3 / v19by3 + pfe6v25by3 / v25by3 + pfe6v28by3 / v28by3 
		    + pfe6v31by3 / v31by3 + pfe6v37by3 / v37by3 
		    + pfl6v19by3 / v19by3 * logv) * v6; /* 3PN e0^2 corrections */
        phasing += (pfa5 + pfl5 * logv) * v5;
	phasing += (pfe5v19by3 / v19by3 + pfe5v25by3 / v25by3
		    + pfe5v28by3 / v28by3 + pfe5v34by3 / v34by3) * v5; /* 2.5PN e0^2 corrections */
        phasing += pfa4 * v4;
	phasing += (pfe4v19by3 / v19by3 + pfe4v25by3 / v25by3 
		    + pfe4v31by3 / v31by3) * v4; /* 2PN e0^2 corrections */
        phasing += pfa3 * v3;
	phasing += (pfe3v19by3 / v19by3 + pfe3v28by3 / v28by3) * v3; /* 1.5PN e0^2 corrections */
        phasing += pfa2 * v2;
	phasing += (pfe2v19by3 / v19by3 + pfe2v25by3 / v25by3) * v2; /* 1PN e0^2 corrections */
        phasing += pfaN;
	phasing += pfeNv19by3 / v19by3; /* Newt. e0^2 corrections */

        /* Tidal terms in phasing */
        phasing += pft12 * v12;
        phasing += pft10 * v10;


    /* WARNING! Amplitude orders beyond 0 have NOT been reviewed!
     * Use at your own risk. The default is to turn them off.
     * These do not currently include spin corrections.
     * Note that these are not higher PN corrections to the amplitude.
     * They are the corrections to the leading-order amplitude arising
     * from the stationary phase approximation. See for instance
     * Eq 6.9 of arXiv:0810.5336
     */
	switch (amplitudeO)
        {
            case 7:
                flux += FTa7 * v7;
            case 6:
                flux += (FTa6 + FTl6*logv) * v6;
                dEnergy += dETa3 * v6;
            case 5:
                flux += FTa5 * v5;
            case 4:
                flux += FTa4 * v4;
                dEnergy += dETa2 * v4;
            case 3:
                flux += FTa3 * v3;
            case 2:
                flux += FTa2 * v2;
                dEnergy += dETa1 * v2;
            case -1: /* Default to no SPA amplitude corrections */
            case 0:
                flux += 1.;
                dEnergy += 1.;
        }

        phasing /= v5;
        flux *= FTaN * v10;
        dEnergy *= dETaN * v;
        // Note the factor of 2 b/c phi_ref is orbital phase
        phasing += shft * f - 2.*phi_ref - ref_phasing;
        amp = amp0 * sqrt(-dEnergy/flux) * v;
        data[i+iStart] = amp * cos(phasing - LAL_PI_4)
                - amp * sin(phasing - LAL_PI_4) * 1.0j;
    }

    *htilde_out = htilde;
    return XLAL_SUCCESS;
}

/**
 * Computes the stationary phase approximation to the Fourier transform of
 * a chirp waveform. The amplitude is given by expanding \f$1/\sqrt{\dot{F}}\f$.
 * If the PN order is set to -1, then the highest implemented order is used.
 *
 * @note f_ref is the GW frequency at which phi_ref is defined. The most common
 * choice in the literature is to choose the reference point as "coalescence",
 * when the frequency becomes infinite. This is the behavior of the code when
 * f_ref==0. If f_ref > 0, phi_ref sets the orbital phase at that GW frequency.
 *
 * See arXiv:0810.5336 and arXiv:astro-ph/0504538 for spin corrections
 * to the phasing.
 * See arXiv:1303.7412 for spin-orbit phasing corrections at 3 and 3.5PN order
 *
 * See LIGO-P1500148-v1 for eccentricity corrections 
 *
 * The spin and tidal order enums are defined in LALSimInspiralWaveformFlags.h
 */
int XLALSimInspiralEccTF2(
        COMPLEX16FrequencySeries **htilde_out, /**< FD waveform */
        const REAL8 phi_ref,                   /**< reference orbital phase (rad) */
        const REAL8 deltaF,                    /**< frequency resolution */
        const REAL8 m1_SI,                     /**< mass of companion 1 (kg) */
        const REAL8 m2_SI,                     /**< mass of companion 2 (kg) */
        const REAL8 S1z,                       /**<  z component of the spin of companion 1 */
        const REAL8 S2z,                       /**<  z component of the spin of companion 2  */
        const REAL8 fStart,                    /**< start GW frequency (Hz) */
        const REAL8 fEnd,                      /**< highest GW frequency (Hz) of waveform generation - if 0, end at Schwarzschild ISCO */
        const REAL8 f_ref,                     /**< Reference GW frequency (Hz) - if 0 reference point is coalescence */
        const REAL8 r,                         /**< distance of source (m) */
        const REAL8 quadparam1,                /**< quadrupole deformation parameter of body 1 (dimensionless, 1 for BH) */
        const REAL8 quadparam2,                /**< quadrupole deformation parameter of body 2 (dimensionless, 1 for BH) */
        const REAL8 lambda1,                   /**< (tidal deformation of body 1)/(mass of body 1)^5 */
        const REAL8 lambda2,                   /**< (tidal deformation of body 2)/(mass of body 2)^5 */
        const REAL8 e_min,                     /**< initial eccentricity at f_min */
        const LALSimInspiralSpinOrder spinO,  /**< twice PN order of spin effects */
        const LALSimInspiralTidalOrder tideO,  /**< flag to control tidal effects */
        const INT4 phaseO,                     /**< twice PN phase order */
        const INT4 amplitudeO                  /**< twice PN amplitude order */
        )
{
    /* external: SI; internal: solar masses */
    const REAL8 m1 = m1_SI / LAL_MSUN_SI;
    const REAL8 m2 = m2_SI / LAL_MSUN_SI;
    const REAL8 m = m1 + m2;
    const REAL8 m_sec = m * LAL_MTSUN_SI;  /* total mass in seconds */
    // const REAL8 eta = m1 * m2 / (m * m);
    const REAL8 piM = LAL_PI * m_sec;
    const REAL8 vISCO = 1. / sqrt(6.);
    const REAL8 fISCO = vISCO * vISCO * vISCO / piM;
    //const REAL8 m1OverM = m1 / m;
    // const REAL8 m2OverM = m2 / m;
    REAL8 shft, f_max;
    size_t i, n;
    INT4 iStart;
    REAL8Sequence *freqs = NULL;
    LIGOTimeGPS tC = {0, 0};
    int ret;

    COMPLEX16FrequencySeries *htilde = NULL;

    /* Perform some initial checks */
    if (!htilde_out) XLAL_ERROR(XLAL_EFAULT);
    if (*htilde_out) XLAL_ERROR(XLAL_EFAULT);
    if (m1_SI <= 0) XLAL_ERROR(XLAL_EDOM);
    if (m2_SI <= 0) XLAL_ERROR(XLAL_EDOM);
    if (fStart <= 0) XLAL_ERROR(XLAL_EDOM);
    if (f_ref < 0) XLAL_ERROR(XLAL_EDOM);
    if (r <= 0) XLAL_ERROR(XLAL_EDOM);

    /* allocate htilde */
    if ( fEnd == 0. ) // End at ISCO
        f_max = fISCO;
    else // End at user-specified freq.
        f_max = fEnd;
    if (f_max <= fStart) XLAL_ERROR(XLAL_EDOM);

    n = (size_t) (f_max / deltaF + 1);
    XLALGPSAdd(&tC, -1 / deltaF);  /* coalesce at t=0 */
    htilde = XLALCreateCOMPLEX16FrequencySeries("htilde: FD waveform", &tC, 0.0, deltaF, &lalStrainUnit, n);
    if (!htilde) XLAL_ERROR(XLAL_EFUNC);
    memset(htilde->data->data, 0, n * sizeof(COMPLEX16));
    XLALUnitMultiply(&htilde->sampleUnits, &htilde->sampleUnits, &lalSecondUnit);

    /* Fill with non-zero vals from fStart to f_max */
    iStart = (INT4) ceil(fStart / deltaF);

    /* Sequence of frequencies where waveform model is to be evaluated */
    freqs = XLALCreateREAL8Sequence(n - iStart);

    /* extrinsic parameters */
    shft = LAL_TWOPI * (tC.gpsSeconds + 1e-9 * tC.gpsNanoSeconds);

    #pragma omp parallel for
    for (i = iStart; i < n; i++) {
        freqs->data[i-iStart] = i * deltaF;
    }
    ret = XLALSimInspiralEccTF2Core(&htilde, freqs, phi_ref, m1_SI, m2_SI,
                                      S1z, S2z, f_ref, shft, r, quadparam1, quadparam2,
                                      lambda1, lambda2, fStart, e_min, spinO, tideO, phaseO, amplitudeO);

    XLALDestroyREAL8Sequence(freqs);

    *htilde_out = htilde;

    return ret;
}

/** @} */
/** @} */
