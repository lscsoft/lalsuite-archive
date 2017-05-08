/*
 *  Copyright (C) 2017 Walter Del Pozzo
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
#include <lal/LALDatatypes.h>
#include <lal/LALSimInspiral.h>
#include <lal/LALSimInspiralTestingGRCorrections.h>
#include <lal/Sequence.h>
#include "LALSimInspiralPNCoefficients.c"

UNUSED static inline REAL8
GetNRSpinPeakOmegaV4 (INT4 UNUSED l, INT4 UNUSED m, REAL8 UNUSED eta, REAL8 a)
{
    REAL8 chi = a;
    REAL8 res;
    res = 0.5626787200433265 + (-0.08706198756945482 +
                                0.0017434519312586804 * chi) *
    log (10.26207326082448 -
         chi * (7.629921628648589 -
                72.75949266353584 * (-0.25 + eta)) -
         62.353217004599784 * (-0.25 + eta));
    //    printf("w %.16e\n", res);
    return res;
}

int XLALSimInspiralTestingGRCorrections(COMPLEX16FrequencySeries *htilde,       /**< input htilde, will be modified in place */
                                        const REAL8 m1_SI,
                                        const REAL8 m2_SI,
                                        const REAL8 f_low,
                                        const REAL8 deltaF,
                                        const LALSimInspiralTestGRParam *pnCorrections    /**< input linked list of testing gr parameters */
)
{
    /* external: SI; internal: solar masses */
    const REAL8 m1 = m1_SI / LAL_MSUN_SI;
    const REAL8 m2 = m2_SI / LAL_MSUN_SI;
    const REAL8 m = m1 + m2;
    const REAL8 m_sec = m * LAL_MTSUN_SI;  /* total mass in seconds */
    const REAL8 eta = m1 * m2 / (m * m);
    const REAL8 piM = LAL_PI * m_sec;
    const REAL8 vISCO = 1. / sqrt(6.);
    const REAL8 fISCO = vISCO * vISCO * vISCO / piM;
    const REAL8 fPeak = GetNRSpinPeakOmegaV4(2, 2, eta, 0.) / (2. * LAL_PI * m_sec);
    
    INT4 i;
    INT4 n = (INT4) htilde->data->length;
    
    INT4 iStart, iEnd;
    /* Fill with non-zero vals from fStart to fEnd */
    iStart = (INT4) ceil(f_low / deltaF);
    iEnd = (INT4) ceil(fISCO / deltaF);
    iEnd = (INT4) ceil(fPeak / deltaF);
    /* Sequence of frequencies where corrections to the model need to be evaluated */
    REAL8Sequence *freqs =NULL;
    freqs = XLALCreateREAL8Sequence(n - iStart);
    
    for (i = iStart; i < n; i++)
    {
        freqs->data[i-iStart] = i * deltaF;
    }
    PNPhasingSeries pfa;
    XLALSimInspiralNonSpinningPNCorrections(&pfa, eta, pnCorrections);
    XLALSimInspiralPhaseCorrectionsPhasing(htilde,freqs,iStart,iEnd,pfa,m_sec);
    XLALDestroyREAL8Sequence(freqs);
    return 0;
}

void XLALSimInspiralNonSpinningPNCorrections(PNPhasingSeries *pfa, const REAL8 eta, const LALSimInspiralTestGRParam *pnCorrections)
{
    const REAL8 pfaN = 3.L/(128.L * eta);
    /* initialise the PN correction  coefficients to 0 identically */
    memset(pfa, 0, sizeof(PNPhasingSeries));
    
    /* Non-spinning corrections to phasing terms - see arXiv:0907.0700, Eq. 3.18 */
    if (XLALSimInspiralTestGRParamExists(pnCorrections,"dchi0")) pfa->v[0] = XLALSimInspiralGetTestGRParam(pnCorrections,"dchi0");
    if (XLALSimInspiralTestGRParamExists(pnCorrections,"dchi1")) pfa->v[1] = XLALSimInspiralGetTestGRParam(pnCorrections,"dchi1");
    
    if (XLALSimInspiralTestGRParamExists(pnCorrections,"dchi2"))
    {
        pfa->v[2] = 5.L*(743.L/84.L + 11.L * eta)/9.L;
        pfa->v[2] *= XLALSimInspiralGetTestGRParam(pnCorrections,"dchi2");
    }
    if (XLALSimInspiralTestGRParamExists(pnCorrections,"dchi3"))
    {
        pfa->v[3] = -16.L*LAL_PI;
        pfa->v[3] *= XLALSimInspiralGetTestGRParam(pnCorrections,"dchi3");
    }
    if (XLALSimInspiralTestGRParamExists(pnCorrections,"dchi4"))
    {
        pfa->v[4] = 5.L*(3058.673L/7.056L + 5429.L/7.L * eta
                        + 617.L * eta*eta)/72.L;
        pfa->v[4] *= XLALSimInspiralGetTestGRParam(pnCorrections,"dchi4");
    }
    if (XLALSimInspiralTestGRParamExists(pnCorrections,"dchi5"))
    {
        pfa->v[5] = 5.L/9.L * (7729.L/84.L - 13.L * eta) * LAL_PI;
        pfa->v[5] *= XLALSimInspiralGetTestGRParam(pnCorrections,"dchi5");
    }
    if (XLALSimInspiralTestGRParamExists(pnCorrections,"dchi5l"))
    {
        pfa->vlogv[5] = 5.L/3.L * (7729.L/84.L - 13.L * eta) * LAL_PI;
        pfa->vlogv[5] *= XLALSimInspiralGetTestGRParam(pnCorrections,"dchi5l");
    }
    if (XLALSimInspiralTestGRParamExists(pnCorrections,"dchi6"))
    {
        pfa->v[6] = (11583.231236531L/4.694215680L
                     - 640.L/3.L * LAL_PI * LAL_PI - 6848.L/21.L*LAL_GAMMA)
        + eta * (-15737.765635L/3.048192L
                 + 2255./12. * LAL_PI * LAL_PI)
        + eta*eta * 76055.L/1728.L
        - eta*eta*eta * 127825.L/1296.L;
        pfa->v[6] += (-6848.L/21.L)*log(4.);
        pfa->v[6] *= XLALSimInspiralGetTestGRParam(pnCorrections,"dchi6");
    }
    if (XLALSimInspiralTestGRParamExists(pnCorrections,"dchi6l"))
    {
        pfa->vlogv[6] = -6848.L/21.L;
        pfa->vlogv[6] *= XLALSimInspiralGetTestGRParam(pnCorrections,"dchi6l");
    }
    if (XLALSimInspiralTestGRParamExists(pnCorrections,"dchi7"))
    {
        pfa->v[7] = LAL_PI * ( 77096675.L/254016.L
                              + 378515.L/1512.L * eta - 74045.L/756.L * eta*eta);
        pfa->v[7] *= XLALSimInspiralGetTestGRParam(pnCorrections,"dchi7");
    }
    
    /* At the very end, multiply everything in the series by the newtonian term */
    for(int ii = 0; ii <= PN_PHASING_SERIES_MAX_ORDER; ii++)
    {
        pfa->v[ii] *= pfaN;
        pfa->vlogv[ii] *= pfaN;
        pfa->vlogvsq[ii] *= pfaN;
    }
}

int XLALSimInspiralPhaseCorrectionsPhasing(COMPLEX16FrequencySeries *htilde,       /**< input htilde, will be modified in place */
                                           const REAL8Sequence *freqs,
                                           const INT4 iStart,
                                           const INT4 iEnd,
                                           PNPhasingSeries pfa,
                                           const REAL8 mtot) /** this must be in seconds **/
{
    const REAL8 piM = LAL_PI * mtot;
    const REAL8 pfa7 = pfa.v[7];
    const REAL8 pfa6 = pfa.v[6];
    const REAL8 pfl6 = pfa.vlogv[6];
    const REAL8 pfa5 = pfa.v[5];
    const REAL8 pfl5 = pfa.vlogv[5];
    const REAL8 pfa4 = pfa.v[4];
    const REAL8 pfa3 = pfa.v[3];
    const REAL8 pfa2 = pfa.v[2];
    const REAL8 pfa1 = pfa.v[1];
    const REAL8 pfaN = pfa.v[0];
    
    COMPLEX16 *data = NULL;
    /* Compute the SPA phase at the reference point
     * N.B. f_ref == 0 means we define the reference time/phase at "coalescence"
     * when the frequency approaches infinity. In that case,
     * the integrals Eq. 3.15 of arXiv:0907.0700 vanish when evaluated at
     * f_ref == infinity. If f_ref is finite, we must compute the SPA phase
     * evaluated at f_ref, store it as ref_phasing and subtract it off.
     */
    REAL8 ref_phasing = 0.;
    data = htilde->data->data;
    /*
    if( f_ref != 0. ) {
        const REAL8 vref = cbrt(piM*f_ref);
        const REAL8 logvref = log(vref);
        const REAL8 v2ref = vref * vref;
        const REAL8 v3ref = vref * v2ref;
        const REAL8 v4ref = vref * v3ref;
        const REAL8 v5ref = vref * v4ref;
        const REAL8 v6ref = vref * v5ref;
        const REAL8 v7ref = vref * v6ref;

        ref_phasing += pfa7 * v7ref;
        ref_phasing += (pfa6 + pfl6 * logvref) * v6ref;
        ref_phasing += (pfa5 + pfl5 * logvref) * v5ref;
        ref_phasing += pfa4 * v4ref;
        ref_phasing += pfa3 * v3ref;
        ref_phasing += pfa2 * v2ref;
        ref_phasing += pfa1 * vref;
        ref_phasing += pfaN;
        ref_phasing /= v5ref;
    }
    */

    const REAL8 vPeak = cbrt(piM*freqs->data[iEnd-iStart]); //Center of the tapering step function in v-space
    const INT4 NCyclesStep = 5; //Choose no. of GW cycles over which to taper the non-GR phase correction
    const REAL8 width = NCyclesStep * ( sqrt( 2./3.*piM*LAL_PI/vPeak + 0.25*vPeak*vPeak ) - 0.5*vPeak );
    for (INT4 i = iStart; i < (INT4) htilde->data->length; i++)
    {
        const REAL8 f = freqs->data[i-iStart];
        if (f>0)
        {
            const REAL8 v = cbrt(piM*f);
            const REAL8 logv = log(v);
            const REAL8 v2 = v * v;
            const REAL8 v3 = v * v2;
            const REAL8 v4 = v * v3;
            const REAL8 v5 = v * v4;
            const REAL8 v6 = v * v5;
            const REAL8 v7 = v * v6;
            REAL8 phasing = 0.0;
            
            phasing += pfa7 * v7;
            phasing += (pfa6 + pfl6 * logv) * v6;
            phasing += (pfa5 + pfl5 * logv) * v5;
            phasing += pfa4 * v4;
            phasing += pfa3 * v3;
            phasing += pfa2 * v2;
            phasing += pfa1 * v;
            phasing += pfaN;
            phasing /= v5;
            
            phasing -= ref_phasing;
            phasing *= 1. / (1. + exp( (v - vPeak) / width ) );
            data[i] *= cos(phasing)- sin(phasing) * 1.0j;
//            printf("i %d freq: %.2f data %e + %e i\n", i, f, creal(data[i+iStart]), cimag(data[i+iStart]));
        }
    }
    return 0;
}
