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

#ifndef _LALSIMINSPIRALTESTINGGRCORRECTIONS_H
#define _LALSIMINSPIRALTESTINGGRCORRECTIONS_H

#if defined(__cplusplus)
extern "C" {
#elif 0
} /* so that editors will match preceding brace */
#endif

#include <stdlib.h>
#include <math.h>
#include <lal/Date.h>
#include <lal/FrequencySeries.h>
#include <lal/LALConstants.h>
#include <lal/LALDatatypes.h>
#include <lal/LALSimInspiral.h>

int XLALSimInspiralTestingGRCorrections(COMPLEX16FrequencySeries *htilde,       /**< input htilde, will be modified in place */
                                        const REAL8 distance,
                                        const REAL8 m1_SI,
                                        const REAL8 m2_SI,
                                        const REAL8 f_low,
                                        const REAL8 f_ref,
                                        const LALSimInspiralTestGRParam *pnCorrections    /**< input linked list of testing gr parameters */
);

void XLALSimInspiralNonSpinningPNCorrections(PNPhasingSeries *pfa, const REAL8 eta, const LALSimInspiralTestGRParam *pnCorrections);

int XLALSimInspiralPhaseCorrectionsPhasing(COMPLEX16FrequencySeries *htilde,       /**< input htilde, will be modified in place */
                                           const REAL8 distance,
                                           const REAL8Sequence *freqs,
                                           const UINT4 iStart,
                                           const UINT4 iEnd,
                                           PNPhasingSeries pfa,
                                           const REAL8 mtot,
                                           const REAL8 eta,
                                           const REAL8 f_ref); /** this must be in seconds **/

#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif /* _LALSIMINSPIRALTESTINGGRCORRECTIONS_H */
