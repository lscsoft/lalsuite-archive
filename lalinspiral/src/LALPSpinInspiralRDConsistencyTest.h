/*
*  Copyright (C) 2011 Michalis Agathos
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

/* <lalVerbatim> */
#ifndef _LALPSPININSPIRALRDTEST_H  /* Double-include protection. */
#define _LALPSPININSPIRALRDTEST_H

#ifdef  __cplusplus   /* C++ protection. */
extern "C" {
#endif


NRCSID( LALPSPININSPIRALRDTESTH, "$Id$" );

    void LALPSpinInspiralRDConsistencyTest (
                                           LALStatus        *status,
                                           REAL4Vector      *signalvec,
                                           InspiralTemplate *params,
                                           REAL8 dphi[10]);

    void LALPSpinInspiralRDEngineCT (
                                    LALStatus	*status,
                                    REAL4Vector	*signalvec1,
                                    REAL4Vector * signalvec2,
                                    REAL4Vector * hh,
                                    REAL4Vector * ff,
                                    REAL8Vector * phi,
                                    REAL4Vector * shift,
                                    UINT4 * countback,
                                    InspiralTemplate * params,
                                    InspiralInit * paramsInit, 
                                    REAL8 dphi[10]);

    void LALPSpinInspiralRDForInjectionCT(LALStatus * status,
				    CoherentGW * waveform,
				    InspiralTemplate * params,
				    PPNParamStruc * ppnParams,
				    REAL8 dphi[10]);

    void LALPSpinInspiralRDFreqDomCT(LALStatus * status,
			       REAL4Vector * signalvec,
			       InspiralTemplate * params,
			       REAL8 dphi[10]);

    void LALPSpinInspiralRDTemplatesCT(LALStatus * status,
				 REAL4Vector * signalvec1,
				 REAL4Vector * signalvec2,
				 InspiralTemplate * params,
				 REAL8 dphi[10]);

    void LALPSpinInspiralRDderivativesCT(REAL8Vector * values,
				   REAL8Vector * dvalues, void *mparams);


#ifdef  __cplusplus
}
#endif  /* C++ protection. */
#endif  /* Double-include protection. */
/* </lalVerbatim> */
